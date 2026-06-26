import pandas as pd
import numpy as np
import re
from utils.constants import CONSEQUENCE_IMPACT, IMPACT_NUMERIC

def _numeric_column(out: pd.DataFrame, column: str, default=np.nan) -> pd.Series:
    if column not in out.columns:
        return pd.Series(default, index=out.index, dtype="float64")
    return pd.to_numeric(out[column], errors="coerce")


def _truncation_norm(hgvsp, consequence) -> float:
    if pd.isna(hgvsp) or consequence != "stop_gained":
        return 0.0

    match = re.search(r"p\.[A-Za-z*]+(\d+)(?:Ter|\*)", str(hgvsp))
    if not match:
        return 0.0

    aa_pos = float(match.group(1))
    return float(np.clip(1.0 - (aa_pos / 1200.0), 0.0, 1.0))


def rank_variants(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    popfreq_cols = [c for c in out.columns if c.startswith("dbnsfp_popfreq_")]
    numeric_inputs = {
        "hpo_match_count": _numeric_column(out, "hpo_match_count").fillna(0.0),
        "phenotype_score": _numeric_column(out, "phenotype_score").fillna(0.0),
        "hpo_overlap_ratio": _numeric_column(out, "hpo_overlap_ratio").fillna(0.0),
        "am_score": _numeric_column(out, "am_score"),
        "cadd_from_vep": _numeric_column(out, "cadd_from_vep"),
        "dbnsfp_revel_score": _numeric_column(out, "dbnsfp_revel_score"),
        "dbnsfp_cadd_phred": _numeric_column(out, "dbnsfp_cadd_phred"),
        "dbnsfp_alphamissense_score": _numeric_column(out, "dbnsfp_alphamissense_score"),
        "dbnsfp_sift_score": _numeric_column(out, "dbnsfp_sift_score"),
        "dbnsfp_polyphen2_hdiv_score": _numeric_column(out, "dbnsfp_polyphen2_hdiv_score"),
        "dbnsfp_gnomad_af": _numeric_column(out, "dbnsfp_gnomad_af"),
        "dbnsfp_gerp_rs": _numeric_column(out, "dbnsfp_gerp_rs"),
        "dbnsfp_aapos": _numeric_column(out, "dbnsfp_aapos"),
        "dbnsfp_metarnn_score": _numeric_column(out, "dbnsfp_metarnn_score"),
    }
    numeric_inputs.update({col: _numeric_column(out, col) for col in popfreq_cols})
    out = pd.concat([out.drop(columns=list(numeric_inputs.keys()), errors="ignore"), pd.DataFrame(numeric_inputs)], axis=1).copy()

    out["consequence_score"] = out["consequence"].map(CONSEQUENCE_IMPACT).fillna(0.0)
    out["impact_score"] = out["impact"].map(IMPACT_NUMERIC).fillna(0.0)

    out["hpo_norm"] = (out["hpo_match_count"] / 5.0).clip(0, 1)
    out["cons_norm"] = (out["consequence_score"] / 8.0).clip(0, 1)
    out["impact_norm"] = out["impact_score"].clip(0, 1)

    out["am_effective"] = out["dbnsfp_alphamissense_score"].combine_first(out["am_score"])
    out["cadd_effective"] = out["dbnsfp_cadd_phred"].combine_first(out["cadd_from_vep"])
    out["am_norm"] = out["am_effective"].fillna(0.5).clip(0, 1)
    out["cadd_norm"] = (out["cadd_effective"].fillna(20) / 60.0).clip(0, 1)
    out["revel_norm"] = out["dbnsfp_revel_score"].fillna(0.5).clip(0, 1)

    # SIFT is inverted: lower score means stronger predicted damage.
    out["sift_norm"] = (1 - out["dbnsfp_sift_score"]).fillna(0.5).clip(0, 1)
    out["polyphen_norm"] = out["dbnsfp_polyphen2_hdiv_score"].fillna(0.5).clip(0, 1)
    out["gerp_norm"] = ((out["dbnsfp_gerp_rs"] + 12.3) / 18.5).fillna(0.5).clip(0, 1)
    out["metarnn_norm"] = out["dbnsfp_metarnn_score"].fillna(0.5).clip(0, 1)

    af_sources = []
    if "dbnsfp_gnomad_af" in out.columns:
        af_sources.append(out["dbnsfp_gnomad_af"])
    af_sources.extend([out[col] for col in popfreq_cols])
    out["population_max_af"] = pd.NA
    if af_sources:
        out["population_max_af"] = pd.concat(af_sources, axis=1).max(axis=1, skipna=True)

    out["af_penalty"] = pd.to_numeric(out["population_max_af"], errors="coerce").fillna(0.0).clip(0, 0.05) / 0.05
    out["truncation_norm"] = out.apply(
        lambda row: _truncation_norm(
            row.get("hgvsp") if pd.notna(row.get("hgvsp")) else row.get("dbnsfp_hgvsp"),
            row.get("consequence"),
        ),
        axis=1,
    )

    if "gene_symbol" in out.columns:
        out["gene_hash"] = out["gene_symbol"].apply(
            lambda x: (sum(ord(c) for c in str(x)) % 1000) / 1000 if pd.notna(x) else 0.5
        )
    else:
        out["gene_hash"] = 0.5

    if "chrom" in out.columns and "pos" in out.columns:
        def chrom_to_num(chrom):
            chrom_str = str(chrom).replace("chr", "").upper()
            if chrom_str == "X":
                return 23
            elif chrom_str == "Y":
                return 24
            elif chrom_str in ["MT", "M"]:
                return 25
            else:
                try:
                    return float(chrom_str)
                except:
                    return 0
        
        out["chrom_num"] = out["chrom"].apply(chrom_to_num)
        out["pos_norm"] = (
            (out["chrom_num"].fillna(0) / 25.0 + 
            (out["pos"].astype(float).fillna(0) % 1000000) / 1000000.0)
        ).clip(0, 1) * 0.02  # Small weight for tie-breaking
    else:
        out["pos_norm"] = 0.0

    phenotype_block = (
        0.65 * out["phenotype_score"] +
        0.35 * out["hpo_norm"]
    )
    molecular_block = (
        0.26 * out["revel_norm"] +
        0.16 * out["am_norm"] +
        0.16 * out["cadd_norm"] +
        0.09 * out["sift_norm"] +
        0.07 * out["polyphen_norm"] +
        0.08 * out["gerp_norm"] +
        0.08 * out["metarnn_norm"] +
        0.10 * out["truncation_norm"]
    )
    consequence_block = (
        0.65 * out["cons_norm"] +
        0.35 * out["impact_norm"]
    )
    population_block = 1.0 - out["af_penalty"]

    out["phenotype_block"] = phenotype_block.clip(0, 1)
    out["molecular_block"] = molecular_block.clip(0, 1)
    out["consequence_block"] = consequence_block.clip(0, 1)
    out["population_block"] = population_block.clip(0, 1)

    out["pathogenicity_score"] = (
        0.34 * out["phenotype_block"] +
        0.28 * out["molecular_block"] +
        0.24 * out["consequence_block"] +
        0.14 * out["population_block"] +
        0.01 * out["gene_hash"] +
        0.01 * out["pos_norm"]
    ).clip(0, 1)

    return out.sort_values("pathogenicity_score", ascending=False).reset_index(drop=True)

def score_label(score: float):
    if score >= 0.70:
        return "🔴", "#CA1818"
    if score >= 0.40:
        return "🟡", "#BB6D15"
    return "⚪", "#444441"

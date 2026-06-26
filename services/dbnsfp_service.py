from io import BytesIO

import pandas as pd


KEY_ALIASES = {
    "chrom": ["chrom", "chr", "#chr", "chromosome"],
    "pos": ["pos", "pos(1-based)", "position", "start", "hg38_pos(1-based)", "hg19_pos(1-based)"],
    "ref": ["ref", "reference", "reference_allele"],
    "alt": ["alt", "alternative", "alternate", "alternative_allele"],
}

ANNOTATION_ALIASES = {
    "dbnsfp_cadd_phred": ["cadd_phred", "cadd", "dbnsfp_cadd_phred"],
    "dbnsfp_revel_score": ["revel_score", "revel", "dbnsfp_revel_score"],
    "dbnsfp_alphamissense_score": [
        "alphamissense_score",
        "alphamissense_pathogenicity",
        "am_score",
        "dbnsfp_alphamissense_score",
    ],
    "dbnsfp_alphamissense_pred": [
        "alphamissense_pred",
        "alphamissense_class",
        "am_class",
        "dbnsfp_alphamissense_pred",
    ],
    "dbnsfp_metarnn_score": ["metarnn_score", "dbnsfp_metarnn_score"],
    "dbnsfp_metarnn_pred": ["metarnn_pred", "dbnsfp_metarnn_pred"],
    "dbnsfp_sift_score": ["sift_score", "sift", "dbnsfp_sift_score"],
    "dbnsfp_sift_pred": ["sift_pred", "dbnsfp_sift_pred"],
    "dbnsfp_polyphen2_hdiv_score": [
        "polyphen2_hdiv_score",
        "polyphen2_hdiv",
        "polyphen_score",
        "dbnsfp_polyphen2_hdiv_score",
    ],
    "dbnsfp_polyphen2_hdiv_pred": [
        "polyphen2_hdiv_pred",
        "polyphen_pred",
        "dbnsfp_polyphen2_hdiv_pred",
    ],
    "dbnsfp_gnomad_af": [
        "gnomad_genomes_af",
        "gnomad_exomes_af",
        "gnomad4.1_joint_popmax_af",
        "gnomad_af",
        "popmax_af",
        "dbnsfp_popmax_af",
        "dbnsfp_gnomad_genomes_af",
        "dbnsfp_gnomad_exomes_af",
    ],
    "dbnsfp_gerp_rs": ["gerp++_rs", "gerp_rs", "dbnsfp_gerp_rs"],
    "dbnsfp_aapos": ["aapos", "aa_pos", "dbnsfp_aapos"],
    "dbnsfp_hgvsp": ["hgvsp_vep", "hgvsp_snpeff", "hgvsp", "dbnsfp_hgvsp"],
}

def _clean_name(value: str) -> str:
    return (
        str(value)
        .strip()
        .replace("dbNSFP_", "")
        .replace("dbnsfp_", "")
        .replace("#", "")
        .lower()
    )


def _find_column(columns, aliases):
    cleaned = {_clean_name(c): c for c in columns}
    for alias in aliases:
        hit = cleaned.get(_clean_name(alias))
        if hit is not None:
            return hit
    return None


def _is_population_column(column_name: str) -> bool:
    name = _clean_name(column_name)
    return (
        name.endswith("_af")
        or name.endswith("af")
        or "popmax_af" in name
        or "frequency" in name
        or name.endswith("_freq")
    )


def _prediction_to_score(value):
    if pd.isna(value):
        return pd.NA

    raw = str(value).strip().lower()
    if not raw or raw == ".":
        return pd.NA

    scores = []
    for token in raw.replace(",", ";").split(";"):
        token = token.strip()
        if token in {"p", "pathogenic", "d", "damaging", "deleterious", "disease_causing"}:
            scores.append(1.0)
        elif token in {"possibly_damaging", "possibly_pathogenic"}:
            scores.append(0.75)
        elif token in {"b", "benign", "t", "tolerated", "n", "neutral"}:
            scores.append(0.0)
        elif token in {"u", "unknown", "uncertain", "ambiguous"}:
            scores.append(0.5)

    return max(scores) if scores else pd.NA


def _split_multi_value(value, numeric=False):
    if pd.isna(value):
        return pd.NA

    parts = [part for part in str(value).replace(",", ";").split(";") if part and part != "."]
    if not parts:
        return pd.NA

    if not numeric:
        return parts[0]

    vals = pd.to_numeric(pd.Series(parts), errors="coerce").dropna()
    if vals.empty:
        return pd.NA
    return float(vals.max())


def load_dbnsfp_annotations(uploaded_file) -> pd.DataFrame:
    data = uploaded_file.getvalue()
    sep = "," if uploaded_file.name.lower().endswith(".csv") else "\t"
    raw = pd.read_csv(BytesIO(data), sep=sep, low_memory=False)

    rename = {}
    for canonical, aliases in KEY_ALIASES.items():
        column = _find_column(raw.columns, aliases)
        if column is None:
            raise ValueError(
                "Fisierul dbNSFP trebuie sa contina coloanele pentru chrom/pos/ref/alt "
            )
        rename[column] = canonical

    for canonical, aliases in ANNOTATION_ALIASES.items():
        column = _find_column(raw.columns, aliases)
        if column is not None:
            rename[column] = canonical

    for column in raw.columns:
        if column not in rename and _is_population_column(column):
            rename[column] = f"dbnsfp_popfreq_{_clean_name(column)}"

    out = raw[list(rename.keys())].rename(columns=rename).copy()
    out["chrom"] = out["chrom"].astype(str).str.replace("chr", "", case=False, regex=False)
    out["pos"] = pd.to_numeric(out["pos"], errors="coerce").astype("Int64")
    out["ref"] = out["ref"].astype(str).str.upper()
    out["alt"] = out["alt"].astype(str).str.upper()
    out = out.dropna(subset=["chrom", "pos", "ref", "alt"])

    numeric_cols = [
        "dbnsfp_cadd_phred",
        "dbnsfp_revel_score",
        "dbnsfp_alphamissense_score",
        "dbnsfp_sift_score",
        "dbnsfp_polyphen2_hdiv_score",
        "dbnsfp_gnomad_af",
        "dbnsfp_gerp_rs",
        "dbnsfp_aapos",
        "dbnsfp_metarnn_score",
    ]
    numeric_cols.extend([col for col in out.columns if col.startswith("dbnsfp_popfreq_")])
    for col in numeric_cols:
        if col in out.columns:
            out[col] = out[col].apply(lambda v: _split_multi_value(v, numeric=True))

    text_cols = [
        "dbnsfp_alphamissense_pred",
        "dbnsfp_sift_pred",
        "dbnsfp_polyphen2_hdiv_pred",
        "dbnsfp_hgvsp",
        "dbnsfp_metarnn_pred",
    ]
    for col in text_cols:
        if col in out.columns:
            out[col] = out[col].apply(_split_multi_value)

    if "dbnsfp_alphamissense_pred" in out.columns:
        pred_score = out["dbnsfp_alphamissense_pred"].apply(_prediction_to_score)
        if "dbnsfp_alphamissense_score" in out.columns:
            out["dbnsfp_alphamissense_score"] = out["dbnsfp_alphamissense_score"].combine_first(pred_score)
        else:
            out["dbnsfp_alphamissense_score"] = pred_score

    if "dbnsfp_metarnn_pred" in out.columns:
        pred_score = out["dbnsfp_metarnn_pred"].apply(_prediction_to_score)
        if "dbnsfp_metarnn_score" in out.columns:
            out["dbnsfp_metarnn_score"] = out["dbnsfp_metarnn_score"].combine_first(pred_score)
        else:
            out["dbnsfp_metarnn_score"] = pred_score

    if "dbnsfp_metarnn_score" in out.columns and "dbnsfp_metarnn_score" not in numeric_cols:
        numeric_cols.append("dbnsfp_metarnn_score")
    if "dbnsfp_alphamissense_score" in out.columns and "dbnsfp_alphamissense_score" not in numeric_cols:
        numeric_cols.append("dbnsfp_alphamissense_score")

    aggregations = {
        col: "max" for col in numeric_cols if col in out.columns
    }
    aggregations.update({col: "first" for col in text_cols if col in out.columns})

    return (
        out.groupby(["chrom", "pos", "ref", "alt"], as_index=False)
        .agg(aggregations)
    )


def merge_dbnsfp_annotations(variants_df: pd.DataFrame, dbnsfp_df: pd.DataFrame | None) -> pd.DataFrame:
    if dbnsfp_df is None or dbnsfp_df.empty:
        return variants_df.copy()

    left = variants_df.copy()
    left["chrom"] = left["chrom"].astype(str).str.replace("chr", "", case=False, regex=False)
    left["pos"] = pd.to_numeric(left["pos"], errors="coerce").astype("Int64")
    left["ref"] = left["ref"].astype(str).str.upper()
    left["alt"] = left["alt"].astype(str).str.upper()

    merged = left.merge(dbnsfp_df, on=["chrom", "pos", "ref", "alt"], how="left", suffixes=("", "_dbnsfp"))

    dup_cols = [c for c in merged.columns if c.endswith("_dbnsfp")]
    for dup_col in dup_cols:
        base_col = dup_col[:-8]
        if base_col in merged.columns:
            merged[base_col] = merged[base_col].combine_first(merged[dup_col])
            merged = merged.drop(columns=[dup_col])
        else:
            merged = merged.rename(columns={dup_col: base_col})

    return merged

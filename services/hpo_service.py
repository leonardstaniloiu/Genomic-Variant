import pandas as pd
import numpy as np
import streamlit as st
from pathlib import Path

@st.cache_resource(show_spinner="Se incarca baza de date HPO!")
def load_hpo_index(hpo_file: Path):
    if not hpo_file.exists():
        raise FileNotFoundError(f"fisier inexistent: {hpo_file}")

    df = pd.read_csv(hpo_file, sep="\t", low_memory=False)
    df = df[df["gene_symbol"] != "-"]

    gene_to_hpo = df.groupby("gene_symbol")["hpo_id"].apply(set).to_dict()

    hpo_terms = (
        df.groupby(["hpo_id", "hpo_name"])["gene_symbol"]
        .nunique()
        .reset_index()
        .rename(columns={"gene_symbol": "gene_count"})
        .sort_values("gene_count", ascending=False)
    )

    exclude = {
        "HP:0000007", "HP:0000006", "HP:0001466", "HP:0001423",
        "HP:0001417", "HP:0001419", "HP:0001425", "HP:0001426", "HP:0001427"
    }
    hpo_terms = hpo_terms[~hpo_terms["hpo_id"].isin(exclude)]

    hpo_weight_map = (
        hpo_terms.assign(
            hpo_weight=lambda d: 1.0 / (1.0 + np.log1p(d["gene_count"]))
        )
        .set_index("hpo_id")["hpo_weight"]
        .to_dict()
    )

    return gene_to_hpo, hpo_terms, hpo_weight_map


def hpo_match(df: pd.DataFrame, patient_hpo: list, gene_to_hpo: dict, hpo_weight_map: dict | None = None) -> pd.DataFrame:
    patient_set = set(patient_hpo)
    hpo_weight_map = hpo_weight_map or {}

    hpo_match_col = []
    hpo_count_col = []
    matched_col = []
    phenotype_score_col = []
    overlap_ratio_col = []

    for _, row in df.iterrows():
        gene = row.get("gene_symbol")
        gene_terms = gene_to_hpo.get(gene, set()) if pd.notna(gene) else set()
        matched = gene_terms.intersection(patient_set)
        matched_weight = sum(hpo_weight_map.get(term, 0.0) for term in matched)
        overlap_ratio = (len(matched) / len(gene_terms)) if gene_terms else 0.0

        hpo_match_col.append(1 if matched else 0)
        hpo_count_col.append(len(matched))
        matched_col.append("|".join(sorted(matched)) if matched else "")
        phenotype_score_col.append(float(min(1.0, matched_weight / 3.0 + overlap_ratio * 0.4)))
        overlap_ratio_col.append(overlap_ratio)

    out = df.copy()
    out["hpo_match"] = hpo_match_col
    out["hpo_match_count"] = hpo_count_col
    out["matched_hpo_terms"] = matched_col
    out["phenotype_score"] = phenotype_score_col
    out["hpo_overlap_ratio"] = overlap_ratio_col

    return out

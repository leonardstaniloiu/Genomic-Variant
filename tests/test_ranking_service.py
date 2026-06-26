import pandas as pd

from services.ranking_service import rank_variants


def test_ranking_varints():
    variants = pd.DataFrame(
        [
            {
                "variant_id": "weak",
                "gene_symbol": "GENE2",
                "consequence": "synonymous_variant",
                "impact": "LOW",
                "hpo_match_count": 0,
                "phenotype_score": 0.0,
                "dbnsfp_revel_score": 0.1,
                "dbnsfp_gnomad_af": 0.04,
            },
            {
                "variant_id": "strong",
                "gene_symbol": "GENE1",
                "consequence": "stop_gained",
                "impact": "HIGH",
                "hpo_match_count": 4,
                "phenotype_score": 0.9,
                "dbnsfp_revel_score": 0.95,
                "dbnsfp_gnomad_af": 0.0,
                "hgvsp": "p.Arg100Ter",
            },
        ]
    )

    ranked = rank_variants(variants)

    assert ranked.loc[0, "variant_id"] == "strong"
    assert ranked["pathogenicity_score"].between(0, 1).all()
    assert ranked.loc[0, "pathogenicity_score"] > ranked.loc[1, "pathogenicity_score"]

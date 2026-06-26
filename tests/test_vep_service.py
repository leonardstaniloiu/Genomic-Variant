import json

import pandas as pd

from services.vep_service import annotate_vep, format_region


class FakeResponse:
    def __init__(self, status_code, payload=None, text=""):
        self.status_code = status_code
        self._payload = payload or []
        self.text = text

    def json(self):
        return self._payload


def test_format_region():
    row = {"chrom": "chr1", "pos": 12345, "alt": "T"}

    assert format_region(row) == "1:12345-12345:1/T"


def test_annotate_vep(monkeypatch):
    variants = pd.DataFrame(
        [{"chrom": "chr1", "pos": 12345, "ref": "C", "alt": "T", "qual": 90.0}]
    )

    def fake_get(url, headers, timeout):
        assert url == "https://rest.ensembl.org/info/ping"
        return FakeResponse(200, {"ping": 1})

    def fake_post(url, headers, data, timeout):
        payload = json.loads(data)
        assert payload == {"variants": ["1:12345-12345:1/T"]}
        return FakeResponse(
            200,
            [
                {
                    "input": "1:12345-12345:1/T",
                    "most_severe_consequence": "missense_variant",
                    "transcript_consequences": [
                        {
                            "gene_symbol": "GENE1",
                            "gene_id": "ENSG000001",
                            "impact": "MODERATE",
                            "canonical": 1,
                            "hgvsp": "p.Arg123Trp",
                            "cadd_phred": 24.1,
                            "alphamissense": {
                                "am_pathogenicity": 0.82,
                                "am_class": "likely_pathogenic",
                            },
                        }
                    ],
                }
            ],
        )

    monkeypatch.setattr("services.vep_service.requests.get", fake_get)
    monkeypatch.setattr("services.vep_service.requests.post", fake_post)
    monkeypatch.setattr("services.vep_service.time.sleep", lambda seconds: None)

    annotated = annotate_vep(variants)

    assert annotated.loc[0, "variant_id"] == "1_12345_C/T"
    assert annotated.loc[0, "gene_symbol"] == "GENE1"
    assert annotated.loc[0, "consequence"] == "missense_variant"
    assert annotated.loc[0, "impact"] == "MODERATE"
    assert annotated.loc[0, "am_score"] == 0.82
    assert annotated.loc[0, "cadd_from_vep"] == 24.1

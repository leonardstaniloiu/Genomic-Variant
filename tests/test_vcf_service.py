import time

import pandas as pd
import pytest

from services.vcf_service import parse_vcf


def test_parse_vcf():
    content = "\n".join(
        [
            "##fileformat=VCFv4.2",
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
            "1\t100\t.\tA\tG\t50\tPASS\t.",
            "chr2\t200\t.\tC\tT\t.\t.\t.",
            "3\t300\t.\tG\tA\t10\tPASS\t.",
            "4\t400\t.\tT\tC\t60\tLowQual\t.",
            "GL0001\t500\t.\tA\tC\t60\tPASS\t.",
            "5\t600\t.\tA\tC,T\t60\tPASS\t.",
        ]
    )

    result = parse_vcf(content, min_qual=20.0)

    assert list(result["chrom"]) == ["1", "chr2"]
    assert list(result["pos"]) == [100, 200]
    assert list(result["alt"]) == ["G", "T"]
    assert result.loc[0, "qual"] == 50.0
    assert pd.isna(result.loc[1, "qual"])


@pytest.mark.performance
def test_parse_vcf_5000_variants():
    lines = ["#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"]
    lines.extend(f"1\t{i}\t.\tA\tG\t60\tPASS\t." for i in range(1, 5001))

    start = time.perf_counter()
    result = parse_vcf("\n".join(lines))
    elapsed = time.perf_counter() - start

    assert len(result) == 5000
    assert elapsed < 2.0

import pandas as pd

def parse_vcf(content: str, min_qual: float = 20.0) -> pd.DataFrame:
    rows = []

    valid_chroms = (
        {str(i) for i in range(1, 23)}
        | {"X", "Y", "MT", "M"}
        | {f"chr{i}" for i in range(1, 23)}
        | {"chrX", "chrY", "chrM", "chrMT"}
    )

    for line in content.splitlines():
        if not line or line.startswith("#"):
            continue

        fields = line.strip().split("\t")
        if len(fields) < 8:
            continue

        chrom, pos, id_, ref, alt, qual, filter_, info = fields[:8]
        chrom = str(chrom)

        if chrom not in valid_chroms:
            continue
        if "," in alt or alt == ".":
            continue

        if qual == ".":
            qual_val = None
        else:
            try:
                qual_val = float(qual)
            except ValueError:
                qual_val = None

        if filter_ not in ("PASS", "."):
            continue

        if qual_val is not None and qual_val < min_qual:
            continue

        try:
            pos_val = int(pos)
        except ValueError:
            continue

        rows.append({
            "chrom": chrom,
            "pos": pos_val,
            "ref": str(ref),
            "alt": str(alt),
            "qual": qual_val,
            "filter": filter_,
            "info": info,
        })

    return pd.DataFrame(rows) if rows else pd.DataFrame()
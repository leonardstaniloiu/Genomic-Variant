from pathlib import Path

INPUT_VCF = Path("data/clinvar.vcf")
OUTPUT_VCF = Path("data/retinitis_demo.vcf")

TARGET_GENES = {
    "RPE65",
    "CRX",
    "AIPL1",
    "LRAT",
    "IMPDH1"
}

TARGET_SIG = {
    "Pathogenic",
    "Likely_pathogenic",
    "Pathogenic/Likely_pathogenic",
}

VALID_CHROMS = {str(i) for i in range(1, 23)} | {"X", "Y", "MT", "M", "chrX", "chrY", "chrM", "chrMT"} | {f"chr{i}" for i in range(1, 23)}

def parse_info(info_str: str) -> dict:
    info = {}
    for item in info_str.split(";"):
        if "=" in item:
            k, v = item.split("=", 1)
            info[k] = v
        else:
            info[item] = True
    return info


def gene_matches(geneinfo: str) -> bool:
    if not geneinfo:
        return False
    for entry in geneinfo.split("|"):
        gene = entry.split(":")[0].strip()
        if gene in TARGET_GENES:
            return True
    return False


def sig_matches(clnsig: str) -> bool:
    if not clnsig:
        return False
    values = [x.strip() for x in clnsig.split("|")]
    return any(v in TARGET_SIG for v in values)


def is_simple_variant(ref: str, alt: str) -> bool:
    if alt == "." or "," in alt:
        return False
    return True


def main():
    kept = 0

    with INPUT_VCF.open("r", encoding="utf-8") as fin, OUTPUT_VCF.open("w", encoding="utf-8") as fout:
        for line in fin:
            if line.startswith("##"):
                fout.write(line)
                continue

            if line.startswith("#CHROM"):
                fout.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                continue

            chrom, pos, var_id, ref, alt, qual, filt, info_str = fields[:8]

            if chrom not in VALID_CHROMS:
                continue

            if not is_simple_variant(ref, alt):
                continue

            info = parse_info(info_str)
            geneinfo = info.get("GENEINFO", "")
            clnsig = info.get("CLNSIG", "")

            if not gene_matches(geneinfo):
                continue

            if not sig_matches(clnsig):
                continue

            qual_out = qual if qual != "." else "100"
            filt_out = filt if filt != "." else "PASS"

            fout.write(
                "\t".join([chrom, pos, var_id, ref, alt, qual_out, filt_out, info_str]) + "\n"
            )
            kept += 1

    print(f"Done. Extracted {kept} variants to {OUTPUT_VCF}")


if __name__ == "__main__":
    main()
import json
import time
import requests
import pandas as pd

def format_region(row) -> str:
    chrom = str(row["chrom"]).replace("chr", "")
    pos = int(row["pos"])
    alt = str(row["alt"])
    return f"{chrom}:{pos}-{pos}:1/{alt}"


def annotate_vep(variants_df: pd.DataFrame, progress_cb=None) -> pd.DataFrame:
    VEP_URL = "https://rest.ensembl.org/vep/human/region"
    PARAMS = "canonical=1&hgvs=1&AlphaMissense=1&CADD=1"
    HEADERS = {"Content-Type": "application/json", "Accept": "application/json"}

    ping = requests.get(
        "https://rest.ensembl.org/info/ping",
        headers={"Accept": "application/json"},
        timeout=10
    )
    if ping.status_code != 200:
        raise RuntimeError("Ensembl VEP API nu este disponibil.")

    regions = [format_region(r) for _, r in variants_df.iterrows()]
    id_map = {format_region(r): r.to_dict() for _, r in variants_df.iterrows()}
    annotated = {}

    batch_size = 100

    for i in range(0, len(regions), batch_size):
        batch = regions[i:i + batch_size]

        if progress_cb:
            progress_cb(i / max(len(regions), 1), f"VEP batch {i // batch_size + 1}")

        resp = requests.post(
            f"{VEP_URL}?{PARAMS}",
            headers=HEADERS,
            data=json.dumps({"variants": batch}),
            timeout=120
        )
        if resp.status_code != 200:
            raise RuntimeError(f"Raspuns vep - {resp.status_code}: {resp.text[:300]}")

        for res in resp.json():
            tcs = res.get("transcript_consequences", [])
            best = None
            best_score = -1

            for t in tcs:
                impact = t.get("impact", "MODIFIER")
                canonical = t.get("canonical", 0)
                score = {"HIGH": 4, "MODERATE": 3, "LOW": 2, "MODIFIER": 1}.get(impact, 1) * 10 + canonical
                if score > best_score:
                    best_score = score
                    best = t

            am = (best or {}).get("alphamissense", {}) or {}

            annotated[res.get("input", "")] = {
                "gene_symbol": (best or {}).get("gene_symbol"),
                "gene_id": (best or {}).get("gene_id"),
                "consequence": res.get("most_severe_consequence"),
                "impact": (best or {}).get("impact"),
                "hgvsp": (best or {}).get("hgvsp"),
                "am_score": am.get("am_pathogenicity"),
                "am_class": am.get("am_class"),
                "cadd_from_vep": (best or {}).get("cadd_phred"),
            }

        time.sleep(0.25)

    rows = []
    for region in regions:
        orig = id_map[region]
        chrom_c = str(orig["chrom"]).replace("chr", "")
        variant_id = f"{chrom_c}_{orig['pos']}_{orig['ref']}/{orig['alt']}"
        ann = annotated.get(region)

        rows.append({
            **orig,
            "variant_id": variant_id,
            "gene_symbol": ann.get("gene_symbol") if ann else None,
            "gene_id": ann.get("gene_id") if ann else None,
            "consequence": ann.get("consequence") if ann else None,
            "impact": ann.get("impact") if ann else None,
            "hgvsp": ann.get("hgvsp") if ann else None,
            "am_score": ann.get("am_score") if ann else None,
            "am_class": ann.get("am_class") if ann else None,
            "cadd_from_vep": ann.get("cadd_from_vep") if ann else None,
        })

    return pd.DataFrame(rows)
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent.parent
HPO_FILE = BASE_DIR / "data" / "genes_to_phenotype.txt"

CONSEQUENCE_IMPACT = {
    "transcript_ablation": 8,
    "stop_gained": 7,
    "frameshift_variant": 7,
    "splice_acceptor_variant": 6,
    "splice_donor_variant": 6,
    "start_lost": 6,
    "stop_lost": 5,
    "missense_variant": 4,
    "inframe_insertion": 3,
    "inframe_deletion": 3,
    "splice_region_variant": 2,
    "synonymous_variant": 1,
    "3_prime_UTR_variant": 1,
    "5_prime_UTR_variant": 1,
    "intron_variant": 0,
    "intergenic_variant": 0,
}

IMPACT_NUMERIC = {
    "HIGH": 1.0,
    "MODERATE": 0.7,
    "LOW": 0.35,
    "MODIFIER": 0.1,
}
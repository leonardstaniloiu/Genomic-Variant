# Teste

Rulare completa:

```powershell
pytest
```

Rulare fara testele simple de performanta:

```powershell
pytest -m "not performance"
```

Ce acopera testele:

- `test_vcf_service.py`: parsare VCF si un control simplu de performanta pentru 5.000 variante.
- `test_ranking_service.py`: ordonarea variantelor dupa semnal fenotipic, molecular si consecinta.
- `test_vep_service.py`: formatul regiunii Ensembl si adnotarea VEP cu API mock-uit, fara apeluri reale la internet.

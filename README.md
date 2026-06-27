# Genomic Variant Prioritizer [Live demo](https://genomic-ranking.streamlit.app/)

A web application for prioritizing genetic variants based on patient phenotype and genomic annotations.

## Overview

This application helps prioritize candidate genetic variants by integrating:

- Human Phenotype Ontology (HPO)
- Ensembl Variant Effect Predictor (VEP)
- dbNSFP pathogenicity scores
- Variant Call Format (VCF) files

The platform was developed as part of a Master's dissertation in Bioinformatics.

## Features

- Upload VCF files
- Import dbNSFP annotation files
- Select Human Phenotype Ontology (HPO) terms
- Automatic variant prioritization
- Interactive visualization of results
- Export prioritized variants

## Technologies

- Python
- Streamlit
- Pandas
- Requests
- Ensembl VEP REST API

## Project Structure

```
app.py
components/
services/
models/
utils/
data/
tests/
```

## Installation

Clone the repository:

```bash
git clone https://github.com/username/genomic-variant.git
cd genomic-variant
```

Install dependencies:

```bash
pip install -r requirements.txt
```

Run the application:

```bash
streamlit run app.py
```

OR

```bash
python -m streamlit run app.py
```

## Input

The application accepts:

- VCF files
- dbNSFP annotation files
- Patient information
- HPO terms

## Output

The application generates:

- Ranked genetic variants
- Pathogenicity scores
- Functional annotations
- Interactive tables and charts

## Future Improvements

- Machine Learning based prioritization
- Explainable AI (SHAP/LIME)
- GPU accelerated annotation (G-VEP)

## Author

Leonard Stăniloiu

Leonard Stăniloiu

Master's Dissertation – 2026

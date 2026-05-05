# CpG Methylation ML Pipeline

> Machine learning analysis of DNA methylation patterns across a shared genomic region in multiple independent cohorts.

---

## Overview

This project applies unsupervised and supervised machine learning to CpG-level methylation data from a shared ~7–8kb genomic region (CpG island + gene body) across four independent datasets 

---

## Repository Structure

```
CpG_Methylation_ML/
├── data/
│   ├── raw/                  # Raw beta value matrices per dataset
│   └── processed/            # QC-filtered, aligned matrices
├── notebooks/                # Exploratory Jupyter notebooks
├── src/
│   ├── preprocessing/        # Data loading, QC, alignment
│   ├── unsupervised/         # PCA, UMAP, clustering
│   ├── supervised/           # ElasticNet, Random Forest
│   ├── integration/          # MOFA+, cross-dataset analysis
│   └── visualization/        # Plotting utilities
├── results/
│   ├── figures/              # Output plots
│   ├── tables/               # Output tables
│   └── models/               # Saved model objects
├── docs/                     # Methods documentation
├── tests/                    # Unit tests
├── requirements.txt
├── environment.yml
└── README.md
```

---

## Methods


### Machine Learning Options

#### 1. Preprocessing
- Coverage filtering (minimum reads per CpG per sample)
- Variance filtering (remove near-zero variance CpGs)
- Alignment of CpG positions across datasets

#### 2. Unsupervised Analysis
- PCA — variance structure, sample separation
- UMAP — non-linear dimensionality reduction
- Hierarchical clustering heatmap (samples × CpGs)

#### 3. Supervised Analysis
- ElasticNet logistic regression (primary — regularized for small N)
- Random Forest (largest N)
- Cross-validation: LOOCV for small datasets, 5-fold for stress cohort

#### 4. Multi-Dataset Integration (MOFA+)
- Latent factor decomposition across all datasets simultaneously
- Identifies shared vs. exposure-specific methylation factors
- Implemented via `mofapy2`

---

## Installation

```bash
git clone https://github.com/Ensiuc/CpG_Methylation_ML.git
cd CpG_Methylation_ML
conda env create -f environment.yml
conda activate cpg_ml
```

Or with pip:
```bash
pip install -r requirements.txt
```

---

## Usage

```bash
# 1. Preprocess all datasets
python src/preprocessing/preprocess.py --config configs/config.yaml

# 2. Run unsupervised analysis
python src/unsupervised/run_unsupervised.py

# 3. Run supervised analysis
python src/supervised/run_supervised.py

# 4. Run MOFA+ integration
python src/integration/run_mofa.py
```

---

## Key Dependencies

| Package | Version | Purpose |
|---|---|---|
| numpy | ≥1.24 | Numerical computing |
| pandas | ≥2.0 | Data manipulation |
| scikit-learn | ≥1.3 | ML models, CV |
| umap-learn | ≥0.5 | UMAP |
| mofapy2 | ≥0.7 | Multi-dataset integration |
| matplotlib | ≥3.7 | Plotting |
| seaborn | ≥0.12 | Statistical visualization |
| scipy | ≥1.11 | Statistics |

---

## Citation

If you use this pipeline, please cite the relevant tools:

- MOFA+: Argelaguet et al. (2020) *Genome Biology*
- UMAP: McInnes et al. (2018) *JOSS*

---

## Author

**Ensieh Habibi** | GitHub: [@Ensiuc](https://github.com/Ensiuc)

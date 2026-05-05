# CpG Methylation ML Pipeline

> Machine learning analysis of DNA methylation patterns across a shared genomic region in multiple independent cohorts exposed to distinct early-life stressors.

---

## Overview

This project applies unsupervised and supervised machine learning to CpG-level methylation data from a shared ~7–8kb genomic region (CpG island + gene body) across four independent datasets representing different early-life exposures:

| Dataset | Exposure | Tissue | N (Control / Exposed) |
|---|---|---|---|
| Wildfire Smoke | PM2.5 exposure (2008 CA wildfires) | Nasal epithelial | 14 / 8 |
| Obesity | Maternal obesity | Brain regions | 6 / 7 |
| Preconception Stress | Maternal stress | TBD | ~45 / ~25 |
| cfDNA (Longitudinal) | — | Cell-free fetal DNA | 4–7 / trimester |

The core biological question: **do different early-life stressors converge on shared methylation patterns within this region?**

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

### Prior Differential Methylation Analysis
Before ML, differential methylation was assessed using:
- **dmrseq** — region-based DMR detection with GLS + permutation testing
- **Feature-level t-tests** — sample-level mean methylation per feature (CpG island, gene body) with global FDR correction
- Unit of analysis: biological sample (not individual CpG positions)

### Machine Learning Pipeline

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
- Random Forest (preconception stress dataset only, largest N)
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
- dmrseq: Korthauer et al. (2019) *Biostatistics*
- MOFA+: Argelaguet et al. (2020) *Genome Biology*
- UMAP: McInnes et al. (2018) *JOSS*

---

## Author

**Ensiuc** | GitHub: [@Ensiuc](https://github.com/Ensiuc)

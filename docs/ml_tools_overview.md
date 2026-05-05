# ML Tools Overview for CpG Methylation Analysis

> This document covers all machine learning tools considered for analyzing CpG-level methylation data across a shared genomic region (~7–8kb) in four independent cohorts. For each tool, we describe what it does, its pros and cons in the context of this project, and why we selected or excluded it.

---

## Table of Contents

1. [Dimensionality Reduction & Unsupervised](#1-dimensionality-reduction--unsupervised)
   - [PCA](#11-pca-principal-component-analysis)
   - [UMAP](#12-umap-uniform-manifold-approximation-and-projection)
   - [t-SNE](#13-t-sne-t-distributed-stochastic-neighbor-embedding)
   - [Hierarchical Clustering](#14-hierarchical-clustering)
   - [k-Means Clustering](#15-k-means-clustering)
2. [Supervised Classification](#2-supervised-classification)
   - [ElasticNet Logistic Regression](#21-elasticnet-logistic-regression)
   - [LASSO Logistic Regression](#22-lasso-logistic-regression)
   - [Random Forest](#23-random-forest)
   - [Support Vector Machine (SVM)](#24-support-vector-machine-svm)
   - [Gradient Boosting (XGBoost / LightGBM)](#25-gradient-boosting-xgboost--lightgbm)
   - [Neural Networks / Deep Learning](#26-neural-networks--deep-learning)
3. [Multi-Dataset Integration](#3-multi-dataset-integration)
   - [MOFA+](#31-mofa-multi-omics-factor-analysis)
   - [Harmony](#32-harmony)
   - [Combat / Combat-seq](#33-combat--combat-seq)
   - [Multi-Task Learning](#34-multi-task-learning)
   - [Meta-Analysis (Statistical)](#35-meta-analysis-statistical)
4. [Selected Tools & Rationale](#4-selected-tools--rationale)
5. [Excluded Tools & Why](#5-excluded-tools--why)

---

## 1. Dimensionality Reduction & Unsupervised

These methods explore structure in the data without using condition labels. They are used to ask: *do samples naturally separate by condition, tissue, or dataset?*

---

### 1.1 PCA (Principal Component Analysis)

**What it does:**
PCA is a linear method that finds the directions (principal components) of maximum variance in the data. Each PC is a linear combination of CpG features. The first PC captures the most variance, the second captures the most remaining variance orthogonal to the first, and so on.

**In this context:**
- Rows = biological samples
- Columns = CpG beta values within the region
- Output = low-dimensional embedding (PC1, PC2, PC3...) showing how samples relate

**Pros:**
- Fully interpretable — loadings tell you exactly which CpGs drive each component
- Fast and deterministic (same result every run)
- Works well even with small sample sizes
- Variance explained per PC is quantifiable (scree plot)
- No hyperparameters to tune
- Gold standard first step in any omics analysis

**Cons:**
- Linear only — cannot capture non-linear methylation patterns
- Sensitive to outliers
- PCs may reflect batch effects or technical variance rather than biology
- With many CpGs relative to samples (high-dimensional), interpretation can be noisy

**Learn more:**
- [scikit-learn PCA documentation](https://scikit-learn.org/stable/modules/generated/sklearn.decomposition.PCA.html)
- [StatQuest PCA explanation (YouTube)](https://www.youtube.com/watch?v=FgakZw6K1QQ)
- [Bioconductor PCA for methylation](https://bioconductor.org/packages/release/bioc/vignettes/minfi/inst/doc/minfi.html)

**Status: ✅ SELECTED** — See [Section 4](#4-selected-tools--rationale)

---

### 1.2 UMAP (Uniform Manifold Approximation and Projection)

**What it does:**
UMAP is a non-linear dimensionality reduction method that preserves both local and global structure in high-dimensional data. It constructs a high-dimensional graph of the data and then optimizes a low-dimensional layout to match it.

**In this context:**
- Captures complex, non-linear methylation patterns across CpGs
- Reveals clusters that PCA might miss
- Often used after PCA as a follow-up

**Pros:**
- Captures non-linear structure — better than PCA for complex biological patterns
- Preserves global structure better than t-SNE
- Faster than t-SNE
- Handles high-dimensional data well
- Visually intuitive — clusters are usually clear

**Cons:**
- Results vary with random seed and hyperparameters (n_neighbors, min_dist)
- Less interpretable than PCA — no direct loadings
- Can be misleading with very small datasets (n < 20)
- Distances between clusters in UMAP space are not directly meaningful
- Not deterministic by default

**Learn more:**
- [UMAP documentation](https://umap-learn.readthedocs.io/en/latest/)
- [Original UMAP paper — McInnes et al. 2018](https://arxiv.org/abs/1802.03426)
- [Understanding UMAP (interactive)](https://pair-code.github.io/understanding-umap/)

**Status: ✅ SELECTED** — See [Section 4](#4-selected-tools--rationale)

---

### 1.3 t-SNE (t-Distributed Stochastic Neighbor Embedding)

**What it does:**
t-SNE is a non-linear dimensionality reduction method that focuses on preserving local neighborhood structure. It maps high-dimensional points to 2D by minimizing the divergence between probability distributions over pairwise distances.

**Pros:**
- Excellent at revealing local cluster structure
- Widely used in single-cell and methylation literature
- Visually striking separation of clusters

**Cons:**
- Does NOT preserve global structure — distances between clusters are meaningless
- Very slow on large datasets
- Highly sensitive to perplexity hyperparameter
- Results are not reproducible without fixed seed
- Cannot embed new data points (no transform function)
- Generally superseded by UMAP for most use cases

**Learn more:**
- [scikit-learn t-SNE](https://scikit-learn.org/stable/modules/generated/sklearn.manifold.TSNE.html)
- [How to use t-SNE effectively](https://distill.pub/2016/misread-tsne/)

**Status: ❌ NOT SELECTED** — See [Section 5](#5-excluded-tools--why)

---

### 1.4 Hierarchical Clustering (Heatmap)

**What it does:**
Hierarchical clustering builds a tree (dendrogram) of samples and/or features based on a distance metric (e.g., Euclidean, correlation). When visualized as a heatmap with dendrograms on both axes, it reveals methylation patterns across the region and how samples group.

**In this context:**
- Rows = samples, columns = CpG positions in the region
- Color = beta methylation value (0–1)
- Dendrograms show which samples and which CpGs cluster together

**Pros:**
- Highly interpretable — you see the actual methylation values
- No need to specify number of clusters in advance
- Reveals spatial patterns within the region (e.g., a sub-region that is consistently hypermethylated in exposed samples)
- Works with any sample size
- Standard in epigenomics publications

**Cons:**
- Sensitive to distance metric and linkage method choice
- Does not scale well to very large numbers of CpGs (but fine for one region)
- Clusters are not statistically validated by default
- Color scale choice affects visual interpretation

**Learn more:**
- [seaborn clustermap](https://seaborn.pydata.org/generated/seaborn.clustermap.html)
- [scipy hierarchical clustering](https://docs.scipy.org/doc/scipy/reference/cluster.hierarchy.html)

**Status: ✅ SELECTED** — See [Section 4](#4-selected-tools--rationale)

---

### 1.5 k-Means Clustering

**What it does:**
k-Means partitions samples into k clusters by minimizing within-cluster variance. It iteratively assigns samples to the nearest centroid and updates centroids until convergence.

**Pros:**
- Simple and fast
- Works well when clusters are spherical and similar in size

**Cons:**
- Requires specifying k in advance — not obvious in methylation data
- Assumes spherical clusters — methylation patterns are rarely this clean
- Sensitive to outliers and initialization
- Not appropriate for small, unbalanced datasets (your case)
- Deterministic only with fixed seed

**Learn more:**
- [scikit-learn k-Means](https://scikit-learn.org/stable/modules/generated/sklearn.cluster.KMeans.html)

**Status: ❌ NOT SELECTED** — See [Section 5](#5-excluded-tools--why)

---

## 2. Supervised Classification

These methods use condition labels (exposed vs. control) to train a model that distinguishes groups based on CpG methylation patterns. They also provide feature importance — identifying which CpGs within the region best discriminate conditions.

---

### 2.1 ElasticNet Logistic Regression

**What it does:**
ElasticNet combines L1 (LASSO) and L2 (Ridge) regularization penalties in logistic regression. The L1 penalty drives some coefficients exactly to zero (feature selection), while L2 handles correlated features by distributing coefficients across them rather than zeroing all but one.

**Model:**
```
penalty = α * L1 + (1-α) * L2
```
where α (l1_ratio) controls the mix.

**In this context:**
- Input: sample × CpG beta matrix
- Output: binary classification (exposed vs. control) + non-zero CpG coefficients
- Non-zero coefficients = CpGs that best predict condition

**Pros:**
- Handles correlated features well (nearby CpGs are correlated — this is critical for methylation)
- Performs feature selection automatically (L1 component)
- Works with small sample sizes — regularization prevents overfitting
- Coefficients are directly interpretable (direction + magnitude of effect)
- Fast to train and cross-validate
- Well-validated in methylation and genomics literature

**Cons:**
- Assumes linear relationship between methylation and outcome
- Sensitive to feature scaling (beta values 0–1 helps, but coverage variation can still matter)
- Optimal hyperparameters (alpha, l1_ratio) require cross-validation
- With very small N (<10), any supervised method is unreliable

**Learn more:**
- [scikit-learn ElasticNet](https://scikit-learn.org/stable/modules/generated/sklearn.linear_model.ElasticNet.html)
- [Zou & Hastie (2005) — Original ElasticNet paper](https://rss.onlinelibrary.wiley.com/doi/10.1111/j.1467-9868.2005.00503.x)
- [Regularization in genomics (tutorial)](https://scikit-learn.org/stable/modules/linear_model.html#elastic-net)

**Status: ✅ SELECTED** — See [Section 4](#4-selected-tools--rationale)

---

### 2.2 LASSO Logistic Regression

**What it does:**
LASSO (Least Absolute Shrinkage and Selection Operator) applies only an L1 penalty, which forces sparse solutions — many coefficients become exactly zero. It is a special case of ElasticNet with l1_ratio = 1.

**Pros:**
- Strong feature selection — produces very sparse models
- Simple and interpretable
- Well-established in genomics

**Cons:**
- When CpGs are correlated (which they are within a region), LASSO arbitrarily selects one and zeros the rest — this loses biologically meaningful spatial information
- ElasticNet handles this better by distributing weight across correlated features

**Learn more:**
- [scikit-learn LASSO](https://scikit-learn.org/stable/modules/generated/sklearn.linear_model.Lasso.html)
- [Tibshirani (1996) — Original LASSO paper](https://www.jstor.org/stable/2346178)

**Status: ⚠️ CONSIDERED but superseded by ElasticNet** — See [Section 5](#5-excluded-tools--why)

---

### 2.3 Random Forest

**What it does:**
Random Forest builds an ensemble of decision trees, each trained on a random subset of samples (bootstrap) and features (CpGs). Final prediction is by majority vote. Feature importance is derived from how much each CpG reduces impurity across all trees.

**In this context:**
- Can capture non-linear relationships between CpG methylation and condition
- Feature importance gives a ranked list of most informative CpGs
- Does not assume any particular distribution of methylation values

**Pros:**
- Handles non-linear patterns and interactions between CpGs
- Robust to outliers
- No need for feature scaling
- Built-in feature importance (mean decrease in impurity or permutation importance)
- Works with correlated features
- Out-of-bag (OOB) error provides internal validation

**Cons:**
- Needs sufficient sample size — generally N > 20 per class recommended
- Black box — less interpretable than ElasticNet coefficients
- Feature importance can be biased toward high-cardinality or continuous features
- Overfits easily with very small N
- Not appropriate for obesity (n=13) or cfDNA (n=8–14) datasets

**Learn more:**
- [scikit-learn Random Forest](https://scikit-learn.org/stable/modules/generated/sklearn.ensemble.RandomForestClassifier.html)
- [Breiman (2001) — Original Random Forest paper](https://link.springer.com/article/10.1023/A:1010933404324)
- [Understanding Random Forest feature importance](https://scikit-learn.org/stable/auto_examples/inspection/plot_permutation_importance.html)

**Status: ✅ SELECTED (preconception stress dataset only)** — See [Section 4](#4-selected-tools--rationale)

---

### 2.4 Support Vector Machine (SVM)

**What it does:**
SVM finds the hyperplane that maximally separates classes in feature space. With a kernel (e.g., RBF), it can find non-linear decision boundaries. The support vectors are the samples closest to the decision boundary.

**Pros:**
- Works well in high-dimensional spaces (many CpGs relative to samples)
- Effective when classes are clearly separable
- Kernel trick allows non-linear classification

**Cons:**
- Does not provide feature importance natively (requires SHAP or permutation)
- Slow to tune (kernel, C, gamma all need cross-validation)
- Not interpretable compared to ElasticNet
- No probability output by default (requires Platt scaling)
- With small N, performance is unreliable and cross-validation folds are too small

**Learn more:**
- [scikit-learn SVM](https://scikit-learn.org/stable/modules/svm.html)
- [Understanding SVMs](https://www.youtube.com/watch?v=efR1C6CvhmE)

**Status: ❌ NOT SELECTED** — See [Section 5](#5-excluded-tools--why)

---

### 2.5 Gradient Boosting (XGBoost / LightGBM)

**What it does:**
Gradient boosting builds an ensemble of weak learners (shallow trees) sequentially, where each tree corrects the errors of the previous. XGBoost and LightGBM are highly optimized implementations.

**Pros:**
- State-of-the-art performance on tabular data
- Handles missing values natively
- Built-in feature importance
- Fast training with LightGBM

**Cons:**
- Requires large sample sizes to perform well — poorly calibrated with small N
- Many hyperparameters to tune (n_estimators, max_depth, learning_rate, etc.)
- High risk of overfitting with small datasets
- Overkill for a focused single-region analysis

**Learn more:**
- [XGBoost documentation](https://xgboost.readthedocs.io/en/stable/)
- [LightGBM documentation](https://lightgbm.readthedocs.io/en/stable/)

**Status: ❌ NOT SELECTED** — See [Section 5](#5-excluded-tools--why)

---

### 2.6 Neural Networks / Deep Learning

**What it does:**
Neural networks learn hierarchical representations through multiple layers of non-linear transformations. Deep learning models (CNNs, Transformers) have been applied to methylation arrays (450K, EPIC) for disease classification.

**Pros:**
- Can learn complex non-linear patterns
- Can model spatial relationships between CpGs (with CNNs)
- State of the art for large methylation datasets

**Cons:**
- Requires large datasets — typically N > 100 per class minimum
- Completely inappropriate for small N studies (your largest dataset has ~70 samples total)
- Black box — very difficult to interpret which CpGs matter
- Computationally expensive
- High risk of overfitting

**Learn more:**
- [PyTorch for genomics](https://pytorch.org/)
- [DeepMethyl paper](https://academic.oup.com/bioinformatics)

**Status: ❌ NOT SELECTED** — See [Section 5](#5-excluded-tools--why)

---

## 3. Multi-Dataset Integration

These methods combine data from multiple datasets to find shared and dataset-specific methylation patterns.

---

### 3.1 MOFA+ (Multi-Omics Factor Analysis)

**What it does:**
MOFA+ is a probabilistic latent factor model designed for multi-group, multi-omics data. It decomposes the data into a set of latent factors, each of which captures a source of variation. Crucially, it models which factors are shared across groups (datasets) and which are specific to one group.

**In this context:**
- One data type: CpG methylation
- Multiple groups: obesity, wildfire, stress, cfDNA trimesters
- Output: latent factors with loadings on specific CpGs
- Factor 1 might be "shared stress response across all exposures"
- Factor 2 might be "obesity-specific"
- Factor 3 might be "trimester trajectory"

**Pros:**
- Explicitly designed for multiple groups with the same features — a perfect fit
- Handles unequal sample sizes across groups naturally
- No batch correction needed — dataset is modeled as a grouping variable
- Separates shared from exposure-specific signals
- Latent factors are biologically interpretable
- CpG loadings on each factor tell you which positions drive the shared signal
- Handles missing CpGs per dataset (not all CpGs need to be present in all groups)
- Has a Python API (mofapy2)

**Cons:**
- Requires careful selection of number of factors (K)
- Assumes Gaussian or other parametric distributions for methylation values
- Convergence depends on initialization — should run multiple times
- Primarily developed in R (MOFA2); Python API (mofapy2) is functional but less documented
- Harder to interpret than simple PCA

**Learn more:**
- [MOFA+ paper — Argelaguet et al. 2020, Genome Biology](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02015-1)
- [mofapy2 Python package](https://github.com/bioFAM/mofapy2)
- [MOFA+ tutorials](https://biofam.github.io/MOFA2/tutorials.html)
- [MOFA+ vignette for multi-group](https://raw.githack.com/bioFAM/MOFA2_tutorials/master/R_tutorials/getting_started_R.html)

**Status: ✅ SELECTED** — See [Section 4](#4-selected-tools--rationale)

---

### 3.2 Harmony

**What it does:**
Harmony is a batch correction method originally developed for single-cell RNA-seq. It iteratively adjusts sample embeddings (e.g., PCA coordinates) to remove batch effects while preserving biological variation.

**Pros:**
- Very effective at removing technical batch effects
- Works on PCA embeddings — fast and easy to use
- Widely used in single-cell and bulk omics

**Cons:**
- Designed for batch correction, not for modeling multiple biological groups simultaneously
- Can remove real biological signal if exposure effects correlate with dataset
- In this project, tissue type and dataset are confounded with each other — Harmony would struggle to separate technical from biological variation
- Does not identify shared factors — only removes unwanted variation

**Learn more:**
- [Harmony paper — Korsunsky et al. 2019](https://www.nature.com/articles/s41592-019-0619-0)
- [harmonypy Python package](https://github.com/slowkow/harmonypy)

**Status: ❌ NOT SELECTED** — See [Section 5](#5-excluded-tools--why)

---

### 3.3 Combat / Combat-seq

**What it does:**
ComBat is an empirical Bayes method for batch correction. It models and removes batch effects from a data matrix while preserving group-level (biological) differences. ComBat-seq is adapted for count data.

**Pros:**
- Well-validated for methylation array data (450K/EPIC)
- Removes technical batch effects effectively
- Preserves biological group differences when properly specified

**Cons:**
- Requires that batch and biological group are not completely confounded — in this project, each dataset has its own tissue and exposure, making it difficult to separate batch from biology
- Overcorrection is a real risk — can remove the signal of interest
- Not designed for combining datasets with fundamentally different tissues
- Designed for correction before analysis, not for integration itself

**Learn more:**
- [ComBat paper — Johnson et al. 2007](https://academic.oup.com/biostatistics/article/8/1/118/252073)
- [sva package (R)](https://bioconductor.org/packages/release/bioc/html/sva.html)
- [pyComBat (Python)](https://github.com/epigenelabs/pyComBat)

**Status: ❌ NOT SELECTED** — See [Section 5](#5-excluded-tools--why)

---

### 3.4 Multi-Task Learning

**What it does:**
Multi-task learning trains a single model on multiple related tasks simultaneously, sharing a common representation while allowing task-specific components. In this context, the shared task would be predicting "exposed vs. control" across all datasets, with dataset-specific output layers.

**Pros:**
- Explicitly models shared structure across exposures
- Can improve performance on small datasets by leveraging shared signal
- Identifies CpGs important across all tasks (shared features)

**Cons:**
- Requires custom implementation (no off-the-shelf tool for this exact use case)
- Assumes the tasks (predicting exposure in each dataset) are related — may not hold if tissues are too different
- Needs careful architecture design
- With small N in most datasets, even the shared representation may overfit
- More complex to validate and interpret than MOFA+

**Learn more:**
- [Multi-task learning overview](https://ruder.io/multi-task/)
- [PyTorch multi-task tutorial](https://pytorch.org/tutorials/)

**Status: ⚠️ CONSIDERED as secondary approach** — superseded by MOFA+ for this project

---

### 3.5 Meta-Analysis (Statistical)

**What it does:**
Meta-analysis combines statistical results (effect sizes, p-values) from separate analyses of each dataset using fixed-effects or random-effects models. It does not combine raw data — it combines results.

**In this context:**
- Run ElasticNet or t-test separately on each dataset
- Combine effect sizes (Cohen's d or beta coefficients) per CpG across datasets
- Fixed-effects model: assumes true effect is the same in all datasets
- Random-effects model: allows effect to vary across datasets

**Pros:**
- Statistically rigorous and well-validated
- Does not require raw data alignment across datasets
- Handles different sample sizes naturally (weights by precision)
- Can be applied CpG-by-CpG across the region
- Already used in your previous analysis (metafor in R)

**Cons:**
- Operates on summary statistics, not raw data — loses power compared to joint modeling
- Requires same CpGs to be tested in each dataset
- Fixed-effects assumption may not hold if tissues are very different
- Does not identify latent shared factors (MOFA+ does this better)

**Learn more:**
- [metafor R package](https://www.metafor-project.org/)
- [statsmodels meta-analysis (Python)](https://www.statsmodels.org/)
- [Borenstein et al. (2009) — Introduction to Meta-Analysis](https://www.wiley.com/en-us/Introduction+to+Meta+Analysis-p-9780470057247)

**Status: ✅ SELECTED as validation step** — See [Section 4](#4-selected-tools--rationale)

---

## 4. Selected Tools & Rationale

The following tools were selected for this project based on sample sizes, data structure, biological question, and interpretability requirements.

---

### Primary Pipeline

| Step | Tool | Dataset(s) | Reason |
|---|---|---|---|
| Unsupervised exploration | PCA | All | Linear, interpretable, CpG loadings, works with any N |
| Unsupervised exploration | UMAP | All (use cautiously for small N) | Captures non-linear structure, visual separation |
| Spatial pattern visualization | Hierarchical clustering heatmap | All | Reveals spatial methylation patterns within region |
| Supervised classification | ElasticNet (LOOCV) | Wildfire (n=22), Obesity (n=13) | Regularized for small N, handles correlated CpGs, interpretable |
| Supervised classification | ElasticNet (5-fold CV) + Random Forest | Preconception stress (n=70) | Largest dataset — can support RF; ElasticNet for comparison |
| Multi-dataset integration | MOFA+ | All datasets combined | Explicitly models multiple groups, shared vs. specific factors |
| Cross-dataset validation | Meta-analysis (per CpG) | All | Validates MOFA findings with statistical rigor |

---

### Why ElasticNet over LASSO for small N:
CpG positions within a genomic region are spatially correlated — nearby sites tend to have similar methylation levels. LASSO, when faced with correlated features, arbitrarily selects one and discards the rest. This loses spatial information that is biologically meaningful (e.g., a sub-region of the CpG island that is consistently differentially methylated). ElasticNet's L2 component distributes coefficients across correlated CpGs, preserving spatial resolution.

### Why MOFA+ over Harmony/ComBat for integration:
Harmony and ComBat are batch correction tools — they remove unwanted variation. Our goal is the opposite: we want to *find* shared biological variation across datasets while acknowledging that each dataset has its own structure. MOFA+ models this explicitly through latent factor decomposition. It treats each dataset as a "group" and asks: which factors explain variation across all groups? Which are specific to one group? This matches our biological question precisely.

### Why not deep learning:
The total number of samples across all datasets is approximately 120. Deep learning models require hundreds to thousands of samples per class to generalize. Applying neural networks here would produce models that are severely overfit, unvalidatable, and uninterpretable — the opposite of what this project needs.

---

## 5. Excluded Tools & Why

| Tool | Reason for Exclusion |
|---|---|
| t-SNE | Superseded by UMAP; does not preserve global structure; slow |
| k-Means | Requires specifying k; assumes spherical clusters; not suitable for small unbalanced datasets |
| SVM | No native feature importance; requires extensive tuning; not interpretable |
| XGBoost / LightGBM | Requires large N; many hyperparameters; overkill for single-region analysis |
| Neural Networks | Completely inappropriate for N~120 total; uninterpretable; guaranteed to overfit |
| LASSO | Arbitrarily selects one of correlated CpGs — loses spatial information; ElasticNet is strictly better here |
| Harmony | Batch correction tool, not integration tool; risk of removing real biological signal |
| ComBat | Same issue as Harmony; tissue and dataset are confounded, making correction unreliable |
| Multi-Task Learning | MOFA+ achieves the same goal more rigorously with less implementation complexity |

---

*Document version: 1.0 | Last updated: May 2026*
*Author: Ensiuc | Project: CpG_Methylation_ML*

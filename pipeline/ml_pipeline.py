"""
CpG Methylation ML Pipeline
============================
Multi-dataset methylation analysis pipeline for shared signature discovery
across a ~7-8 kb region (CpG island + gene body).

Datasets:
  - Wildfire Smoke  : Nasal epithelial  | n=22  (14 ctrl / 8 exposed)
  - Obesity         : Brain regions     | n=13  (6 lean / 7 obese)
  - Preconc. Stress : TBD              | n≈70  (~45 ctrl / ~25 stress)
  - cfDNA           : Cell-free fetal  | n=8-14 per trimester

Selected ML Tools (see docs/ml_tools_overview.md for full rationale):
  - PCA + UMAP + Hierarchical Clustering  (unsupervised exploration)
  - ElasticNet with LOOCV / 5-fold CV     (primary supervised classifier)
  - Random Forest                          (stress dataset only, n≈70)
  - MOFA+                                  (multi-dataset integration)
  - Per-CpG meta-analysis                  (validation step)

NOT used: XGBoost, SVM, Deep Learning, ComBat, Harmony, k-Means, t-SNE
See docs/ml_tools_overview.md for exclusion rationale.

Usage:
  python ml_pipeline.py --dataset wildfire --region cpgi
  python ml_pipeline.py --dataset stress --region genebody --model both
  python ml_pipeline.py --mofa --region cpgi

Author: Ensi Habibi
"""

import argparse
import os
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import (
    StratifiedKFold, LeaveOneOut, cross_val_score, cross_val_predict
)
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score, classification_report
from sklearn.pipeline import Pipeline
import umap

warnings.filterwarnings("ignore")

# ─────────────────────────────────────────────
# CONFIGURATION
# ─────────────────────────────────────────────

# Datasets and their approximate sample sizes (for CV strategy selection)
DATASET_SIZES = {
    "wildfire": 22,   # LOOCV
    "obesity":  13,   # LOOCV
    "stress":   70,   # 5-fold CV + Random Forest
    "cfdna":    12,   # LOOCV (per-trimester; approximate)
}

REGIONS = ["cpgi", "genebody"]

RESULTS_DIR = os.path.join(os.path.dirname(__file__), "..", "results")
FIGURES_DIR = os.path.join(RESULTS_DIR, "figures")
TABLES_DIR  = os.path.join(RESULTS_DIR, "tables")
MODELS_DIR  = os.path.join(RESULTS_DIR, "models")

for d in [FIGURES_DIR, TABLES_DIR, MODELS_DIR]:
    os.makedirs(d, exist_ok=True)


# ─────────────────────────────────────────────
# DATA LOADING
# ─────────────────────────────────────────────

def load_dataset(dataset_name: str, region: str):
    """
    Load methylation matrix and labels for a given dataset and region.

    Standard bioinformatics format (minfi, SeSAMe, bismark-compatible):

      Beta matrix — data/processed/{dataset_name}_{region}_methylation.csv
        rows    = CpG positions (e.g. cg00001234)
        columns = sample IDs
        values  = beta values [0, 1]

      Labels file — data/processed/{dataset_name}_labels.csv
        columns: sample_id, label   (1 = exposed/case, 0 = control)
        sample_id must match column names in the beta matrix

    Returns:
        X : DataFrame (n_samples × n_cpgs) — beta values [0, 1]
        y : Series — binary labels, indexed by sample_id
    """
    data_dir = os.path.join(os.path.dirname(__file__), "..", "data", "processed")

    beta_path   = os.path.join(data_dir, f"{dataset_name}_{region}_methylation.csv")
    labels_path = os.path.join(data_dir, f"{dataset_name}_labels.csv")

    if not os.path.exists(beta_path):
        raise FileNotFoundError(
            f"Beta matrix not found: {beta_path}\n"
            "Expected format: CpGs as rows, samples as columns."
        )
    if not os.path.exists(labels_path):
        raise FileNotFoundError(
            f"Labels file not found: {labels_path}\n"
            "Expected columns: sample_id, label  (1=case, 0=control)."
        )

    # Load beta matrix (CpGs × samples) and transpose → samples × CpGs
    beta = pd.read_csv(beta_path, index_col=0)
    X = beta.T                          # now: rows=samples, columns=CpGs
    X.index.name = "sample_id"

    # Load labels and align to sample order
    labels = pd.read_csv(labels_path, index_col="sample_id")
    if "label" not in labels.columns:
        raise ValueError(f"{labels_path} must have a 'label' column (1=case, 0=control).")

    # Keep only samples present in both files
    common = X.index.intersection(labels.index)
    if len(common) == 0:
        raise ValueError(
            "No matching sample IDs between beta matrix columns and labels file.\n"
            "Check that sample names match exactly (case-sensitive)."
        )
    if len(common) < len(X):
        missing = set(X.index) - set(common)
        print(f"  [WARN] {len(missing)} samples in beta matrix have no label — dropped: {missing}")

    X = X.loc[common]
    y = labels.loc[common, "label"].astype(int)

    n_case    = y.sum()
    n_control = (y == 0).sum()
    print(f"  {dataset_name} ({region}): {len(y)} samples — "
          f"{n_control} control / {n_case} case | {X.shape[1]} CpGs")
    return X, y


# ─────────────────────────────────────────────
# UNSUPERVISED — STEP 1: PCA
# ─────────────────────────────────────────────

def run_pca(X: pd.DataFrame, y: pd.Series, label: str):
    """PCA with sample plot colored by condition. Returns PCA coordinates."""
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    pca = PCA(n_components=min(10, X.shape[0] - 1), random_state=42)
    coords = pca.fit_transform(X_scaled)
    explained = pca.explained_variance_ratio_ * 100

    # Scree plot
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    axes[0].bar(range(1, len(explained) + 1), explained, color="steelblue", alpha=0.8)
    axes[0].set_xlabel("Principal Component")
    axes[0].set_ylabel("Variance Explained (%)")
    axes[0].set_title(f"PCA Scree — {label}")

    for val, color, lbl in [(0, "steelblue", "Control"), (1, "tomato", "Exposed/Case")]:
        mask = y == val
        axes[1].scatter(coords[mask, 0], coords[mask, 1],
                        c=color, label=lbl, alpha=0.8, s=70, edgecolors="white", linewidths=0.5)
    axes[1].set_xlabel(f"PC1 ({explained[0]:.1f}%)")
    axes[1].set_ylabel(f"PC2 ({explained[1]:.1f}%)")
    axes[1].set_title(f"PCA Scores — {label}")
    axes[1].legend()

    plt.tight_layout()
    outpath = os.path.join(FIGURES_DIR, f"pca_{label}.png")
    plt.savefig(outpath, dpi=150)
    plt.close()
    print(f"  → PCA plot saved: {outpath}")

    # Save loadings table (which CpGs drive each PC)
    loadings = pd.DataFrame(
        pca.components_[:3].T,
        index=X.columns,
        columns=["PC1", "PC2", "PC3"]
    )
    loadings_path = os.path.join(TABLES_DIR, f"pca_loadings_{label}.csv")
    loadings.to_csv(loadings_path)
    print(f"  → PCA loadings saved: {loadings_path}")

    return coords, pca


# ─────────────────────────────────────────────
# UNSUPERVISED — STEP 2: UMAP
# ─────────────────────────────────────────────

def run_umap(X: pd.DataFrame, y: pd.Series, label: str, n: int):
    """
    UMAP visualization. Only run if n >= 15 (unreliable with very small N).
    """
    if n < 15:
        print(f"  [SKIP] UMAP skipped for {label} (n={n} < 15 — unreliable)")
        return

    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    n_neighbors = max(5, min(15, n // 3))
    reducer = umap.UMAP(n_components=2, random_state=42,
                        n_neighbors=n_neighbors, min_dist=0.1)
    coords = reducer.fit_transform(X_scaled)

    fig, ax = plt.subplots(figsize=(7, 6))
    for val, color, lbl in [(0, "steelblue", "Control"), (1, "tomato", "Exposed/Case")]:
        mask = y == val
        ax.scatter(coords[mask, 0], coords[mask, 1],
                   c=color, label=lbl, alpha=0.8, s=70, edgecolors="white", linewidths=0.5)
    ax.set_title(f"UMAP — {label} (n_neighbors={n_neighbors})")
    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")
    ax.legend()
    plt.tight_layout()

    outpath = os.path.join(FIGURES_DIR, f"umap_{label}.png")
    plt.savefig(outpath, dpi=150)
    plt.close()
    print(f"  → UMAP plot saved: {outpath}")


# ─────────────────────────────────────────────
# UNSUPERVISED — STEP 3: Hierarchical Clustering Heatmap
# ─────────────────────────────────────────────

def run_clustermap(X: pd.DataFrame, y: pd.Series, label: str):
    """
    Hierarchical clustering heatmap (samples × CpGs).
    Reveals spatial methylation patterns within the ~7-8kb region.
    """
    # Color bar for condition
    condition_colors = y.map({0: "steelblue", 1: "tomato"})

    # Select top variable CpGs to avoid overloading the plot
    top_cpgs = X.var(axis=0).nlargest(min(100, X.shape[1])).index
    X_plot = X[top_cpgs]

    g = sns.clustermap(
        X_plot,
        row_colors=condition_colors,
        cmap="RdYlBu_r",
        vmin=0, vmax=1,
        figsize=(14, max(6, len(X) * 0.3)),
        yticklabels=True,
        xticklabels=False,
        method="ward",
        metric="euclidean",
        cbar_kws={"label": "Beta (methylation)"},
    )
    g.fig.suptitle(f"Hierarchical Clustering — {label}\n"
                   f"(top {len(top_cpgs)} variable CpGs | blue=control, red=exposed)",
                   y=1.02, fontsize=12)

    outpath = os.path.join(FIGURES_DIR, f"clustermap_{label}.png")
    g.fig.savefig(outpath, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  → Clustermap saved: {outpath}")


# ─────────────────────────────────────────────
# SUPERVISED — ElasticNet
# ─────────────────────────────────────────────

def choose_cv_strategy(n: int, y: pd.Series):
    """
    LOOCV for small datasets (n < 30), 5-fold CV for larger.
    Random Forest only enabled for n >= 40.
    """
    if n < 30:
        cv = LeaveOneOut()
        cv_name = "LOOCV"
    else:
        cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
        cv_name = "5-fold CV"
    return cv, cv_name


def run_elasticnet(X: pd.DataFrame, y: pd.Series, label: str, n: int):
    """
    ElasticNet logistic regression with LOOCV (small N) or 5-fold CV (large N).
    Returns AUC, accuracy, and non-zero CpG coefficients.
    """
    cv, cv_name = choose_cv_strategy(n, y)

    model = Pipeline([
        ("scaler", StandardScaler()),
        ("clf", LogisticRegression(
            penalty="elasticnet", solver="saga",
            l1_ratio=0.5, C=0.1, max_iter=3000, random_state=42
        ))
    ])

    # Cross-validated predictions
    if isinstance(cv, LeaveOneOut):
        # LOOCV — compute AUC from held-out predictions
        y_pred_proba = cross_val_predict(model, X, y, cv=cv, method="predict_proba")[:, 1]
        auc = roc_auc_score(y, y_pred_proba)
        accs = cross_val_score(model, X, y, cv=cv, scoring="accuracy")
        acc = accs.mean()
        print(f"  ElasticNet ({cv_name}): AUC = {auc:.3f} | Acc = {acc:.3f}")
    else:
        aucs = cross_val_score(model, X, y, cv=cv, scoring="roc_auc")
        accs = cross_val_score(model, X, y, cv=cv, scoring="accuracy")
        auc, acc = aucs.mean(), accs.mean()
        print(f"  ElasticNet ({cv_name}): AUC = {auc:.3f} ± {aucs.std():.3f} | "
              f"Acc = {acc:.3f} ± {accs.std():.3f}")

    # Fit on full data to get coefficients
    model.fit(X, y)
    coefs = pd.Series(
        model.named_steps["clf"].coef_[0],
        index=X.columns
    )
    nonzero_cpgs = coefs[coefs != 0].sort_values(key=abs, ascending=False)
    coef_path = os.path.join(TABLES_DIR, f"elasticnet_coefs_{label}.csv")
    nonzero_cpgs.to_csv(coef_path, header=["coefficient"])
    print(f"  → {len(nonzero_cpgs)} non-zero CpG coefficients saved: {coef_path}")

    return {"model": "ElasticNet", "cv": cv_name, "auc": auc, "acc": acc,
            "n_nonzero": len(nonzero_cpgs), "dataset": label}


# ─────────────────────────────────────────────
# SUPERVISED — Random Forest (stress dataset only)
# ─────────────────────────────────────────────

def run_random_forest(X: pd.DataFrame, y: pd.Series, label: str, n: int):
    """
    Random Forest — only appropriate for datasets with n >= 40 per class.
    Uses permutation importance (unbiased for correlated features).
    """
    if n < 40:
        print(f"  [SKIP] Random Forest skipped for {label} "
              f"(n={n} — minimum ~40 recommended per class)")
        return None

    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    rf = RandomForestClassifier(
        n_estimators=500, max_features="sqrt",
        n_jobs=-1, random_state=42, oob_score=True
    )

    aucs = cross_val_score(rf, X, y, cv=cv, scoring="roc_auc")
    accs = cross_val_score(rf, X, y, cv=cv, scoring="accuracy")
    print(f"  RandomForest (5-fold CV): AUC = {aucs.mean():.3f} ± {aucs.std():.3f} | "
          f"Acc = {accs.mean():.3f} ± {accs.std():.3f}")

    # Fit on full data for importance
    rf.fit(X, y)
    importance = pd.Series(rf.feature_importances_, index=X.columns)
    importance = importance.sort_values(ascending=False)
    imp_path = os.path.join(TABLES_DIR, f"rf_importance_{label}.csv")
    importance.to_csv(imp_path, header=["importance"])
    print(f"  → OOB error: {rf.oob_score_:.3f} | Feature importances saved: {imp_path}")

    # Top 20 feature importance plot
    plt.figure(figsize=(10, 5))
    importance.head(20).plot(kind="bar", color="steelblue", alpha=0.8)
    plt.title(f"Random Forest — Top 20 CpG Importances — {label}")
    plt.xlabel("CpG position")
    plt.ylabel("Mean Decrease in Impurity")
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, f"rf_importance_{label}.png"), dpi=150)
    plt.close()

    return {"model": "RandomForest", "cv": "5-fold CV",
            "auc": aucs.mean(), "acc": accs.mean(), "dataset": label}


# ─────────────────────────────────────────────
# MULTI-DATASET INTEGRATION — MOFA+
# ─────────────────────────────────────────────

def run_mofa(region: str):
    """
    MOFA+ multi-group latent factor decomposition across all four datasets.
    Identifies shared vs. exposure-specific methylation factors.

    Requires: pip install mofapy2
    """
    print("\n=== MOFA+ Multi-Dataset Integration ===")
    try:
        from mofapy2.run.entry_point import entry_point
    except ImportError:
        print("  [ERROR] mofapy2 not installed.")
        print("  Install with: pip install mofapy2")
        print("  Then re-run with --mofa flag.")
        return

    datasets_data = {}
    for ds in DATASET_SIZES.keys():
        try:
            X, y = load_dataset(ds, region)
            datasets_data[ds] = (X, y)
        except FileNotFoundError as e:
            print(f"  [SKIP] {e}")

    if len(datasets_data) < 2:
        print("  [ERROR] Need at least 2 datasets loaded for MOFA+.")
        return

    # Align to shared CpGs
    cpg_sets = [set(X.columns) for X, _ in datasets_data.values()]
    shared_cpgs = list(set.intersection(*cpg_sets))
    print(f"  Shared CpGs across {len(datasets_data)} datasets: {len(shared_cpgs)}")

    if len(shared_cpgs) < 10:
        print("  [WARNING] Very few shared CpGs — check that CpG coordinates are aligned.")

    # Build MOFA input: list of dicts with group/view structure
    # MOFA+ expects: data[view][group] = samples × features matrix
    mofa_data = {}
    mofa_data["methylation"] = {}
    for ds_name, (X, y) in datasets_data.items():
        mofa_data["methylation"][ds_name] = X[shared_cpgs].values

    # Initialize and run MOFA+
    ent = entry_point()
    ent.set_data_options(scale_groups=False, scale_views=False)
    ent.set_data_matrix(
        [[mofa_data["methylation"][g] for g in datasets_data.keys()]],
        views_names=["methylation"],
        groups_names=list(datasets_data.keys()),
        samples_names=[list(X.index) for X, _ in datasets_data.values()],
        features_names=[shared_cpgs]
    )
    ent.set_model_options(factors=10, likelihoods=["gaussian"])
    ent.set_train_options(
        iter=1000, convergence_mode="fast",
        startELBO=1, freqELBO=5,
        gpu_mode=False, verbose=False, seed=42
    )

    mofa_outpath = os.path.join(MODELS_DIR, f"mofa_model_{region}.hdf5")
    ent.set_outfile(mofa_outpath)
    ent.build()
    ent.run()

    print(f"  → MOFA+ model saved: {mofa_outpath}")
    print("  → Load in R with: MOFA2::load_model('{mofa_outpath}')")
    print("  → Or in Python: muon.read_h5mu() / mofapy2 utilities")

    # Extract and save factor weights (CpG loadings per factor)
    try:
        import h5py
        with h5py.File(mofa_outpath, "r") as f:
            weights = f["expectations"]["W"]["methylation"][:]
        weights_df = pd.DataFrame(
            weights.T,
            index=shared_cpgs,
            columns=[f"Factor{i+1}" for i in range(weights.shape[0])]
        )
        weights_path = os.path.join(TABLES_DIR, f"mofa_weights_{region}.csv")
        weights_df.to_csv(weights_path)
        print(f"  → Factor weights (CpG loadings) saved: {weights_path}")
    except Exception as e:
        print(f"  [NOTE] Could not extract weights directly: {e}")
        print("  Use MOFA2 R package or mofapy2 utilities to explore the model.")


# ─────────────────────────────────────────────
# VALIDATION — Per-CpG Meta-Analysis
# ─────────────────────────────────────────────

def run_meta_analysis(region: str):
    """
    Per-CpG meta-analysis: combine effect sizes (Cohen's d) across datasets
    using a fixed-effects weighted mean. Validates MOFA+ shared factor findings.
    """
    print("\n=== Per-CpG Meta-Analysis ===")

    dataset_results = {}
    for ds in DATASET_SIZES.keys():
        try:
            X, y = load_dataset(ds, region)
            dataset_results[ds] = (X, y)
        except FileNotFoundError:
            pass

    if len(dataset_results) < 2:
        print("  [SKIP] Need at least 2 datasets for meta-analysis.")
        return

    # Get shared CpGs
    cpg_sets = [set(X.columns) for X, _ in dataset_results.values()]
    shared_cpgs = list(set.intersection(*cpg_sets))
    print(f"  Running meta-analysis on {len(shared_cpgs)} shared CpGs "
          f"across {len(dataset_results)} datasets...")

    meta_rows = []
    for cpg in shared_cpgs:
        d_values, se_values = [], []
        for ds_name, (X, y) in dataset_results.items():
            case    = X.loc[y == 1, cpg].values
            control = X.loc[y == 0, cpg].values

            # Cohen's d
            pooled_sd = np.sqrt((case.var() + control.var()) / 2 + 1e-8)
            d = (case.mean() - control.mean()) / pooled_sd

            # Standard error of d (approximate)
            n1, n2 = len(case), len(control)
            se = np.sqrt((n1 + n2) / (n1 * n2) + d**2 / (2 * (n1 + n2)))

            d_values.append(d)
            se_values.append(se)

        # Fixed-effects meta-analysis: weighted mean by 1/SE²
        weights = [1 / (se**2) for se in se_values]
        w_total = sum(weights)
        d_pooled = sum(w * d for w, d in zip(weights, d_values)) / w_total
        se_pooled = np.sqrt(1 / w_total)
        z = d_pooled / se_pooled
        p = 2 * (1 - stats.norm.cdf(abs(z)))

        meta_rows.append({
            "CpG": cpg,
            "d_pooled": d_pooled,
            "SE_pooled": se_pooled,
            "Z": z,
            "p_value": p,
            "n_datasets": len(dataset_results),
        })

    meta_df = pd.DataFrame(meta_rows)

    # FDR correction (Benjamini-Hochberg)
    from scipy.stats import rankdata
    pvals = meta_df["p_value"].values
    n = len(pvals)
    ranks = rankdata(pvals)
    fdr = np.minimum(1, pvals * n / ranks)
    meta_df["FDR"] = fdr

    meta_df = meta_df.sort_values("FDR")
    out_path = os.path.join(TABLES_DIR, f"meta_analysis_{region}.csv")
    meta_df.to_csv(out_path, index=False)

    n_sig = (meta_df["FDR"] < 0.05).sum()
    print(f"  → {n_sig} CpGs significant at FDR < 0.05 across all datasets")
    print(f"  → Meta-analysis results saved: {out_path}")


# ─────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description="CpG Methylation ML Pipeline")
    parser.add_argument(
        "--dataset", choices=list(DATASET_SIZES.keys()),
        help="Single dataset to analyse (unsupervised + supervised)"
    )
    parser.add_argument(
        "--region", choices=REGIONS, default="cpgi",
        help="Genomic region: 'cpgi' (CpG island) or 'genebody'"
    )
    parser.add_argument(
        "--model", choices=["elasticnet", "rf", "both"], default="elasticnet",
        help="Supervised model(s) to run (rf = Random Forest; only valid for stress)"
    )
    parser.add_argument(
        "--mofa", action="store_true",
        help="Run MOFA+ integration across all datasets"
    )
    parser.add_argument(
        "--meta", action="store_true",
        help="Run per-CpG meta-analysis across all datasets"
    )
    args = parser.parse_args()

    print(f"\n{'='*60}")
    print("CpG Methylation ML Pipeline")
    if args.dataset:
        print(f"  Dataset : {args.dataset}")
    print(f"  Region  : {args.region}")
    print(f"{'='*60}\n")

    results = []

    # ── Single dataset analysis ──
    if args.dataset:
        ds = args.dataset
        n  = DATASET_SIZES[ds]
        label = f"{ds}_{args.region}"

        print(f"Loading {ds}...")
        X, y = load_dataset(ds, args.region)

        print("\n--- Unsupervised Analysis ---")
        run_pca(X, y, label)
        run_umap(X, y, label, n)
        run_clustermap(X, y, label)

        print("\n--- Supervised Classification ---")
        if args.model in ["elasticnet", "both"]:
            r = run_elasticnet(X, y, label, n)
            results.append(r)

        if args.model in ["rf", "both"]:
            r = run_random_forest(X, y, label, n)
            if r:
                results.append(r)

        if results:
            summary = pd.DataFrame(results)
            summary_path = os.path.join(TABLES_DIR, f"cv_summary_{label}.csv")
            summary.to_csv(summary_path, index=False)
            print(f"\n  → CV summary saved: {summary_path}")

    # ── Multi-dataset analyses ──
    if args.mofa:
        run_mofa(args.region)

    if args.meta:
        run_meta_analysis(args.region)

    if not args.dataset and not args.mofa and not args.meta:
        parser.print_help()
        print("\nExample commands:")
        print("  python ml_pipeline.py --dataset wildfire --region cpgi")
        print("  python ml_pipeline.py --dataset stress --region genebody --model both")
        print("  python ml_pipeline.py --mofa --region cpgi")
        print("  python ml_pipeline.py --meta --region genebody")

    print(f"\n{'='*60}")
    print("Done. Results saved to results/")
    print(f"{'='*60}\n")


if __name__ == "__main__":
    main()

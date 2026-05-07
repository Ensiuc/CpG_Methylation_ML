"""
preprocess_wgbs.py
==================
Converts per-sample Bismark ROI TSV files into ML-ready beta matrices.

What this script does, step by step:
  1. Loops over all four datasets
  2. For each dataset, finds all *.roi1kb.tsv.gz sample files
  3. Loads each sample, filters to region of interest, computes beta values
  4. Combines all samples into one CpG × samples matrix
  5. Drops CpGs that are missing in too many samples
  6. Splits into CpG island and gene body sub-matrices
  7. Saves three files per dataset:
       {dataset}_cpgi_methylation.csv      (CpGs × samples, beta values)
       {dataset}_genebody_methylation.csv  (CpGs × samples, beta values)
       {dataset}_labels.csv               (sample_id, label: 1=case, 0=control)

Output files go to: CpG_Methylation_ML/data/processed/
These are exactly what ml_pipeline.py expects as input.

Author: Ensi Habibi
"""

import os
import glob
import numpy as np
import pandas as pd
from pathlib import Path

# ============================================================
# SECTION 1: CONFIGURATION
# ============================================================
# Think of this block as the "settings panel" — the only part
# you ever need to change if paths or parameters shift.

# Parent directory that contains all four dataset folders
PARENT = Path("~/project/nhip_macaque").expanduser()
# expanduser() converts ~ to your actual home directory path
# e.g. /home/ensi/project/nhip_macaque

# Where to save the processed output files
# __file__ is the path of this script itself; we go one level up (..)
# to reach the CpG_Methylation_ML root, then into data/processed/
OUT_DIR = Path(__file__).resolve().parent.parent / "data" / "processed"
OUT_DIR.mkdir(parents=True, exist_ok=True)
# mkdir(parents=True)  → creates all intermediate folders if needed
# mkdir(exist_ok=True) → no error if the folder already exists

# ── Genomic region of interest ──────────────────────────────
CHROM        = "chr10"
REGION_START = 2432261   # full region start (used when loading TSV)
REGION_END   = 2443263   # full region end

# Sub-feature windows (used when splitting the final matrix)
CPG_START      = 2433175
CPG_END        = 2433562
GENEBODY_START = REGION_START + 1250  # 2433511
GENEBODY_END   = REGION_END - 1000   # 2442263

# ── Filtering parameters ─────────────────────────────────────
MIN_COV     = 1     # minimum read coverage per CpG per sample
MAX_MISSING = 0.20  # drop CpGs missing in more than 20% of samples
# Why 20%? With small N (e.g. n=13), losing even one sample at a CpG
# is already ~8% missing. Being too strict would remove most CpGs.
# Adjust lower (e.g. 0.0) if you want only fully covered CpGs.

# ── Dataset configuration ────────────────────────────────────
# Keys (e.g. "stress") become the prefix in output filenames:
#   stress_cpgi_methylation.csv, stress_labels.csv, etc.
DATASETS = {
    "stress": {
        "project":    PARENT / "nhip_cumulus",
        "case_label": "Stressed",   # this group gets label = 1
        "ctrl_label": "Control",    # this group gets label = 0
    },
    "obesity": {
        "project":    PARENT / "fromBen_MacaqueObese_brain",
        "case_label": "Obese",
        "ctrl_label": "Control",
    },
    "cfdna": {
        "project":    PARENT / "fromBen_MacaqueObese_cfDna",
        "case_label": "Obese",
        "ctrl_label": "Control",
    },
    "wildfire": {
        "project":    PARENT / "macaque_nasal_HongJi",
        "case_label": "Exposed",
        "ctrl_label": "Control",
    },
}


# ============================================================
# SECTION 2: HELPER FUNCTIONS
# ============================================================

def parse_sample_name(fpath: str) -> str:
    """
    Extract the sample ID from a Bismark ROI filename.

    Filenames look like:
      SAMPLE_merged_name_sorted.deduplicated.bismark.cov.gz.CpG_report...roi1kb.tsv.gz

    Strategy: strip the known suffix, then take everything before
    the first underscore — that's the sample ID.

    Example:
      "NHIP01_merged_...roi1kb.tsv.gz"  →  "NHIP01"
    """
    base = os.path.basename(fpath)                  # strip directory path
    base = base.replace(".roi1kb.tsv.gz", "")       # remove ROI suffix
    base = base.replace(".cov.gz", "")              # or plain .cov.gz suffix
    return base.split("_", 1)[0]                    # take first token only
    # split("_", 1) splits on the FIRST underscore only (maxsplit=1)
    # [0] takes the left side → sample ID


def load_sample(fpath: str) -> pd.DataFrame:
    """
    Load one sample's ROI TSV file and return a clean beta-value table.

    Input file columns (no header):
      chrom | start | end | pct_methylation | meth_reads | unmeth_reads

    What we do:
      1. Read the file
      2. Filter to our chromosome and genomic region
      3. Compute coverage = meth + unmeth
      4. Compute beta = meth / coverage  (value between 0 and 1)
      5. Drop CpGs below minimum coverage
      6. Return a DataFrame with columns: start, beta

    Returns an EMPTY DataFrame if the sample has no usable data in region.
    """
    cols = ["chrom", "start", "end", "pct", "meth", "unmeth"]
    df = pd.read_csv(
        fpath,
        sep="\t",           # tab-separated
        header=None,        # no column names in file
        names=cols,         # we assign names ourselves
        compression="gzip"  # file is gzip-compressed (.tsv.gz)
    )

    if df.empty:
        return df

    # Filter to the right chromosome
    df = df[df["chrom"] == CHROM].copy()
    # .copy() prevents a SettingWithCopyWarning — always use it after filtering

    if df.empty:
        return df

    # Filter to the genomic region of interest
    # We use 'start' position as the CpG identifier (standard in WGBS)
    df = df[
        (df["start"] >= REGION_START) &
        (df["start"] <= REGION_END)
    ].copy()

    if df.empty:
        return df

    # Compute coverage and beta value
    df["cov"]  = df["meth"].astype(float) + df["unmeth"].astype(float)
    df["beta"] = np.where(
        df["cov"] > 0,                              # condition: covered
        df["meth"].astype(float) / df["cov"],       # if yes: compute beta
        np.nan                                       # if no: set to NaN
    )
    # np.where(condition, value_if_true, value_if_false)
    # This is safer than simple division — avoids dividing by zero

    # Apply minimum coverage filter
    df = df[df["cov"] >= MIN_COV].copy()

    # Drop CpGs where beta couldn't be computed
    df = df.dropna(subset=["beta"])

    # Return only what we need: position and beta value
    return df[["start", "beta"]].reset_index(drop=True)


def build_matrix(files: list, name_to_group: dict,
                 case_label: str, ctrl_label: str):
    """
    Load all samples and combine into one CpG × samples matrix.

    Logic:
      - For each file: parse sample name, load data, store as a Series
        (index = CpG start position, value = beta)
      - Combine all Series into a DataFrame
        → rows = CpG positions, columns = sample names
      - Also build a labels dict: sample_name → 0 or 1

    Returns:
      beta_matrix : DataFrame (CpGs × samples), may contain NaN
      labels      : dict {sample_name: 0 or 1}
    """
    sample_series = {}   # will hold one pd.Series per sample
    labels        = {}   # will hold label (0 or 1) per sample

    for fpath in files:
        sample_name = parse_sample_name(fpath)

        # Check sample is in metadata
        if sample_name not in name_to_group:
            print(f"  [SKIP] {sample_name} — not found in metadata")
            continue

        group = name_to_group[sample_name]

        # Check group is one we expect
        if group not in (case_label, ctrl_label):
            print(f"  [SKIP] {sample_name} — unexpected group: {group}")
            continue

        # Load the sample
        df = load_sample(fpath)
        if df.empty:
            print(f"  [SKIP] {sample_name} — no data in region")
            continue

        # Convert to Series: index = CpG position, value = beta
        s = df.set_index("start")["beta"]
        s.name = sample_name
        # set_index("start") makes the CpG position the row label
        # ["beta"] selects just the beta column → becomes a Series

        sample_series[sample_name] = s
        labels[sample_name] = 1 if group == case_label else 0
        # case = 1 (exposed/disease), control = 0 (standard encoding in ML)

    if not sample_series:
        raise ValueError("No samples were successfully loaded. Check paths and metadata.")

    # Combine all samples into one DataFrame
    # pd.DataFrame(dict_of_series) aligns on index automatically
    # Missing CpG positions (not covered in that sample) become NaN
    beta_matrix = pd.DataFrame(sample_series)
    # Result: rows = CpG start positions, columns = sample names

    print(f"  Raw matrix: {beta_matrix.shape[0]} CpGs × {beta_matrix.shape[1]} samples")

    return beta_matrix, labels


def filter_missing(beta_matrix: pd.DataFrame) -> pd.DataFrame:
    """
    Drop CpGs that are missing (NaN) in more than MAX_MISSING fraction of samples.

    Example with MAX_MISSING=0.20 and 13 samples:
      A CpG is kept only if it has valid beta values in at least 11 samples.

    Why do we need this?
      Not every sample will have coverage at every CpG.
      If a CpG is missing in most samples, it's not useful for ML.
      But being too strict would leave us with very few CpGs.
    """
    n_samples    = beta_matrix.shape[1]
    missing_frac = beta_matrix.isna().sum(axis=1) / n_samples
    # isna()         → True/False matrix of missing values
    # .sum(axis=1)   → count missing per ROW (per CpG)
    # / n_samples    → fraction missing

    kept = beta_matrix[missing_frac <= MAX_MISSING].copy()
    dropped = beta_matrix.shape[0] - kept.shape[0]
    print(f"  Coverage filter: dropped {dropped} CpGs "
          f"(>{int(MAX_MISSING*100)}% missing) → {kept.shape[0]} CpGs retained")
    return kept


def subset_to_window(beta_matrix: pd.DataFrame,
                     start: int, end: int, name: str) -> pd.DataFrame:
    """
    Subset the full region matrix to a specific genomic window.

    Uses the row index (CpG start positions) to filter.
    """
    mask   = (beta_matrix.index >= start) & (beta_matrix.index <= end)
    subset = beta_matrix.loc[mask].copy()
    print(f"  {name} window ({start}–{end}): {subset.shape[0]} CpGs")
    return subset


# ============================================================
# SECTION 3: MAIN LOOP — process each dataset
# ============================================================

def process_dataset(dataset_name: str, config: dict):
    """
    Full preprocessing pipeline for one dataset.

    Steps:
      1. Locate all ROI TSV files
      2. Load metadata (sample name → group)
      3. Build beta matrix
      4. Filter missing CpGs
      5. Split into CpG island and gene body sub-matrices
      6. Save all output files
    """
    print(f"\n{'='*60}")
    print(f"Processing: {dataset_name}")
    print(f"{'='*60}")

    project    = config["project"]
    case_label = config["case_label"]
    ctrl_label = config["ctrl_label"]

    # ── Find input files ─────────────────────────────────────
    roi_dir = project / "roi_tsvfiles"
    files   = sorted(glob.glob(str(roi_dir / "*.roi1kb.tsv.gz")))
    if not files:
        files = sorted(glob.glob(str(roi_dir / "*.cov.gz")))
    if not files:
        print(f"  [ERROR] No TSV files found in {roi_dir}. Skipping.")
        return
    print(f"  Found {len(files)} sample files")

    # ── Load metadata ─────────────────────────────────────────
    meta_path = project / "sample_info.csv"
    meta      = pd.read_csv(meta_path, sep="\t")
    # We read with tab separator — adjust to "," if your file is comma-separated

    meta["Name"]  = meta["Name"].astype(str)
    meta["Group"] = meta["Group"].astype(str)
    name_to_group = dict(zip(meta["Name"], meta["Group"]))
    # zip() pairs up two lists: [name1, name2, ...] + [group1, group2, ...]
    # dict() converts those pairs into a lookup dictionary

    n_case = (meta["Group"] == case_label).sum()
    n_ctrl = (meta["Group"] == ctrl_label).sum()
    print(f"  Metadata: {n_ctrl} {ctrl_label} / {n_case} {case_label}")

    # ── Build beta matrix ─────────────────────────────────────
    beta_matrix, labels = build_matrix(files, name_to_group, case_label, ctrl_label)

    # ── Filter missing CpGs ───────────────────────────────────
    beta_matrix = filter_missing(beta_matrix)

    # ── Split into sub-features ───────────────────────────────
    cpgi_matrix     = subset_to_window(beta_matrix, CPG_START,      CPG_END,      "CpG island")
    genebody_matrix = subset_to_window(beta_matrix, GENEBODY_START, GENEBODY_END, "Gene body")

    # ── Save output files ─────────────────────────────────────
    # Beta matrices: CpGs as rows, samples as columns
    cpgi_path     = OUT_DIR / f"{dataset_name}_cpgi_methylation.csv"
    genebody_path = OUT_DIR / f"{dataset_name}_genebody_methylation.csv"
    labels_path   = OUT_DIR / f"{dataset_name}_labels.csv"

    cpgi_matrix.to_csv(cpgi_path)
    genebody_matrix.to_csv(genebody_path)
    # .to_csv() writes the DataFrame index (CpG positions) as the first column
    # and column names (sample IDs) as the header row

    # Labels file: one row per sample
    labels_df = pd.DataFrame(
        [(name, lbl) for name, lbl in labels.items()],
        columns=["sample_id", "label"]
    )
    labels_df.to_csv(labels_path, index=False)
    # index=False → don't write the row numbers (0,1,2,...) to the file

    print(f"\n  Saved:")
    print(f"    {cpgi_path}")
    print(f"    {genebody_path}")
    print(f"    {labels_path}")


# ============================================================
# SECTION 4: ENTRY POINT
# ============================================================
# This block only runs when you execute the script directly:
#   python preprocess_wgbs.py
# It does NOT run if you import this file from another script.
# This is a standard Python pattern you'll see everywhere.

if __name__ == "__main__":
    for dataset_name, config in DATASETS.items():
        process_dataset(dataset_name, config)

    print(f"\n{'='*60}")
    print("All datasets processed.")
    print(f"Output files in: {OUT_DIR}")
    print(f"{'='*60}\n")

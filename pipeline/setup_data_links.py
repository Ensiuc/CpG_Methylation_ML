"""
setup_data_links.py
===================
Creates the standardized data/raw/ directory structure for the ML pipeline.

What this script does:
  - Creates one folder per dataset under data/raw/
  - Symlinks all TSV files from their original project locations
  - Copies (small) metadata CSV files into each dataset folder

Why symlinks and not copies?
  WGBS coverage files are large (hundreds of MB each).
  A symlink is just a pointer — it looks like the file is in data/raw/
  but the actual bytes stay in the original project directory.
  This is standard practice in bioinformatics to avoid duplicating data.

Run once before any analysis:
  python pipeline/setup_data_links.py

After running, every other script just looks in data/raw/{dataset}/
and never needs to know about the original project paths again.

Author: Ensi Habibi
"""

import os
import glob
import shutil
from pathlib import Path

# ============================================================
# SECTION 1: PATHS
# ============================================================

# Root of the CpG_Methylation_ML project
# __file__ = this script = pipeline/setup_data_links.py
# .parent   = pipeline/
# .parent   = CpG_Methylation_ML/   ← project root
PROJECT_ROOT = Path(__file__).resolve().parent.parent

# Where raw data folders will be created
RAW_DIR = PROJECT_ROOT / "data" / "raw"

# Where the original project data lives on the server
NHIP = Path("/project/nhip_macaque")


# ============================================================
# SECTION 2: DATASET CONFIGURATION
# ============================================================
# Each entry defines one dataset:
#   "key"        → becomes the folder name in data/raw/
#   "tsv_dir"    → where to find the *.tsv.gz files to symlink
#   "meta_src"   → original metadata file to copy in
#   "meta_name"  → what to call the metadata file inside data/raw/{key}/
#   "pattern"    → glob pattern to find TSV files

DATASETS = {

    # ── Stress (NHIP preconception stress) ──────────────────
    "stress": {
        "tsv_dir":   NHIP / "nhip_cumulus" / "roi_tsvfiles",
        "meta_src":  NHIP / "nhip_cumulus" / "sample_info.csv",
        "meta_name": "sample_info.csv",
        "pattern":   "*.roi1kb.tsv.gz",
    },

    # ── Wildfire (nasal epithelial) ──────────────────────────
    "wildfire": {
        "tsv_dir":   NHIP / "macaque_nasal_HongJi" / "roi_tsvfiles",
        "meta_src":  NHIP / "macaque_nasal_HongJi" / "sample_info.csv",
        "meta_name": "sample_info.csv",
        "pattern":   "*.roi1kb.tsv.gz",
    },

    # ── Obesity – brain regions (three separate datasets) ────
    # All three share the same metadata file (sample_info_master.csv)
    # but their TSV files live in different OvC subdirectories
    "obesity_hippocampus": {
        "tsv_dir":   NHIP / "fromBen_MacaqueObese_brain" / "DMRs" / "Hippocampus_OvC",
        "meta_src":  NHIP / "fromBen_MacaqueObese_brain" / "DMRs" / "sample_info_master.csv",
        "meta_name": "sample_info.csv",
        "pattern":   "*.roi.tsv.gz",
    },
    "obesity_hypothalamus": {
        "tsv_dir":   NHIP / "fromBen_MacaqueObese_brain" / "DMRs" / "Hypothalamus_OvC",
        "meta_src":  NHIP / "fromBen_MacaqueObese_brain" / "DMRs" / "sample_info_master.csv",
        "meta_name": "sample_info.csv",
        "pattern":   "*.roi.tsv.gz",
    },
    "obesity_prefrontalcortex": {
        "tsv_dir":   NHIP / "fromBen_MacaqueObese_brain" / "DMRs" / "PrefrontalCortex_OvC",
        "meta_src":  NHIP / "fromBen_MacaqueObese_brain" / "DMRs" / "sample_info_master.csv",
        "meta_name": "sample_info.csv",
        "pattern":   "*.roi.tsv.gz",
    },

    # ── cfDNA – gestational days (three timepoints) ──────────
    # Each gestational day is a separate OvC subdirectory
    "cfdna_GD90": {
        "tsv_dir":   NHIP / "fromBen_MacaqueObese_cfDna" / "DMRs" / "GD90_OvC",
        "meta_src":  NHIP / "fromBen_MacaqueObese_cfDna" / "DMRs" / "master_sample_info_cfDNA.csv",
        "meta_name": "sample_info.csv",
        "pattern":   "*.roi.tsv.gz",
    },
    "cfdna_GD120": {
        "tsv_dir":   NHIP / "fromBen_MacaqueObese_cfDna" / "DMRs" / "GD120_OvC",
        "meta_src":  NHIP / "fromBen_MacaqueObese_cfDna" / "DMRs" / "master_sample_info_cfDNA.csv",
        "meta_name": "sample_info.csv",
        "pattern":   "*.roi.tsv.gz",
    },
    "cfdna_GD150": {
        "tsv_dir":   NHIP / "fromBen_MacaqueObese_cfDna" / "DMRs" / "GD150_OvC",
        "meta_src":  NHIP / "fromBen_MacaqueObese_cfDna" / "DMRs" / "master_sample_info_cfDNA.csv",
        "meta_name": "sample_info.csv",
        "pattern":   "*.roi.tsv.gz",
    },
}


# ============================================================
# SECTION 3: SETUP FUNCTION
# ============================================================

def setup_dataset(name: str, config: dict):
    """
    Create data/raw/{name}/ and populate it with:
      - Symlinks to all TSV files
      - A copy of the metadata file

    Key concepts used here:
      os.symlink(src, dst)  → create a symlink at dst pointing to src
      shutil.copy2(src, dst)→ copy a file (preserving metadata)
      path.exists()         → check if a path already exists
    """
    print(f"\n── {name} ──")

    # Create the dataset folder inside data/raw/
    dest_dir = RAW_DIR / name
    dest_dir.mkdir(parents=True, exist_ok=True)
    # parents=True  → create data/ and data/raw/ if they don't exist yet
    # exist_ok=True → no error if folder already exists (safe to re-run)

    tsv_dir = config["tsv_dir"]
    pattern = config["pattern"]
    meta_src = config["meta_src"]

    # ── Symlink TSV files ────────────────────────────────────
    if not tsv_dir.exists():
        print(f"  [WARN] TSV directory not found: {tsv_dir}")
        print(f"         Skipping file linking for this dataset.")
    else:
        files = sorted(glob.glob(str(tsv_dir / pattern)))
        if not files:
            print(f"  [WARN] No files matching '{pattern}' in {tsv_dir}")
        else:
            n_linked = 0
            n_skipped = 0
            for fpath in files:
                src  = Path(fpath).resolve()
                # resolve() converts to absolute path — required for symlinks
                # Symlinks must use absolute paths, otherwise they break
                # when accessed from a different working directory

                dst  = dest_dir / src.name
                # src.name = just the filename (no directory)

                if dst.exists() or dst.is_symlink():
                    n_skipped += 1
                    continue
                    # Skip if symlink already exists (safe to re-run)

                os.symlink(src, dst)
                # Creates a symlink: data/raw/{name}/{filename} → original path
                n_linked += 1

            print(f"  TSV files: {n_linked} linked, {n_skipped} already existed")

    # ── Copy metadata file ───────────────────────────────────
    # We copy (not symlink) the metadata because:
    # 1. It's a small CSV — no disk space concern
    # 2. You may want to edit it locally without affecting the original
    meta_dst = dest_dir / config["meta_name"]
    if meta_dst.exists():
        print(f"  Metadata: already exists, skipped")
    elif not meta_src.exists():
        print(f"  [WARN] Metadata file not found: {meta_src}")
    else:
        shutil.copy2(str(meta_src), str(meta_dst))
        # copy2 preserves timestamps and permissions
        print(f"  Metadata: copied → {meta_dst.name}")


# ============================================================
# SECTION 4: MAIN
# ============================================================

if __name__ == "__main__":
    print(f"Setting up data/raw/ under: {RAW_DIR}")
    print(f"Datasets to configure: {len(DATASETS)}")

    RAW_DIR.mkdir(parents=True, exist_ok=True)

    for name, config in DATASETS.items():
        setup_dataset(name, config)

    print(f"\n{'='*50}")
    print("Done. Directory structure:")
    print(f"{'='*50}")

    # Print a summary of what was created
    for name in DATASETS:
        d = RAW_DIR / name
        if d.exists():
            n_tsv  = len(list(d.glob("*.tsv.gz")))
            has_meta = (d / "sample_info.csv").exists()
            print(f"  data/raw/{name}/")
            print(f"    {n_tsv} TSV files | metadata: {'✓' if has_meta else '✗'}")

    print(f"\nNext step: python pipeline/preprocess_wgbs.py")

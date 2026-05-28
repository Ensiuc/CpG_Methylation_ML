"""
check_hotspot_in_factors.py
──────────────────────────────────────────────────────────────────────
Check whether the Phase 1 hotspot CpGs (chr10:2,435,505–2,435,579)
appear in any MOFA+ factor in BOTH the Phase 1 and block analysis runs.

Run from the project root:
    python check_hotspot_in_factors.py
"""
import pandas as pd
import os
from pathlib import Path

PROJECT_ROOT = Path(__file__).parent
DATA_DIR     = PROJECT_ROOT / 'data' / 'processed'
TABLES_DIR   = PROJECT_ROOT / 'results' / 'tables'

# ── Hotspot definition ─────────────────────────────────────────────────
HOTSPOT_S = 2435505
HOTSPOT_E = 2435579
WINDOW    = 500   # also report CpGs within ±500 bp of hotspot

N_TOP = 10        # how many top CpGs to print per factor for context

print('=' * 65)
print('Checking Phase 1 hotspot in MOFA+ factors')
print(f'Hotspot: chr10:{HOTSPOT_S:,}–{HOTSPOT_E:,} (74 bp)')
print('=' * 65)

# ══════════════════════════════════════════════════════════════════════
# Helper: map featureN_viewM → real CpG position
# ══════════════════════════════════════════════════════════════════════
def build_cpg_list(region, subdir=''):
    """Reconstruct the sorted, coverage-filtered CpG list used by MOFA+."""
    dsets = ['wildfire','stress','obesity_hippocampus','obesity_hypothalamus',
             'obesity_prefrontalcortex','cfdna_GD45','cfdna_GD90',
             'cfdna_GD120','cfdna_GD150']
    d = DATA_DIR / subdir if subdir else DATA_DIR
    mats = {}
    for ds in dsets:
        p = d / f'{ds}_{region}_methylation.csv'
        if p.exists():
            mats[ds] = pd.read_csv(p, index_col=0)
    if not mats:
        return None
    all_cpgs = sorted(set().union(*[set(m.index) for m in mats.values()]))
    cov = pd.DataFrame({ds: (~mats[ds].isna().all(axis=1)).astype(int)
                        for ds in mats})
    frac = cov.sum(axis=1) / len(mats)
    return sorted(frac[frac >= 0.5].index)


def map_feature(name, keep_cpgs):
    try:
        idx = int(str(name).split('_')[0].replace('feature', ''))
        return int(keep_cpgs[idx])
    except Exception:
        return name   # already real position


def resolve_weights(wpath, region, subdir=''):
    """Load weights CSV, map featureN → real position if needed."""
    if not wpath.exists():
        print(f'  [MISSING] {wpath.relative_to(PROJECT_ROOT)}')
        return None
    w = pd.read_csv(wpath, index_col=0)
    if str(w.index[0]).startswith('feature'):
        keep = build_cpg_list(region, subdir)
        if keep is None:
            print(f'  [WARN] Could not build CpG list for {region} — skipping mapping')
            return w
        w.index = [map_feature(f, keep) for f in w.index]
    w.index = w.index.astype(int)
    return w


def check_hotspot(w, label):
    print(f'\n── {label} ──')
    print(f'   Weight matrix: {w.shape[0]} CpGs × {w.shape[1]} factors')

    hotspot_cpgs = [c for c in w.index if HOTSPOT_S <= c <= HOTSPOT_E]
    window_cpgs  = [c for c in w.index if HOTSPOT_S - WINDOW <= c <= HOTSPOT_E + WINDOW]

    print(f'   CpGs in hotspot ({HOTSPOT_S:,}–{HOTSPOT_E:,}): {len(hotspot_cpgs)}')
    if hotspot_cpgs:
        print(f'   Positions: {hotspot_cpgs}')
    print(f'   CpGs within ±{WINDOW} bp of hotspot: {len(window_cpgs)}')
    if window_cpgs:
        print(f'   Positions: {window_cpgs}')

    target = hotspot_cpgs if hotspot_cpgs else window_cpgs
    if not target:
        print('   → No CpGs found near the hotspot in this matrix.')
        return

    print()
    print('   Weights and rank for hotspot/nearby CpGs:')
    n_total = len(w)
    for cpg in sorted(target):
        row = []
        for factor in w.columns:
            weight = w.loc[cpg, factor]
            rank   = int(w[factor].abs().rank(ascending=False)[cpg])
            row.append(f'{factor}: w={weight:+.4f} rank={rank}/{n_total}')
        print(f'   CpG {cpg}: {" | ".join(row)}')

    print()
    print('   Top 5 CpGs per factor (for context):')
    for factor in w.columns:
        top5 = w[factor].abs().sort_values(ascending=False).head(5).index.tolist()
        print(f'   {factor}: {top5}')


# ══════════════════════════════════════════════════════════════════════
# 1. Phase 1 — genebody (hotspot lives here)
# ══════════════════════════════════════════════════════════════════════
print('\n>>> PHASE 1 ANALYSIS')
w_gen = resolve_weights(TABLES_DIR / 'mofa_weights_genebody.csv', 'genebody')
if w_gen is not None:
    check_hotspot(w_gen, 'Phase 1 — genebody')

w_cgi = resolve_weights(TABLES_DIR / 'mofa_weights_cpgi.csv', 'cpgi')
if w_cgi is not None:
    check_hotspot(w_cgi, 'Phase 1 — cpgi')

# ══════════════════════════════════════════════════════════════════════
# 2. Block analysis — fullblock
# ══════════════════════════════════════════════════════════════════════
print('\n>>> BLOCK ANALYSIS')
ba_tables = TABLES_DIR / 'block_analysis'
w_full = resolve_weights(
    ba_tables / 'mofa_weights_fullblock.csv',
    'fullblock', 'block_analysis'
)
if w_full is not None:
    check_hotspot(w_full, 'Block analysis — fullblock')

print('\n' + '=' * 65)
print('Done.')

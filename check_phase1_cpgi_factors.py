"""
check_phase1_cpgi_factors.py
──────────────────────────────────────────────────────────────────────
Print top CpGs for every factor in the Phase 1 cpgi and genebody
MOFA+ weight matrices, with their positions and weights.

Run from the project root:
    python check_phase1_cpgi_factors.py
"""
import pandas as pd
from pathlib import Path

PROJECT_ROOT = Path(__file__).parent
DATA_DIR     = PROJECT_ROOT / 'data' / 'processed'
TABLES_DIR   = PROJECT_ROOT / 'results' / 'tables'

N_TOP = 20   # top CpGs to show per factor

DATASETS = [
    'wildfire', 'stress', 'obesity_hippocampus', 'obesity_hypothalamus',
    'obesity_prefrontalcortex', 'cfdna_GD45', 'cfdna_GD90',
    'cfdna_GD120', 'cfdna_GD150'
]

# ── Helper: reconstruct sorted coverage-filtered CpG list ─────────────
def build_cpg_list(region):
    mats = {}
    for ds in DATASETS:
        p = DATA_DIR / f'{ds}_{region}_methylation.csv'
        if p.exists():
            mats[ds] = pd.read_csv(p, index_col=0)
    if not mats:
        return None
    all_cpgs = sorted(set().union(*[set(m.index) for m in mats.values()]))
    cov = pd.DataFrame({ds: (~mats[ds].isna().all(axis=1)).astype(int)
                        for ds in mats}, index=all_cpgs).fillna(0)
    frac = cov.sum(axis=1) / len(mats)
    return sorted(frac[frac >= 0.5].index)


def map_feature(name, keep_cpgs):
    try:
        idx = int(str(name).split('_')[0].replace('feature', ''))
        return int(keep_cpgs[idx])
    except Exception:
        return name


def load_weights(region):
    wpath = TABLES_DIR / f'mofa_weights_{region}.csv'
    if not wpath.exists():
        print(f'[MISSING] {wpath.name}')
        return None
    w = pd.read_csv(wpath, index_col=0)
    if str(w.index[0]).startswith('feature'):
        keep = build_cpg_list(region)
        if keep is None:
            print(f'[WARN] Could not build CpG list for {region}')
            return w
        w.index = [map_feature(f, keep) for f in w.index]
    w.index = w.index.astype(int)
    return w


def report(w, label):
    print(f'\n{"="*60}')
    print(f'{label}  ({len(w)} CpGs × {len(w.columns)} factors)')
    print(f'{"="*60}')
    for factor in w.columns:
        top = w[factor].abs().sort_values(ascending=False).head(N_TOP)
        print(f'\n  {factor} — top {N_TOP}:')
        print(f'  {"CpG position":<14} {"weight":>10}  {"rank":>12}')
        print(f'  {"-"*40}')
        n = len(w)
        for rank_i, (cpg, _) in enumerate(top.items(), 1):
            weight = w.loc[cpg, factor]
            print(f'  {cpg:<14} {weight:>+10.4f}  {rank_i:>4}/{n}')


# ── Phase 1 cpgi ──────────────────────────────────────────────────────
w_cgi = load_weights('cpgi')
if w_cgi is not None:
    report(w_cgi, 'Phase 1 — cpgi  (chr10:2,433,175–2,433,562)')

# ── Phase 1 genebody ─────────────────────────────────────────────────
w_gen = load_weights('genebody')
if w_gen is not None:
    report(w_gen, 'Phase 1 — genebody  (chr10:2,433,511–2,441,516)')

print('\nDone.')

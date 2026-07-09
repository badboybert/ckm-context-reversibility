#!/usr/bin/env python
"""downsampling_sensitivity.py — TASK B6 power-proxy / detectability robustness (review sec 4f).

STANDALONE belt-and-suspenders that the transcriptome reversibility RANKING is not a
power/detectability artifact. Does NOT touch build_determinants.py or the canonical
determinant_*.tsv, and does NOT overwrite score_reversibility outputs. Writes ONE new file:
    results/downsampling_sensitivity.tsv   (columns: metric, panel, rho, n, ...extras)

CANONICAL SCORE (the thing under test) = Z_responsiveness in
results/reversibility_score_transcriptome.tsv = | HK-RE-meta pooled blom-rank-INT(rev_beta) |
(rank-normalized, direction-aware, inverse-variance weighted; see score_reversibility.py).

Two sensitivities:
  (1) RAW-EFFECT RE-RANKING. Re-rank marks by the RAW reversal effect instead of the
      rank-INT/variance-weighted score, and Spearman-correlate the two rankings.
        - per-panel  : within each score-contributing panel, Spearman(|rev_beta|, Z_responsiveness)
        - pooled     : two raw pooled magnitudes vs Z_responsiveness
            raw_absbeta  = mean over panels of within-panel |rev_beta| percentile (direction-AGNOSTIC)
            raw_signed   = | mean over panels of within-panel SD-standardized signed rev_beta |
                           (direction-AWARE, LINEAR scale -> isolates rank-INT vs linear choice)
      High rho => the score tracks real effect size, not a detectability/normalization artifact.
      A gap between raw_absbeta and raw_signed localizes any divergence to DIRECTION-consistency
      (which the score rewards and raw magnitude does not) rather than to power.

  (2) DOWNSAMPLING to reduced effective n. Panel data are SUMMARY-level (per-mark beta+SE), so
      subject-level resampling is impossible; we simulate n' = f*n by SE-scaled Gaussian noise
      (SE ~ 1/sqrt(n)  =>  Var scales by 1/f  =>  extra_var = se^2 * (1/f - 1)), recompute the
      per-panel score z = blom(noisy_beta), re-rank by |z|, and Spearman vs the full-panel |z|
      ranking. Run on the two largest panels (most subjects = GSE141221 n=220; most marks =
      GSE199063). Mean rho over R replicates; stable rho under reduced power = robust.

Reported honestly: if raw-vs-score is only moderate, we say so and attribute it (direction).
"""
import os, numpy as np, pandas as pd
from scipy import stats

SIG = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))   # signature-pivot/
PROJ = os.path.dirname(SIG)                                         # weight-loss.omics/ (shared root)
def P(*a): return os.path.join(PROJ, *a)
def S(*a): return os.path.join(SIG, *a)

RNG = np.random.default_rng(20260708)
R_REPLICATES = 50

def blom(x):
    """Blom rank-INT, IDENTICAL to score_reversibility.py."""
    x = np.asarray(x, float); r = stats.rankdata(x); n = len(x)
    return stats.norm.ppf((r - 0.5) / n)

def spearman(a, b):
    a = np.asarray(a, float); b = np.asarray(b, float)
    m = np.isfinite(a) & np.isfinite(b)
    if m.sum() < 4: return np.nan, int(m.sum())
    rho = stats.spearmanr(a[m], b[m]).correlation
    return float(rho), int(m.sum())

# ---- canonical per-mark score --------------------------------------------------------------
SC = pd.read_csv(S('results/reversibility_score_transcriptome.tsv'), sep='\t')
Zmap = SC.set_index('mark_key')['Z_responsiveness']          # canonical score, per mark
Kmap = SC.set_index('mark_key')['k_panels']

# ---- score-contributing transcriptome panels (match the canonical build) -------------------
MAN = pd.read_csv(S('results/panel_manifest.tsv'), sep='\t')
PANELS = MAN[(MAN.layer == 'transcriptome') & (MAN.contributes_to_score)].copy()

def load_panel(row):
    """Return per-mark rev_beta / rev_se for one panel, deduped exactly like the score build."""
    df = pd.read_csv(P(row.path), sep='\t')
    out = pd.DataFrame({'mark_key': df['mark_id'].astype(str),
                        'beta': pd.to_numeric(df['rev_beta'], errors='coerce'),
                        'se':   pd.to_numeric(df['rev_se'],   errors='coerce')})
    out = out.dropna(subset=['beta']).drop_duplicates('mark_key').reset_index(drop=True)
    return out

rows_out = []

# ============================================================================================
# SENSITIVITY 1 : RAW-EFFECT RE-RANKING
# ============================================================================================
# Collect per-panel percentile / standardized contributions for the pooled raw metrics.
pool_abs = {}      # mark -> list of within-panel |beta| percentile (0..1)
pool_sgn = {}      # mark -> list of within-panel SD-standardized signed beta

for r in PANELS.itertuples(index=False):
    d = load_panel(r)
    absb = d['beta'].abs().values
    # within-panel percentile of |beta| (direction-agnostic raw magnitude, common 0..1 scale)
    pct = stats.rankdata(absb) / len(absb)
    # within-panel SD-standardized signed beta (direction-aware, LINEAR raw scale)
    sd = d['beta'].std(ddof=1)
    sgn = d['beta'].values / sd if sd and np.isfinite(sd) and sd > 0 else np.zeros(len(d))
    d['_pct'] = pct; d['_sgn'] = sgn

    # per-panel: raw |beta| vs canonical pooled score, on this panel's marks
    z_here = d['mark_key'].map(Zmap).values
    rho, n = spearman(absb, z_here)
    rows_out.append(dict(metric='raw_absbeta_vs_score', panel=r.panel, rho=rho, n=n,
                         frac=np.nan, n_replicates=np.nan,
                         note='within-panel raw |rev_beta| vs canonical Z_responsiveness'))

    for k, pv, sv in zip(d['mark_key'].values, d['_pct'].values, d['_sgn'].values):
        pool_abs.setdefault(k, []).append(pv)
        pool_sgn.setdefault(k, []).append(sv)

# pooled raw metrics, per mark, vs canonical score
marks = [m for m in pool_abs.keys() if m in Zmap.index]
raw_abs = np.array([np.mean(pool_abs[m]) for m in marks])                 # mean |beta| percentile
raw_sgn = np.array([abs(np.mean(pool_sgn[m])) for m in marks])            # |mean signed z|
zc = np.array([Zmap[m] for m in marks])
kc = np.array([Kmap.get(m, 1) for m in marks])

rho, n = spearman(raw_abs, zc)
rows_out.append(dict(metric='raw_absbeta_vs_score', panel='POOLED', rho=rho, n=n,
                     frac=np.nan, n_replicates=np.nan,
                     note='pooled mean-|rev_beta|-percentile (direction-agnostic) vs Z_responsiveness'))
rho, n = spearman(raw_sgn, zc)
rows_out.append(dict(metric='raw_signed_vs_score', panel='POOLED', rho=rho, n=n,
                     frac=np.nan, n_replicates=np.nan,
                     note='pooled |mean SD-standardized signed rev_beta| (direction-aware,linear) vs Z_responsiveness'))
# multi-panel subset (parallels the k>=2 headline classes)
mm = kc >= 2
rho, n = spearman(raw_abs[mm], zc[mm])
rows_out.append(dict(metric='raw_absbeta_vs_score', panel='POOLED_k>=2', rho=rho, n=n,
                     frac=np.nan, n_replicates=np.nan,
                     note='pooled raw |beta| percentile vs score, marks in >=2 panels'))
rho, n = spearman(raw_sgn[mm], zc[mm])
rows_out.append(dict(metric='raw_signed_vs_score', panel='POOLED_k>=2', rho=rho, n=n,
                     frac=np.nan, n_replicates=np.nan,
                     note='pooled raw signed vs score, marks in >=2 panels'))

# ============================================================================================
# SENSITIVITY 2 : DOWNSAMPLING TO REDUCED EFFECTIVE n  (SE-scaled noise injection)
# ============================================================================================
# For a single panel, the per-mark score reduces to |blom(rev_beta)| (HK-meta with k=1 returns
# the value itself), so the full-panel ranking = |blom(beta)|.  We inject SE-scaled noise to
# emulate n' = f*n and re-rank.
FRACS = [0.75, 0.50, 0.25, 0.10]
# largest by SUBJECTS (n_pairs) and largest by MARKS
DS_PANELS = [
    ('gse141221_adipose_table', 'largest by subjects (n_pairs=220)'),
    ('gse199063_adipose_table', 'largest by marks (~81k)'),
]
for pname, why in DS_PANELS:
    prow = PANELS[PANELS.panel == pname].iloc[0]
    d = load_panel(prow)
    beta = d['beta'].values
    se = d['se'].values
    n_pairs = int(prow.n_pairs)
    full_z = np.abs(blom(beta))                       # full-panel score ranking
    # se may have NaN for a few marks; median-impute for the noise draw only
    se_imp = np.where(np.isfinite(se) & (se > 0), se, np.nanmedian(se[np.isfinite(se) & (se > 0)]))
    for f in FRACS:
        extra_var = np.clip(se_imp**2 * (1.0 / f - 1.0), 0, None)
        sd_noise = np.sqrt(extra_var)
        rhos = []
        for _ in range(R_REPLICATES):
            noisy = beta + RNG.normal(0.0, sd_noise)
            z_ds = np.abs(blom(noisy))
            rho = stats.spearmanr(full_z, z_ds).correlation
            rhos.append(rho)
        rhos = np.array(rhos, float)
        n_eff = int(round(f * n_pairs))
        rows_out.append(dict(metric='downsample_effN', panel=pname,
                             rho=float(np.mean(rhos)), n=n_eff,
                             frac=f, n_replicates=R_REPLICATES,
                             note=f'{why}; full-panel |blom(beta)| vs n\'={n_eff} SE-noise re-rank; '
                                  f'rho_sd={np.std(rhos):.4f} over {len(beta)} marks'))

# ---- ANCHOR: within-panel blom magnitude vs raw |beta| (should be ~1) -----------------------
for r in PANELS.itertuples(index=False):
    d = load_panel(r)
    rho, n = spearman(d['beta'].abs().values, np.abs(blom(d['beta'].values)))
    rows_out.append(dict(metric='anchor_blom_vs_rawabs', panel=r.panel, rho=rho, n=n,
                         frac=np.nan, n_replicates=np.nan,
                         note='baseline: within-panel |blom(beta)| ranking vs raw |beta| ranking'))

# ============================================================================================
OUT = pd.DataFrame(rows_out)[['metric', 'panel', 'rho', 'n', 'frac', 'n_replicates', 'note']]
op = S('results/downsampling_sensitivity.tsv')
OUT.to_csv(op, sep='\t', index=False)

# ---- console summary -----------------------------------------------------------------------
pd.set_option('display.width', 200); pd.set_option('display.max_colwidth', 90)
def show(df, cols=('metric', 'panel', 'rho', 'n', 'frac')):
    print(df[list(cols)].to_string(index=False))

print('\n=== SENSITIVITY 1: RAW-EFFECT RE-RANKING vs canonical Z_responsiveness ===')
s1 = OUT[OUT.metric.isin(['raw_absbeta_vs_score', 'raw_signed_vs_score'])]
show(s1[s1.panel.isin(['POOLED', 'POOLED_k>=2'])])
print('\n  per-panel raw |rev_beta| vs pooled score:')
show(s1[(s1.metric == 'raw_absbeta_vs_score') & (~s1.panel.isin(['POOLED', 'POOLED_k>=2']))])

print('\n=== SENSITIVITY 2: DOWNSAMPLING (SE-scaled reduced effective n), mean rho over %d reps ===' % R_REPLICATES)
show(OUT[OUT.metric == 'downsample_effN'])

print('\n=== ANCHOR: within-panel |blom(beta)| vs raw |beta| (expect ~1.0) ===')
show(OUT[OUT.metric == 'anchor_blom_vs_rawabs'], cols=('panel', 'rho', 'n'))

print(f'\nwrote {op}  ({len(OUT)} rows)')

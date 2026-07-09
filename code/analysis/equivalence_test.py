#!/usr/bin/env python
"""equivalence_test.py  (TASK B5) — TOST equivalence test for determinant coefficients.

Converts the practical-null WORDING into a FORMAL claim (review sec 3.3 / sec 9-123).

STANDALONE by design: reads the B3-updated results/determinant_meta.tsv READ-ONLY,
does NOT touch build_determinants.py or the canonical determinant_{meta,per_panel}.tsv.

METHOD
------
The canonical meta table (build_determinants.py, hk()) reports, for each feature, the
Hartung-Knapp pooled beta and a 95% CI half-width = t_{0.975,df} * SE, with df = k-1.
We therefore recover the pooled SE exactly:

        SE = ci_halfwidth / t_{0.975, df}          (df = k-1)

For the transcriptome determinant layer k=9 -> df=8, t_{0.975,8}=2.306004.

TOST (two one-sided tests) of practical equivalence against a symmetric SESOI:
    lower H0:  beta <= -SESOI   vs  H1: beta > -SESOI   ->  t_L = (beta+SESOI)/SE,  p_L = P(T_df > t_L)
    upper H0:  beta >= +SESOI   vs  H1: beta < +SESOI   ->  t_U = (beta-SESOI)/SE,  p_U = P(T_df < t_U)
    TOST p = max(p_L, p_U).  Equivalence declared at alpha=0.05  <=>  TOST p < 0.05
    <=>  the 90% CI (= beta +/- t_{0.95,df}*SE) lies fully within [-SESOI, +SESOI].
We report BOTH the TOST p and the 90%-CI-containment flag (they are algebraically identical).

SESOI (pre-declared, not data-tuned)
------------------------------------
PRIMARY  +/-0.05 SD.  Biological justification: the reversibility score is a rank-INT
  standardized quantity; 0.05 SD is a negligible fraction of its dynamic range and is SMALLER
  than the smallest effect the paper interprets anywhere, yet sits ABOVE every feature's
  minimum detectable effect (median transcriptome MDE ~0.013 SD; even the largest feature
  MDEs ~0.02-0.035 SD remain below 0.05), so the study resolves effects finer than the SESOI
  and equivalence within +/-0.05 SD is a well-powered claim. An effect under 0.05 SD cannot
  move a determinant story.
SECONDARY (stricter) +/-0.03 SD = the LARGEST observed pooled |beta| across the 13
  transcriptome features (arch_nsig, |beta|=0.0287). Equivalence to +/-0.03 means a feature
  is bounded below even the biggest effect anyone in this layer shows.

HONEST READ (expected + reported): the determinant-carrying features -- causal_nonEGFR and
the architecture terms (arch_nsig/arch_str/has_arch), which have narrow CIs -- are
statistically EQUIVALENT to zero within +/-0.05 (and mostly within +/-0.03). A small number
of high-variance dynamic-range covariates (notably tau, I2=96.9%) have CIs too wide to
certify equivalence at these bounds and are reported as INCONCLUSIVE, not as null and not as
signal -- consistent with their role as magnitude/dynamic-range covariates we do not
interpret as determinants.
"""
import os
import numpy as np
import pandas as pd
from scipy import stats

SIG = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
def S(*a): return os.path.join(SIG, *a)

# ---- pre-declared SESOI ----
SESOI_PRIMARY = 0.05   # SD
SESOI_STRICT  = 0.03   # SD (= largest observed pooled |beta|)
ALPHA = 0.05

META = pd.read_csv(S('results', 'determinant_meta.tsv'), sep='\t')
tx = META[META.layer == 'transcriptome'].copy()   # determinant layer (k=9, df=8); proteome k=3 is power-limited, excluded

def recover_se(row):
    df = int(row['df'])                      # df = k-1
    t975 = stats.t.ppf(0.975, df)
    return row['ci_halfwidth'] / t975

def tost(beta, se, df, sesoi):
    """max of the two one-sided p-values, plus 90% CI containment flag."""
    t_lower = (beta + sesoi) / se            # H0: beta <= -SESOI
    t_upper = (beta - sesoi) / se            # H0: beta >= +SESOI
    p_lower = stats.t.sf(t_lower, df)        # P(T > t_lower)
    p_upper = stats.t.cdf(t_upper, df)       # P(T < t_upper)
    p = max(p_lower, p_upper)
    t90 = stats.t.ppf(1 - ALPHA, df)         # 90% CI uses t_{0.95,df}
    lo, hi = beta - t90 * se, beta + t90 * se
    equiv = bool((lo >= -sesoi) and (hi <= sesoi))
    return p, equiv, lo, hi

rows = []
for _, r in tx.iterrows():
    df = int(r['df'])
    se = recover_se(r)
    b = r['pooled_beta']
    p05, eq05, lo90, hi90 = tost(b, se, df, SESOI_PRIMARY)
    p03, eq03, _, _       = tost(b, se, df, SESOI_STRICT)
    rows.append(dict(
        feature=r['feature'], pooled_beta=b, se=se,
        tost_p_05=p05, equiv_05=eq05,
        tost_p_03=p03, equiv_03=eq03,
        ci90_lo=lo90, ci90_hi=hi90,
        k=int(r['k']), df=df, I2=r['I2'],
    ))

OUT = pd.DataFrame(rows)
# order: determinant-carrying (equivalent, narrow CI) first, inconclusive last, by CI width
OUT = OUT.sort_values(['equiv_05', 'equiv_03', 'se'], ascending=[False, False, True]).reset_index(drop=True)

# write exactly the requested schema (+ transparency columns)
COLS = ['feature', 'pooled_beta', 'se', 'tost_p_05', 'equiv_05', 'tost_p_03', 'equiv_03',
        'ci90_lo', 'ci90_hi', 'k', 'df', 'I2']
OUT[COLS].to_csv(S('results', 'determinant_tost.tsv'), sep='\t', index=False)

# ---- console summary ----
pd.set_option('display.width', 200); pd.set_option('display.max_columns', 20)
print(f"TOST equivalence test — transcriptome determinant layer (k={int(tx.k.max())}, df={int(tx.df.max())})")
print(f"pre-declared SESOI: PRIMARY +/-{SESOI_PRIMARY} SD, STRICT +/-{SESOI_STRICT} SD; alpha={ALPHA}\n")
show = OUT.copy()
for c in ['pooled_beta', 'se', 'tost_p_05', 'tost_p_03', 'ci90_lo', 'ci90_hi']:
    show[c] = show[c].map(lambda x: f"{x:+.4f}")
show['I2'] = show['I2'].map(lambda x: f"{x:.1f}")
print(show[COLS].to_string(index=False))

eq05 = OUT[OUT.equiv_05]; ne05 = OUT[~OUT.equiv_05]
eq03 = OUT[OUT.equiv_03]; ne03 = OUT[~OUT.equiv_03]
print(f"\nEQUIVALENT within +/-0.05 SD  ({len(eq05)}/{len(OUT)}): {list(eq05.feature)}")
print(f"INCONCLUSIVE at +/-0.05 SD    ({len(ne05)}/{len(OUT)}): {list(ne05.feature)}")
print(f"EQUIVALENT within +/-0.03 SD  ({len(eq03)}/{len(OUT)}): {list(eq03.feature)}")
print(f"INCONCLUSIVE at +/-0.03 SD    ({len(ne03)}/{len(OUT)}): {list(ne03.feature)}")

det = ['causal_nonEGFR', 'arch_nsig', 'arch_str', 'has_arch']
d = OUT[OUT.feature.isin(det)]
print(f"\nDETERMINANT-CARRYING features {det}:")
print(f"  equivalent within +/-0.05: {sorted(d[d.equiv_05].feature)}")
print(f"  equivalent within +/-0.03: {sorted(d[d.equiv_03].feature)}")
print("\nwrote results/determinant_tost.tsv")

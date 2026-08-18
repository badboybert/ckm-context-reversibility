#!/usr/bin/env python
"""reconcile_drug_class.py — join the two drug-class metrics per intervention pair.

Companion to build_drug_class.py (effect-size Spearman) and drug_class_overlap.py (sig-set overlap).
Reconciles the two Fig 5 metrics: the effect-size rho and the
power-normalised sig-set overlap fold are DIFFERENT tests and must not be conflated.

Metric 1 (effect-size correlation): results/drug_class_corr.tsv
   rho = Spearman of rev_beta over shared UniProt (min overlap 30). Tracks detectability/power.
Metric 2 (sig-set overlap, power-normalised): results/drug_class_overlap.tsv
   fold = obs/exp of jointly-significant proteins (hypergeometric); samedir = same-direction fraction.

A pair is "supported shared-signature" ONLY where BOTH metrics are significant + concordant-positive.

Output:
   results/drug_class_reconciliation.tsv  (wide 28-pair join + derived `agree` flag)
   results/drug_class_fig5_panel.tsv      (tidy long, powered-on-both pairs, for side-by-side Fig 5)
"""
import os, pandas as pd, numpy as np
SIG = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
def S(*a): return os.path.join(SIG, *a)

corr = pd.read_csv(S('results/drug_class_corr.tsv'), sep='\t').rename(columns={'p': 'rho_p'})
ovl = pd.read_csv(S('results/drug_class_overlap.tsv'), sep='\t').rename(
    columns={'U': 'overlap_U', 'fold': 'overlap_fold', 'p': 'overlap_p', 'samedir': 'samedir_frac',
             'na': 'sig_a', 'nb': 'sig_b', 'obs': 'overlap_obs', 'exp': 'overlap_exp'})
# both engines emit pairs in identical itertools.combinations order -> (a,b) keys align exactly
J = corr.merge(ovl[['a', 'b', 'overlap_U', 'sig_a', 'sig_b', 'overlap_obs', 'overlap_exp',
                    'overlap_fold', 'overlap_p', 'samedir_frac']], on=['a', 'b'], how='outer')

MECH = {'metformin': 'biguanide', 'dapagliflozin': 'SGLT2i', 'empagliflozin': 'SGLT2i',
        'liraglutide': 'GLP1RA', 'semaglutide': 'GLP1RA', 'BBS': 'bariatric',
        'MS_bariatric': 'bariatric', 'DiRECT': 'diet'}
CLS = {'metformin': 'DRUG', 'dapagliflozin': 'DRUG', 'empagliflozin': 'DRUG', 'liraglutide': 'DRUG',
       'semaglutide': 'DRUG', 'BBS': 'SURGERY', 'MS_bariatric': 'SURGERY', 'DiRECT': 'DIET'}
J['class_a'] = J['a'].map(CLS); J['class_b'] = J['b'].map(CLS)
J['mech_a'] = J['a'].map(MECH); J['mech_b'] = J['b'].map(MECH)
SIGONLY = {'dapagliflozin', 'liraglutide'}   # 38/38, 151/151 pre-filtered sig-only -> degenerate overlap universe

def rho_state(r):
    if pd.isna(r.rho): return 'underpowered'            # n<30
    if r.rho_p < 0.05 and r.rho > 0: return 'pos'
    if r.rho_p < 0.05 and r.rho < 0: return 'neg'
    return 'null'

def ovl_state(r):
    if pd.isna(r.overlap_fold): return 'underpowered'   # U<50 skipped
    if (r.a in SIGONLY and r.sig_a == r.overlap_U) or (r.b in SIGONLY and r.sig_b == r.overlap_U):
        return 'degenerate'                             # one panel entirely-sig within U -> fold forced ~1
    if r.overlap_p < 0.05 and r.overlap_fold > 1: return 'enriched'
    if r.overlap_p < 0.05 and r.overlap_fold < 1: return 'depleted'
    return 'null'

def agree_flag(rs, os_):
    if rs == 'underpowered' or os_ in ('underpowered', 'degenerate'): return 'underpowered'
    if rs == 'pos' and os_ == 'enriched': return 'concordant-positive'
    if rs == 'null' and os_ == 'null': return 'concordant-null'
    if rs == 'neg' and os_ == 'depleted': return 'concordant-negative'
    return 'discordant'

J['rho_state'] = J.apply(rho_state, axis=1)
J['overlap_state'] = J.apply(ovl_state, axis=1)
J['agree'] = J.apply(lambda r: agree_flag(r.rho_state, r.overlap_state), axis=1)
J['supported_shared'] = (J['agree'] == 'concordant-positive')
J['powered_both'] = (J['rho_state'] != 'underpowered') & (~J['overlap_state'].isin(['underpowered', 'degenerate']))
J['class_pair'] = J.apply(lambda r: 'within-' + r.class_a if r.class_a == r.class_b
                          else ' x '.join(sorted([r.class_a, r.class_b])), axis=1)
for c in ('rho', 'overlap_fold', 'samedir_frac', 'overlap_exp'):
    J[c] = J[c].astype(float).round(3)

cols = ['a', 'b', 'n', 'rho', 'rho_p', 'overlap_U', 'sig_a', 'sig_b', 'overlap_obs', 'overlap_exp',
        'overlap_fold', 'overlap_p', 'samedir_frac', 'class_a', 'class_b', 'mech_a', 'mech_b',
        'rho_state', 'overlap_state', 'agree', 'supported_shared', 'powered_both', 'class_pair']
J[cols].rename(columns={'a': 'pair_a', 'b': 'pair_b'}).to_csv(
    S('results/drug_class_reconciliation.tsv'), sep='\t', index=False)

# tidy long file for a side-by-side Fig 5 panel: powered-on-both pairs only
rows = []
for _, r in J[J.powered_both].iterrows():
    lab = f'{r.a} x {r.b}'
    rows.append(dict(pair_label=lab, class_pair=r.class_pair, metric='effect_size_rho',
                     value=r.rho, p=r.rho_p, significant=bool(r.rho_p < 0.05),
                     n_or_U=int(r.n), samedir_frac=r.samedir_frac, agree=r.agree))
    rows.append(dict(pair_label=lab, class_pair=r.class_pair, metric='overlap_fold',
                     value=r.overlap_fold, p=r.overlap_p,
                     significant=bool(r.overlap_p < 0.05 and r.overlap_fold > 1),
                     n_or_U=int(r.overlap_U), samedir_frac=r.samedir_frac, agree=r.agree))
pd.DataFrame(rows).to_csv(S('results/drug_class_fig5_panel.tsv'), sep='\t', index=False)
print('wrote results/drug_class_reconciliation.tsv + results/drug_class_fig5_panel.tsv')

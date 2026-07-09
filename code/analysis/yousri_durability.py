#!/usr/bin/env python
"""yousri_durability.py — 12-YEAR PROTEOME durability + regain (Yousri 2022, PMID 34796696, Obesity).
Provenance + dedup for the acquired Yousri RYGB SomaLogic tables (acquisition wf_50f0e0b1).
Dedups multi-aptamer reagents to UniProt (min-p; APOE x4 etc.), computes the 2yr->12yr durability
concordance (the proteome complement to the GSE199063 5yr transcriptome durability), and reports the
regain set. rev_beta = RAW post-minus-pre RFU (UP if positive); durability is convention-agnostic
(both arms same sign), so no toward-lean re-sign is needed for the 2yr<->12yr Spearman.
Output: results/yousri_durability.tsv
"""
import os, numpy as np, pandas as pd
from scipy.stats import spearmanr
SIG=os.path.dirname(os.path.dirname(os.path.abspath(__file__))); PROJ=os.path.dirname(SIG)
def P(*a): return os.path.join(PROJ,*a)
def S(*a): return os.path.join(SIG,*a)
def load(f):
    d=pd.read_csv(P('data/proteome',f),sep='\t')
    d['UniProt']=d['UniProt'].astype(str).str.split(r'[;,_]').str[0]
    return d.sort_values('rev_p',na_position='last').drop_duplicates('UniProt').set_index('UniProt')  # min-p per UniProt
y2=load('yousri_rygb_2yr_reversal.tsv'); y12=load('yousri_rygb_12yr.tsv')
ov=y2.index.intersection(y12.index)
rho,p=spearmanr(y2.loc[ov,'rev_beta'],y12.loc[ov,'rev_beta'])
qcol='rev_q' if 'rev_q' in y2 else 'rev_p'; thr=0.10 if qcol=='rev_q' else 0.05
sig2=y2.index[pd.to_numeric(y2[qcol],errors='coerce')<thr].intersection(y12.index)
same=np.mean(np.sign(y2.loc[sig2,'rev_beta'])==np.sign(y12.loc[sig2,'rev_beta']))
# magnitude retention: |12yr|/|2yr| among 2yr-sig (capped) = how much of the effect size remains
ret=np.median(np.clip((y12.loc[sig2,'rev_beta'].abs()/y2.loc[sig2,'rev_beta'].abs().replace(0,np.nan)),0,2).dropna())
rows=[dict(metric='durability_spearman_2yr_12yr',value=round(rho,4)),
      dict(metric='durability_spearman_p',value=p),
      dict(metric='n_overlap_proteins',value=len(ov)),
      dict(metric='n_2yr_sig',value=len(sig2)),
      dict(metric='frac_same_direction_at_12yr',value=round(float(same),4)),
      dict(metric='median_magnitude_retained_12yr',value=round(float(ret),4))]
pd.DataFrame(rows).to_csv(S('results/yousri_durability.tsv'),sep='\t',index=False)
print(f'Yousri RYGB PROTEOME durability (n={len(ov)} proteins, dedup-to-UniProt):')
print(f'  Spearman(2yr,12yr rev_beta) = {rho:.3f} (p={p:.1e})')
print(f'  among {len(sig2)} 2yr-significant proteins: {same*100:.0f}% same direction at 12yr; median |12yr|/|2yr| = {ret:.2f}')
gc=[c for c in y2.columns if c.lower()=="gene"][0]
print('  canonical (2yr -> 12yr RFU change): ', end='')
print(', '.join(f'{g} {y2[y2[gc]==g].rev_beta.iloc[0]:+.0f}->{y12[y12[gc]==g].rev_beta.iloc[0]:+.0f}'
      for g in ['LEP','IGFBP1','IGFBP2'] if len(y2[y2[gc]==g]) and len(y12[y12[gc]==g])))
print('wrote results/yousri_durability.tsv')

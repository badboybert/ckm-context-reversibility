#!/usr/bin/env python
"""prep_ED3_data.py — Extended Data 3: the determinant null is POWERED, not under-powered.
(a) per-panel coefficient forest (9 transcriptome panels) for the anchor determinants + pooled meta + MDE band.
(b) tau dynamic-range diagnostic (per-panel sign-flip).  (c) causal-status per-panel (Bonferroni-NS in every panel).
Sources: determinant_per_panel.tsv + determinant_meta.tsv."""
import os, numpy as np, pandas as pd
HERE=os.path.dirname(os.path.abspath(__file__)); SIG=os.path.dirname(os.path.dirname(HERE))
RES=os.path.join(SIG,'results'); OUT=os.path.join(HERE,'source_data'); os.makedirs(OUT,exist_ok=True)
R=lambda p: pd.read_csv(os.path.join(RES,p),sep='\t')
pp=R('determinant_per_panel.tsv'); meta=R('determinant_meta.tsv')
tx=pp[pp.layer=='transcriptome'].copy()
NAME={'loeuf':'LOEUF','n_drug_log':'druggability','n_gwas_log':'GWAS burden','tau':'tissue-spec. (τ)',
      'is_secreted':'secreted','is_enzyme':'enzyme','is_membrane':'membrane','has_arch':'has cis-QTL','causal_nonEGFR':'causal status'}
ANCHOR=list(NAME)

# a) per-panel betas (long) + pooled meta + MDE for the anchor features
a=tx[tx.feature.isin(ANCHOR)][['panel','feature','beta','p_bonf']].copy()
a['feat']=a.feature.map(NAME)
a.to_csv(os.path.join(OUT,'ED3a_perpanel.csv'),index=False)
mt=meta[(meta.layer=='transcriptome')&(meta.feature.isin(ANCHOR))][['feature','pooled_beta','ci_halfwidth','mde']].copy()
mt['feat']=mt.feature.map(NAME)
mt.to_csv(os.path.join(OUT,'ED3a_meta.csv'),index=False)

# b) tau per-panel (sign-flip diagnostic)
tau=tx[tx.feature=='tau'][['panel','beta','p','p_bonf']].copy()
tau['sig_bonf']=tau.p_bonf<0.05
tau.to_csv(os.path.join(OUT,'ED3b_tau.csv'),index=False)

# c) causal-status per-panel (NS everywhere)
ca=tx[tx.feature=='causal_nonEGFR'][['panel','beta','se','p','p_bonf']].copy()
ca['sig_bonf']=ca.p_bonf<0.05
ca.to_csv(os.path.join(OUT,'ED3c_causal.csv'),index=False)

mde_med=meta[meta.layer=='transcriptome'].mde.median()
print('ED3a: %d anchor features × %d panels; max|pooled|=%.4f; median MDE=%.4f'%(
    len(ANCHOR), tx.panel.nunique(), mt.pooled_beta.abs().max(), mde_med))
print('ED3b tau: %d pos / %d neg; Bonferroni-sig %d/%d panels'%((tau.beta>0).sum(),(tau.beta<0).sum(),tau.sig_bonf.sum(),len(tau)))
print('ED3c causal: pooled-ish mean=%.4f; Bonferroni-sig %d/%d (should be 0)'%(ca.beta.mean(),ca.sig_bonf.sum(),len(ca)))
print('per-panel Bonferroni-sig counts by feature:', tx[tx.p_bonf<0.05].feature.value_counts().to_dict())

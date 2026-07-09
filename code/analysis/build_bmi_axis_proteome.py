#!/usr/bin/env python
"""Proteome BMI-restoration axis from Goudswaard 2023 ST5 (one-sample MR: effect of BMI
on plasma protein, beta = mean protein-level diff per SD higher BMI). This is the external
lean axis for the restoration score (SCORE_SPEC_FROZEN.md s1-2, PI decision #1).
A protein's restoration = intervention moves it OPPOSITE to its BMI (higher-adiposity) sign.
Aggregates the multiple SomaLogic aptamers per UniProt to one axis (median beta, min p)."""
import os, numpy as np, pandas as pd
SIG=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))   # signature-pivot/
PROJ=os.path.dirname(SIG)                                          # weight-loss.omics/ (shared root)
def P(*a): return os.path.join(PROJ,*a)   # SHARED inputs (Goudswaard suppl, reversal tables)
def S(*a): return os.path.join(SIG,*a)    # signature-pivot outputs

df=pd.read_excel(P('data/proteome/goudswaard2023_supp.xlsx'),sheet_name='ST5',header=2)
df=df.rename(columns={'UniProtID':'UniProt','Beta coefficient':'bmi_beta','P_val':'bmi_p','SE':'bmi_se'})
df=df[['UniProt','bmi_beta','bmi_p','bmi_se']].copy()
for c in ['bmi_beta','bmi_p','bmi_se']: df[c]=pd.to_numeric(df[c],errors='coerce')
df=df.dropna(subset=['UniProt','bmi_beta'])

# aggregate aptamers -> one axis per UniProt; require sign-consistency across aptamers (else flag)
def agg(g):
    return pd.Series(dict(bmi_beta=g['bmi_beta'].median(), bmi_p=g['bmi_p'].min(),
        n_aptamers=len(g), sign_consistent=bool((g['bmi_beta']>0).all() or (g['bmi_beta']<0).all())))
ax=df.groupby('UniProt').apply(agg).reset_index()
ax['bmi_sign']=np.sign(ax['bmi_beta']).astype(int)
op=S('results/bmi_axis_proteome.tsv'); ax.to_csv(op,sep='\t',index=False)

print(f'proteins with BMI axis: {len(ax)}  (sign-consistent across aptamers: {int(ax.sign_consistent.sum())}, '
      f'{ax.sign_consistent.mean()*100:.0f}%)')
print('controls (BMI should RAISE these):')
for u,name in [('P41159','LEP/leptin'),('P15090','FABP4'),('P05231','IL6')]:
    r=ax[ax.UniProt==u]
    if len(r): print(f'  {name:14s} {u}: bmi_beta={r.bmi_beta.iloc[0]:+.3f} p={r.bmi_p.iloc[0]:.1e} n_apt={int(r.n_aptamers.iloc[0])}')
# overlap with the proteome reversal panels
for f in ['reversal_olink_bbs','reversal_somascan_direct','reversal_emperor']:
    rv=pd.read_csv(P(f'data/proteome/{f}.tsv'),sep='\t')
    ov=rv['UniProt'].isin(set(ax.UniProt)).sum()
    print(f'  {f:26s}: {ov}/{len(rv)} proteins have a BMI axis')
print(f'wrote {op}')

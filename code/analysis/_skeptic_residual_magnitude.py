#!/usr/bin/env python
"""Quantify the MAGNITUDE of the residual within-decile gradient and finer-bin sensitivity.
The logit p<<0 on methylation is N-driven; question is whether the residual is MATERIAL."""
import os, numpy as np, pandas as pd
from scipy import stats
SIG=os.path.dirname(os.path.dirname(os.path.abspath(__file__))); PROJ=os.path.dirname(SIG)
def P(*a): return os.path.join(PROJ,*a)
def S(*a): return os.path.join(SIG,*a)
def nm(s): return str(s).strip().rstrip('*').strip().lower()
def pick(c,cands):
    for x in cands:
        if x in c: return x
MAN=pd.read_csv(S('results/panel_manifest.tsv'),sep='\t'); ELIG=MAN[MAN.persistence_eligible]
EFF=['rev_beta','effect','delta_beta']
bsd_tx=pd.read_csv(S('results/base_sd_transcriptome.tsv'),sep='\t')
bsd_mb=pd.read_csv(S('results/base_sd_metabolome.tsv'),sep='\t')
def base_sd_for(row,d):
    if row.layer=='transcriptome':
        m=bsd_tx[bsd_tx.panel==row.panel].set_index('mark_id')['base_sd']; return d['mark_key'].map(m)
    if row.layer=='metabolome':
        m=bsd_mb.set_index('mark_key')['base_sd']; return d['mark_key'].map(m)
def load(row):
    df=pd.read_csv(P(row.path),sep='\t'); cols=list(df.columns); ec=pick(cols,EFF)
    key=(df[pick(cols,['mark_id','metabolite'])].map(nm) if row.layer=='metabolome'
         else df[{'transcriptome':'mark_id','methylation':'cpg'}[row.layer]].astype(str))
    eff=pd.to_numeric(df[ec],errors='coerce'); d=pd.DataFrame({'mark_key':key,'eff':eff})
    if row.layer=='methylation': d['base_sd']=pd.to_numeric(df['base_sd'],errors='coerce')
    d=d.dropna(subset=['eff']).drop_duplicates('mark_key')
    if row.layer!='methylation': d['base_sd']=base_sd_for(row,d).values
    d=d.dropna(subset=['base_sd']); d['za']=np.abs(d['eff'].values); return d

def classify(d,nbin):
    try: d['dec']=pd.qcut(d['base_sd'],nbin,labels=False,duplicates='drop')
    except ValueError: d['dec']=0
    pers=pd.Series(False,index=d.index); rev=pd.Series(False,index=d.index)
    for _,g in d.groupby('dec'):
        lo=g['za'].quantile(1/3); hi=g['za'].quantile(2/3)
        pers.loc[g.index[g['za']<=lo]]=True; rev.loc[g.index[g['za']>=hi]]=True
    return pers,rev

for row in ELIG.itertuples(index=False):
    d=load(row)
    print(f"\n[{row.panel}] n={len(d)} base_sd range {d.base_sd.min():.3g}-{d.base_sd.max():.3g} (span x{d.base_sd.max()/max(d.base_sd.min(),1e-12):.0f})")
    for nbin in [10,20,50,100]:
        if len(d)<nbin*6: continue
        pers,rev=classify(d.copy(),nbin)
        bp=d.loc[pers,'base_sd'].median(); br=d.loc[rev,'base_sd'].median()
        # AUC-style: fraction of (persistent,reversible) pairs where persistent has LOWER base_sd
        # use rank-biserial via Mann-Whitney
        mw=stats.mannwhitneyu(d.loc[pers,'base_sd'],d.loc[rev,'base_sd'],alternative='two-sided')
        n1=pers.sum(); n2=rev.sum(); auc=mw.statistic/(n1*n2)  # P(persistent>reversible)
        print(f"   {nbin:3d} bins: ratio={bp/br:.3f}  P(base_sd_pers>base_sd_rev)={auc:.3f} (0.5=clean) p={mw.pvalue:.1e}")

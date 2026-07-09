#!/usr/bin/env python
"""prep_ED4_data.py — Extended Data 4: persistence is a variance-matched measurability/noise floor, NOT memory.
(a) base_sd vs Z_responsiveness by class (transcriptome) — persistent = low-variance, low-responsiveness floor.
(b) variance-conditional class counts (transcriptome + methylation).
(c) THE FALSIFICATION: raw cross-tissue methylation persistent overlap COLLAPSES to chance on base-mean conditioning,
    while the reversible pole holds. Sources: marks_classified_* + cross_tissue_methylation_persistence{,_edgecond}.tsv."""
import os, numpy as np, pandas as pd
HERE=os.path.dirname(os.path.abspath(__file__)); SIG=os.path.dirname(os.path.dirname(HERE))
RES=os.path.join(SIG,'results'); OUT=os.path.join(HERE,'source_data'); os.makedirs(OUT,exist_ok=True)
R=lambda p: pd.read_csv(os.path.join(RES,p),sep='\t')

# a) base_sd vs Z_responsiveness by class (transcriptome), persistent + reversible poles
tx=R('marks_classified_transcriptome.tsv')
pa=tx[tx['class'].isin(['persistent','consistently-reversible'])][['mark_key','class','base_sd_mean','Z_responsiveness']].copy()
pa['class']=pa['class'].map({'persistent':'persistent','consistently-reversible':'reversible'})
pa.to_csv(os.path.join(OUT,'ED4a_class_scatter.csv'),index=False)
med=pa.groupby('class').median(numeric_only=True).round(3)
ratio=med.loc['persistent','base_sd_mean']/med.loc['reversible','base_sd_mean']

# b) class counts (transcriptome + methylation)
rows=[]
for lay,f in [('transcriptome','marks_classified_transcriptome.tsv'),('methylation','marks_classified_methylation.tsv')]:
    vc=R(f)['class'].value_counts()
    for cls in ['persistent','consistently-reversible','intermediate']:
        rows.append({'layer':lay,'class':{'consistently-reversible':'reversible'}.get(cls,cls),'n':int(vc.get(cls,0))})
pd.DataFrame(rows).to_csv(os.path.join(OUT,'ED4b_class_counts.csv'),index=False)

# c) falsification: raw vs edge-conditioned cross-tissue methylation overlap, persistent vs reversible
raw=R('cross_tissue_methylation_persistence.tsv'); ec=R('cross_tissue_methylation_persistence_edgecond.tsv')
m=raw.merge(ec,on='pair',suffixes=('','_ec'))
crows=[]
for _,x in m.iterrows():
    pr=x['pair'].replace('x',' × ')
    crows.append({'pair':pr,'pole':'persistent','cond':'raw','fold':x['pers_top5_fold']})
    crows.append({'pair':pr,'pole':'persistent','cond':'edge-conditioned','fold':x['pers_ec_fold']})
    crows.append({'pair':pr,'pole':'reversible','cond':'raw','fold':x['rev_fold']})
    crows.append({'pair':pr,'pole':'reversible','cond':'edge-conditioned','fold':x['rev_ec_fold']})
pd.DataFrame(crows).to_csv(os.path.join(OUT,'ED4c_falsification.csv'),index=False)
print('ED4a: persistent base_sd med=%.3f vs reversible %.3f (ratio %.3f); Z med persist=%.3f rev=%.3f'%(
    med.loc['persistent','base_sd_mean'],med.loc['reversible','base_sd_mean'],ratio,
    med.loc['persistent','Z_responsiveness'],med.loc['reversible','Z_responsiveness']))
print('ED4c persistent raw->ec folds:', [(r['pair'],round(r['fold'],2)) for r in crows if r['pole']=='persistent'])

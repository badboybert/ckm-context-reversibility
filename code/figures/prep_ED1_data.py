#!/usr/bin/env python
"""prep_ED1_data.py — Extended Data 1: reversibility-score validity.
(a) per-layer Z_responsiveness distributions (folded magnitude); (b) biology face-validity by SIGNED Z
(SAA1/SAA2/CRP/LEP down, GDF15 up — recovers known weight-loss biology); (c) score-vs-measurability (corr~0,
a NEW axis vs F1c's log-N). Read Z by column NAME (methylation has a different column order)."""
import os, numpy as np, pandas as pd
from scipy.stats import spearmanr
HERE=os.path.dirname(os.path.abspath(__file__)); SIG=os.path.dirname(os.path.dirname(HERE))
RES=os.path.join(SIG,'results'); OUT=os.path.join(HERE,'source_data'); os.makedirs(OUT,exist_ok=True)
R=lambda p: pd.read_csv(os.path.join(RES,p),sep='\t')

# a) per-layer Z distributions (read Z_responsiveness by NAME)
parts=[]
for lay,f in [('proteome','reversibility_score_proteome.tsv'),('transcriptome','reversibility_score_transcriptome.tsv'),
              ('metabolome','reversibility_score_metabolome.tsv'),('methylation','reversibility_score_methylation.tsv')]:
    d=R(f); z=pd.to_numeric(d['Z_responsiveness'],errors='coerce').dropna()
    parts.append(pd.DataFrame({'layer':lay,'Z':z.values}))
dist=pd.concat(parts)
# subsample methylation (866k) for plotting; keep all others
dist=pd.concat([dist[dist.layer!='methylation'], dist[dist.layer=='methylation'].sample(8000,random_state=1)])
dist.to_csv(os.path.join(OUT,'ED1a_z_dist.csv'),index=False)
kcov={lay:int(R(f)['k_panels'].median()) for lay,f in [('proteome','reversibility_score_proteome.tsv'),
      ('transcriptome','reversibility_score_transcriptome.tsv'),('metabolome','reversibility_score_metabolome.tsv'),
      ('methylation','reversibility_score_methylation.tsv')]}

# b) biology face-validity (transcriptome signed Z; gene symbols)
tx=R('reversibility_score_transcriptome.tsv')
bio=['SAA1','SAA2','CRP','LEP','IGFBP1','GDF15']
b=tx[tx.mark_key.isin(bio)][['mark_key','pooled_signed','Z_responsiveness','sign_concord']].copy()
b.columns=['gene','signed_Z','Z','sign_concord']
b.to_csv(os.path.join(OUT,'ED1b_biology.csv'),index=False)

# c) score vs measurability (transcriptome): join Z with meas_mean from marks_classified
mc=R('marks_classified_transcriptome.tsv')[['mark_key','meas_mean']]
j=tx[['mark_key','Z_responsiveness']].merge(mc,on='mark_key',how='inner').dropna()
rho=spearmanr(j.Z_responsiveness,j.meas_mean).correlation
j.sample(min(4000,len(j)),random_state=2).to_csv(os.path.join(OUT,'ED1c_measurability.csv'),index=False)
pd.DataFrame([{'corr':round(rho,3),'n':len(j)}]).to_csv(os.path.join(OUT,'ED1c_corr.csv'),index=False)
print('k coverage (median):',kcov)
print('biology signed-Z:', dict(zip(b.gene,b.signed_Z.round(2))))
print('score-vs-measurability (txn): corr=%.3f n=%d'%(rho,len(j)))

#!/usr/bin/env python
"""prep_F1ef_data.py — source data for the two NEW Figure 1 panels.
F1e = within-cohort split-half reproducibility (internal stability; distinct from F1d cross-cohort transfer).
F1f = cross-ancestry generalizability + its limits (directional-only, underpowered n=12)."""
import os, pandas as pd
HERE=os.path.dirname(os.path.abspath(__file__)); SIG=os.path.dirname(os.path.dirname(HERE))
RES=os.path.join(SIG,'results'); OUT=os.path.join(HERE,'source_data'); os.makedirs(OUT,exist_ok=True)
def metr(p):
    d=pd.read_csv(os.path.join(RES,p),sep='\t'); return dict(zip(d.metric,d.value)), dict(zip(d.metric,d.get('sd',pd.Series())))

# e) split-half: signed Spearman + reversible/persistent recurrence, per layer
ta,tasd=metr('orphan/split_half_adipose.tsv'); ma,masd=metr('orphan/split_half_methylation.tsv')
pd.DataFrame([
 {'layer':'transcriptome','sp_signed':round(ta['sp_signed'],3),'sd':round(tasd.get('sp_signed',0),3),'rev_enr':round(ta['rev_enr'],2),'persist_enr':round(ta['persist_enr'],2)},
 {'layer':'methylome','sp_signed':round(ma['sp_signed'],3),'sd':round(masd.get('sp_signed',0),3),'rev_enr':round(ma['rev_enr'],2),'persist_enr':round(ma['persist_enr'],2)},
]).to_csv(os.path.join(OUT,'F1e_split_half.csv'),index=False)

# f) ancestry: directional concordance + genome-wide rho (with EUR-vs-EUR ceiling)
a=pd.read_csv(os.path.join(RES,'ancestry_replication_gse161643.tsv'),sep='\t')
conc=float(a[a.metric=='directional_concordance']['value'].iloc[0])           # 0.6471 (binom p=0.0088, n=85 EUR-consensus genes)
rho_be=float(a[a.comparison.str.contains('consensus\\(mean',na=False)]['rho'].iloc[0])  # BLACK vs EUR consensus genome-wide
rho_ee=float(a[a.metric=='spearman_rev_beta'][a.note.str.contains('ceiling',na=False)]['rho'].iloc[0])  # EUR-vs-EUR ceiling
pd.DataFrame([
 {'kind':'directional concordance','metric':'Black cohort vs European (co-significant genes)','value':round(conc,3),'ref':0.5,'panel':'concordance'},
 {'kind':'genome-wide ρ','metric':'Black cohort vs European','value':round(rho_be,3),'ref':0.0,'panel':'rho'},
 {'kind':'genome-wide ρ','metric':'EUR vs EUR (ceiling)','value':round(rho_ee,3),'ref':0.0,'panel':'rho'},
]).to_csv(os.path.join(OUT,'F1f_ancestry.csv'),index=False)
print('F1e: txn sp_signed=%.2f rev_enr=%.1f | meth sp_signed=%.2f rev_enr=%.1f'%(ta['sp_signed'],ta['rev_enr'],ma['sp_signed'],ma['rev_enr']))
print('F1f: concordance=%.3f (chance .5) | BLACK-EUR rho=%.3f | EUR-EUR ceiling=%.3f'%(conc,rho_be,rho_ee))

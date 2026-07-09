#!/usr/bin/env python
"""prep_ED2_ED6_data.py — ED2 (cross-ancestry detail) + ED6 (drug-class signature overlap)."""
import os, numpy as np, pandas as pd
HERE=os.path.dirname(os.path.abspath(__file__)); SIG=os.path.dirname(os.path.dirname(HERE))
RES=os.path.join(SIG,'results'); OUT=os.path.join(HERE,'source_data'); os.makedirs(OUT,exist_ok=True)
R=lambda p: pd.read_csv(os.path.join(RES,p),sep='\t')

# ===== ED2: cross-ancestry detail =====
a=R('ancestry_replication_gse161643.tsv')
sp=a[a.metric=='spearman_rev_beta'].copy()
def lab(c):
    s=str(c).lower()
    if 'consensus' in s: return 'BLACK vs EUR consensus'
    if 'masters' in s and 'strride' in s: return 'EUR vs EUR (ceiling)'
    if 'masters' in s: return 'BLACK vs EUR (Masters)'
    if 'strride' in s: return 'BLACK vs EUR (STRRIDE)'
    return str(c)[:30]
sp['comparison']=sp['comparison'].map(lab); sp['rho']=pd.to_numeric(sp['rho'],errors='coerce')
sp=sp.dropna(subset=['rho'])[['comparison','rho','n_shared']]
sp['kind']=np.where(sp.comparison.str.contains('ceiling'),'EUR ceiling','BLACK vs EUR')
sp.to_csv(os.path.join(OUT,'ED2a_genomewide.csv'),index=False)
# canonical genes (CANON:: rows): BLACK beta + EUR betas + sign-match
can=a[a.comparison.astype(str).str.startswith('CANON::')].copy()
can['gene']=can['comparison'].str.replace('CANON::','',regex=False)
bcol=[c for c in a.columns if 'BLACK' in c and c.endswith('beta')][0]
ecol=[c for c in a.columns if 'GSE157585_EUR' in c and c.endswith('beta')][0]
can['black_beta']=pd.to_numeric(can[bcol],errors='coerce'); can['eur_beta']=pd.to_numeric(can[ecol],errors='coerce')
can['match']=can['black_matches_eur_sign'].astype(str)
can[['gene','black_beta','eur_beta','match']].dropna(subset=['black_beta']).to_csv(os.path.join(OUT,'ED2b_canonical.csv'),index=False)
conc=float(a[a.metric=='directional_concordance']['value'].iloc[0])
pd.DataFrame([{'concordance':round(conc,3),'p':0.0088,'n':85}]).to_csv(os.path.join(OUT,'ED2_summary.csv'),index=False)

# ===== ED6: drug-class signature overlap (fold-enrichment, orthogonal to F5e correlation) =====
ov=R('drug_class_overlap.tsv')
ov=ov[pd.notna(ov['fold'])].copy(); ov['fold']=pd.to_numeric(ov['fold'],errors='coerce')
ov=ov[ov['U']>=50]  # powered pairs only
# class-pair grouping
def cls(ca,cb):
    s={ca,cb}
    if s=={'DRUG'}: return 'within-DRUG'
    if s=={'DRUG','SURGERY'}: return 'DRUG×SURGERY'
    if s=={'SURGERY','DIET'}: return 'SURGERY×DIET'
    if s=={'DRUG','DIET'}: return 'DRUG×DIET'
    if s=={'SURGERY'}: return 'within-SURGERY'
    return '×'.join(sorted(s))
ov['classpair']=[cls(r.ca,r.cb) for r in ov.itertuples()]
ov[['a','b','U','fold','p','samedir','classpair']].to_csv(os.path.join(OUT,'ED6_pairs.csv'),index=False)
grp=ov.groupby('classpair').agg(mean_fold=('fold','mean'),n=('fold','size')).reset_index().sort_values('mean_fold')
grp.to_csv(os.path.join(OUT,'ED6_classmeans.csv'),index=False)
print('ED2 genome-wide:', dict(zip(sp.comparison,sp.rho.round(3))))
print('ED2 canonical genes:', list(can.gene))
print('ED2 concordance=%.3f'%conc)
print('ED6 class means:', dict(zip(grp.classpair,grp.mean_fold.round(3))))
print('ED6 within-DRUG depletion pair: metformin×empagliflozin fold=%.3f'%ov[(ov.a=='metformin')&(ov.b=='empagliflozin')].fold.iloc[0])

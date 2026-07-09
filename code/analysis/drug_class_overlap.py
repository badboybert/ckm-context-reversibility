#!/usr/bin/env python
"""drug_class_overlap.py — companion to build_drug_class.py.
Earlier tissue/intervention work reported SIG-SET OVERLAP fold-enrichment (e.g. "1.68x within-tissue").
Do the same here for the drug/surgery/diet plasma panels: for every pair, is the set of SIGNIFICANTLY
moved proteins shared more than chance?  This is orthogonal to the Spearman (which uses effect sizes):
overlap asks "do the same proteins RESPOND", correlation asks "do they move the same WAY/AMOUNT".

For each pair (shared-measured universe U):
  obs   = |sig_a & sig_b|           (sig = q<0.10 within U, undirected)
  exp   = |sig_a|*|sig_b|/|U|        (chance, hypergeometric mean)
  fold  = obs/exp                    (>1 = enriched overlap)
  p     = hypergeometric survival (>= obs)
  dir   = of the shared-sig, fraction SAME direction (links to the core)
Class-grouped mean fold (within-DRUG / within-SURGERY / DRUG x SURGERY / DRUG x DIET / SURGERY x DIET).
Output: results/drug_class_overlap.tsv + appends to results/drug_class_summary.txt
"""
import os, itertools, numpy as np, pandas as pd
from scipy.stats import hypergeom
SIG=os.path.dirname(os.path.dirname(os.path.abspath(__file__))); PROJ=os.path.dirname(SIG)
def P(*a): return os.path.join(PROJ,*a)
def S(*a): return os.path.join(SIG,*a)
MINU=50      # need a reasonable shared-measured universe
SIGT=0.10
PAN=[('metformin','data/proteome/nowak_metformin_proteome_reversal.tsv','DRUG'),
 ('dapagliflozin','data/proteome/dapagliflozin_proteome_reversal.tsv','DRUG'),
 ('empagliflozin','data/proteome/reversal_emperor.tsv','DRUG'),
 ('liraglutide','data/proteome/liraglutide_proteome_reversal.tsv','DRUG'),
 ('semaglutide','data/proteome/step_semaglutide_proteome_reversal.tsv','DRUG'),
 ('BBS','data/proteome/reversal_olink_bbs.tsv','SURGERY'),
 ('MS_bariatric','data/proteome/bariatric_ms_proteome_reversal.tsv','SURGERY'),
 ('DiRECT','data/proteome/reversal_somascan_direct.tsv','DIET')]
def load(f):
    d=pd.read_csv(P(f),sep='\t'); d=d.rename(columns={c:'UniProt' for c in d.columns if c.lower()=='uniprot'})
    d['UniProt']=d['UniProt'].astype(str).str.split(r'[;,_]').str[0].str.strip()
    d=d[d.UniProt.str.match(r'^[A-Z][A-Z0-9]{5,9}$',na=False)]
    qc=next((c for c in ['rev_q','rev_p','q','p'] if c in d.columns),None)
    d['q']=pd.to_numeric(d[qc],errors='coerce') if qc else np.nan
    d['b']=pd.to_numeric(d['rev_beta'],errors='coerce')
    return d.dropna(subset=['b']).sort_values('q',na_position='last').drop_duplicates('UniProt').set_index('UniProt')
D={l:load(f) for l,f,_ in PAN}; CLS={l:c for l,_,c in PAN}; labels=list(D)
for l in labels:
    nsig=int((D[l]['q']<SIGT).sum())
    print(f'  {l:14s} measured={len(D[l]):5d}  sig(q<{SIGT})={nsig}')
rows=[]; out=['','===== SIG-SET OVERLAP (fold-enrichment over chance; q<%.2f, shared-measured universe) ====='%SIGT]
def pr(s): out.append(s); print(s)
for a,b in itertools.combinations(labels,2):
    U=D[a].index.intersection(D[b].index)
    if len(U)<MINU:
        rows.append(dict(a=a,b=b,U=len(U),obs=np.nan,exp=np.nan,fold=np.nan,p=np.nan,samedir=np.nan)); continue
    sa=set(U[D[a].loc[U,'q']<SIGT]); sb=set(U[D[b].loc[U,'q']<SIGT])
    obs=len(sa&sb); exp=len(sa)*len(sb)/len(U) if len(U) else np.nan
    fold=obs/exp if exp and exp>0 else np.nan
    # hypergeom survival P(X>=obs): sf(obs-1)
    p=float(hypergeom.sf(obs-1,len(U),len(sa),len(sb))) if (len(sa) and len(sb)) else np.nan
    sh=list(sa&sb)
    samedir=np.mean([np.sign(D[a].loc[u,'b'])==np.sign(D[b].loc[u,'b']) for u in sh]) if sh else np.nan
    rows.append(dict(a=a,b=b,ca=CLS[a],cb=CLS[b],U=len(U),na=len(sa),nb=len(sb),obs=obs,exp=round(exp,2),fold=fold,p=p,samedir=samedir))
M=pd.DataFrame(rows); M.to_csv(S('results/drug_class_overlap.tsv'),sep='\t',index=False)
for _,r in M.iterrows():
    if pd.isna(r.fold): pr(f'  {r.a:14s} x {r.b:14s} U={int(r.U):4d}  (universe<%d or no sig — skipped)'%MINU); continue
    pr(f'  {r.a:14s} x {r.b:14s} U={int(r.U):4d} sig {int(r.na):3d}/{int(r.nb):3d}  obs={int(r.obs):3d} exp={r.exp:5.1f}  fold={r.fold:4.2f}x  p={r.p:.1e}  same-dir={r.samedir:.2f}')
def grp(pred):
    v=[r.fold for _,r in M.iterrows() if pd.notna(r.fold) and pred(r)]
    return (np.mean(v),len(v)) if v else (np.nan,0)
pr('\n  -- class-grouped mean fold-enrichment --')
B={'within-DRUG':lambda r:r.ca=='DRUG' and r.cb=='DRUG','within-SURGERY':lambda r:r.ca=='SURGERY' and r.cb=='SURGERY',
   'DRUG x SURGERY':lambda r:{r.ca,r.cb}=={'DRUG','SURGERY'},'DRUG x DIET':lambda r:{r.ca,r.cb}=={'DRUG','DIET'},
   'SURGERY x DIET':lambda r:{r.ca,r.cb}=={'SURGERY','DIET'}}
for nm,pred in B.items():
    m,k=grp(pred); pr(f'    {nm:16s} mean_fold={m:4.2f}x  (k={k})' if k else f'    {nm:16s} (no pair)')
open(S('results/drug_class_summary.txt'),'a',encoding='utf-8').write('\n'.join(out))
print('\nwrote results/drug_class_overlap.tsv + appended drug_class_summary.txt')

#!/usr/bin/env python
"""build_drug_class.py — DRUG-CLASS analysis in plasma proteome.

Question (user): group all DRUGS together (metformin / SGLT2i / GLP-1RA) and test whether they
share a plasma reversal signature DISTINCT from lifestyle (diet) and surgery (bariatric).

Method (direction-only, per weight-loss.omics/CLAUDE.md cross-platform rule):
  - key every reversal table by UniProt; sign = post-minus-pre (toward-lean).
  - pairwise Spearman of rev_beta on shared proteins (min overlap MINOV).
  - class-grouped: within-DRUG vs DRUG-vs-SURGERY vs DRUG-vs-DIET vs SURGERY-vs-DIET.
  - drug-shared core = same-direction reversal in >=2/ k_drug drug panels at sig, NOT shared by surgery/diet.
  - universal core = same-direction across all 3 classes.
  - SGLT2i probe: does SGLT2i (low weight loss) track the drug axis or separate?
Output: results/drug_class_corr.tsv (matrix) + results/drug_class_summary.txt
"""
import os, itertools, numpy as np, pandas as pd
from scipy.stats import spearmanr
SIG=os.path.dirname(os.path.dirname(os.path.abspath(__file__))); PROJ=os.path.dirname(SIG)
def P(*a): return os.path.join(PROJ,*a)
def S(*a): return os.path.join(SIG,*a)
MINOV=30
# (label, file, class, mechanism, paired?)
PAN=[
 ('metformin',   'data/proteome/nowak_metformin_proteome_reversal.tsv', 'DRUG',    'biguanide', True),
 ('dapagliflozin','data/proteome/dapagliflozin_proteome_reversal.tsv',  'DRUG',    'SGLT2i',    True),
 ('empagliflozin','data/proteome/reversal_emperor.tsv',                 'DRUG',    'SGLT2i',    False),
 ('liraglutide', 'data/proteome/liraglutide_proteome_reversal.tsv',     'DRUG',    'GLP1RA',    True),
 ('semaglutide', 'data/proteome/step_semaglutide_proteome_reversal.tsv','DRUG',    'GLP1RA',    False),
 ('BBS',         'data/proteome/reversal_olink_bbs.tsv',                'SURGERY', 'bariatric', True),
 ('MS_bariatric','data/proteome/bariatric_ms_proteome_reversal.tsv',    'SURGERY', 'bariatric', True),
 ('DiRECT',      'data/proteome/reversal_somascan_direct.tsv',          'DIET',    'diet',      True),
]
def load(f):
    d=pd.read_csv(P(f),sep='\t')
    d=d.rename(columns={c:'UniProt' for c in d.columns if c.lower()=='uniprot'})
    d['UniProt']=d['UniProt'].astype(str).str.split(r'[;,_]').str[0].str.strip()
    d=d[d.UniProt.str.match(r'^[A-Z][A-Z0-9]{5,9}$',na=False)]   # UniProt: O/P/Q or A-N,R-Z first char
    qc=next((c for c in ['rev_q','rev_p','q','p'] if c in d.columns),None)
    d['q']=pd.to_numeric(d[qc],errors='coerce') if qc else np.nan
    d['b']=pd.to_numeric(d['rev_beta'],errors='coerce')
    d=d.dropna(subset=['b']).sort_values('q',na_position='last').drop_duplicates('UniProt')
    return d.set_index('UniProt')
D={l:load(f) for l,f,_,_,_ in PAN}
META={l:(cls,mech) for l,_,cls,mech,_ in PAN}
labels=[l for l,_,_,_,_ in PAN]
out=[]
def pr(s): out.append(s); print(s)
pr('# panels: '+', '.join(f'{l}(n={len(D[l])},{META[l][0]}/{META[l][1]})' for l in labels))
# ---- pairwise Spearman matrix ----
rows=[]
for a,b in itertools.combinations(labels,2):
    ov=D[a].index.intersection(D[b].index)
    if len(ov)<MINOV:
        rows.append(dict(a=a,b=b,n=len(ov),rho=np.nan,p=np.nan)); continue
    rho,p=spearmanr(D[a].loc[ov,'b'],D[b].loc[ov,'b'])
    rows.append(dict(a=a,b=b,n=len(ov),rho=rho,p=p,
        ca=META[a][0],cb=META[b][0],ma=META[a][1],mb=META[b][1]))
M=pd.DataFrame(rows); M.to_csv(S('results/drug_class_corr.tsv'),sep='\t',index=False)
pr('\n===== PAIRWISE SPEARMAN (rev_beta, shared UniProt, min overlap %d) ====='%MINOV)
for _,r in M.iterrows():
    rs=f'{r.rho:+.3f}' if pd.notna(r.rho) else '  NA '
    pr(f'  {r.a:14s} x {r.b:14s} n={int(r.n):4d}  rho={rs}'+('' if pd.isna(r.rho) else f'  p={r.p:.1e}'))
# ---- class-grouped within vs between ----
def grp(pred):
    v=[r.rho for _,r in M.iterrows() if pd.notna(r.rho) and pred(r)]
    return (np.mean(v),len(v),v) if v else (np.nan,0,[])
pr('\n===== CLASS-GROUPED MEAN Spearman =====')
buckets={
 'within-DRUG'     : lambda r: r.ca=='DRUG' and r.cb=='DRUG',
 'within-SURGERY'  : lambda r: r.ca=='SURGERY' and r.cb=='SURGERY',
 'DRUG x SURGERY'  : lambda r: {r.ca,r.cb}=={'DRUG','SURGERY'},
 'DRUG x DIET'     : lambda r: {r.ca,r.cb}=={'DRUG','DIET'},
 'SURGERY x DIET'  : lambda r: {r.ca,r.cb}=={'SURGERY','DIET'},
}
for nm,pred in buckets.items():
    m,k,v=grp(pred)
    pr(f'  {nm:16s} mean_rho={m:+.3f}  (k={k} pairs)' if k else f'  {nm:16s} (no pair with enough overlap)')
# within-drug by mechanism pair (is SGLT2i an outlier among drugs?)
pr('\n  -- within-DRUG detail (mechanism pairs) --')
for _,r in M.iterrows():
    if pd.notna(r.rho) and r.ca=='DRUG' and r.cb=='DRUG':
        pr(f'     {r.ma:9s}({r.a}) x {r.mb:9s}({r.b}): rho={r.rho:+.3f}')
# ---- shared cores (direction + significance) ----
SIGT=0.10
def sigdir(l):
    d=D[l]; s=d[(pd.to_numeric(d['q'],errors='coerce')<SIGT)]
    return {u:np.sign(b) for u,b in s['b'].items() if b==b and b!=0}
drug_labels=[l for l in labels if META[l][0]=='DRUG']
sd={l:sigdir(l) for l in labels}
# drug-shared core: same-direction sig in >=2 drug panels
from collections import defaultdict
ddir=defaultdict(list)
for l in drug_labels:
    for u,s in sd[l].items(): ddir[u].append(s)
drug_core={u:ss[0] for u,ss in ddir.items() if len(ss)>=2 and len(set(ss))==1}
surg_diet=set()
for l in labels:
    if META[l][0] in ('SURGERY','DIET'): surg_diet|=set(sd[l].keys())
drug_specific={u for u in drug_core if u not in surg_diet}
# universal core: same-direction sig in >=1 of each class
def classdir(cls):
    dd=defaultdict(list)
    for l in labels:
        if META[l][0]==cls:
            for u,s in sd[l].items(): dd[u].append(s)
    return {u:ss[0] for u,ss in dd.items() if len(set(ss))==1}
cd={c:classdir(c) for c in ['DRUG','SURGERY','DIET']}
universal={u for u in cd['DRUG'] if u in cd['SURGERY'] and u in cd['DIET']
           and cd['DRUG'][u]==cd['SURGERY'][u]==cd['DIET'][u]}
def names(uset,ref):
    g=[]
    for u in uset:
        for l in labels:
            if 'gene' in D[l].columns and u in D[l].index:
                gg=D[l].loc[u,'gene']
                if isinstance(gg,str): g.append(gg); break
        else: g.append(u)
    return g
pr('\n===== SHARED CORES (q<%.2f, direction-consistent) ====='%SIGT)
pr(f'  DRUG-shared core (>=2/{len(drug_labels)} drug panels, same dir): {len(drug_core)} proteins')
pr(f'  DRUG-SPECIFIC (drug core NOT moved by surgery/diet): {len(drug_specific)} -> '+', '.join(names(drug_specific,None)[:15]))
pr(f'  UNIVERSAL core (same dir in DRUG & SURGERY & DIET): {len(universal)} -> '+', '.join(names(universal,None)[:15]))
open(S('results/drug_class_summary.txt'),'w',encoding='utf-8').write('\n'.join(out))
print('\nwrote results/drug_class_corr.tsv + results/drug_class_summary.txt')

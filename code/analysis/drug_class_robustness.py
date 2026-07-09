#!/usr/bin/env python
"""drug_class_robustness.py — adversarial self-checks on the drug-class result before it is trusted.
Tests three refutations:
 (R1) PLATFORM confound: is the metformin "hub" / within-DRUG signal just same-platform pairs?
      -> stratify every pairwise rho by same-platform vs cross-platform.
 (R2) Is the within-DRUG mean (+0.296) load-bearing on ONE pair? -> report k and the actual pairs.
 (R3) Is the "drug axis" just the universal core resurfacing? -> recompute metformin x {empagliflozin,BBS}
      AFTER removing the universal-core proteins; if rho collapses, there is no drug-specific axis.
 (R4) Universal-39 core robustness: leave-one-panel-out direction consistency.
Reads the same panels as build_drug_class.py. Output: results/drug_class_robustness.txt
"""
import os, itertools, numpy as np, pandas as pd
from collections import defaultdict
from scipy.stats import spearmanr
SIG=os.path.dirname(os.path.dirname(os.path.abspath(__file__))); PROJ=os.path.dirname(SIG)
def P(*a): return os.path.join(PROJ,*a)
def S(*a): return os.path.join(SIG,*a)
MINOV=30; SIGT=0.10
PAN=[('metformin','data/proteome/nowak_metformin_proteome_reversal.tsv','DRUG','Olink+SomaScan'),
 ('dapagliflozin','data/proteome/dapagliflozin_proteome_reversal.tsv','DRUG','Olink'),
 ('empagliflozin','data/proteome/reversal_emperor.tsv','DRUG','Olink'),
 ('liraglutide','data/proteome/liraglutide_proteome_reversal.tsv','DRUG','SomaScan'),
 ('semaglutide','data/proteome/step_semaglutide_proteome_reversal.tsv','DRUG','SomaScan'),
 ('BBS','data/proteome/reversal_olink_bbs.tsv','SURGERY','Olink'),
 ('MS_bariatric','data/proteome/bariatric_ms_proteome_reversal.tsv','SURGERY','MS'),
 ('DiRECT','data/proteome/reversal_somascan_direct.tsv','DIET','SomaScan')]
def load(f):
    d=pd.read_csv(P(f),sep='\t'); d=d.rename(columns={c:'UniProt' for c in d.columns if c.lower()=='uniprot'})
    d['UniProt']=d['UniProt'].astype(str).str.split(r'[;,_]').str[0].str.strip()
    d=d[d.UniProt.str.match(r'^[A-Z][A-Z0-9]{5,9}$',na=False)]
    qc=next((c for c in ['rev_q','rev_p','q','p'] if c in d.columns),None)
    d['q']=pd.to_numeric(d[qc],errors='coerce') if qc else np.nan
    d['b']=pd.to_numeric(d['rev_beta'],errors='coerce')
    return d.dropna(subset=['b']).sort_values('q',na_position='last').drop_duplicates('UniProt').set_index('UniProt')
D={l:load(f) for l,f,_,_ in PAN}; CLS={l:c for l,_,c,_ in PAN}; PLAT={l:p for l,_,_,p in PAN}
labels=list(D)
out=[]; pr=lambda s:(out.append(s),print(s))
def rho(a,b,drop=None):
    ov=D[a].index.intersection(D[b].index)
    if drop is not None: ov=ov.difference(drop)
    if len(ov)<MINOV: return np.nan,len(ov)
    return spearmanr(D[a].loc[ov,'b'],D[b].loc[ov,'b'])[0],len(ov)
def plat_match(a,b):
    pa,pb=set(PLAT[a].split('+')),set(PLAT[b].split('+'))
    return 'same' if pa&pb else 'cross'
# universal core (recompute, same logic as main)
def sigdir(l):
    d=D[l]; s=d[pd.to_numeric(d['q'],errors='coerce')<SIGT]
    return {u:np.sign(b) for u,b in s['b'].items() if b==b and b!=0}
sd={l:sigdir(l) for l in labels}
def classdir(c):
    dd=defaultdict(list)
    for l in labels:
        if CLS[l]==c:
            for u,s in sd[l].items(): dd[u].append(s)
    return {u:ss[0] for u,ss in dd.items() if len(set(ss))==1}
cd={c:classdir(c) for c in ['DRUG','SURGERY','DIET']}
universal={u for u in cd['DRUG'] if u in cd['SURGERY'] and u in cd['DIET'] and cd['DRUG'][u]==cd['SURGERY'][u]==cd['DIET'][u]}
pr('===== R1: PLATFORM-STRATIFIED pairwise rho (same vs cross platform) =====')
same=[]; cross=[]
for a,b in itertools.combinations(labels,2):
    r,n=rho(a,b)
    if np.isnan(r): continue
    pm=plat_match(a,b); tag=f'{CLS[a][:3]}x{CLS[b][:3]}'
    (same if pm=='same' else cross).append(r)
    pr(f'  {a:13s}({PLAT[a]:13s}) x {b:13s}({PLAT[b]:13s}) {tag:8s} {pm:5s} rho={r:+.3f} n={n}')
pr(f'\n  SAME-platform mean rho = {np.mean(same):+.3f} (k={len(same)})')
pr(f'  CROSS-platform mean rho = {np.mean(cross):+.3f} (k={len(cross)})')
pr('  -> if same>>cross, the class-correlations are platform-driven, NOT biology.')
pr('\n===== R2: within-DRUG powered pairs =====')
wd=[(a,b,rho(a,b)) for a,b in itertools.combinations([l for l in labels if CLS[l]=="DRUG"],2)]
for a,b,(r,n) in wd:
    pr(f'  {a} x {b}: rho={r:+.3f} n={n}' if not np.isnan(r) else f'  {a} x {b}: n={n} (<{MINOV}, EXCLUDED)')
pr('\n===== R3: is the "drug axis" just the universal core? (drop universal, recompute) =====')
for a,b in [('metformin','empagliflozin'),('metformin','BBS'),('empagliflozin','DiRECT')]:
    r0,n0=rho(a,b); r1,n1=rho(a,b,drop=universal)
    pr(f'  {a} x {b}: full rho={r0:+.3f}(n={n0})  ->  minus-universal rho={r1:+.3f}(n={n1})  Δ={r1-r0:+.3f}')
pr('  -> if rho collapses toward 0 after removing the universal core, there is NO drug-specific axis.')
pr('\n===== R4: universal-39 core leave-one-panel-out robustness =====')
pr(f'  full universal core = {len(universal)} proteins')
for drop_l in labels:
    keep=[l for l in labels if l!=drop_l]
    cdk={c:None for c in ['DRUG','SURGERY','DIET']}
    dd={c:defaultdict(list) for c in ['DRUG','SURGERY','DIET']}
    for l in keep:
        for u,s in sd[l].items(): dd[CLS[l]][u].append(s)
    cdk={c:{u:ss[0] for u,ss in dd[c].items() if len(set(ss))==1} for c in dd}
    uni2={u for u in cdk['DRUG'] if u in cdk['SURGERY'] and u in cdk['DIET'] and cdk['DRUG'][u]==cdk['SURGERY'][u]==cdk['DIET'][u]}
    pr(f'  drop {drop_l:13s} -> universal core = {len(uni2)}')
open(S('results/drug_class_robustness.txt'),'w',encoding='utf-8').write('\n'.join(out))
print('\nwrote results/drug_class_robustness.txt')

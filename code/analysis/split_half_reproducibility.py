#!/usr/bin/env python
"""split_half_reproducibility.py — is PERSISTENCE reproducible within a SINGLE well-powered,
single-context panel? Splits one cohort in half (same tissue + intervention + cohort), recomputes
reversal + persistence in each half, and measures cross-half reproducibility. This isolates
POWER+NOISE from CONTEXT-SPECIFICITY: the cross-PANEL concordance was ~chance (verifier wf_0c3f3bf3),
but those panels differ in tissue/intervention/power. If split-half ALSO ~chance -> persistence is
noise-limited (low-movers unmeasurable). If split-half HIGH but cross-panel low -> persistence is
real but CONTEXT-SPECIFIC (a determinant story, modeled per-context).

Panel = GSE141221 (DiOGenes adipose, LCD diet, n=220 pairs). 10 random splits.
Reports: discrete persistence concordance (obs/expected), reversible concordance, and CONTINUOUS
reproducibility (Spearman of signed rev_beta and of |rev_beta| between halves).
"""
import os, gzip, json, numpy as np, pandas as pd
from scipy import stats
SIG=os.path.dirname(os.path.dirname(os.path.abspath(__file__))); PROJ=os.path.dirname(SIG)
def P(*a): return os.path.join(PROJ,*a)
def S(*a): return os.path.join(SIG,*a)
RAW=P('data/transcriptome/raw')

# --- pairs (replicate builder) ---
gmeta=json.load(open(os.path.join(RAW,'gse141221_meta.json')))
def strip_b(k): return k[:-1] if k.endswith('b') else k
ic={(m['indiv'],m['cid']):strip_b(m['key']) for m in gmeta}
counts=pd.read_csv(gzip.open(os.path.join(RAW,'GSE141221_cleaned_count.csv.gz'),'rt'),index_col=0)
cs=set(counts.columns)
pairs=[(ic[(iv,'CID1')],ic[(iv,'CID2')]) for iv in sorted({m['indiv'] for m in gmeta})
       if ic.get((iv,'CID1')) in cs and ic.get((iv,'CID2')) in cs]
n=len(pairs); print(f'GSE141221 pairs: {n}')

# full-panel logCPM filter + base_sd/abundance for binning (stable baseline property)
allcols=[c for p in pairs for c in p]
lib=counts[allcols].sum(axis=0)
ens2sym=json.load(open(os.path.join(RAW,'gencode_v25_ens2sym.json')))
def logcpm(cols): return np.log2(counts[cols].divide(counts[cols].sum(axis=0),axis=1)*1e6+1.0)
pre_all=[p[0] for p in pairs]
keep=( (counts[allcols].divide(lib,axis=1)*1e6>=1).sum(axis=1) >= 0.2*len(allcols) )
genes=counts.index[keep]
base=logcpm(pre_all).loc[genes]
base_sd=base.std(axis=1,ddof=1); base_mean=base.mean(axis=1)

def rev_beta(half):
    pre=[pairs[i][0] for i in half]; post=[pairs[i][1] for i in half]
    return (logcpm(post).loc[genes].values - logcpm(pre).loc[genes].values).mean(axis=1)

def persistence_flags(rb):
    d=pd.DataFrame({'za':np.abs(rb),'bsd':base_sd.values,'ab':base_mean.values},index=genes).dropna()
    d['sdb']=pd.qcut(d['bsd'].rank(method='first'),20,labels=False,duplicates='drop')
    d['mb']=pd.qcut(d['ab'].rank(method='first'),3,labels=False,duplicates='drop')
    d['cell']=d['sdb'].astype(str)+'|'+d['mb'].astype(str)
    per=pd.Series(False,index=d.index); rev=pd.Series(False,index=d.index)
    for _,g in d.groupby('cell',sort=False):
        if len(g)<6: continue
        lo=g['za'].quantile(1/3); hi=g['za'].quantile(2/3)
        per[g.index[g['za']<=lo]]=True; rev[g.index[g['za']>=hi]]=True
    return per,rev,d['za']

res=[]
for seed in range(10):
    rng=np.random.RandomState(seed); idx=rng.permutation(n); h1=idx[:n//2]; h2=idx[n//2:]
    rb1,rb2=rev_beta(h1),rev_beta(h2)
    p1,r1,za1=persistence_flags(rb1); p2,r2,za2=persistence_flags(rb2)
    common=p1.index.intersection(p2.index)
    p1,p2,r1,r2=p1[common],p2[common],r1[common],r2[common]
    both_p=(p1&p2).mean(); exp_p=p1.mean()*p2.mean()
    both_r=(r1&r2).mean(); exp_r=r1.mean()*r2.mean()
    sp_signed=stats.spearmanr(pd.Series(rb1,index=genes)[common],pd.Series(rb2,index=genes)[common]).correlation
    sp_abs=stats.spearmanr(za1[common],za2[common]).correlation
    res.append(dict(persist_enr=both_p/exp_p,rev_enr=both_r/exp_r,sp_signed=sp_signed,sp_abs=sp_abs))
R=pd.DataFrame(res)
print(f'\nGSE141221 adipose split-half (n=110 vs 110, same tissue+lever+cohort), 10 splits:')
print(f'  DISCRETE persistence concordance enrichment : {R.persist_enr.mean():.2f}x  (sd {R.persist_enr.std():.2f})   [~1=chance, >1=reproducible]')
print(f'  DISCRETE reversible  concordance enrichment : {R.rev_enr.mean():.2f}x  (sd {R.rev_enr.std():.2f})')
print(f'  CONTINUOUS Spearman signed rev_beta (A vs B): {R.sp_signed.mean():.3f}  (sd {R.sp_signed.std():.3f})   [reproducibility of the magnitude/direction]')
print(f'  CONTINUOUS Spearman |rev_beta|       (A vs B): {R.sp_abs.mean():.3f}  (sd {R.sp_abs.std():.3f})')
print('\nINTERPRETATION:')
print('  persist_enr ~1 + sp_signed HIGH  -> magnitude reproducible but DISCRETE low-mover label is noise (persist pole unmodellable as a class; use continuous score)')
print('  persist_enr >>1                  -> persistence reproducible within context -> cross-panel chance = CONTEXT-SPECIFICITY (real biology)')
print('  sp_signed LOW                    -> even continuous reversal not reproducible at n=110 (pure power limit)')

# --- persist orphan figure numbers (output-only; math unchanged) ---
_out=pd.DataFrame([
    dict(metric='persist_enr',        value=R.persist_enr.mean(), sd=R.persist_enr.std()),
    dict(metric='rev_enr',            value=R.rev_enr.mean(),     sd=R.rev_enr.std()),
    dict(metric='sp_signed',          value=R.sp_signed.mean(),   sd=R.sp_signed.std()),
    dict(metric='sp_abs',             value=R.sp_abs.mean(),      sd=R.sp_abs.std()),
])
_op=S('results/orphan/split_half_adipose.tsv'); _out.to_csv(_op,sep='\t',index=False)
print(f'wrote {_op}')

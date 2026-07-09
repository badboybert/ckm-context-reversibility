#!/usr/bin/env python
"""Methylation confirmation of the split-half persistence-reproducibility test.
Panel = DIRECT-PLUS (E-MTAB-12527), whole blood, diet, n=256 paired (the DEFINITIVE powered panel
where the original cross-panel at-chance concordance was measured). Splits 256 individuals 128 vs 128,
recomputes per-CpG delta + persistence in each half, measures cross-half reproducibility.
CpGs SUBSAMPLED (every 8th, ~108k) via chunked streaming read of the 3.5GB beta matrix for tractability.
"""
import os, gzip, numpy as np, pandas as pd
from scipy import stats
SIG=os.path.dirname(os.path.dirname(os.path.abspath(__file__))); PROJ=os.path.dirname(SIG)
def P(*a): return os.path.join(PROJ,*a)
DDIR=P('data/methylation/directplus_e12527')
STEP=8   # keep every 8th CpG

ss=pd.read_csv(os.path.join(DDIR,'blood_samplesheet.tsv'),sep='\t',dtype={'basename':str,'individual':str,'time':str})
pre,post={}, {}
for _,r in ss.iterrows():
    (pre if str(r.time)=='0' else post)[r.individual]=r.basename

print('streaming beta matrix (subsample every %d CpGs)...'%STEP)
chunks=[]
for i,ch in enumerate(pd.read_csv(gzip.open(os.path.join(DDIR,'directplus_blood_betas.tsv.gz'),'rt'),
                                  sep='\t',index_col=0,chunksize=40000)):
    chunks.append(ch.iloc[::STEP].astype('float32'))
beta=pd.concat(chunks); beta.index.name='cpg'
inds=sorted(i for i in pre if i in post and pre[i] in beta.columns and post[i] in beta.columns)
n=len(inds); print(f'beta {beta.shape}; paired individuals {n}')
A=beta[[pre[i] for i in inds]].to_numpy('float32')   # baseline
B=beta[[post[i] for i in inds]].to_numpy('float32')  # post
cpg=beta.index.to_numpy()
base_mean=np.nanmean(A,axis=1); base_sd=np.nanstd(A,axis=1,ddof=1); edge=np.minimum(base_mean,1-base_mean)

def flags(rb):
    d=pd.DataFrame({'za':np.abs(rb),'bsd':base_sd,'ed':edge}).dropna()
    d['sdb']=pd.qcut(d['bsd'].rank(method='first'),50,labels=False,duplicates='drop')
    d['mb']=pd.qcut(d['ed'].rank(method='first'),3,labels=False,duplicates='drop')
    d['cell']=d['sdb'].astype(str)+'|'+d['mb'].astype(str)
    per=pd.Series(False,index=d.index); rev=pd.Series(False,index=d.index)
    for _,g in d.groupby('cell',sort=False):
        if len(g)<6: continue
        lo=g['za'].quantile(1/3); hi=g['za'].quantile(2/3)
        per[g.index[g['za']<=lo]]=True; rev[g.index[g['za']>=hi]]=True
    return per,rev,d['za']

res=[]
for seed in range(8):
    rng=np.random.RandomState(seed); idx=rng.permutation(n); h1=idx[:n//2]; h2=idx[n//2:]
    rb1=np.nanmean(B[:,h1]-A[:,h1],axis=1); rb2=np.nanmean(B[:,h2]-A[:,h2],axis=1)
    p1,r1,za1=flags(rb1); p2,r2,za2=flags(rb2)
    c=p1.index.intersection(p2.index)
    ep=p1[c].mean()*p2[c].mean(); er=r1[c].mean()*r2[c].mean()
    res.append(dict(persist_enr=(p1[c]&p2[c]).mean()/ep, rev_enr=(r1[c]&r2[c]).mean()/er,
        sp_signed=stats.spearmanr(rb1[c],rb2[c]).correlation, sp_abs=stats.spearmanr(za1[c],za2[c]).correlation))
R=pd.DataFrame(res)
print(f'\nDIRECT-PLUS methylation split-half (n=128 vs 128, same tissue+lever+cohort), {len(R)} splits, ~{len(cpg)//1000}k CpGs:')
print(f'  DISCRETE persistence concordance enrichment : {R.persist_enr.mean():.2f}x (sd {R.persist_enr.std():.2f})  [~1=chance]')
print(f'  DISCRETE reversible  concordance enrichment : {R.rev_enr.mean():.2f}x (sd {R.rev_enr.std():.2f})')
print(f'  CONTINUOUS Spearman signed delta (A vs B)   : {R.sp_signed.mean():.3f} (sd {R.sp_signed.std():.3f})')
print(f'  CONTINUOUS Spearman |delta|       (A vs B)   : {R.sp_abs.mean():.3f} (sd {R.sp_abs.std():.3f})')

# --- persist orphan figure numbers (output-only; math unchanged) ---
SIG=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_out=pd.DataFrame([
    dict(metric='persist_enr', value=R.persist_enr.mean(), sd=R.persist_enr.std()),
    dict(metric='rev_enr',     value=R.rev_enr.mean(),     sd=R.rev_enr.std()),
    dict(metric='sp_signed',   value=R.sp_signed.mean(),   sd=R.sp_signed.std()),
    dict(metric='sp_abs',      value=R.sp_abs.mean(),      sd=R.sp_abs.std()),
])
_op=os.path.join(SIG,'results','orphan','split_half_methylation.tsv'); _out.to_csv(_op,sep='\t',index=False)
print(f'wrote {_op}')

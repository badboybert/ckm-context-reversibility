#!/usr/bin/env python
"""responder_stratified.py — does reversibility track CLINICAL RESPONSE? GSE273902 (whole blood,
bariatric surgery) carries T2D-outcome strata. Split paired patients into REMISSION (responders),
NO-REMISSION (non-responders), and NO-DIABETES (reference: surgery without T2D), recompute per-gene
reversal SEPARATELY per group, and ask: (a) does blood reverse MORE in remitters? (b) do the groups
reverse the SAME genes? (c) which genes revert in remitters but PERSIST in non-remitters = residual-risk
candidates? Small n/group -> aggregate measures robust, per-gene differential = hypothesis-generating.
"""
import os, gzip, numpy as np, pandas as pd
from scipy import stats
SIG=os.path.dirname(os.path.dirname(os.path.abspath(__file__))); PROJ=os.path.dirname(SIG)
def P(*a): return os.path.join(PROJ,*a)
def S(*a): return os.path.join(SIG,*a)
RAW=P('data/transcriptome/raw')
si=pd.read_csv(gzip.open(os.path.join(RAW,'GSE273902_Sample_Info.tsv.gz'),'rt'),sep='\t',index_col=0)
cnt=pd.read_csv(gzip.open(os.path.join(RAW,'GSE273902_Counts_Table.tsv.gz'),'rt'),sep='\t',index_col=0)
gi=pd.read_csv(gzip.open(os.path.join(RAW,'GSE273902_Gene_Info.tsv.gz'),'rt'),sep='\t').drop_duplicates('GENEID').set_index('GENEID')['SYMBOL']

def grp(s):
    s=str(s)
    if 'Remission' in s and 'No remission' not in s: return 'remission'
    if 'No remission' in s: return 'no_remission'
    if 'No diabetes' in s: return 'no_diabetes'
    return 'other'
si['grp']=si['Patient subgroup'].map(grp)
# pair each patient: V1 (pre) + latest post visit present in counts
pairs=[]   # (patient, pre, post, grp)
for pt,g in si.groupby('Patient'):
    g=g[g.index.isin(cnt.columns)]
    pre=g[g['Visit']=='V1']; post=g[g['Visit'].isin(['V3','V4'])]
    if len(pre) and len(post):
        po=post.sort_values('Visit').iloc[-1]   # latest post
        pairs.append((pt,pre.index[0],po.name,pre['grp'].iloc[0]))
pr=pd.DataFrame(pairs,columns=['patient','pre','post','grp'])
print('paired patients per stratum:', dict(pr.grp.value_counts()))

def rev(sub):
    pre=list(sub.pre); post=list(sub.post); m=cnt[pre+post].astype(float)
    lib=m.sum(); cpm=m.divide(lib,axis=1)*1e6; keep=(cpm>=1).sum(axis=1)>=0.5*len(m.columns)
    lc=np.log2(cpm[keep]+1.0)
    A=lc[pre].values; B=lc[post].values; d=B-A; n=d.shape[1]
    rb=d.mean(1); t,p=stats.ttest_rel(B,A,axis=1)
    out=pd.DataFrame({'ensg':lc.index,'rev_beta':rb,'p':np.where(np.isfinite(p),p,1.0),'n':n})
    out['gene']=out['ensg'].map(gi); out=out.dropna(subset=['gene']).sort_values('rev_beta',key=abs,ascending=False).drop_duplicates('gene')
    return out.set_index('gene')

R={g:rev(pr[pr.grp==g]) for g in ['remission','no_remission','no_diabetes'] if (pr.grp==g).sum()>=4}
print('\n=== (a) OVERALL reversal magnitude per stratum (median |rev_beta|, # nominal-sig genes) ===')
for g,d in R.items():
    print(f'  {g:13s} (n_pairs={d.n.iloc[0]:2d}): median|rev|={d.rev_beta.abs().median():.4f}  sig(p<.01)={int((d.p<0.01).sum())}  sig(p<.05)={int((d.p<0.05).sum())}')
print('\n=== (b) do strata reverse the SAME genes? (Spearman rev_beta on shared genes) ===')
import itertools
for a,b in itertools.combinations(R,2):
    sh=R[a].index.intersection(R[b].index); rho=stats.spearmanr(R[a].loc[sh,'rev_beta'],R[b].loc[sh,'rev_beta']).correlation
    print(f'  {a:13s} vs {b:13s}: Spearman={rho:+.2f} (n={len(sh)})')
if 'remission' in R and 'no_remission' in R:
    rem,nor=R['remission'],R['no_remission']; sh=rem.index.intersection(nor.index)
    m=pd.DataFrame({'rem':rem.loc[sh,'rev_beta'],'rem_p':rem.loc[sh,'p'],'nor':nor.loc[sh,'rev_beta'],'nor_p':nor.loc[sh,'p']})
    # residual-risk candidates: revert in REMISSION (sig, large) but NOT in NO-REMISSION (small)
    cand=m[(m.rem_p<0.05)&(m.rem.abs()>0.5)&(m.nor.abs()<0.2)].reindex(m.rem.abs().sort_values(ascending=False).index).dropna()
    cand=cand[~cand.index.astype(str).str.contains('-AS|LINC|^AC|^AL|^IG[HKL]|orf|^RNU|^MIR|^HB',regex=True,na=False)]
    print('\n=== (c) RESIDUAL-RISK candidates: revert in REMITTERS but PERSIST (fail to revert) in non-remitters ===')
    print('  '+', '.join(cand.index[:15]) if len(cand) else '  (none pass the strict filter at this n)')
    m.assign(gene=m.index).to_csv(S('results/responder_stratified_GSE273902.tsv'),sep='\t',index=False)
    print('  wrote results/responder_stratified_GSE273902.tsv (per-gene rev_beta by stratum)')
print('\nCAVEAT: small n/stratum -> aggregate measures (a,b) robust; per-gene residual-risk list (c) is hypothesis-generating.')

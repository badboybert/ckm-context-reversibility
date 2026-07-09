#!/usr/bin/env python
"""responder_dissociation_hardened.py — RE-RUN + HARDEN the responder / clinical-dissociation arm.

GSE273902 (whole blood, bariatric surgery) carries T2D-outcome strata: REMISSION (n=7 paired),
NO-REMISSION (n=8), NO-DIABETES (n=8, surgery without T2D = generic-surgical reference).

This script does NOT change the per-gene reversal computation in responder_stratified.py — it reuses
the SAME rev() pipeline (CPM>=1 in >=50% samples; log2; paired delta = post-pre; paired t-test) and
adds hardening:
  (1) per-stratum molecular-reversal MAGNITUDE = median |rev_beta| among SIGNIFICANT genes (p<.05 AND p<.01),
      remitters vs non-remitters (and no-diabetes reference).
  (2) bootstrap + leave-one-subject-out over SUBJECTS within strata -> magnitude-difference 95% CI and the
      fraction of resamples where non-remitter >= remitter -> does the dissociation survive at this n?
  (3) RESEMBLANCE: Spearman(non-remitter delta, no-diabetes delta) and Spearman(remitter delta, no-diabetes delta)
      on shared genes -> corrects the record (+0.72 was nor~no_diabetes, NOT rem-vs-nor = +0.18).

Persists results/responder_dissociation_hardened.tsv + results/responder_dissociation_bootstrap.tsv.
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

pairs=[]
for pt,g in si.groupby('Patient'):
    g=g[g.index.isin(cnt.columns)]
    pre=g[g['Visit']=='V1']; post=g[g['Visit'].isin(['V3','V4'])]
    if len(pre) and len(post):
        po=post.sort_values('Visit').iloc[-1]
        pairs.append((pt,pre.index[0],po.name,pre['grp'].iloc[0]))
pr=pd.DataFrame(pairs,columns=['patient','pre','post','grp'])
print('paired patients per stratum:', dict(pr.grp.value_counts()))

def rev(sub):
    """IDENTICAL to responder_stratified.rev(): per-gene reversal on the given subject pairs."""
    pre=list(sub.pre); post=list(sub.post); m=cnt[pre+post].astype(float)
    lib=m.sum(); cpm=m.divide(lib,axis=1)*1e6; keep=(cpm>=1).sum(axis=1)>=0.5*len(m.columns)
    lc=np.log2(cpm[keep]+1.0)
    A=lc[pre].values; B=lc[post].values; d=B-A; n=d.shape[1]
    rb=d.mean(1); t,p=stats.ttest_rel(B,A,axis=1)
    out=pd.DataFrame({'ensg':lc.index,'rev_beta':rb,'p':np.where(np.isfinite(p),p,1.0),'n':n})
    out['gene']=out['ensg'].map(gi); out=out.dropna(subset=['gene']).sort_values('rev_beta',key=abs,ascending=False).drop_duplicates('gene')
    return out.set_index('gene')

STRATA=['remission','no_remission','no_diabetes']
R={g:rev(pr[pr.grp==g]) for g in STRATA}

# ---------- (1) MAGNITUDE among SIGNIFICANT genes ----------
print('\n=== (1) molecular-reversal MAGNITUDE per stratum (median |rev_beta|) ===')
rows=[]
for g in STRATA:
    d=R[g]; npair=int(d.n.iloc[0])
    med_all=d.rev_beta.abs().median()
    sig05=d[d.p<0.05]; sig01=d[d.p<0.01]
    med05=sig05.rev_beta.abs().median(); med01=sig01.rev_beta.abs().median()
    rows.append(dict(stratum=g,n_pairs=npair,n_genes=len(d),
                     median_abs_all=round(med_all,4),
                     n_sig05=len(sig05),median_abs_sig05=round(med05,4),
                     n_sig01=len(sig01),median_abs_sig01=round(med01,4)))
    print(f'  {g:13s} (n={npair}): all median|rev|={med_all:.4f} | sig(p<.05) n={len(sig05):5d} median|rev|={med05:.4f} | sig(p<.01) n={len(sig01):4d} median|rev|={med01:.4f}')
mag=pd.DataFrame(rows)

# point dissociation: remitters vs non-remitters, among SIG (p<.05) genes
rem_m=mag.set_index('stratum').loc['remission','median_abs_sig05']
nor_m=mag.set_index('stratum').loc['no_remission','median_abs_sig05']
print(f'\n  DISSOCIATION (median|rev| among p<.05 genes): remission={rem_m:.4f}  no_remission={nor_m:.4f}  diff(nor-rem)={nor_m-rem_m:+.4f}')

# ---------- (2) BOOTSTRAP + LOSO over SUBJECTS ----------
def stratum_sig_median(sub,alpha=0.05):
    """recompute rev() on a (possibly resampled) subject set, return median|rev_beta| among p<alpha genes + count."""
    if len(sub)<2: return np.nan,0
    d=rev(sub)
    s=d[d.p<alpha]
    return (s.rev_beta.abs().median() if len(s) else np.nan), len(s)

rng=np.random.default_rng(20260624)
NB=2000
rem_sub=pr[pr.grp=='remission']; nor_sub=pr[pr.grp=='no_remission']
print(f'\n=== (2) BOOTSTRAP (subject resampling within strata, NB={NB}) — magnitude diff nor-rem, sig p<.05 ===')
diffs=[]; rem_b=[]; nor_b=[]
for _ in range(NB):
    rs=rem_sub.iloc[rng.integers(0,len(rem_sub),len(rem_sub))]
    ns=nor_sub.iloc[rng.integers(0,len(nor_sub),len(nor_sub))]
    rm,_=stratum_sig_median(rs); nm,_=stratum_sig_median(ns)
    if np.isfinite(rm) and np.isfinite(nm):
        rem_b.append(rm); nor_b.append(nm); diffs.append(nm-rm)
diffs=np.array(diffs); rem_b=np.array(rem_b); nor_b=np.array(nor_b)
ci=np.percentile(diffs,[2.5,97.5]); frac_nor_ge_rem=np.mean(diffs>=0)
print(f'  valid resamples: {len(diffs)}/{NB}')
print(f'  rem median|rev| (sig): boot mean={rem_b.mean():.4f} [{np.percentile(rem_b,2.5):.4f},{np.percentile(rem_b,97.5):.4f}]')
print(f'  nor median|rev| (sig): boot mean={nor_b.mean():.4f} [{np.percentile(nor_b,2.5):.4f},{np.percentile(nor_b,97.5):.4f}]')
print(f'  diff (nor-rem): mean={diffs.mean():+.4f}  95% CI [{ci[0]:+.4f}, {ci[1]:+.4f}]')
print(f'  fraction of resamples with nor>=rem (dissociation direction holds): {frac_nor_ge_rem:.3f}')
survive_boot = ci[0] > 0
print(f'  -> CI excludes 0 (dissociation survives): {survive_boot}')

# leave-one-subject-out: drop one subject from each stratum at a time, recompute diff
print('\n--- leave-one-subject-out (drop 1 subject, recompute magnitude diff) ---')
loso=[]
# LOSO on remitters (hold nor full), then on non-remitters (hold rem full)
nor_full,_=stratum_sig_median(nor_sub)
rem_full,_=stratum_sig_median(rem_sub)
for i in range(len(rem_sub)):
    rm,_=stratum_sig_median(rem_sub.drop(rem_sub.index[i]))
    loso.append(('drop_rem',rem_sub.iloc[i].patient,nor_full-rm))
for i in range(len(nor_sub)):
    nm,_=stratum_sig_median(nor_sub.drop(nor_sub.index[i]))
    loso.append(('drop_nor',nor_sub.iloc[i].patient,nm-rem_full))
loso=pd.DataFrame(loso,columns=['which','subject','diff_nor_minus_rem'])
print(f'  LOSO diff range: [{loso.diff_nor_minus_rem.min():+.4f}, {loso.diff_nor_minus_rem.max():+.4f}]  all>0: {(loso.diff_nor_minus_rem>0).all()}')

# ---------- (3) RESEMBLANCE ----------
print('\n=== (3) RESEMBLANCE: Spearman of per-gene rev_beta vs no_diabetes (generic-surgical) ===')
nd=R['no_diabetes']
res={}
for g in ['remission','no_remission']:
    sh=R[g].index.intersection(nd.index)
    rho=stats.spearmanr(R[g].loc[sh,'rev_beta'],nd.loc[sh,'rev_beta']).correlation
    res[g]=(rho,len(sh))
    print(f'  {g:13s} vs no_diabetes: Spearman={rho:+.3f} (n={len(sh)})')
# rem vs nor for the record
sh=R['remission'].index.intersection(R['no_remission'].index)
rho_rn=stats.spearmanr(R['remission'].loc[sh,'rev_beta'],R['no_remission'].loc[sh,'rev_beta']).correlation
print(f'  remission     vs no_remission: Spearman={rho_rn:+.3f} (n={len(sh)})  [record-correction: this is the rem-vs-nor value]')

# ---------- PERSIST ----------
mag.to_csv(S('results/responder_dissociation_hardened.tsv'),sep='\t',index=False)
boot=pd.DataFrame([
    dict(metric='median_abs_sig05_remission_point',value=round(rem_m,4)),
    dict(metric='median_abs_sig05_no_remission_point',value=round(nor_m,4)),
    dict(metric='diff_nor_minus_rem_point',value=round(nor_m-rem_m,4)),
    dict(metric='boot_diff_mean',value=round(float(diffs.mean()),4)),
    dict(metric='boot_diff_ci_lo',value=round(float(ci[0]),4)),
    dict(metric='boot_diff_ci_hi',value=round(float(ci[1]),4)),
    dict(metric='boot_frac_nor_ge_rem',value=round(float(frac_nor_ge_rem),4)),
    dict(metric='boot_ci_excludes_0',value=int(survive_boot)),
    dict(metric='loso_min_diff',value=round(float(loso.diff_nor_minus_rem.min()),4)),
    dict(metric='loso_max_diff',value=round(float(loso.diff_nor_minus_rem.max()),4)),
    dict(metric='loso_all_positive',value=int((loso.diff_nor_minus_rem>0).all())),
    dict(metric='spearman_no_remission_vs_no_diabetes',value=round(float(res['no_remission'][0]),3)),
    dict(metric='spearman_remission_vs_no_diabetes',value=round(float(res['remission'][0]),3)),
    dict(metric='spearman_remission_vs_no_remission',value=round(float(rho_rn),3)),
    dict(metric='n_pairs_remission',value=int(R['remission'].n.iloc[0])),
    dict(metric='n_pairs_no_remission',value=int(R['no_remission'].n.iloc[0])),
    dict(metric='n_pairs_no_diabetes',value=int(R['no_diabetes'].n.iloc[0])),
    dict(metric='n_bootstrap',value=NB),
])
boot.to_csv(S('results/responder_dissociation_bootstrap.tsv'),sep='\t',index=False)
loso.to_csv(S('results/responder_dissociation_loso.tsv'),sep='\t',index=False)
print('\nwrote results/responder_dissociation_hardened.tsv, responder_dissociation_bootstrap.tsv, responder_dissociation_loso.tsv')

#!/usr/bin/env python
"""ADVERSARIAL confound check (independent re-derivation).
Reproduces the classifier's per-panel persistence call, then attacks:
  A. recompute median base_sd ratio persistent/reversible (per-panel + per-layer class)
  B. logistic regression persistent_in_panel ~ base_sd (WITHIN panel) -> residual gradient?
  C. within-decile corr(|za|, base_sd) -> is 10 deciles too coarse (low-edge clustering)?
  D. persistent-rate across base_sd quintiles (raw, no within-decile correction) for contrast
  E. tertile-split sensitivity (quartile / median split)
  F. low-edge test: within each decile, is base_sd of persistent < base_sd of reversible?
  G. base_sd sample mismatch: report n_base vs n_pairs(reversal)
"""
import os, numpy as np, pandas as pd
from scipy import stats
SIG=os.path.dirname(os.path.dirname(os.path.abspath(__file__))); PROJ=os.path.dirname(SIG)
def P(*a): return os.path.join(PROJ,*a)
def S(*a): return os.path.join(SIG,*a)
def nm(s): return str(s).strip().rstrip('*').strip().lower()
def pick(c,cands):
    for x in cands:
        if x in c: return x

MAN=pd.read_csv(S('results/panel_manifest.tsv'),sep='\t')
ELIG=MAN[MAN.persistence_eligible]
EFF=['rev_beta','effect','delta_beta']
bsd_tx=pd.read_csv(S('results/base_sd_transcriptome.tsv'),sep='\t')
bsd_mb=pd.read_csv(S('results/base_sd_metabolome.tsv'),sep='\t')

def base_sd_for(row,d):
    if row.layer=='transcriptome':
        # the classifier maps by mark_id but base_sd file uses panel name match
        pname={'liver_table':'liver_table','gse141221_adipose_table':'gse141221_adipose_table'}[row.panel]
        m=bsd_tx[bsd_tx.panel==pname].set_index('mark_id')['base_sd']
        return d['mark_key'].map(m)
    if row.layer=='metabolome':
        m=bsd_mb.set_index('mark_key')['base_sd']
        return d['mark_key'].map(m)

def panel_df(row):
    df=pd.read_csv(P(row.path),sep='\t'); cols=list(df.columns)
    ec=pick(cols,EFF)
    key=(df[pick(cols,['mark_id','metabolite'])].map(nm) if row.layer=='metabolome'
         else df[{'transcriptome':'mark_id','methylation':'cpg'}[row.layer]].astype(str))
    eff=pd.to_numeric(df[ec],errors='coerce')
    d=pd.DataFrame({'mark_key':key,'eff':eff})
    if row.layer=='methylation': d['base_sd']=pd.to_numeric(df['base_sd'],errors='coerce')
    d=d.dropna(subset=['eff']).drop_duplicates('mark_key')
    if row.layer!='methylation': d['base_sd']=base_sd_for(row,d).values
    d=d.dropna(subset=['base_sd'])
    d['za']=np.abs(d['eff'].values)
    try: d['dec']=pd.qcut(d['base_sd'],10,labels=False,duplicates='drop')
    except ValueError: d['dec']=0
    d['persistent']=False; d['reversible']=False
    for _,g in d.groupby('dec'):
        lo=g['za'].quantile(1/3); hi=g['za'].quantile(2/3)
        d.loc[g.index[g['za']<=lo],'persistent']=True
        d.loc[g.index[g['za']>=hi],'reversible']=True
    d['panel']=row.panel; d['layer']=row.layer
    return d

print("="*90)
print("ATTACK A/B/C/F: per-panel confound diagnostics")
print("="*90)
allp={}
for row in ELIG.itertuples(index=False):
    d=panel_df(row)
    allp[row.panel]=d
    # A: within-panel median base_sd ratio persistent/reversible
    bp=d.loc[d.persistent,'base_sd'].median()
    br=d.loc[d.reversible,'base_sd'].median()
    ratio=bp/br
    # B: logistic regression persistent ~ standardized base_sd (within panel)
    x=(np.log(d['base_sd'])-np.log(d['base_sd']).mean())/np.log(d['base_sd']).std()
    y=d['persistent'].astype(int).values
    # simple logit via statsmodels-free Newton; fall back to point-biserial corr
    try:
        import statsmodels.api as sm
        X=sm.add_constant(x.values)
        res=sm.Logit(y,X).fit(disp=0)
        b_logsd=res.params[1]; p_logsd=res.pvalues[1]
        reg=f"logit beta(log base_sd)={b_logsd:+.3f} p={p_logsd:.1e}"
    except Exception as e:
        r=stats.pointbiserialr(y,x.values)
        reg=f"pt-biserial r(persistent,log base_sd)={r.statistic:+.3f} p={r.pvalue:.1e}"
    # C: within-decile spearman corr(|za|, base_sd) pooled across deciles (residual gradient)
    cs=[]
    for _,g in d.groupby('dec'):
        if g['base_sd'].nunique()>3 and len(g)>=10:
            cs.append(stats.spearmanr(g['za'],g['base_sd']).statistic)
    wd_corr=np.nanmean(cs) if cs else float('nan')
    # F: low-edge — within each decile, base_sd(persistent) vs base_sd(reversible)
    edge=[]
    for _,g in d.groupby('dec'):
        gp=g.loc[g.persistent,'base_sd']; gr=g.loc[g.reversible,'base_sd']
        if len(gp)>=3 and len(gr)>=3:
            edge.append((gp.median()-gr.median())/g['base_sd'].median())
    edge_med=np.nanmedian(edge) if edge else float('nan')
    print(f"\n[{row.panel}] layer={row.layer} n_marks={len(d)} n_dec={d['dec'].nunique()}")
    print(f"  A median base_sd persistent={bp:.4g} reversible={br:.4g} RATIO={ratio:.3f}  (want ~1; <<1 = confound)")
    print(f"  B {reg}")
    print(f"  C within-decile mean Spearman(|za|,base_sd) = {wd_corr:+.3f}  (want ~0; >0 = persistent still low-base_sd within decile)")
    print(f"  F within-decile relative median base_sd (persistent-reversible)/med = {edge_med:+.3f}  (want ~0; <0 = persistent at LOW edge)")

print("\n"+"="*90)
print("ATTACK D/E: raw persistent-rate gradient & tertile-split sensitivity (DiRECT-largest panel each layer)")
print("="*90)
for pname in ['directplus','liver_table','gse141221_adipose_table','palau']:
    d=allp[pname]
    print(f"\n[{pname}] n={len(d)}")
    # D raw quintile persistent-rate using SINGLE global tertile (no within-decile)
    lo=d['za'].quantile(1/3); hi=d['za'].quantile(2/3)
    d['_glob_pers']=d['za']<=lo
    d['_q']=pd.qcut(d['base_sd'],5,labels=False,duplicates='drop')
    rate=d.groupby('_q')['_glob_pers'].mean()
    print("  D GLOBAL-tertile persistent rate by base_sd quintile (no within-decile correction):")
    print("    "+"  ".join(f"Q{int(q)+1}={r:.2f}" for q,r in rate.items()))
    # E sensitivity: redo within-decile with quartile & median splits, recompute ratio
    for frac,label in [(0.25,'quartile'),(0.5,'median')]:
        pers=pd.Series(False,index=d.index); rev=pd.Series(False,index=d.index)
        for _,g in d.groupby('dec'):
            if frac==0.5:
                med=g['za'].median(); pers.loc[g.index[g['za']<med]]=True; rev.loc[g.index[g['za']>=med]]=True
            else:
                lo2=g['za'].quantile(frac); hi2=g['za'].quantile(1-frac)
                pers.loc[g.index[g['za']<=lo2]]=True; rev.loc[g.index[g['za']>=hi2]]=True
        bp=d.loc[pers,'base_sd'].median(); br=d.loc[rev,'base_sd'].median()
        print(f"  E {label:8s} split ratio persistent/reversible base_sd = {bp/br:.3f}")

print("\n"+"="*90)
print("ATTACK G: base_sd SAMPLE-SET mismatch (n_base for variance vs n_pairs for reversal)")
print("="*90)
print(bsd_tx.groupby('panel')['n_base'].first())
print("metabolome palau n_base range:", bsd_mb['n_base'].min(), "-", bsd_mb['n_base'].max())
for row in ELIG.itertuples(index=False):
    print(f"  {row.panel}: reversal n_pairs={row.n_pairs}")

print("\n"+"="*90)
print("ATTACK H: LAYER-CLASS confound ratio reproduce (what the manuscript reports, base_sd_mean over panels)")
print("="*90)
for layer in ['methylation','transcriptome','metabolome']:
    rows=ELIG[ELIG.layer==layer]
    pp=pd.concat([panel_df(r) for r in rows.itertuples(index=False)],ignore_index=True)
    npan=pp.groupby('mark_key').size().rename('n_elig')
    agg=pp.groupby('mark_key').agg(n_persistent=('persistent','sum'),n_reversible=('reversible','sum'),
                                   base_sd_mean=('base_sd','mean')).join(npan)
    agg['class']='intermediate'
    agg.loc[agg.n_persistent==agg.n_elig,'class']='persistent'
    agg.loc[agg.n_reversible==agg.n_elig,'class']='consistently-reversible'
    pe=agg[agg['class']=='persistent']; rv=agg[agg['class']=='consistently-reversible']
    bp=pe.base_sd_mean.median(); br=rv.base_sd_mean.median()
    # also Mann-Whitney on full distribution
    mw=stats.mannwhitneyu(pe.base_sd_mean.dropna(), rv.base_sd_mean.dropna(), alternative='two-sided')
    print(f"[{layer}] n_per={len(pe)} n_rev={len(rv)} ratio={bp/br:.3f}  MannWhitney p={mw.pvalue:.2e} "
          f"(full-distribution test; ratio alone hides shape)")

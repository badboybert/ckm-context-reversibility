#!/usr/bin/env python
"""determinant_robustness.py — TASK B7 (review sec 4g).
Is the transcriptome determinant NULL an artifact of the measurability-adjustment
CHOICE?  Refit the SAME per-panel OLS + hk() meta under 3 OUTCOME definitions that
differ ONLY in how mark measurability enters the outcome, holding the gene UNIVERSE
and the feature matrix FIXED (the base_sd-eligible set) so the comparison isolates
the adjustment, not the gene set:

  (1) UNADJUSTED   : y = rank-INT(|rev_beta|) globally, NO measurability info.
  (2) COVARIATE    : y = rank-INT(|rev_beta|) globally, with base_sd + log-abundance
                     (base_mean, already log2-expression scale) added as OLS COVARIATES
                     (nuisance, partialled out) — NOT cell-residualized.
  (3) RESIDUALIZED : within (base_sd x abundance) cell rank-INT (CANONICAL) — must
                     reproduce results/determinant_meta.tsv pooled_beta (baseline anchor).

STANDALONE: fit logic COPIED from src/build_determinants.py (fit_panel, hk, FEATS,
gene_features). Does NOT import/execute build_determinants.py; does NOT touch its
canonical outputs. Writes ONLY results/determinant_modes.tsv (+ prints).
"""
import os, numpy as np, pandas as pd
import statsmodels.api as sm
from scipy import stats
SIG=os.path.dirname(os.path.dirname(os.path.abspath(__file__))); PROJ=os.path.dirname(SIG)
def P(*a): return os.path.join(PROJ,*a)
def S(*a): return os.path.join(SIG,*a)
def rint(x):
    x=np.asarray(x,float); r=stats.rankdata(x); return stats.norm.ppf((r-0.5)/len(x))

# ---- external gene features (COPIED verbatim from build_determinants.py) ----
F=lambda n: pd.read_csv(S('data/features',n),sep='\t')
lo=F('loeuf.tsv'); da=F('disease_assoc.tsv'); ta=F('tissue_tau.tsv')[['gene','tau']]
dr=F('druggability.tsv'); dr=dr[dr.gene.astype(str).str.match(r'^[A-Z][A-Z0-9-]*$')]
fc=F('functional_class.tsv'); fc=fc.assign(gene=fc.gene.astype(str).str.split('; ')).explode('gene')
ext=lo[['gene','loeuf']].merge(dr[['gene','n_drug_interactions']],on='gene',how='outer').merge(
    da[['gene','n_gwas_assoc']],on='gene',how='outer').merge(ta,on='gene',how='outer').merge(
    fc[['gene','is_secreted','is_enzyme','is_membrane']],on='gene',how='outer')
ext=ext.groupby('gene',as_index=False).max(numeric_only=True)
ext['n_drug_log']=np.log1p(ext['n_drug_interactions'].fillna(0)); ext['n_gwas_log']=np.log1p(ext['n_gwas_assoc'].fillna(0))

tx_loc=pd.read_csv(S('data/features/transcriptome_local.tsv'),sep='\t')
bsd_tx=pd.read_csv(S('results/base_sd_transcriptome.tsv'),sep='\t')
MAN=pd.read_csv(S('results/panel_manifest.tsv'),sep='\t'); MAN=MAN[MAN.contributes_to_score]

EXTNUM=['loeuf','tau','n_drug_log','n_gwas_log']; EXTBIN=['is_secreted','is_enzyme','is_membrane']
ARCH=['arch_nsig','arch_str']; OTHER=['causal_nonEGFR','has_arch']
FEATS=EXTNUM+EXTBIN+ARCH+OTHER+[f+'_miss' for f in ['loeuf','tau']]

def gene_features_tx():
    loc=tx_loc.rename(columns={'gene':'mark_key'}).copy(); loc['gene']=loc['mark_key']
    loc['arch_nsig']=loc['eqtl_n_cs'].fillna(0); loc['arch_str']=loc['eqtl_max_pip'].fillna(0); loc['has_arch']=1
    return loc

def build_d(row, loc):
    """Load one transcriptome panel; restrict to base_sd-eligible genes (the CANONICAL
    universe); attach standardized-ready features + measurability columns (bsd, log-abund,
    cell rank-INT outcome). Returns d with 'y_unadj','y_resid','csd','cab' + FEATS present."""
    df=pd.read_csv(P(row.path),sep='\t'); cols=list(df.columns)
    ec='rev_beta' if 'rev_beta' in cols else ('effect' if 'effect' in cols else 'delta_beta')
    key='mark_id'
    d=pd.DataFrame({'mark_key':df[key].astype(str),'rb':pd.to_numeric(df[ec],errors='coerce')}).dropna().drop_duplicates('mark_key')
    d['za']=d['rb'].abs()
    # base_sd-eligible universe (identical restriction to canonical residualize path)
    b=bsd_tx[bsd_tx.panel==row.panel].set_index('mark_id')
    d['bsd']=d['mark_key'].map(b['base_sd']); d['ab']=d['mark_key'].map(b['base_mean']); d=d.dropna(subset=['bsd','ab'])
    # (mode 3) within-cell rank-INT outcome (canonical)
    d['sb']=pd.qcut(d['bsd'].rank(method='first'),20,labels=False,duplicates='drop')
    d['mb']=pd.qcut(d['ab'].rank(method='first'),5,labels=False,duplicates='drop')
    d['y_resid']=d.groupby([d.sb,d.mb])['za'].transform(lambda s: rint(s.values) if len(s)>5 else np.nan)
    # features (COPIED merge/fill logic)
    d=d.merge(loc,on='mark_key',how='left')
    d['gene']=d['mark_key']
    d=d.merge(ext[['gene']+EXTNUM+EXTBIN],on='gene',how='left')
    d['has_arch']=d['has_arch'].fillna(0); d['arch_nsig']=d['arch_nsig'].fillna(0); d['arch_str']=d['arch_str'].fillna(0)
    d['causal_nonEGFR']=d['causal_nonEGFR'].fillna(0)
    for f in EXTBIN: d[f]=d[f].fillna(0)
    for f in ['loeuf','tau']: d[f+'_miss']=d[f].isna().astype(int)   # miss = per-gene isna (row-set independent)
    # NOTE: EXTNUM median-imputation is deferred to fit_mode so each mode imputes on ITS OWN
    # analysis set (exactly as canonical fit_panel does after its y-drop) -> mode 3 == canonical.
    return d

def fit_mode(d, mode):
    """OLS(HC3) of the mode's outcome on the standardized FEATS. For COVARIATE mode,
    base_sd + log-abundance are appended as standardized nuisance covariates and NOT
    reported. Returns per-feature (feature,beta,se) over the FEATS 'used' set."""
    feats=[f for f in FEATS if f in d.columns]
    if mode=='residualized':
        dd=d.dropna(subset=['y_resid']).copy(); y=dd['y_resid'].values
    else:
        dd=d.copy(); y=rint(dd['za'].values)          # global rank-INT, no measurability
    for f in EXTNUM: dd[f]=dd[f].fillna(dd[f].median())   # impute on THIS mode's analysis set
    X=dd[feats].apply(pd.to_numeric,errors='coerce').fillna(0)
    X=(X-X.mean())/X.std(ddof=0); X=X.loc[:,X.std()>0]; used=list(X.columns)
    if mode=='covariate':
        cov=pd.DataFrame({'cov_base_sd':dd['bsd'].values,'cov_log_abund':dd['ab'].values},index=X.index)
        cov=(cov-cov.mean())/cov.std(ddof=0); cov=cov.loc[:,cov.std()>0]
        Xf=pd.concat([X,cov],axis=1)
    else:
        Xf=X
    res=sm.OLS(y,sm.add_constant(Xf)).fit(cov_type='HC3')
    return pd.DataFrame([dict(feature=f,beta=res.params[f],se=res.bse[f],n=len(dd)) for f in used])

def hk(b,s):   # COPIED verbatim from build_determinants.py
    b=np.asarray(b,float); v=np.asarray(s,float)**2; k=len(b)
    if k==1: return b[0],np.sqrt(v[0]),k,0.0,0.0
    w=1/v; sw=w.sum(); bf=(w*b).sum()/sw; Q=(w*(b-bf)**2).sum(); df=k-1
    C=sw-(w**2).sum()/sw; t2=max(0,(Q-df)/C) if C>0 else 0; ws=1/(v+t2); br=(ws*b).sum()/ws.sum()
    se=np.sqrt(max((ws*(b-br)**2).sum()/(df*ws.sum()),1e-12)); return br,se*stats.t.ppf(0.975,df),k,Q,t2

loc=gene_features_tx()
rows=MAN[(MAN.layer=='transcriptome')&(MAN.panel.isin(set(bsd_tx.panel)))]
# pre-build d once per panel (feature matrix identical across modes)
DS={r.panel:build_d(r,loc) for r in rows.itertuples(index=False)}
print('panels (k=%d): '%len(DS)+', '.join(DS.keys()))

MODES=['unadjusted','covariate','residualized']
percoef=[]
for panel,d in DS.items():
    for mode in MODES:
        r=fit_mode(d,mode); r['panel']=panel; r['mode']=mode; percoef.append(r)
PC=pd.concat(percoef,ignore_index=True)

mrows=[]
for (mode,feature),g in PC.groupby(['mode','feature']):
    br,ci,k,Q,t2=hk(g.beta.values,g.se.values)
    df=k-1; I2=max(0.0,(Q-df)/Q)*100 if Q>0 else 0.0
    mrows.append(dict(mode=mode,feature=feature,k=k,pooled_beta=br,ci_halfwidth=ci,I2=I2))
M=pd.DataFrame(mrows)
M['mode']=pd.Categorical(M['mode'],categories=MODES,ordered=True)
M=M.sort_values(['mode','feature']).reset_index(drop=True)
M.to_csv(S('results/determinant_modes.tsv'),sep='\t',index=False)

# ---------- reporting ----------
pd.set_option('display.width',160,'display.max_columns',20)
print('\n=== BASELINE ANCHOR: mode 3 (residualized) vs canonical determinant_meta.tsv (transcriptome) ===')
canon=pd.read_csv(S('results/determinant_meta.tsv'),sep='\t')
canon=canon[canon.layer=='transcriptome'][['feature','pooled_beta','ci_halfwidth','I2']].set_index('feature')
m3=M[M['mode']=='residualized'][['feature','pooled_beta','ci_halfwidth','I2']].set_index('feature')
cmp=m3.join(canon,lsuffix='_mine',rsuffix='_canon')
cmp['dbeta']=(cmp.pooled_beta_mine-cmp.pooled_beta_canon).abs()
print(cmp[['pooled_beta_mine','pooled_beta_canon','dbeta']].round(8).to_string())
print('MAX |pooled_beta mode3 - canonical| = %.2e  (should be ~0)'%cmp['dbeta'].max())

print('\n=== pooled_beta by mode (all 13 features) ===')
piv=M.pivot_table(index='feature',columns='mode',values='pooled_beta',observed=False)[MODES]
print(piv.round(4).to_string())
print('\nmax |pooled_beta| per mode:')
for mode in MODES:
    s=M[M['mode']==mode]; f=s.loc[s.pooled_beta.abs().idxmax()]
    print('  %-12s max|beta|=%.4f  (%s)   any |beta|>=0.03 ? %s'%(
        mode,s.pooled_beta.abs().max(),f.feature,bool((s.pooled_beta.abs()>=0.03).any())))

print('\n=== TASK-REQUESTED features across modes ===')
for feature in ['causal_nonEGFR','arch_nsig']:
    print('  [%s]'%feature)
    for mode in MODES:
        r=M[(M['mode']==mode)&(M.feature==feature)].iloc[0]
        lo_,hi_=r.pooled_beta-r.ci_halfwidth,r.pooled_beta+r.ci_halfwidth
        print('    %-12s beta=%+.4f  CI[%+.4f,%+.4f]  half=%.4f  I2=%.1f  |beta|>=0.03? %s  excl0? %s'%(
            mode,r.pooled_beta,lo_,hi_,r.ci_halfwidth,r.I2,
            abs(r.pooled_beta)>=0.03, (lo_>0 or hi_<0)))

print('\nWROTE', S('results/determinant_modes.tsv'))

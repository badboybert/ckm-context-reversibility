#!/usr/bin/env python
"""determinant_sensitivity.py — TASK B8 robustness for the transcriptome determinant null.

STANDALONE. Copies the fit machinery from build_determinants.py verbatim; it does NOT import
from, mutate, or overwrite any canonical output (determinant_per_panel.tsv / determinant_meta.tsv
are read-only here — never written).

The transcriptome determinant is the POWERED object (k=9 panels, n~14-28k genes/panel); the null is
a PRACTICAL-EQUIVALENCE result (coefficients resolvable at n~150k but negligibly small, |beta|<~0.03
SD). These three sensitivities pressure-test that the null is not an artifact of median-imputation of
missing intrinsic annotations, nor of the not-testable 0-by-construction genes inflating the
causal-status contrast, and they tabulate the missingness floor.

Produces (NEW files only):
 (a) results/determinant_completecase.tsv   — determinant re-fit on genes with COMPLETE intrinsic
     annotations (loeuf_miss==0 AND tau_miss==0), with the *_miss indicator features DROPPED, so no
     loeuf/tau imputation enters the model. Per-panel betas + Hartung-Knapp meta, all features.
 (b) results/causal_nominatable_universe.tsv — causal-status contrast restricted to the
     cis-eQTL-instrumented, MR-TESTED (= nominatable) universe: nominated (causal_nonEGFR>=1) vs
     tested-but-not-nominated, EXCLUDING not-testable genes (which are 0 by construction in the
     canonical full-universe fit). causal_nonEGFR recomputed multivariate + univariate + raw contrast.
 (c) results/intrinsic_missingness.tsv       — per-feature missing-gene counts, per panel + pooled +
     unique-gene union, for the transcriptome determinant fit set.

Testability definition (documented): a gene is TESTABLE / nominatable iff it carries a cis-eQTL
instrument and was MR-tested, i.e. it is present in Paper B's transcriptome_local.tsv with has_eqtl==1
(all 1,432 rows of that file satisfy has_eqtl==1). Genes absent from that file could NOT have been
nominated and enter the canonical causal_nonEGFR contrast only as structural zeros; the nominatable
universe drops them.
"""
import os, numpy as np, pandas as pd
import statsmodels.api as sm
from scipy import stats

SIG=os.path.dirname(os.path.dirname(os.path.abspath(__file__))); PROJ=os.path.dirname(SIG)
def P(*a): return os.path.join(PROJ,*a)
def S(*a): return os.path.join(SIG,*a)
def rint(x):
    x=np.asarray(x,float); r=stats.rankdata(x); return stats.norm.ppf((r-0.5)/len(x))

# ---- external gene features (identical construction to build_determinants.py) ----
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
MISS=['loeuf_miss','tau_miss']
FEATS_FULL=EXTNUM+EXTBIN+ARCH+OTHER+MISS               # canonical model (with imputation + miss flags)
FEATS_CC  =EXTNUM+EXTBIN+ARCH+OTHER                    # complete-case model (no miss flags)

# nominatable / testable universe: cis-eQTL-instrumented, MR-tested genes
TESTABLE=set(tx_loc.loc[tx_loc.has_eqtl==1,'gene'].astype(str))

def gene_features_tx():
    loc=tx_loc.rename(columns={'gene':'mark_key'}).copy(); loc['gene']=loc['mark_key']
    loc['arch_nsig']=loc['eqtl_n_cs'].fillna(0); loc['arch_str']=loc['eqtl_max_pip'].fillna(0); loc['has_arch']=1
    return loc

def build_design(row, loc):
    """Return the per-gene design d (with residualized y + every feature), replicating fit_panel's
    transcriptome branch EXACTLY, plus explicit pre-imputation missingness flags for each feature."""
    df=pd.read_csv(P(row.path),sep='\t'); cols=list(df.columns)
    ec='rev_beta' if 'rev_beta' in cols else ('effect' if 'effect' in cols else 'delta_beta')
    d=pd.DataFrame({'mark_key':df['mark_id'].astype(str),'rb':pd.to_numeric(df[ec],errors='coerce')}).dropna().drop_duplicates('mark_key')
    d['za']=d['rb'].abs()
    # measurability residualization (transcriptome): within base_sd x abundance cells, rank-INT |rev_beta|
    b=bsd_tx[bsd_tx.panel==row.panel].set_index('mark_id')
    d['bsd']=d['mark_key'].map(b['base_sd']); d['ab']=d['mark_key'].map(b['base_mean']); d=d.dropna(subset=['bsd','ab'])
    d['sb']=pd.qcut(d['bsd'].rank(method='first'),20,labels=False,duplicates='drop')
    d['mb']=pd.qcut(d['ab'].rank(method='first'),5,labels=False,duplicates='drop')
    d['y']=d.groupby([d.sb,d.mb])['za'].transform(lambda s: rint(s.values) if len(s)>5 else np.nan); d=d.dropna(subset=['y'])
    # left-join local (architecture + causal tier) + external features onto ALL panel genes
    d=d.merge(loc,on='mark_key',how='left'); d['gene']=d['mark_key']
    d=d.merge(ext[['gene']+EXTNUM+EXTBIN],on='gene',how='left')
    # --- capture missingness BEFORE any imputation ---
    d['miss_loeuf']=d['loeuf'].isna().astype(int)
    d['miss_tau']=d['tau'].isna().astype(int)
    d['miss_secreted']=d['is_secreted'].isna().astype(int)
    d['miss_enzyme']=d['is_enzyme'].isna().astype(int)
    d['miss_membrane']=d['is_membrane'].isna().astype(int)
    d['miss_arch']=d['has_arch'].isna().astype(int)                 # gene absent from transcriptome_local
    d['in_nominatable']=d['mark_key'].isin(TESTABLE).astype(int)     # testable == cis-eQTL MR-tested
    # --- canonical imputation (identical to fit_panel) ---
    d['has_arch']=d['has_arch'].fillna(0); d['arch_nsig']=d['arch_nsig'].fillna(0); d['arch_str']=d['arch_str'].fillna(0)
    d['causal_nonEGFR']=d['causal_nonEGFR'].fillna(0)
    for f in EXTBIN: d[f]=d[f].fillna(0)
    d['loeuf_miss']=d['loeuf'].isna().astype(int); d['tau_miss']=d['tau'].isna().astype(int)
    for f in EXTNUM: d[f]=d[f].fillna(d[f].median())
    return d

def fit(d, feats):
    """OLS(HC3) of y on the given standardized features; drop zero-variance cols. Returns DataFrame."""
    X=d[[f for f in feats if f in d.columns]].apply(pd.to_numeric,errors='coerce').fillna(0)
    X=(X-X.mean())/X.std(ddof=0); X=X.loc[:,X.std()>0]; used=list(X.columns)
    res=sm.OLS(d['y'].values,sm.add_constant(X)).fit(cov_type='HC3')
    return res, used

def hk(b,s):
    """Hartung-Knapp random-effects meta (identical to build_determinants.py)."""
    b=np.asarray(b,float); v=np.asarray(s,float)**2; k=len(b)
    if k==1: return b[0],np.sqrt(v[0]),k,0.0,0.0
    w=1/v; sw=w.sum(); bf=(w*b).sum()/sw; Q=(w*(b-bf)**2).sum(); df=k-1
    C=sw-(w**2).sum()/sw; t2=max(0,(Q-df)/C) if C>0 else 0; ws=1/(v+t2); br=(ws*b).sum()/ws.sum()
    se=np.sqrt(max((ws*(b-br)**2).sum()/(df*ws.sum()),1e-12)); return br,se*stats.t.ppf(0.975,df),k,Q,t2

def meta_row(sub, layer, feature, tag):
    br,ci,k,Q,t2=hk(sub['beta'].values,sub['se'].values); df=k-1
    I2=max(0.0,(Q-df)/Q)*100 if Q>0 else 0.0
    return dict(scope='META_'+tag,layer=layer,feature=feature,beta=br,se=np.nan,p=np.nan,
                ci_halfwidth=ci,k=k,Q=Q,I2=I2,tau2=t2,n=int(sub['n'].median()),n_min=int(sub['n'].min()))

# =====================================================================================
loc=gene_features_tx()
rows=[r for r in MAN[MAN.layer=='transcriptome'].itertuples(index=False) if r.panel in set(bsd_tx.panel)]
designs={r.panel:build_design(r,loc) for r in rows}
print(f'transcriptome panels fitted: {len(designs)}  ({sorted(designs)})')
print(f'nominatable (testable) universe size: {len(TESTABLE)} genes; nominated (causal_nonEGFR>=1): '
      f'{int((tx_loc.causal_nonEGFR>=1).sum())}')

# ---- sanity: reproduce canonical full-universe causal_nonEGFR (baseline anchor) ----
canon=[]
for panel,d in designs.items():
    res,used=fit(d,FEATS_FULL)
    canon.append(dict(panel=panel,beta=res.params.get('causal_nonEGFR',np.nan),
                      se=res.bse.get('causal_nonEGFR',np.nan),n=len(d)))
canon=pd.DataFrame(canon); cb,cci,ck,cQ,ct2=hk(canon.beta.values,canon.se.values)
print(f'\n[ANCHOR] canonical full-universe causal_nonEGFR meta (this script) = {cb:+.5f} (CI half {cci:.5f}); '
      f'canonical determinant_meta.tsv = -0.00362 (CI half 0.00338).')

# =====================================================================================
# (a) COMPLETE-CASE determinant: genes complete on loeuf AND tau; *_miss features dropped.
# =====================================================================================
cc_pp=[]
for panel,d in designs.items():
    dcc=d[(d['loeuf_miss']==0)&(d['tau_miss']==0)].copy()
    res,used=fit(dcc,FEATS_CC)
    for f in used:
        cc_pp.append(dict(scope=panel,layer='transcriptome',feature=f,beta=res.params[f],se=res.bse[f],
                          p=res.pvalues[f],ci_halfwidth=np.nan,k=1,Q=np.nan,I2=np.nan,tau2=np.nan,
                          n=len(dcc),n_min=len(dcc)))
cc_pp=pd.DataFrame(cc_pp)
cc_meta=[meta_row(g,'transcriptome',f,'completecase') for f,g in cc_pp.groupby('feature') if len(g)>=2]
cc_out=pd.concat([cc_pp,pd.DataFrame(cc_meta)],ignore_index=True)
cc_out.to_csv(S('results/determinant_completecase.tsv'),sep='\t',index=False)
ccm=pd.DataFrame(cc_meta).set_index('feature')
print('\n=== (a) COMPLETE-CASE determinant meta (loeuf & tau observed; no imputation; k=9) ===')
print(f"  n/panel (complete-case): "+", ".join(f"{p}={cc_pp[cc_pp.scope==p].n.iloc[0]}" for p in designs))
for f in ['loeuf','tau','causal_nonEGFR','arch_nsig','arch_str','has_arch','is_secreted','is_enzyme','is_membrane','n_drug_log','n_gwas_log']:
    if f in ccm.index:
        r=ccm.loc[f]; print(f"  {f:14s} pooled_beta={r.beta:+.5f}  CI_half={r.ci_halfwidth:.5f}  I2={r.I2:4.1f}  {'CI-excl-0' if r.ci_halfwidth<abs(r.beta) else 'incl-0'}")

# =====================================================================================
# (b) NOMINATABLE-UNIVERSE causal-status: nominated vs tested-but-not-nominated only.
# =====================================================================================
nom_pp=[]
for panel,d in designs.items():
    du=d[d['in_nominatable']==1].copy()
    nn=int((du['causal_nonEGFR']>=1).sum()); ntot=len(du)
    # multivariate within the tested universe (same canonical feature spec, degenerate cols dropped)
    res_mv,used_mv=fit(du,FEATS_FULL)
    bmv=res_mv.params.get('causal_nonEGFR',np.nan); smv=res_mv.bse.get('causal_nonEGFR',np.nan); pmv=res_mv.pvalues.get('causal_nonEGFR',np.nan)
    # univariate contrast within the tested universe
    Xc=du[['causal_nonEGFR']].astype(float); Xc=(Xc-Xc.mean())/Xc.std(ddof=0)
    if Xc['causal_nonEGFR'].std()>0:
        ru=sm.OLS(du['y'].values,sm.add_constant(Xc)).fit(cov_type='HC3')
        buni=ru.params['causal_nonEGFR']; suni=ru.bse['causal_nonEGFR']; puni=ru.pvalues['causal_nonEGFR']
    else:
        buni=suni=puni=np.nan
    my1=du.loc[du['causal_nonEGFR']>=1,'y'].mean(); my0=du.loc[du['causal_nonEGFR']<1,'y'].mean()
    nom_pp.append(dict(scope=panel,layer='transcriptome',n_universe=ntot,n_nominated=nn,
                       beta_mv=bmv,se_mv=smv,p_mv=pmv,beta_uni=buni,se_uni=suni,p_uni=puni,
                       mean_y_nominated=my1,mean_y_tested_notnom=my0,mean_y_diff=my1-my0))
nom_pp=pd.DataFrame(nom_pp)
# meta over per-panel causal betas (multivariate + univariate)
mvm=hk(nom_pp.beta_mv.values,nom_pp.se_mv.values)
unm=hk(nom_pp.dropna(subset=['se_uni']).beta_uni.values,nom_pp.dropna(subset=['se_uni']).se_uni.values)
def _i2(Q,k): df=k-1; return max(0.0,(Q-df)/Q)*100 if Q>0 else 0.0
nom_meta=pd.DataFrame([
    dict(scope='META_nominatable_mv',layer='transcriptome',n_universe=int(nom_pp.n_universe.median()),
         n_nominated=int(nom_pp.n_nominated.median()),beta_mv=mvm[0],se_mv=np.nan,p_mv=np.nan,
         beta_uni=np.nan,se_uni=np.nan,p_uni=np.nan,ci_halfwidth=mvm[1],k=mvm[2],Q=mvm[3],
         I2=_i2(mvm[3],mvm[2]),tau2=mvm[4],mean_y_nominated=np.nan,mean_y_tested_notnom=np.nan,mean_y_diff=np.nan),
    dict(scope='META_nominatable_uni',layer='transcriptome',n_universe=int(nom_pp.n_universe.median()),
         n_nominated=int(nom_pp.n_nominated.median()),beta_mv=np.nan,se_mv=np.nan,p_mv=np.nan,
         beta_uni=unm[0],se_uni=np.nan,p_uni=np.nan,ci_halfwidth=unm[1],k=unm[2],Q=unm[3],
         I2=_i2(unm[3],unm[2]),tau2=unm[4],mean_y_nominated=np.nan,mean_y_tested_notnom=np.nan,mean_y_diff=np.nan),
])
nom_out=pd.concat([nom_pp,nom_meta],ignore_index=True)
nom_out.to_csv(S('results/causal_nominatable_universe.tsv'),sep='\t',index=False)
print('\n=== (b) NOMINATABLE-UNIVERSE causal_nonEGFR (nominated vs tested-but-not-nominated; not-testable EXCLUDED) ===')
print(f"  universe n/panel ~{int(nom_pp.n_universe.median())}, nominated n/panel median={int(nom_pp.n_nominated.median())} "
      f"(vs canonical full-universe n~14-28k with the same {int((tx_loc.causal_nonEGFR>=1).sum())} nominated)")
print(f"  META multivariate  pooled_beta={mvm[0]:+.5f}  CI_half={mvm[1]:.5f}  I2={_i2(mvm[3],mvm[2]):.1f}  ({'CI-excl-0' if mvm[1]<abs(mvm[0]) else 'incl-0'})")
print(f"  META univariate    pooled_beta={unm[0]:+.5f}  CI_half={unm[1]:.5f}  I2={_i2(unm[3],unm[2]):.1f}  ({'CI-excl-0' if unm[1]<abs(unm[0]) else 'incl-0'})")
print(f"  canonical full-universe (for reference) = -0.00362")

# =====================================================================================
# (c) PER-FEATURE MISSINGNESS: per panel + pooled-sum + unique-gene union.
# =====================================================================================
MISS_MAP={'loeuf':'miss_loeuf','tau':'miss_tau','architecture':'miss_arch',
          'is_secreted':'miss_secreted','is_enzyme':'miss_enzyme','is_membrane':'miss_membrane',
          'causal_status':'miss_arch'}  # causal status testable iff gene in transcriptome_local (== arch present)
mrows=[]
for panel,d in designs.items():
    n=len(d)
    for feat,col in MISS_MAP.items():
        nm=int(d[col].sum())
        mrows.append(dict(scope=panel,feature=feat,n_genes=n,n_missing=nm,pct_missing=round(100*nm/n,2)))
# pooled sum across panels (gene-instances)
mm=pd.DataFrame(mrows)
for feat in MISS_MAP:
    sub=mm[mm.feature==feat]
    mrows.append(dict(scope='POOLED_SUM',feature=feat,n_genes=int(sub.n_genes.sum()),
                      n_missing=int(sub.n_missing.sum()),pct_missing=round(100*sub.n_missing.sum()/sub.n_genes.sum(),2)))
# unique-gene union across all transcriptome fit sets (feature missingness is gene-level)
uni=pd.concat([d[['mark_key']+list(dict.fromkeys(MISS_MAP.values()))] for d in designs.values()],ignore_index=True)
uni=uni.groupby('mark_key',as_index=False).max()  # a gene counts as missing a feature if missing wherever it appears
nu=len(uni)
for feat,col in MISS_MAP.items():
    nm=int(uni[col].sum())
    mrows.append(dict(scope='UNIQUE_UNION',feature=feat,n_genes=nu,n_missing=nm,pct_missing=round(100*nm/nu,2)))
miss_out=pd.DataFrame(mrows)
miss_out.to_csv(S('results/intrinsic_missingness.tsv'),sep='\t',index=False)
print('\n=== (c) PER-FEATURE MISSINGNESS (transcriptome determinant fit set) ===')
piv=miss_out[miss_out.scope=='UNIQUE_UNION'].set_index('feature')
print(f"  unique genes across 9 panels = {nu}")
for feat in MISS_MAP:
    r=piv.loc[feat]; print(f"  {feat:14s} missing {int(r.n_missing):6d}/{nu} = {r.pct_missing:5.1f}%")

print('\nWROTE:')
for f in ['determinant_completecase.tsv','causal_nominatable_universe.tsv','intrinsic_missingness.tsv']:
    print('  results/'+f)

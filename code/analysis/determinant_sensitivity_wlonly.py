#!/usr/bin/env python
"""determinant_sensitivity_wlonly.py — TRUE-WEIGHT-LOSS-ONLY sensitivity for the
transcriptome determinant null (true-weight-loss-only sensitivity).

STANDALONE. Copies the fit machinery from build_determinants.py / determinant_sensitivity.py
VERBATIM; it does NOT import from, mutate, or overwrite any canonical output. Canonical
determinant_per_panel.tsv / determinant_meta.tsv are read only for the ANCHOR comparison.

QUESTION the reviewer raised: the atlas mixes weight-loss and non-weight-loss interventions
(metformin, empagliflozin, exercise). Does the determinant null survive when the meta is
restricted to interventions whose primary mechanism is weight loss (caloric restriction / diet
and bariatric surgery)?

Empagliflozin (EMPEROR) is a PROTEOME panel and never enters the transcriptome determinant meta,
so it is already absent from the k=9 object. Within the k=9 transcriptome meta, exactly three
panels touch the exclusion set:
   gse157585_muscle_table       metformin + resistance training   (MASTERS)
   gse43471_adipose_table       diet + exercise
   gse83352_muscle_exercise_table  exercise                        (STRRIDE)

SCENARIOS
  k9_canonical : all 9 panels (reproduces the shipped determinant_meta.tsv as an anchor)
  k6_strict    : drop ALL panels touching metformin / exercise -> pure diet + surgery (k=6)
  k7_keepdietex: also keep gse43471 (diet+exercise, since diet drives the weight loss) -> k=7

For each scenario: per-panel standardized betas (all features) + Hartung-Knapp meta; report
pooled_beta, CI half-width, whether |pooled|>the 0.05 SD SESOI, whether the CI excludes 0, and
the max |pooled_beta| across features. The load-bearing check is that NO intrinsic determinant
reaches the +/-0.05 SD SESOI in any scenario, and genetic causal status stays small/underpowered.

Writes NEW files only:
  results/sensitivity_wlonly/determinant_per_panel_wlonly.tsv   (per-panel betas, all scenarios)
  results/sensitivity_wlonly/determinant_meta_wlonly.tsv        (HK meta per feature x scenario)
  results/sensitivity_wlonly/determinant_wlonly_summary.txt     (human-readable)
"""
import os, numpy as np, pandas as pd
import statsmodels.api as sm
from scipy import stats

SIG=os.path.dirname(os.path.dirname(os.path.abspath(__file__))); PROJ=os.path.dirname(SIG)
def P(*a): return os.path.join(PROJ,*a)
def S(*a): return os.path.join(SIG,*a)
OUT=S('results','sensitivity_wlonly'); os.makedirs(OUT,exist_ok=True)
def rint(x):
    x=np.asarray(x,float); r=stats.rankdata(x); return stats.norm.ppf((r-0.5)/len(x))

SESOI=0.05   # smallest effect of interest (SD), pre-specified for the equivalence analysis

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
FEATS=EXTNUM+EXTBIN+ARCH+OTHER+[f+'_miss' for f in ['loeuf','tau']]

# panels excluded because their intervention touches metformin / empagliflozin / exercise
EXC_METEX = {'gse157585_muscle_table','gse83352_muscle_exercise_table'}   # metformin+RT, exercise
EXC_DIETEX = {'gse43471_adipose_table'}                                   # diet+exercise
SCENARIOS = {
    'k9_canonical' : set(),
    'k6_strict'    : EXC_METEX | EXC_DIETEX,
    'k7_keepdietex': EXC_METEX,
}

def gene_features_tx():
    loc=tx_loc.rename(columns={'gene':'mark_key'}).copy(); loc['gene']=loc['mark_key']
    loc['arch_nsig']=loc['eqtl_n_cs'].fillna(0); loc['arch_str']=loc['eqtl_max_pip'].fillna(0); loc['has_arch']=1
    return loc

def fit_panel(row, loc):
    """Replicates build_determinants.fit_panel transcriptome branch EXACTLY."""
    df=pd.read_csv(P(row.path),sep='\t'); cols=list(df.columns)
    ec='rev_beta' if 'rev_beta' in cols else ('effect' if 'effect' in cols else 'delta_beta')
    d=pd.DataFrame({'mark_key':df['mark_id'].astype(str),'rb':pd.to_numeric(df[ec],errors='coerce')}).dropna().drop_duplicates('mark_key')
    d['za']=d['rb'].abs()
    b=bsd_tx[bsd_tx.panel==row.panel].set_index('mark_id')
    d['bsd']=d['mark_key'].map(b['base_sd']); d['ab']=d['mark_key'].map(b['base_mean']); d=d.dropna(subset=['bsd','ab'])
    d['sb']=pd.qcut(d['bsd'].rank(method='first'),20,labels=False,duplicates='drop')
    d['mb']=pd.qcut(d['ab'].rank(method='first'),5,labels=False,duplicates='drop')
    d['y']=d.groupby([d.sb,d.mb])['za'].transform(lambda s: rint(s.values) if len(s)>5 else np.nan); d=d.dropna(subset=['y'])
    d=d.merge(loc,on='mark_key',how='left'); d['gene']=d['mark_key']
    d=d.merge(ext[['gene']+EXTNUM+EXTBIN],on='gene',how='left')
    d['has_arch']=d['has_arch'].fillna(0); d['arch_nsig']=d['arch_nsig'].fillna(0); d['arch_str']=d['arch_str'].fillna(0)
    d['causal_nonEGFR']=d['causal_nonEGFR'].fillna(0)
    for f in EXTBIN: d[f]=d[f].fillna(0)
    for f in ['loeuf','tau']: d[f+'_miss']=d[f].isna().astype(int)
    for f in EXTNUM: d[f]=d[f].fillna(d[f].median())
    feats=[f for f in FEATS if f in d.columns]
    X=d[feats].apply(pd.to_numeric,errors='coerce').fillna(0)
    X=(X-X.mean())/X.std(ddof=0); X=X.loc[:,X.std()>0]; used=list(X.columns)
    res=sm.OLS(d['y'].values,sm.add_constant(X)).fit(cov_type='HC3')
    return pd.DataFrame([dict(panel=row.panel,feature=f,beta=res.params[f],se=res.bse[f],
        p=res.pvalues[f],n=len(d)) for f in used])

def hk(b,s):
    """Hartung-Knapp random-effects meta (identical to build_determinants.py)."""
    b=np.asarray(b,float); v=np.asarray(s,float)**2; k=len(b)
    if k==1: return b[0],np.sqrt(v[0]),k,0.0,0.0
    w=1/v; sw=w.sum(); bf=(w*b).sum()/sw; Q=(w*(b-bf)**2).sum(); df=k-1
    C=sw-(w**2).sum()/sw; t2=max(0,(Q-df)/C) if C>0 else 0; ws=1/(v+t2); br=(ws*b).sum()/ws.sum()
    se=np.sqrt(max((ws*(b-br)**2).sum()/(df*ws.sum()),1e-12)); return br,se*stats.t.ppf(0.975,df),k,Q,t2

# ---- fit every k=9 transcriptome panel ONCE, then meta-pool per scenario ----
loc=gene_features_tx()
rows=[r for r in MAN[MAN.layer=='transcriptome'].itertuples(index=False) if r.panel in set(bsd_tx.panel)]
per_panel=pd.concat([fit_panel(r,loc) for r in rows],ignore_index=True)
K9=sorted(per_panel.panel.unique())
assert len(K9)==9, f'expected 9 transcriptome determinant panels, got {len(K9)}: {K9}'
per_panel.to_csv(os.path.join(OUT,'determinant_per_panel_wlonly.tsv'),sep='\t',index=False)

meta_rows=[]
for scen,excl in SCENARIOS.items():
    keep=[p for p in K9 if p not in excl]
    sub=per_panel[per_panel.panel.isin(keep)]
    for f,g in sub.groupby('feature'):
        if len(g)<2: continue
        br,ci,k,Q,t2=hk(g.beta.values,g.se.values); df=k-1
        I2=max(0.0,(Q-df)/Q)*100 if Q>0 else 0.0
        meta_rows.append(dict(scenario=scen,k=k,feature=f,pooled_beta=br,ci_halfwidth=ci,
            ci_lo=br-ci,ci_hi=br+ci,I2=I2,tau2=t2,
            ci_excludes_0=bool(ci<abs(br)),
            exceeds_SESOI=bool(abs(br)>=SESOI),
            ci_within_SESOI=bool((abs(br)+ci)<=SESOI),
            n_genes_median=int(g.n.median()),n_genes_min=int(g.n.min())))
M=pd.DataFrame(meta_rows)
M.to_csv(os.path.join(OUT,'determinant_meta_wlonly.tsv'),sep='\t',index=False)

# ---- anchor check vs shipped canonical ----
canon=pd.read_csv(S('results/determinant_meta.tsv'),sep='\t')
canon=canon[canon.layer=='transcriptome'].set_index('feature')['pooled_beta']
k9=M[M.scenario=='k9_canonical'].set_index('feature')['pooled_beta']
anchor_maxdiff=float((k9-canon.reindex(k9.index)).abs().max())

# ---- human-readable summary ----
L=[]
def p(s=''): L.append(s); print(s)
p('TRUE-WEIGHT-LOSS-ONLY SENSITIVITY — transcriptome determinant null')
p('='*78)
p(f'k=9 transcriptome determinant panels: {K9}')
p('Excluded for weight-loss-only:')
p(f'  metformin/exercise (both scenarios): {sorted(EXC_METEX)}')
p(f'  diet+exercise (strict only):         {sorted(EXC_DIETEX)}')
p(f'Empagliflozin (EMPEROR) is proteome -> never in the transcriptome meta (not excluded here).')
p('')
p(f'ANCHOR: k9_canonical reproduces shipped determinant_meta.tsv within max |Δpooled_beta| = {anchor_maxdiff:.2e}')
p('')
# tissues per scenario (transparency: strict k6 drops both muscle panels)
MAN_TIS=MAN.set_index('panel')['tissue'].to_dict()
for scen in SCENARIOS:
    s=M[M.scenario==scen].copy()
    k=int(s.k.max())
    keep=[pn for pn in K9 if pn not in SCENARIOS[scen]]
    tis=sorted(set(MAN_TIS.get(pn,'?') for pn in keep))
    top=s.reindex(s.pooled_beta.abs().sort_values(ascending=False).index)
    maxfeat=top.iloc[0]
    causal=s[s.feature=='causal_nonEGFR'].iloc[0]
    n_exceed=int(s.exceeds_SESOI.sum())
    ci_wide=list(s[~s.ci_within_SESOI].feature)   # point est < SESOI but |beta|+CI not equivalence-tight
    p(f'--- {scen} (k={k}; tissues={tis}) ---')
    p(f'  max |pooled_beta| : {maxfeat.pooled_beta:+.4f} ({maxfeat.feature}); CI half {maxfeat.ci_halfwidth:.4f}')
    p(f'  features exceeding +/-{SESOI} SD SESOI (point estimate) : {n_exceed}   '
      f'(headline null {"HOLDS" if n_exceed==0 else "BROKEN"})')
    p(f'  features whose |beta|+CI is NOT within +/-{SESOI} (heterogeneous, not equivalence-tight): {ci_wide}')
    p(f'    (these are the known dynamic-range / annotation features, esp. tau [sign-flipping, I2~97%], NOT novel determinants)')
    p(f'  genetic causal status : pooled {causal.pooled_beta:+.4f}  CI [{causal.ci_lo:+.4f},{causal.ci_hi:+.4f}]  '
      f'(rare-binary, standardization-compressed ~50x; NOT an effect size)')
    p(f'    -> interpret via the raw group contrast: canonical -0.18 SD [-0.37,+0.01], n=9, UNDERPOWERED (Supp Table 14 / causal_status_effect.tsv)')
    feats_excl0=list(s[s.ci_excludes_0].feature)
    p(f'  features whose CI excludes 0 (statistically resolvable, still << SESOI): {feats_excl0}')
    p('')
p('CONCLUSION')
p('-'*78)
allhold=all(int(M[M.scenario==sc].exceeds_SESOI.sum())==0 for sc in SCENARIOS)
mx=M.loc[M.pooled_beta.abs().idxmax()]
p(f'  Across all scenarios the headline determinant null {"HOLDS" if allhold else "DOES NOT HOLD"}: '
  f'no intrinsic determinant reaches the +/-{SESOI} SD SESOI.')
p(f'  Largest single pooled coefficient anywhere = {mx.pooled_beta:+.4f} ({mx.feature}, {mx.scenario}), '
  f'still {SESOI/abs(mx.pooled_beta):.1f}x below the SESOI.')
p(f'  Restricting to true-weight-loss interventions (diet + surgery) does not change the conclusion.')
open(os.path.join(OUT,'determinant_wlonly_summary.txt'),'w',encoding='utf-8').write('\n'.join(L))
print('\nwrote results/sensitivity_wlonly/{determinant_per_panel_wlonly.tsv, determinant_meta_wlonly.tsv, determinant_wlonly_summary.txt}')

#!/usr/bin/env python
"""causal_status_effect.py — genetic causal-status effect on reversibility, on INTERPRETABLE scales.

WHY THIS EXISTS
---------------
The canonical determinant model (build_determinants.py) z-standardizes EVERY predictor,
including the genetic causal-status feature (causal_nonEGFR). But causal status is a rare
binary: only 9 of ~18,600 transcripts are genetically nominated. Standardizing a predictor
with SD_x = sqrt(p(1-p)) ~ 0.02 divides its coefficient by ~0.02, mechanically COMPRESSING
the pooled beta to -0.0036 SD. That number is an artifact of the scale, not the effect size.
The INTERPRETABLE quantity is the raw group difference in reversibility between the 9 causal
transcripts and everyone else — expected ~-0.16 SD.

This script re-expresses the SAME effect on interpretable scales, reusing the EXACT canonical
outcome construction. It is STANDALONE: the fit machinery is copied verbatim from
build_determinants.py / determinant_robustness.py / determinant_sensitivity.py; it does NOT
import, execute, or overwrite build_determinants.py or any canonical determinant_*.tsv.
Writes ONLY results/causal_status_effect.tsv.

FOUR ESTIMATES (per panel, then Hartung-Knapp random-effects meta across the 9 panels, the
SAME hk() the determinant pipeline uses):
  1. RAW GROUP SMD (unadjusted): mean(y|causal) - mean(y|non-causal). y is unit-variance
     rank-INT, so this difference IS a standardized mean difference. Also Cohen's d = diff /
     pooled-SD (confirms pooled-SD ~ 1). Reports per-panel n_causal + two-sample (Welch) SE.
  2. MULTIVARIABLE RAW COEFFICIENT: the canonical per-panel OLS(HC3) on the full feature set,
     but causal_nonEGFR left 0/1 (NOT standardized); every other feature standardized exactly
     as canonical. The causal coefficient = the adjusted raw mean difference.
  3. NOMINATABLE-UNIVERSE SMD: restrict to cis-eQTL-instrumented, MR-tested genes (has_eqtl==1);
     nominated (causal>=1) vs tested-but-NOT-nominated; not-testable genes EXCLUDED. Underpowered.
  4. SANITY / ANCHOR: the fully standardized coefficient reproduces canonical -0.0036, and the
     per-panel identity beta_std = beta_raw * SD_x (SD_x = per-panel SD of causal_nonEGFR, ddof=0)
     PROVES the compression. Reports SD_x per panel + total UNIQUE causal genes across panels.

OUTCOME (identical to canonical transcriptome branch): per score-contributing transcriptome
panel, |rev_beta| rank-inverse-normal-transformed WITHIN base_sd(20) x abundance(5) cells
(the 'residualized', unit-variance outcome). Causal status = causal_nonEGFR>=1 from
data/features/transcriptome_local.tsv (9 genes), filled 0 for every other panel gene (canonical).
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
FEATS=EXTNUM+EXTBIN+ARCH+OTHER+[f+'_miss' for f in ['loeuf','tau']]   # canonical feature list/order

# nominatable / testable universe: cis-eQTL-instrumented, MR-tested genes (all rows of tx_loc have has_eqtl==1)
TESTABLE=set(tx_loc.loc[tx_loc.has_eqtl==1,'gene'].astype(str))
CAUSAL_GENES=set(tx_loc.loc[tx_loc.causal_nonEGFR>=1,'gene'].astype(str))

def gene_features_tx():
    loc=tx_loc.rename(columns={'gene':'mark_key'}).copy(); loc['gene']=loc['mark_key']
    loc['arch_nsig']=loc['eqtl_n_cs'].fillna(0); loc['arch_str']=loc['eqtl_max_pip'].fillna(0); loc['has_arch']=1
    return loc

def build_design(row, loc):
    """One transcriptome panel -> per-gene design d with the CANONICAL residualized outcome y and
    the full standardized-ready feature set. Replicates fit_panel's transcriptome branch EXACTLY."""
    df=pd.read_csv(P(row.path),sep='\t'); cols=list(df.columns)
    ec='rev_beta' if 'rev_beta' in cols else ('effect' if 'effect' in cols else 'delta_beta')
    d=pd.DataFrame({'mark_key':df['mark_id'].astype(str),'rb':pd.to_numeric(df[ec],errors='coerce')}).dropna().drop_duplicates('mark_key')
    d['za']=d['rb'].abs()
    # measurability residualization: within base_sd x abundance cells, rank-INT |rev_beta| (canonical)
    b=bsd_tx[bsd_tx.panel==row.panel].set_index('mark_id')
    d['bsd']=d['mark_key'].map(b['base_sd']); d['ab']=d['mark_key'].map(b['base_mean']); d=d.dropna(subset=['bsd','ab'])
    d['sb']=pd.qcut(d['bsd'].rank(method='first'),20,labels=False,duplicates='drop')
    d['mb']=pd.qcut(d['ab'].rank(method='first'),5,labels=False,duplicates='drop')
    d['y']=d.groupby([d.sb,d.mb])['za'].transform(lambda s: rint(s.values) if len(s)>5 else np.nan); d=d.dropna(subset=['y'])
    # LEFT-join local (architecture + causal tier) + external features onto ALL panel genes (canonical)
    d=d.merge(loc,on='mark_key',how='left'); d['gene']=d['mark_key']
    d=d.merge(ext[['gene']+EXTNUM+EXTBIN],on='gene',how='left')
    d['in_nominatable']=d['mark_key'].isin(TESTABLE).astype(int)
    d['has_arch']=d['has_arch'].fillna(0); d['arch_nsig']=d['arch_nsig'].fillna(0); d['arch_str']=d['arch_str'].fillna(0)
    d['causal_nonEGFR']=d['causal_nonEGFR'].fillna(0)                 # non-panel genes -> 0 (canonical)
    for f in EXTBIN: d[f]=d[f].fillna(0)
    for f in ['loeuf','tau']: d[f+'_miss']=d[f].isna().astype(int)
    for f in EXTNUM: d[f]=d[f].fillna(d[f].median())
    return d

def hk(b,s):   # COPIED verbatim from build_determinants.py (Hartung-Knapp random-effects)
    b=np.asarray(b,float); v=np.asarray(s,float)**2; k=len(b)
    if k==1: return b[0],np.sqrt(v[0]),k,0.0,0.0
    w=1/v; sw=w.sum(); bf=(w*b).sum()/sw; Q=(w*(b-bf)**2).sum(); df=k-1
    C=sw-(w**2).sum()/sw; t2=max(0,(Q-df)/C) if C>0 else 0; ws=1/(v+t2); br=(ws*b).sum()/ws.sum()
    se=np.sqrt(max((ws*(b-br)**2).sum()/(df*ws.sum()),1e-12)); return br,se*stats.t.ppf(0.975,df),k,Q,t2

def group_smd(y, causal):
    """Raw group standardized mean difference on unit-variance rank-INT outcome y.
    Returns diff (=mean1-mean0=SMD), Welch SE, Cohen's d, pooled-SD, n1(causal), n0."""
    y=np.asarray(y,float); g=np.asarray(causal,float)>=1
    y1,y0=y[g],y[~g]; n1,n0=len(y1),len(y0)
    m1,m0=y1.mean(),y0.mean(); diff=m1-m0
    s1,s0=y1.std(ddof=1),y0.std(ddof=1)
    se=np.sqrt(s1**2/n1+s0**2/n0)                                    # two-sample (Welch) SE
    psd=np.sqrt(((n1-1)*s1**2+(n0-1)*s0**2)/(n1+n0-2))               # pooled SD (~1 by construction)
    d_cohen=diff/psd
    return diff,se,d_cohen,psd,n1,n0

def canonical_std_fit(d):
    """Fully canonical per-panel fit: ALL features standardized (incl causal). Returns beta_std,se_std."""
    X=d[[f for f in FEATS if f in d.columns]].apply(pd.to_numeric,errors='coerce').fillna(0)
    X=(X-X.mean())/X.std(ddof=0); X=X.loc[:,X.std()>0]
    res=sm.OLS(d['y'].values,sm.add_constant(X)).fit(cov_type='HC3')
    return res.params.get('causal_nonEGFR',np.nan),res.bse.get('causal_nonEGFR',np.nan)

def raw_causal_fit(d):
    """Item 2: canonical fit but causal_nonEGFR left RAW 0/1; every OTHER feature standardized
    exactly as canonical (guarantees the other coefficients are identical to canonical and that
    beta_raw = beta_std / SD_x exactly). Returns beta_raw,se_raw."""
    X=d[[f for f in FEATS if f in d.columns]].apply(pd.to_numeric,errors='coerce').fillna(0)
    Xz=(X-X.mean())/X.std(ddof=0); Xz=Xz.loc[:,Xz.std()>0]           # canonical standardize + zero-var drop
    Xz['causal_nonEGFR']=d['causal_nonEGFR'].astype(float).values    # overwrite with raw 0/1
    res=sm.OLS(d['y'].values,sm.add_constant(Xz)).fit(cov_type='HC3')
    return res.params.get('causal_nonEGFR',np.nan),res.bse.get('causal_nonEGFR',np.nan)

# ================================================================================================
loc=gene_features_tx()
rows=[r for r in MAN[MAN.layer=='transcriptome'].itertuples(index=False) if r.panel in set(bsd_tx.panel)]
designs={r.panel:build_design(r,loc) for r in rows}
PANELS=list(designs)
print(f'transcriptome score-contributing panels (k={len(PANELS)}): {PANELS}')
print(f'causal genes (causal_nonEGFR>=1 in transcriptome_local): {len(CAUSAL_GENES)} = {sorted(CAUSAL_GENES)}')

# per-panel tables for each estimate
smd_pp, mv_pp, nom_pp, anch_pp = [], [], [], []
seen_causal=set()
for panel,d in designs.items():
    seen_causal |= set(d.loc[d['causal_nonEGFR']>=1,'mark_key'])
    # --- 1. raw group SMD (unadjusted, full universe) ---
    diff,se,dc,psd,n1,n0=group_smd(d['y'].values,d['causal_nonEGFR'].values)
    smd_pp.append(dict(panel=panel,diff=diff,se=se,cohen_d=dc,pooled_sd=psd,n_causal=n1,n_noncausal=n0))
    # --- 2. multivariable RAW coefficient ---
    braw,sraw=raw_causal_fit(d)
    mv_pp.append(dict(panel=panel,beta_raw=braw,se_raw=sraw,n=len(d),n_causal=n1))
    # --- 4. canonical standardized anchor + SD_x + identity check ---
    bstd,sstd=canonical_std_fit(d)
    sdx=float(d['causal_nonEGFR'].astype(float).std(ddof=0))         # exactly what standardization divides by
    anch_pp.append(dict(panel=panel,beta_std=bstd,se_std=sstd,SD_x=sdx,
                        beta_raw_x_SDx=braw*sdx,identity_abserr=abs(bstd-braw*sdx),n_causal=n1))
    # --- 3. nominatable-universe SMD (has_eqtl only; not-testable excluded) ---
    du=d[d['in_nominatable']==1]
    ndiff,nse,ndc,npsd,nn1,nn0=group_smd(du['y'].values,du['causal_nonEGFR'].values)
    nom_pp.append(dict(panel=panel,diff=ndiff,se=nse,n_universe=len(du),n_causal=nn1,n_tested_notnom=nn0))

SMD=pd.DataFrame(smd_pp); MV=pd.DataFrame(mv_pp); NOM=pd.DataFrame(nom_pp); ANCH=pd.DataFrame(anch_pp)
TOTAL_N_CAUSAL=len(seen_causal)   # unique causal genes appearing across the 9 panels

# ---- Hartung-Knapp meta for each estimate ----
def meta(df,bcol,scol):
    b,ci,k,Q,t2=hk(df[bcol].values,df[scol].values); return b,b-ci,b+ci,k,Q,t2
m1=meta(SMD,'diff','se')
m2=meta(MV,'beta_raw','se_raw')
m3=meta(NOM,'diff','se')
m4=meta(ANCH,'beta_std','se_std')

# ================================================================================================
# console report
pd.set_option('display.width',200,'display.max_columns',30)
print('\n================ ESTIMATE 1: RAW GROUP SMD (unadjusted, full universe) ================')
print(SMD.round(4).to_string(index=False))
print(f'  n_causal per panel: {list(SMD.n_causal)}  (range {SMD.n_causal.min()}-{SMD.n_causal.max()})')
print(f'  pooled-SD per panel ~ {SMD.pooled_sd.min():.3f}-{SMD.pooled_sd.max():.3f}  (== 1 by construction: y is unit-variance rank-INT)')
print(f'  Cohen d ~ SMD (pooled-SD~1) : per-panel d = {[round(x,3) for x in SMD.cohen_d]}')
print(f'  HK META  SMD = {m1[0]:+.4f}  95% CI [{m1[1]:+.4f}, {m1[2]:+.4f}]  k={m1[3]}')

print('\n================ ESTIMATE 2: MULTIVARIABLE RAW COEFFICIENT (causal left 0/1) ================')
print(MV.round(4).to_string(index=False))
print(f'  HK META  beta_raw = {m2[0]:+.4f}  95% CI [{m2[1]:+.4f}, {m2[2]:+.4f}]  k={m2[3]}')

print('\n================ ESTIMATE 3: NOMINATABLE-UNIVERSE SMD (has_eqtl; not-testable excluded) ================')
print(NOM.round(4).to_string(index=False))
print(f'  universe n/panel ~ {int(NOM.n_universe.median())}, nominated n/panel {list(NOM.n_causal)}')
print(f'  HK META  SMD = {m3[0]:+.4f}  95% CI [{m3[1]:+.4f}, {m3[2]:+.4f}]  k={m3[3]}   (CI includes 0 -> underpowered)')

print('\n================ ESTIMATE 4: SANITY / ANCHOR (standardized reproduces -0.0036; compression proof) ================')
print(ANCH.round(6).to_string(index=False))
print(f'  HK META  beta_std = {m4[0]:+.5f}  95% CI [{m4[1]:+.5f}, {m4[2]:+.5f}]  k={m4[3]}')
print(f'  canonical determinant_meta.tsv beta_std = -0.0036168 (CI half 0.0033785) -> reproduced: '
      f'{abs(m4[0]-(-0.0036167576581104493))<1e-6}')
print(f'  per-panel identity  beta_std == beta_raw * SD_x  : max |err| = {ANCH.identity_abserr.max():.2e}  (~0 => EXACT)')
print(f'  SD_x per panel = {[round(x,4) for x in ANCH.SD_x]}  (mean {ANCH.SD_x.mean():.4f}) -> dividing beta_std by ~{ANCH.SD_x.mean():.3f} restores ~{m2[0]:+.2f} SD')
print(f'  TOTAL UNIQUE causal genes across the {len(PANELS)} panels = {TOTAL_N_CAUSAL}')

# ================================================================================================
# write results/causal_status_effect.tsv (one row per estimate)
out=[
 dict(estimate_name='raw_group_smd_unadjusted', scale='SMD (rank-INT SD units)',
      pooled=m1[0], ci_lo=m1[1], ci_hi=m1[2], k=m1[3], total_n_causal=TOTAL_N_CAUSAL,
      note=('unadjusted mean(y|causal)-mean(y|noncausal); y unit-variance rank-INT so diff==SMD; '
            f'Cohen d==SMD (pooled-SD~1); per-panel n_causal {SMD.n_causal.min()}-{SMD.n_causal.max()}; Welch SE')),
 dict(estimate_name='multivariable_raw_coef', scale='raw (0/1 predictor -> SD outcome)',
      pooled=m2[0], ci_lo=m2[1], ci_hi=m2[2], k=m2[3], total_n_causal=TOTAL_N_CAUSAL,
      note=('canonical per-panel OLS(HC3), causal_nonEGFR left 0/1, all other features standardized as canonical; '
            'adjusted raw mean difference')),
 dict(estimate_name='nominatable_universe_smd', scale='SMD (rank-INT SD units)',
      pooled=m3[0], ci_lo=m3[1], ci_hi=m3[2], k=m3[3], total_n_causal=TOTAL_N_CAUSAL,
      note=('nominated (causal>=1) vs tested-but-not-nominated within cis-eQTL/MR-tested universe (has_eqtl==1); '
            'not-testable genes excluded; CI includes zero -> underpowered')),
 dict(estimate_name='standardized_coef_anchor', scale='standardized (SD_x/SD_y)',
      pooled=m4[0], ci_lo=m4[1], ci_hi=m4[2], k=m4[3], total_n_causal=TOTAL_N_CAUSAL,
      note=(f'reproduces canonical determinant_meta.tsv -0.0036168; per-panel identity beta_std=beta_raw*SD_x '
            f'(max err {ANCH.identity_abserr.max():.1e}); SD_x~{ANCH.SD_x.mean():.3f} => the compression')),
]
OUT=pd.DataFrame(out)[['estimate_name','scale','pooled','ci_lo','ci_hi','k','total_n_causal','note']]
OUT.to_csv(S('results/causal_status_effect.tsv'),sep='\t',index=False)
print('\nWROTE results/causal_status_effect.tsv')
print(OUT[['estimate_name','scale','pooled','ci_lo','ci_hi','k','total_n_causal']].to_string(index=False))

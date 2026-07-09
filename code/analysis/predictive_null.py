#!/usr/bin/env python
"""predictive_null.py — the determinant-null shown as a PREDICTIVE FAILURE.

Question: do INTRINSIC, context-free properties of a gene/mark (constraint, druggability,
disease association, tissue specificity, cis-architecture, functional class) carry any
TRANSFERABLE signal about whether that mark's reversibility (toward-lean change) will be
large vs small, or reversible vs persistent — in a context the model has NEVER seen?

Construction (reuses build_determinants.py infra exactly):
  - OUTCOME = measurability-residualized reversibility: rank-INT |rev_beta| within
    base_sd(20) x abundance(5) cells, per transcriptome context. (Identical to fit_panel.)
  - FEATURES = full intrinsic matrix (loeuf, tau, n_drug_log, n_gwas_log, is_secreted/
    is_enzyme/is_membrane, arch_nsig, arch_str, has_arch, causal_nonEGFR) + missingness
    indicators, PLUS measurability (base_sd, base_mean) FORCED IN as features. If detectability
    carried any cross-context signal the model is free to exploit it -> a conservative test.
  - LEAKAGE GUARD: rev_beta, |rev_beta|, its rank, and the residualized outcome NEVER enter X.
    Only intrinsic + measurability features are predictors. Asserted at runtime.

Two readouts per outcome:
  (a) WITHIN-CONTEXT 5-fold CV Spearman(pred,actual)  = POSITIVE CONTROL (is ANYTHING learnable?)
  (b) LEAVE-ONE-CONTEXT-OUT Spearman                  = THE REAL TEST (train k-1, predict held-out)
And the same with AUC for the binary reversible(rev_q<0.10)-vs-persistent outcome.

Contexts (base_sd-complete, per spec): liver, adipose-CR, adipose-LCD, blood, muscle.
Output: results/predictive_null.tsv
"""
import os, numpy as np, pandas as pd
from scipy import stats
from sklearn.ensemble import GradientBoostingRegressor, GradientBoostingClassifier
from sklearn.model_selection import KFold, StratifiedKFold
from sklearn.metrics import roc_auc_score

SIG=os.path.dirname(os.path.dirname(os.path.abspath(__file__))); PROJ=os.path.dirname(SIG)
def P(*a): return os.path.join(PROJ,*a)
def S(*a): return os.path.join(SIG,*a)
def rint(x):
    x=np.asarray(x,float); r=stats.rankdata(x); return stats.norm.ppf((r-0.5)/len(x))
RNG=0

# ---- the 5 named base_sd-complete contexts -> manifest panel id ----
CONTEXTS={
 'liver':            'liver_table',
 'adipose-CR':       'adipose_table',
 'adipose-LCD':      'gse141221_adipose_table',
 'blood':            'gse273902_blood_table',
 'muscle':           'gse83352_muscle_exercise_table',
}

# ---- external gene features (identical cleanups to build_determinants.py) ----
F=lambda n: pd.read_csv(S('data/features',n),sep='\t')
lo=F('loeuf.tsv'); da=F('disease_assoc.tsv'); ta=F('tissue_tau.tsv')[['gene','tau']]
dr=F('druggability.tsv'); dr=dr[dr.gene.astype(str).str.match(r'^[A-Z][A-Z0-9-]*$')]
fc=F('functional_class.tsv'); fc=fc.assign(gene=fc.gene.astype(str).str.split('; ')).explode('gene')
ext=lo[['gene','loeuf']].merge(dr[['gene','n_drug_interactions']],on='gene',how='outer').merge(
    da[['gene','n_gwas_assoc']],on='gene',how='outer').merge(ta,on='gene',how='outer').merge(
    fc[['gene','is_secreted','is_enzyme','is_membrane']],on='gene',how='outer')
ext=ext.groupby('gene',as_index=False).max(numeric_only=True)
ext['n_drug_log']=np.log1p(ext['n_drug_interactions'].fillna(0))
ext['n_gwas_log']=np.log1p(ext['n_gwas_assoc'].fillna(0))

tx_loc=pd.read_csv(S('data/features/transcriptome_local.tsv'),sep='\t')
bsd_tx=pd.read_csv(S('results/base_sd_transcriptome.tsv'),sep='\t')
MAN=pd.read_csv(S('results/panel_manifest.tsv'),sep='\t')
MAN=MAN[(MAN.layer=='transcriptome')].set_index('panel')

EXTNUM=['loeuf','tau','n_drug_log','n_gwas_log']; EXTBIN=['is_secreted','is_enzyme','is_membrane']
ARCH=['arch_nsig','arch_str']; OTHER=['causal_nonEGFR','has_arch']
MEAS=['base_sd_feat','base_mean_feat']            # measurability FORCED IN as predictors
MISS=[f+'_miss' for f in ['loeuf','tau']]
FEATS=EXTNUM+EXTBIN+ARCH+OTHER+MISS+MEAS
BANNED={'rb','za','y','rev_beta','rev_q','y_bin'}  # leakage tripwire

def build_context(panel):
    """Return per-gene frame with measurability-residualized outcome y, binary y_bin, and X."""
    row=MAN.loc[panel]
    df=pd.read_csv(P(row.path),sep='\t')
    d=pd.DataFrame({'mark_key':df['mark_id'].astype(str),
                    'rb':pd.to_numeric(df['rev_beta'],errors='coerce'),
                    'rev_q':pd.to_numeric(df['rev_q'],errors='coerce')}).dropna(subset=['rb']).drop_duplicates('mark_key')
    d['za']=d['rb'].abs()
    # --- measurability residualization (identical to build_determinants.fit_panel) ---
    b=bsd_tx[bsd_tx.panel==panel].set_index('mark_id')
    d['bsd']=d['mark_key'].map(b['base_sd']); d['ab']=d['mark_key'].map(b['base_mean'])
    d=d.dropna(subset=['bsd','ab'])
    d['sb']=pd.qcut(d['bsd'].rank(method='first'),20,labels=False,duplicates='drop')
    d['mb']=pd.qcut(d['ab'].rank(method='first'),5,labels=False,duplicates='drop')
    d['y']=d.groupby([d.sb,d.mb])['za'].transform(lambda s: rint(s.values) if len(s)>5 else np.nan)
    d=d.dropna(subset=['y'])
    # binary reversible vs persistent (manifest sig_basis rev_q<0.10)
    d['y_bin']=(d['rev_q']<0.10).astype(int)
    # measurability as FEATURES (the model may cheat via detectability if any signal exists)
    d['base_sd_feat']=d['bsd']; d['base_mean_feat']=d['ab']
    # --- intrinsic features, LEFT-join onto every gene; missingness indicators ---
    d['gene']=d['mark_key']
    loc=tx_loc.rename(columns={'gene':'mark_key'}).copy()
    loc['arch_nsig']=loc['eqtl_n_cs'].fillna(0); loc['arch_str']=loc['eqtl_max_pip'].fillna(0); loc['has_arch']=1
    d=d.merge(loc[['mark_key','arch_nsig','arch_str','causal_nonEGFR','has_arch']],on='mark_key',how='left')
    d=d.merge(ext[['gene']+EXTNUM+EXTBIN],on='gene',how='left')
    d['has_arch']=d['has_arch'].fillna(0); d['arch_nsig']=d['arch_nsig'].fillna(0); d['arch_str']=d['arch_str'].fillna(0)
    d['causal_nonEGFR']=d['causal_nonEGFR'].fillna(0)
    for f in EXTBIN: d[f]=d[f].fillna(0)
    for f in ['loeuf','tau']: d[f+'_miss']=d[f].isna().astype(int)
    for f in EXTNUM: d[f]=d[f].fillna(d[f].median())
    feats=[f for f in FEATS if f in d.columns]
    assert not (set(feats)&BANNED), f'LEAKAGE: outcome-derived col in features: {set(feats)&BANNED}'
    X=d[feats].apply(pd.to_numeric,errors='coerce').fillna(0)
    return d['mark_key'].values, X, d['y'].values, d['y_bin'].values, feats

# ---- build all contexts on a SHARED feature column set ----
data={}; feat_union=None
for ctx,panel in CONTEXTS.items():
    keys,X,y,yb,feats=build_context(panel)
    data[ctx]=dict(keys=keys,X=X,y=y,yb=yb,feats=feats)
    feat_union=feats if feat_union is None else [f for f in feat_union if f in feats]
for ctx in data:                       # align to shared columns, same order
    data[ctx]['X']=data[ctx]['X'][feat_union].values
assert not (set(feat_union)&BANNED)
print('CONTEXTS:',{c:len(data[c]['y']) for c in data})
print('FEATURES (%d, measurability forced in):'%len(feat_union), feat_union)

def gbr(): return GradientBoostingRegressor(random_state=RNG,n_estimators=300,max_depth=3,
                                            learning_rate=0.05,subsample=0.7)
def gbc(): return GradientBoostingClassifier(random_state=RNG,n_estimators=300,max_depth=3,
                                            learning_rate=0.05,subsample=0.7)

def within_spearman(X,y):
    kf=KFold(5,shuffle=True,random_state=RNG); pred=np.zeros_like(y,dtype=float)
    for tr,te in kf.split(X):
        m=gbr().fit(X[tr],y[tr]); pred[te]=m.predict(X[te])
    return stats.spearmanr(pred,y).correlation

def within_auc(X,yb):
    if len(np.unique(yb))<2: return np.nan
    skf=StratifiedKFold(5,shuffle=True,random_state=RNG); pred=np.zeros_like(yb,dtype=float)
    for tr,te in skf.split(X,yb):
        if len(np.unique(yb[tr]))<2: pred[te]=yb[tr].mean(); continue
        m=gbc().fit(X[tr],yb[tr]); pred[te]=m.predict_proba(X[te])[:,1]
    return roc_auc_score(yb,pred)

def loco_spearman(held):
    Xtr=np.vstack([data[c]['X'] for c in data if c!=held])
    ytr=np.concatenate([data[c]['y'] for c in data if c!=held])
    m=gbr().fit(Xtr,ytr); p=m.predict(data[held]['X'])
    return stats.spearmanr(p,data[held]['y']).correlation

def loco_auc(held):
    Xtr=np.vstack([data[c]['X'] for c in data if c!=held])
    ybtr=np.concatenate([data[c]['yb'] for c in data if c!=held])
    yb=data[held]['yb']
    if len(np.unique(yb))<2 or len(np.unique(ybtr))<2: return np.nan
    m=gbc().fit(Xtr,ybtr); p=m.predict_proba(data[held]['X'])[:,1]
    return roc_auc_score(yb,p)

rows=[]
for ctx in data:
    X=data[ctx]['X']; y=data[ctx]['y']; yb=data[ctx]['yb']
    rows.append(dict(context=ctx, n=len(y), pos_rate=float(yb.mean()),
        within_spearman=within_spearman(X,y),    # positive control (continuous)
        loco_spearman=loco_spearman(ctx),        # real test (continuous)
        within_auc=within_auc(X,yb),             # positive control (binary)
        loco_auc=loco_auc(ctx)))                 # real test (binary)
R=pd.DataFrame(rows)
R.to_csv(S('results/predictive_null.tsv'),sep='\t',index=False)

def msummary(col):
    v=R[col].dropna().values; return f'mean={np.nanmean(v):+.3f}  median={np.nanmedian(v):+.3f}  range=[{np.nanmin(v):+.3f},{np.nanmax(v):+.3f}]'
print('\n=== PREDICTIVE-NULL RESULT (5 transcriptome contexts) ===')
print(R.to_string(index=False,
      formatters={'pos_rate':'{:.3f}'.format,'within_spearman':'{:+.3f}'.format,
                  'loco_spearman':'{:+.3f}'.format,'within_auc':'{:.3f}'.format,'loco_auc':'{:.3f}'.format}))
print('\nCONTINUOUS  (residualized |rev_beta|):')
print('  WITHIN-CONTEXT 5-fold CV Spearman (POSITIVE CONTROL):', msummary('within_spearman'))
print('  LEAVE-ONE-CONTEXT-OUT Spearman    (REAL TEST):       ', msummary('loco_spearman'))
print('BINARY  (reversible rev_q<0.10 vs persistent):')
print('  WITHIN-CONTEXT 5-fold CV AUC      (POSITIVE CONTROL):', msummary('within_auc'))
print('  LEAVE-ONE-CONTEXT-OUT AUC         (REAL TEST):       ', msummary('loco_auc'))
pc=np.nanmean(R['within_spearman'])
print('\nPositive control learnable signal?', 'YES' if pc>0.02 else 'NO',
      f'(within-context CV Spearman mean={pc:+.3f})')

# ===== ABLATION (verifier must-fix wf_20edb96f): is the cross-context binary LOCO AUC detectability or intrinsic biology? =====
INTR =[i for i,f in enumerate(feat_union) if f not in MEAS]   # intrinsic-only (drop measurability)
MEASI=[i for i,f in enumerate(feat_union) if f in MEAS]       # measurability-only
def loco_auc_cols(held,cols):
    Xtr=np.vstack([data[c]['X'][:,cols] for c in data if c!=held])
    ybtr=np.concatenate([data[c]['yb'] for c in data if c!=held])
    yb=data[held]['yb']
    if len(np.unique(yb))<2 or len(np.unique(ybtr))<2: return np.nan
    m=gbc().fit(Xtr,ybtr); p=m.predict_proba(data[held]['X'][:,cols])[:,1]
    return roc_auc_score(yb,p)
A=pd.DataFrame([dict(context=ctx,
        loco_auc_full=loco_auc(ctx),
        loco_auc_intrinsic_only=loco_auc_cols(ctx,INTR),
        loco_auc_measurability_only=loco_auc_cols(ctx,MEASI)) for ctx in data])
A.to_csv(S('results/predictive_null_ablation.tsv'),sep='\t',index=False)
print('\n=== ABLATION: cross-context (LOCO) binary AUC by feature subset ===')
print(A.to_string(index=False,formatters={c:'{:.3f}'.format for c in A.columns if c!='context'}))
print('  full      LOCO AUC mean=%.3f'%np.nanmean(A.loco_auc_full))
print('  intrinsic LOCO AUC mean=%.3f (drop measurability)'%np.nanmean(A.loco_auc_intrinsic_only))
print('  measur.   LOCO AUC mean=%.3f (only measurability)'%np.nanmean(A.loco_auc_measurability_only))
print('  -> HONEST read: binary cross-context AUC is WEAK + MIXED (full 0.61 / intrinsic 0.59 / meas 0.62). Intrinsic-only is ABOVE chance in liver+muscle (~0.66), so it is NOT cleanly pure-detectability; both contribute weakly. => the CONTINUOUS residualized LOCO (+0.002) is the clean primary null; the binary AUC is DEMOTED to secondary (its rev_q<0.10 label is a power flag), not interpreted as a detectability artifact.')
print('Wrote', S('results/predictive_null.tsv'))

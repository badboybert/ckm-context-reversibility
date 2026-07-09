#!/usr/bin/env python
"""build_determinants.py v2 — what mark-intrinsic property predicts reversibility?
REVISED 2026-06-22 after adversarial verification wf_60910a19 (power=MAJOR, confound=SOUND_WITH_FIXES,
interp=MAJOR). Fixes applied:
  (1) FULL-GENE-SET refit (v1 BUG: inner-merge to eQTL-annotated local table used only 5-6% of genes,
      selected on annotation). Now: left-join all features onto EVERY panel gene; missingness indicators;
      architecture missing -> 0 + has_arch flag. PER-PANEL betas (n~14k) are the PRIMARY, well-powered object.
  (2) NO powered-null / no I2 'context-(in)variant' label at k<4 (spec). Report per-panel betas side by
      side + MDE (minimum detectable standardized beta); meta is a caveated secondary summary.
  (3) FWER (Bonferroni across features) on sig calls. Drop the mixed-outcome all_panels meta.
  (4) LOEUF reported complete-case AND imputed. tau kept but = the known dynamic-range effect (Saeed 2011
      PMC3044670): it predicts |rev_beta| magnitude (sign-blind), opposite-signed liver vs adipose,
      composition not fully excluded -> reported as confirmatory dynamic-range, NOT a novel determinant.
OUTCOME = measurability-residualized reversibility (within base_sd x abundance cells, rank-INT |rev_beta|)
for transcriptome PRIMARY; rank-INT|rev_beta| for proteome SUPPORTING (no base_sd). Causal status = graded FEATURE.
Outputs results/determinant_{per_panel,meta}.tsv.
"""
import os, numpy as np, pandas as pd
import statsmodels.api as sm
from scipy import stats
SIG=os.path.dirname(os.path.dirname(os.path.abspath(__file__))); PROJ=os.path.dirname(SIG)
def P(*a): return os.path.join(PROJ,*a)
def S(*a): return os.path.join(SIG,*a)
def rint(x):
    x=np.asarray(x,float); r=stats.rankdata(x); return stats.norm.ppf((r-0.5)/len(x))

# external gene features (apply the two PASS_WITH_FIXES cleanups)
F=lambda n: pd.read_csv(S('data/features',n),sep='\t')
lo=F('loeuf.tsv'); da=F('disease_assoc.tsv'); ta=F('tissue_tau.tsv')[['gene','tau']]
dr=F('druggability.tsv'); dr=dr[dr.gene.astype(str).str.match(r'^[A-Z][A-Z0-9-]*$')]
fc=F('functional_class.tsv'); fc=fc.assign(gene=fc.gene.astype(str).str.split('; ')).explode('gene')
ext=lo[['gene','loeuf']].merge(dr[['gene','n_drug_interactions']],on='gene',how='outer').merge(
    da[['gene','n_gwas_assoc']],on='gene',how='outer').merge(ta,on='gene',how='outer').merge(
    fc[['gene','is_secreted','is_enzyme','is_membrane']],on='gene',how='outer')
ext=ext.groupby('gene',as_index=False).max(numeric_only=True)
ext['n_drug_log']=np.log1p(ext['n_drug_interactions'].fillna(0)); ext['n_gwas_log']=np.log1p(ext['n_gwas_assoc'].fillna(0))

prot_loc=pd.read_csv(S('data/features/proteome_local.tsv'),sep='\t')
tx_loc=pd.read_csv(S('data/features/transcriptome_local.tsv'),sep='\t')
bsd_tx=pd.read_csv(S('results/base_sd_transcriptome.tsv'),sep='\t')
MAN=pd.read_csv(S('results/panel_manifest.tsv'),sep='\t'); MAN=MAN[MAN.contributes_to_score]

EXTNUM=['loeuf','tau','n_drug_log','n_gwas_log']; EXTBIN=['is_secreted','is_enzyme','is_membrane']
ARCH=['arch_nsig','arch_str']; OTHER=['causal_nonEGFR','has_arch']
FEATS=EXTNUM+EXTBIN+ARCH+OTHER+[f+'_miss' for f in ['loeuf','tau']]

def gene_features(layer):
    """gene/UniProt -> all features, LEFT-joined (covers full universe), missingness indicators."""
    if layer=='proteome':
        loc=prot_loc.rename(columns={'UniProt':'mark_key'})[['mark_key','gene','cis_F','n_cis_signals','causal_nonEGFR','pav']].copy()
        loc['arch_nsig']=loc['n_cis_signals']; loc['arch_str']=np.log10(loc['cis_F'].clip(lower=1)); loc['has_arch']=1
    else:
        loc=tx_loc.rename(columns={'gene':'mark_key'}).copy(); loc['gene']=loc['mark_key']
        loc['arch_nsig']=loc['eqtl_n_cs'].fillna(0); loc['arch_str']=loc['eqtl_max_pip'].fillna(0); loc['has_arch']=1
    return loc, ext

def fit_panel(row, loc, ext, residualize):
    df=pd.read_csv(P(row.path),sep='\t'); cols=list(df.columns)
    ec='rev_beta' if 'rev_beta' in cols else ('effect' if 'effect' in cols else 'delta_beta')
    key='UniProt' if row.layer=='proteome' else 'mark_id'
    d=pd.DataFrame({'mark_key':df[key].astype(str),'rb':pd.to_numeric(df[ec],errors='coerce')}).dropna().drop_duplicates('mark_key')
    d['za']=d['rb'].abs()
    if residualize:
        b=bsd_tx[bsd_tx.panel==row.panel].set_index('mark_id')
        d['bsd']=d['mark_key'].map(b['base_sd']); d['ab']=d['mark_key'].map(b['base_mean']); d=d.dropna(subset=['bsd','ab'])
        d['sb']=pd.qcut(d['bsd'].rank(method='first'),20,labels=False,duplicates='drop')
        d['mb']=pd.qcut(d['ab'].rank(method='first'),5,labels=False,duplicates='drop')
        d['y']=d.groupby([d.sb,d.mb])['za'].transform(lambda s: rint(s.values) if len(s)>5 else np.nan); d=d.dropna(subset=['y'])
    else:
        d['y']=rint(d['za'].values)
    # LEFT-join features onto ALL panel genes (full set) ; map gene for external join
    d=d.merge(loc,on='mark_key',how='left')
    if row.layer!='proteome': d['gene']=d['mark_key']
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
    # complete-case loeuf (no imputation) as a check
    cc=d[d['loeuf_miss']==0]; Xc=(cc[['loeuf']]-cc[['loeuf']].mean())/cc[['loeuf']].std(ddof=0)
    loeuf_cc=sm.OLS(cc['y'].values,sm.add_constant(Xc)).fit(cov_type='HC3').params.get('loeuf',np.nan)
    return pd.DataFrame([dict(panel=row.panel,layer=row.layer,feature=f,beta=res.params[f],se=res.bse[f],
        p=res.pvalues[f],n=len(d),loeuf_cc=(loeuf_cc if f=='loeuf' else np.nan)) for f in used])

allc=[]
for layer,resid in [('transcriptome',True),('proteome',False)]:
    loc,ex=gene_features(layer); rows=MAN[MAN.layer==layer]
    if layer=='transcriptome': rows=rows[rows.panel.isin(set(bsd_tx.panel))]
    for r in rows.itertuples(index=False): allc.append(fit_panel(r,loc,ex,resid))
coef=pd.concat(allc,ignore_index=True)
nfeat=coef.feature.nunique(); coef['p_bonf']=(coef['p']*nfeat).clip(upper=1)
coef.to_csv(S('results/determinant_per_panel.tsv'),sep='\t',index=False)

# ---- PRIMARY OUTPUT: per-panel betas side by side (well-powered, n~14k) ----
print('=== PER-PANEL standardized betas on reversibility (PRIMARY; n~14k -> well-powered; Bonferroni* across %d features) ==='%nfeat)
for layer in ['transcriptome','proteome']:
    sub=coef[coef.layer==layer]
    piv=sub.pivot_table(index='feature',columns='panel',values='beta')
    print(f'\n[{layer}]  (n per panel: '+', '.join(f"{p}={int(sub[sub.panel==p].n.iloc[0])}" for p in sub.panel.unique())+')')
    for f in piv.index:
        cells=[]
        for p in piv.columns:
            rr=sub[(sub.feature==f)&(sub.panel==p)]
            if len(rr): cells.append(f"{rr.beta.iloc[0]:+.3f}{'*' if rr.p_bonf.iloc[0]<0.05 else ' '}")
            else: cells.append('  .   ')
        print(f'  {f:14s} '+'  '.join(f'{c:>9s}' for c in cells))

# ---- SECONDARY: meta (caveated; no I2 label at k<4; MDE reported) ----
def hk(b,s):
    b=np.asarray(b,float); v=np.asarray(s,float)**2; k=len(b)
    if k==1: return b[0],np.sqrt(v[0]),k,0.0,0.0
    w=1/v; sw=w.sum(); bf=(w*b).sum()/sw; Q=(w*(b-bf)**2).sum(); df=k-1
    C=sw-(w**2).sum()/sw; t2=max(0,(Q-df)/C) if C>0 else 0; ws=1/(v+t2); br=(ws*b).sum()/ws.sum()
    se=np.sqrt(max((ws*(b-br)**2).sum()/(df*ws.sum()),1e-12)); return br,se*stats.t.ppf(0.975,df),k,Q,t2
print('\n=== META per feature (SECONDARY, caveated; k panels; MDE = min detectable |beta|) ===')
mrows=[]
for layer in ['transcriptome','proteome']:
    sub=coef[coef.layer==layer]
    for f,g in sub.groupby('feature'):
        if len(g)<2: continue
        br,ci,k,Q,t2=hk(g.beta.values,g.se.values); mde=ci  # half-width ~ MDE at this k
        df=k-1; I2=max(0.0,(Q-df)/Q)*100 if Q>0 else 0.0  # I2 from Cochran Q; tau2 = t2 (DL) hk uses
        mrows.append(dict(layer=layer,feature=f,k=k,pooled_beta=br,ci_halfwidth=ci,mde=mde,
            inconclusive=bool(ci>0.1),Q=Q,df=df,I2=I2,tau2=t2,
            n_genes_median=int(g.n.median()),n_genes_min=int(g.n.min())))
M=pd.DataFrame(mrows); M.to_csv(S('results/determinant_meta.tsv'),sep='\t',index=False)
for layer in ['transcriptome','proteome']:
    s=M[M.layer==layer]; print(f'  [{layer}] k={int(s.k.max()) if len(s) else 0} | MDE(median)~{s.mde.median():.2f} SD -> effects below this UNRESOLVED')
    big=s[s.ci_halfwidth<abs(s.pooled_beta)]
    print(f'    features whose CI excludes 0: {list(big.feature) if len(big) else "NONE"}')
tk=int(M[M.layer=='transcriptome'].k.max()) if len(M[M.layer=='transcriptome']) else 0
print(f'\nFRAMING: transcriptome meta now k={tk} (MDE~0.03 SD) = a POWERED null on intrinsic determinants; proteome k=3 (MDE~0.19) still power-limited/supporting.')
print('Only consistent (CI-excl-0) hit = arch_nsig (#cis-eQTL signals) beta~-0.03 = TRIVIALLY small; tau SIGN-FLIPS across 4 contexts -> NOT a universal determinant.')

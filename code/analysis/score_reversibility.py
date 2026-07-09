#!/usr/bin/env python
"""score_reversibility.py — per-mark RESPONSIVENESS score for the signature pivot.
Implements SCORE_SPEC_FROZEN.md sections 2-4 (responsiveness arm only; restoration +
persistence are separate downstream steps that need BMI axes + base_sd).

Per layer, over its SCORE-contributing panels (panel_manifest.tsv, contributes_to_score):
  1. within-platform blom rank-INT of signed effect  -> z_mp           (locked decision #1)
  2. has-SE layers  -> Hartung-Knapp RE-meta of signed z_mp            (NOT DerSimonian-Laird)
                       responsiveness Z = |pooled|, I2 = consistency axis
     methylation    -> vote-count V = fraction of panels sig(q<0.10)   (no SE on disk)
                       + unweighted mean z as a magnitude companion
  3. degeneracy guard: |z|<0.1 (or SE missing) -> se_z = 1/sqrt(n-2), down-weighted
  4. mandatory QC: corr(Z, sum n_pairs) and corr(Z, k) -> is the score a sample-size proxy?
Outputs data/signature/reversibility_score_{layer}.tsv (+ a stacked file for the determinant step).
"""
import os, json, numpy as np, pandas as pd
from scipy import stats

SIG=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))   # signature-pivot/
PROJ=os.path.dirname(SIG)                                          # weight-loss.omics/ (shared root)
def P(*a): return os.path.join(PROJ,*a)   # SHARED inputs (reversal tables etc., referenced)
def S(*a): return os.path.join(SIG,*a)    # signature-pivot outputs

MAN=pd.read_csv(S('results/panel_manifest.tsv'),sep='\t')
MAN=MAN[MAN.contributes_to_score].copy()

KEYCOL={'proteome':'UniProt','transcriptome':'mark_id','methylation':'cpg'}  # metabolome handled by normalization
EFF=['rev_beta','effect','delta_beta']; SE=['rev_se','se']
Q=['q','rev_q','fdr']; PV=['rev_p','p','pval']

def pick(cols,c):
    for x in c:
        if x in cols: return x
    return None
def norm_metab(s):  # within-metabolome key: strip Metabolon confidence asterisks + case/space
    return str(s).strip().rstrip('*').strip().lower()
def blom(x):
    x=np.asarray(x,float); r=stats.rankdata(x); n=len(x)
    return stats.norm.ppf((r-0.5)/n)

def hk_meta(y,v):
    """Hartung-Knapp RE-meta. y=effects, v=variances. -> pooled,se,lo,hi,I2,k."""
    y=np.asarray(y,float); v=np.asarray(v,float); k=len(y)
    if k==1:
        se=np.sqrt(v[0]); return y[0],se,y[0]-1.96*se,y[0]+1.96*se,0.0,1
    w=1.0/v; sw=w.sum(); yfe=(w*y).sum()/sw
    Q_=(w*(y-yfe)**2).sum(); df=k-1
    C=sw-(w**2).sum()/sw; tau2=max(0.0,(Q_-df)/C) if C>0 else 0.0
    ws=1.0/(v+tau2); yre=(ws*y).sum()/ws.sum()
    q_hk=(ws*(y-yre)**2).sum()/(df*ws.sum())          # HK small-sample variance
    se=np.sqrt(max(q_hk,1e-12)); tcrit=stats.t.ppf(0.975,df)
    i2=max(0.0,(Q_-df)/Q_) if Q_>0 else 0.0
    return yre,se,yre-tcrit*se,yre+tcrit*se,i2,k

def load_panel(row):
    df=pd.read_csv(P(row.path),sep='\t')
    cols=list(df.columns)
    ec=pick(cols,EFF); sc=pick(cols,SE); qc=pick(cols,Q); pc=pick(cols,PV)
    if row.layer=='metabolome':
        kc=pick(cols,['mark_id','metabolite']); key=df[kc].map(norm_metab)
    else:
        kc=KEYCOL[row.layer]; key=df[kc].astype(str)
    eff=pd.to_numeric(df[ec],errors='coerce')
    sig=(pd.to_numeric(df[qc],errors='coerce')<0.10) if qc else ((pd.to_numeric(df[pc],errors='coerce')<0.05) if pc else False)
    out=pd.DataFrame({'mark_key':key,'eff':eff,'sig':sig.astype(float) if hasattr(sig,'astype') else 0.0})
    out['se']=pd.to_numeric(df[sc],errors='coerce') if sc else np.nan
    out=out.dropna(subset=['eff']).drop_duplicates('mark_key')
    out['z']=blom(out['eff'].values)
    # delta-method z-scale SE with degeneracy guard
    n=int(row.n_pairs)
    if sc and out['se'].notna().any():
        t=out['eff']/out['se']; sez=(out['z'].abs()/t.abs())
        floor=1.0/np.sqrt(max(n-2,1))
        sez=sez.where((out['z'].abs()>=0.1)&np.isfinite(sez), floor)
    else:
        sez=pd.Series(1.0/np.sqrt(max(n-2,1)),index=out.index)
    out['vz']=sez**2; out['panel']=row.panel; out['n_pairs']=n; out['has_se']=bool(row.has_se)
    return out[['mark_key','panel','z','vz','sig','n_pairs','has_se']]

stacked=[]
for layer in ['proteome','transcriptome','methylation','metabolome']:
    rows=MAN[MAN.layer==layer]
    long=pd.concat([load_panel(r) for r in rows.itertuples(index=False)],ignore_index=True)
    if layer=='methylation':                            # vote-count PRIMARY (no SE) -> fully vectorized (866k CpGs)
        long['zpos']=(long['z']>0).astype(float)
        gb=long.groupby('mark_key',sort=False)
        sc=gb.agg(k_panels=('panel','size'),n_total_pairs=('n_pairs','sum'),
                  V_vote=('sig','mean'),mean_z=('z','mean'),posfrac=('zpos','mean')).reset_index()
        sc['layer']=layer; sc['pooled_signed']=sc['mean_z']; sc['Z_responsiveness']=sc['mean_z'].abs()
        sc['Z_ci_lo']=np.nan; sc['Z_ci_hi']=np.nan; sc['I2']=np.nan; sc['panels']=''
        sc['sign_concord']=np.maximum(sc['posfrac'],1-sc['posfrac']); sc=sc.drop(columns='posfrac')
    else:                                                # HK-RE PRIMARY -> per-mark loop (small N)
        recs=[]
        for key,g in long.groupby('mark_key',sort=False):
            pooled,se,lo,hi,i2,k=hk_meta(g['z'].values,g['vz'].values)
            recs.append(dict(mark_key=key,layer=layer,k_panels=k,n_total_pairs=int(g['n_pairs'].sum()),
                panels=','.join(sorted(g['panel'])),Z_responsiveness=abs(pooled),pooled_signed=pooled,
                Z_ci_lo=lo,Z_ci_hi=hi,I2=i2,V_vote=float(g['sig'].mean()),mean_z=float(g['z'].mean()),
                sign_concord=float(max((g['z']>0).mean(),(g['z']<0).mean()))))
        sc=pd.DataFrame(recs)
    op=S(f'results/reversibility_score_{layer}.tsv'); sc.to_csv(op,sep='\t',index=False)
    stacked.append(sc)
    # QC: is responsiveness a sample-size proxy?
    multi=sc[sc.k_panels>=2]
    def corr(a,b):
        m=np.isfinite(a)&np.isfinite(b);
        return np.corrcoef(a[m],b[m])[0,1] if m.sum()>3 else np.nan
    r_n=corr(sc.Z_responsiveness.values,np.log(sc.n_total_pairs.values))
    r_k=corr(sc.Z_responsiveness.values,sc.k_panels.values.astype(float))
    print(f'[{layer:13s}] marks={len(sc):6d} k>=2={len(multi):6d} k=1={len(sc)-len(multi):6d} '
          f'| QC corr(Z,log Σn)={r_n:+.3f} corr(Z,k)={r_k:+.3f} '
          f'| medZ={sc.Z_responsiveness.median():.3f} '+('medI2=%.2f'%multi.I2.median() if layer!='methylation' and len(multi) else ''))
    print(f'               wrote {op}')

allsc=pd.concat(stacked,ignore_index=True)
allsc.to_csv(S('results/reversibility_score_all_layers.tsv'),sep='\t',index=False)
print(f'\nstacked -> signature-pivot/results/reversibility_score_all_layers.tsv  (n={len(allsc)})  [determinant-step input ONLY; Z never meta-analyzed across layers]')
print('\nNOTE: responsiveness arm only. Restoration (BMI-axis) + persistence (variance-conditional) + class = pending downstream steps.')
print('QC reminder (spec s4): if |corr(Z,log Σn)|>0.3 -> residualize Z on log(Σn) before classification + publish (class x power-tertile).')

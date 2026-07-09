#!/usr/bin/env python
"""build_restoration.py — RESTORATION arm of the signature pivot (SCORE_SPEC_FROZEN.md s1-2;
PI decision #1). Restoration = does the intervention move a mark TOWARD the lean reference,
oriented by an external BMI-association axis (NOT by the mark's own sign = circular).

Direction matters here, so (unlike responsiveness which uses median-centered blom-z for
RELATIVE magnitude) restoration uses a SCALE-ONLY within-platform standardization that
PRESERVES zero + direction:  zr_mp = rev_beta_mp / MAD0(panel),  MAD0 = median(|rev_beta|).
Then restoration r_mp = -sign(bmi_axis_beta_m) * zr_mp   ( >0 = moved toward lean ).
Pool across panels with Hartung-Knapp RE-meta. R>0 & CI_lo>0 => 'restorative'.
Available now for layers with a BMI axis on disk: PROTEOME (Goudswaard ST5) + METABOLOME (BBS bmi_mr).
methylation/transcriptome restoration await BMI-EWAS (Wahl 2017) + a BMI-transcriptome signature.
"""
import os, numpy as np, pandas as pd
from scipy import stats
SIG=os.path.dirname(os.path.dirname(os.path.abspath(__file__))); PROJ=os.path.dirname(SIG)
def P(*a): return os.path.join(PROJ,*a)
def S(*a): return os.path.join(SIG,*a)

MAN=pd.read_csv(S('results/panel_manifest.tsv'),sep='\t'); MAN=MAN[MAN.contributes_to_score]
EFF=['rev_beta','effect','delta_beta']; SE=['rev_se','se']
def pick(c,cands):
    for x in cands:
        if x in c: return x
def nm(s): return str(s).strip().rstrip('*').strip().lower()

def hk_meta(y,v):
    y=np.asarray(y,float); v=np.asarray(v,float); k=len(y)
    if k==1:
        se=np.sqrt(v[0]); return y[0],se,y[0]-1.96*se,y[0]+1.96*se,0.0,1
    w=1/v; sw=w.sum(); yfe=(w*y).sum()/sw; Q=(w*(y-yfe)**2).sum(); df=k-1
    C=sw-(w**2).sum()/sw; tau2=max(0.0,(Q-df)/C) if C>0 else 0.0
    ws=1/(v+tau2); yre=(ws*y).sum()/ws.sum()
    se=np.sqrt(max((ws*(y-yre)**2).sum()/(df*ws.sum()),1e-12)); tc=stats.t.ppf(0.975,df)
    i2=max(0.0,(Q-df)/Q) if Q>0 else 0.0
    return yre,se,yre-tc*se,yre+tc*se,i2,k

# --- BMI axes (lean orientation) ---
prot_ax=pd.read_csv(S('results/bmi_axis_proteome.tsv'),sep='\t').set_index('UniProt')['bmi_beta'].to_dict()
bbs=pd.read_csv(P('data/metabolome/bbs_reversal_metabolon.tsv'),sep='\t')
metab_ax={nm(r['mark_id']):r['bmi_mr_beta'] for _,r in bbs.iterrows() if pd.notna(r.get('bmi_mr_beta'))}
AX={'proteome':('UniProt',prot_ax),'metabolome':(None,metab_ax)}

def panel_long(row,keymode,axis):
    df=pd.read_csv(P(row.path),sep='\t'); cols=list(df.columns)
    ec=pick(cols,EFF); sc=pick(cols,SE)
    key=(df[pick(cols,['mark_id','metabolite'])].map(nm) if row.layer=='metabolome'
         else df[keymode].astype(str))
    eff=pd.to_numeric(df[ec],errors='coerce')
    d=pd.DataFrame({'mark_key':key,'eff':eff}); d['se']=pd.to_numeric(df[sc],errors='coerce') if sc else np.nan
    d=d.dropna(subset=['eff']).drop_duplicates('mark_key')
    bmi=d['mark_key'].map(axis); d=d[bmi.notna()].copy(); bmi=bmi[d.index]
    scale0=np.median(np.abs(d['eff'])) or 1.0
    zr=d['eff']/scale0
    vzr=(d['se']/scale0)**2 if sc and d['se'].notna().any() else pd.Series(np.nan,index=d.index)
    vzr=vzr.fillna(np.nanvar(zr) if np.isfinite(np.nanvar(zr)) else 1.0)
    lean=-np.sign(bmi)
    return pd.DataFrame({'mark_key':d['mark_key'],'r':lean*zr,'v':vzr,'panel':row.panel,'n_pairs':int(row.n_pairs)})

for layer,(keymode,axis) in AX.items():
    rows=MAN[MAN.layer==layer]
    long=pd.concat([panel_long(r,keymode,axis) for r in rows.itertuples(index=False)],ignore_index=True)
    recs=[]
    for key,g in long.groupby('mark_key',sort=False):
        R,se,lo,hi,i2,k=hk_meta(g['r'].values,g['v'].values)
        recs.append(dict(mark_key=key,layer=layer,k_panels=k,n_total_pairs=int(g['n_pairs'].sum()),
            R_restoration=R,R_ci_lo=lo,R_ci_hi=hi,I2=i2,restorative=bool(lo>0),anti_restorative=bool(hi<0)))
    sc=pd.DataFrame(recs); op=S(f'results/restoration_score_{layer}.tsv'); sc.to_csv(op,sep='\t',index=False)
    nrest=int(sc.restorative.sum()); nanti=int(sc.anti_restorative.sum())
    print(f'[{layer:11s}] marks_with_axis={len(sc):5d} k>=2={int((sc.k_panels>=2).sum()):5d} '
          f'| restorative(CI>0)={nrest} ({nrest/len(sc)*100:.0f}%)  anti={nanti} ({nanti/len(sc)*100:.0f}%)  medR={sc.R_restoration.median():+.2f}')
    print(f'              wrote {op}')

# ---- METHYLATION restoration = CONFOUNDED by long-interval AGE-DRIFT (NOT cleanly measurable) ----
# Diagnostic: the 18mo (DIRECT-PLUS) / 5yr (PREDIMED) methylation change correlates POSITIVELY with the
# cross-sectional BMI-methylation axis -> methylation moves in the high-adiposity direction, not toward lean.
# Conventions verified (delta=post-pre; bmi+ = higher with BMI). This is age-related methylation drift
# (epigenetic aging tracks the BMI signature) dominating the tiny weight-loss signal over long intervals.
from scipy import stats as _st
mbmi=pd.read_csv(S('data/features/bmi_axis_methylation.tsv'),sep='\t').drop_duplicates('cpg').set_index('cpg')['bmi_beta']
mpan=MAN[(MAN.layer=='methylation')&(MAN.persistence_eligible)]
print('[methylation ] restoration CONFOUNDED by long-interval age-drift -> NOT reported as a clean restoration score:')
_drift=[]
for r in mpan.itertuples(index=False):
    d=pd.read_csv(P(r.path),sep='\t')[['cpg','delta_beta','q']].copy(); d['bmi']=d['cpg'].map(mbmi); d=d.dropna(subset=['bmi','delta_beta'])
    sig=d[pd.to_numeric(d['q'],errors='coerce')<0.10]
    rho=_st.spearmanr(d.delta_beta,d.bmi_beta).correlation if 'bmi_beta' in d else _st.spearmanr(d.delta_beta,d.bmi).correlation
    tl=(np.sign(sig.delta_beta)==np.sign(-sig.bmi)).mean()
    print(f'    {r.panel:12s}: Spearman(delta, BMI-axis)={rho:+.2f} (>0 = moves WITH adiposity = age-drift); toward-lean among sig={tl:.2f}')
    _drift.append(dict(panel=r.panel,spearman_delta_bmi=rho,toward_lean_sig=tl,n=len(d),n_sig=len(sig)))
_op=S('results/orphan/methylation_age_drift.tsv'); pd.DataFrame(_drift).to_csv(_op,sep='\t',index=False)
print(f'wrote {_op}')
print('    -> methylation RESPONSIVENESS + PERSISTENCE are valid (within-interval); RESTORATION is age-drift-confounded over 18mo-5yr -> reported as a LIMITATION, not a result.')

# sanity: leptin (BMI raises it; weight loss lowers it -> RESTORATIVE)
pr=pd.read_csv(S('results/restoration_score_proteome.tsv'),sep='\t')
for u,n in [('P41159','LEP'),('P08833','IGFBP1'),('P15090','FABP4')]:
    r=pr[pr.mark_key==u]
    if len(r): print(f'  sanity prot {n:7s}: R={r.R_restoration.iloc[0]:+.2f} CI[{r.R_ci_lo.iloc[0]:+.2f},{r.R_ci_hi.iloc[0]:+.2f}] restorative={bool(r.restorative.iloc[0])}')
mb=pd.read_csv(S('results/restoration_score_metabolome.tsv'),sep='\t')
for n in ['valine','leucine','isoleucine','glutamate']:
    r=mb[mb.mark_key==n]
    if len(r): print(f'  sanity metab {n:10s}: R={r.R_restoration.iloc[0]:+.2f} restorative={bool(r.restorative.iloc[0])}')

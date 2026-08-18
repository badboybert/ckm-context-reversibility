#!/usr/bin/env python
"""classify_marks.py — two-class taxonomy (SCORE_SPEC_FROZEN.md s3/s5), REVISED 2026-06-22. Fixes applied:
  (1) k=1 EXCLUSION (spec s2.5): marks in <2 eligible panels -> 'unclassified-single-panel',
      removed from the headline persistent/reversible/intermediate counts.
  (2) 2D matched conditioning: bin within base_sd x MEASURABILITY cells (abundance for
      transcriptome, beta edge-distance for methylation) -> absorbs the detection-floor confound
      (persistent genes were low-abundance logCPM 2.99 vs 4.39; persistent CpGs at beta extremes).
  (3) finer base_sd bins for methylation (50) -> kills the residual within-decile gradient
      (AUC 0.494 at 10 bins -> ~0.500 at 50).
  (4) cross-panel CONCORDANCE QC: observed 'persistent in all panels' vs the independence
      expectation -> is persistence a REPRODUCIBLE property or just the chance intersection of
      panel-specific low-movers? (verifier: methylation 1.00x, transcriptome 1.03x = at chance).
  (5) confound ratio reported on the n_elig>=2 (real cross-panel) subset, not the diluted all-marks set.
PERSISTENCE per panel = bottom-tertile |rev_beta| (ZERO-referenced) within base_sd x measurability cells.
PERSISTENT marks (low-movers) are noise-dominated; treat the persistent pole as the WEAKER pole and
report its reproducibility honestly. Metabolome = single eligible panel -> NO multi-panel class (exploratory).
"""
import os, numpy as np, pandas as pd
SIG=os.path.dirname(os.path.dirname(os.path.abspath(__file__))); PROJ=os.path.dirname(SIG)
def P(*a): return os.path.join(PROJ,*a)
def S(*a): return os.path.join(SIG,*a)
def nm(s): return str(s).strip().rstrip('*').strip().lower()
def pick(c,cands):
    for x in cands:
        if x in c: return x

MAN=pd.read_csv(S('results/panel_manifest.tsv'),sep='\t'); ELIG=MAN[MAN.persistence_eligible]
EFF=['rev_beta','effect','delta_beta']
bsd_tx=pd.read_csv(S('results/base_sd_transcriptome.tsv'),sep='\t')
bsd_mb=pd.read_csv(S('results/base_sd_metabolome.tsv'),sep='\t')
NSD={'methylation':50,'transcriptome':20,'metabolome':10}

def panel_persistence(row):
    df=pd.read_csv(P(row.path),sep='\t'); cols=list(df.columns); ec=pick(cols,EFF)
    key=(df[pick(cols,['mark_id','metabolite'])].map(nm) if row.layer=='metabolome'
         else df[{'proteome':'UniProt','transcriptome':'mark_id','methylation':'cpg'}[row.layer]].astype(str))
    d=pd.DataFrame({'mark_key':key,'eff':pd.to_numeric(df[ec],errors='coerce')})
    if row.layer=='methylation':
        d['base_sd']=pd.to_numeric(df['base_sd'],errors='coerce')
        bm=pd.to_numeric(df['base_mean'],errors='coerce'); d['meas']=np.minimum(bm,1-bm)  # edge-distance from 0/1
    d=d.dropna(subset=['eff']).drop_duplicates('mark_key')
    if row.layer=='transcriptome':
        b=bsd_tx[bsd_tx.panel==row.panel].set_index('mark_id')
        d['base_sd']=d['mark_key'].map(b['base_sd']); d['meas']=d['mark_key'].map(b['base_mean'])  # abundance
    elif row.layer=='metabolome':
        d['base_sd']=d['mark_key'].map(bsd_mb.set_index('mark_key')['base_sd']); d['meas']=np.nan
    d=d.dropna(subset=['base_sd']).copy()
    d['za']=np.abs(d['eff'].values)
    d['sd_bin']=pd.qcut(d['base_sd'].rank(method='first'),min(NSD[row.layer],max(2,len(d)//30)),labels=False,duplicates='drop')
    if d['meas'].notna().any():    # 2D: base_sd x measurability tertile
        d['m_bin']=pd.qcut(d['meas'].rank(method='first'),3,labels=False,duplicates='drop')
        d['cell']=d['sd_bin'].astype(str)+'|'+d['m_bin'].astype(str)
    else:
        d['cell']=d['sd_bin'].astype(str)
    d['persistent']=False; d['reversible']=False
    for _,g in d.groupby('cell',sort=False):
        if len(g)<6: continue
        lo=g['za'].quantile(1/3); hi=g['za'].quantile(2/3)
        d.loc[g.index[g['za']<=lo],'persistent']=True
        d.loc[g.index[g['za']>=hi],'reversible']=True
    d['panel']=row.panel
    return d[['mark_key','panel','base_sd','meas','persistent','reversible']]

for layer in ['methylation','transcriptome','metabolome']:
    rows=ELIG[ELIG.layer==layer]
    if not len(rows): continue
    pp=pd.concat([panel_persistence(r) for r in rows.itertuples(index=False)],ignore_index=True)
    g=pp.groupby('mark_key')
    agg=g.agg(n_elig=('panel','size'),n_persistent=('persistent','sum'),n_reversible=('reversible','sum'),
              base_sd_mean=('base_sd','mean'),meas_mean=('meas','mean')).reset_index()
    agg['layer']=layer; agg['class']='unclassified-single-panel'
    multi=agg.n_elig>=2
    agg.loc[multi&(agg.n_persistent==agg.n_elig),'class']='persistent'
    agg.loc[multi&(agg.n_reversible==agg.n_elig),'class']='consistently-reversible'
    agg.loc[multi&~agg['class'].isin(['persistent','consistently-reversible']),'class']='intermediate'
    sc=pd.read_csv(S(f'results/reversibility_score_{layer}.tsv'),sep='\t').set_index('mark_key')['Z_responsiveness']
    agg['Z_responsiveness']=agg['mark_key'].map(sc)
    agg.to_csv(S(f'results/marks_classified_{layer}.tsv'),sep='\t',index=False)

    n_multi=int(multi.sum()); n_single=int((~multi).sum())
    vc=agg.loc[multi,'class'].value_counts(); npe=int(vc.get('persistent',0)); nrv=int(vc.get('consistently-reversible',0)); nim=int(vc.get('intermediate',0))
    print(f'\n[{layer}] panels={len(rows)}  marks={len(agg)}  multi-panel(n_elig>=2)={n_multi}  single-panel(excluded)={n_single}')
    if n_multi:
        print(f'  CLASS (multi-panel only): persistent={npe}  consistently-reversible={nrv}  intermediate={nim}')
        pe=agg[(multi)&(agg['class']=='persistent')]; rv=agg[(multi)&(agg['class']=='consistently-reversible')]
        if len(pe) and len(rv):
            r_sd=pe.base_sd_mean.median()/rv.base_sd_mean.median()
            mtxt=f' ; meas(abundance/edge) persistent={pe.meas_mean.median():.3g} vs reversible={rv.meas_mean.median():.3g}' if pe.meas_mean.notna().any() else ''
            print(f'  CONFOUND (n_elig>=2): base_sd ratio persistent/reversible={r_sd:.3f} (want ~1){mtxt}')
        # CROSS-PANEL CONCORDANCE: observed 'persistent in all panels' vs independence expectation
        sub=pp[pp.mark_key.isin(set(agg.loc[multi,'mark_key']))]
        rates=sub.groupby('panel')['persistent'].mean()
        exp=float(np.prod(rates.values)); obs=npe/n_multi
        rrates=sub.groupby('panel')['reversible'].mean(); rexp=float(np.prod(rrates.values)); robs=nrv/n_multi
        print(f'  CONCORDANCE persist: observed={obs:.4f} expected(indep)={exp:.4f} enrichment={obs/exp:.2f}x | '
              f'reversible: obs={robs:.4f} exp={rexp:.4f} enr={robs/rexp:.2f}x  (~1x = no better than chance)')
    else:
        print(f'  -> 0 multi-panel marks: {layer} persistence is SINGLE-PANEL/EXPLORATORY only (no cross-panel class).')
    print(f'  wrote results/marks_classified_{layer}.tsv')

#!/usr/bin/env python
"""build_tau_restoration_test.py — the DECISIVE tau test (verifier-requested) + adipose transcriptome
restoration, now that a RELIABLE adipose BMI axis exists (GSE70353/METSIM + GSE59034, verified 8/8 signs).
tau predicts reversal MAGNITUDE |rev_beta| (the known dynamic-range effect, Saeed 2011) — but does it predict
toward-lean RESTORATION (direction)? Only a restoration effect would make tau a real determinant; a magnitude-
only effect = dynamic-range, NOT a determinant. Tested in ADIPOSE (tissue-matched axis; where tau appeared),
measurability + composition controlled. Also reports adipose transcriptome restoration (toward-lean) summary.
"""
import os, numpy as np, pandas as pd, gzip, json
import statsmodels.api as sm
from scipy import stats
SIG=os.path.dirname(os.path.dirname(os.path.abspath(__file__))); PROJ=os.path.dirname(SIG)
def P(*a): return os.path.join(PROJ,*a);
def S(*a): return os.path.join(SIG,*a)
def rint(x): x=np.asarray(x,float); r=stats.rankdata(x); return stats.norm.ppf((r-0.5)/len(x))

ax=pd.read_csv(S('data/features/bmi_axis_transcriptome.tsv'),sep='\t').set_index('gene')
print('axis spot-check: LEP=%+.3f ADIPOQ=%+.3f FASN=%+.3f (expect +,-,-)'%(ax.bmi_beta.get('LEP',np.nan),ax.bmi_beta.get('ADIPOQ',np.nan),ax.bmi_beta.get('FASN',np.nan)))
tau=pd.read_csv(S('data/features/tissue_tau.tsv'),sep='\t').set_index('gene')['tau']
bsd=pd.read_csv(S('results/base_sd_transcriptome.tsv'),sep='\t')

# composition-shift proxy (adipocyte+macrophage fraction-Delta) for GSE141221, reused from the composition check
RAW=P('data/transcriptome/raw'); gm=json.load(open(os.path.join(RAW,'gse141221_meta.json')))
ic={(m['indiv'],m['cid']):(m['key'][:-1] if m['key'].endswith('b') else m['key']) for m in gm}
cn=pd.read_csv(gzip.open(os.path.join(RAW,'GSE141221_cleaned_count.csv.gz'),'rt'),index_col=0); csz=set(cn.columns)
pairs=[(ic[(iv,'CID1')],ic[(iv,'CID2')]) for iv in sorted({m['indiv'] for m in gm}) if ic.get((iv,'CID1')) in csz and ic.get((iv,'CID2')) in csz]
pre=[p[0] for p in pairs]; post=[p[1] for p in pairs]; allc=pre+post
lc=np.log2(cn[allc].divide(cn[allc].sum(),axis=1)*1e6+1.0); e2s=json.load(open(os.path.join(RAW,'gencode_v25_ens2sym.json')))
lc.index=[e2s.get(i.split('.')[0],i) for i in lc.index]; lc=lc.groupby(level=0).mean()
def sc(gs): g=[x for x in gs if x in lc.index]; z=lc.loc[g].apply(lambda r:(r-r.mean())/r.std(),axis=1); return z.mean(0)
fa=sc(['ADIPOQ','LEP','PLIN1','FABP4','CFD','LPL']); fm=sc(['CD68','CD14','LYZ','CD163'])
dA=fa[post].values-fa[pre].values; dM=fm[post].values-fm[pre].values
Dm=lc[post].values-lc[pre].values; X=np.column_stack([np.ones(len(dA)),dA,dM]); slp=np.linalg.lstsq(X,Dm.T,rcond=None)[0]
comp_adj=pd.Series(Dm.mean(1)-(slp[1]*dA.mean()+slp[2]*dM.mean()),index=lc.index)   # composition-adjusted rev_beta

def adipose_panel(panel,path):
    d=pd.read_csv(P(path),sep='\t')[['mark_id','rev_beta']].dropna().drop_duplicates('mark_id')
    b=bsd[bsd.panel==panel].set_index('mark_id'); d['bsd']=d.mark_id.map(b['base_sd']); d['ab']=d.mark_id.map(b['base_mean'])
    d['bmi']=d.mark_id.map(ax['bmi_beta']); d['tau']=d.mark_id.map(tau); d=d.dropna(subset=['bsd','ab','bmi','tau'])
    sc0=np.median(np.abs(d.rev_beta)) or 1.0
    d['restore']= -np.sign(d.bmi)*(d.rev_beta/sc0)          # >0 = toward lean
    d['mag']=np.abs(d.rev_beta)
    return d

print('\n=== ADIPOSE transcriptome RESTORATION (toward-lean; tissue-matched adipose BMI axis) ===')
for panel,path,nm in [('adipose_table','data/transcriptome/adipose_reversal_table.tsv','adipose-CR'),
                      ('gse141221_adipose_table','data/transcriptome/gse141221_adipose_reversal_table.tsv','adipose-LCD')]:
    d=adipose_panel(panel,path); frac=(d.restore>0).mean()
    print(f'  {nm:11s}: n={len(d):6d}  median restoration={d.restore.median():+.3f}  fraction toward-lean={frac:.2f} (>0.5=net restorative)')

print('\n=== DECISIVE tau TEST (adipose-LCD GSE141221): does tau predict toward-lean DIRECTION, or only MAGNITUDE? ===')
d=adipose_panel('gse141221_adipose_table','data/transcriptome/gse141221_adipose_reversal_table.tsv')
d['comp_adj']=d.mark_id.map(comp_adj)
# measurability-residualize each outcome within base_sd x abundance cells
d['sb']=pd.qcut(d.bsd.rank(method='first'),20,labels=False,duplicates='drop'); d['mb']=pd.qcut(d.ab.rank(method='first'),5,labels=False,duplicates='drop')
def cellresid(col):
    return d.groupby([d.sb,d.mb])[col].transform(lambda s: rint(s.values) if len(s)>5 else np.nan)
d['y_mag']=cellresid('mag')      # reversal MAGNITUDE (dynamic-range outcome)
d['y_dir']=cellresid('restore')  # toward-lean DIRECTION (restoration outcome)
# composition-adjusted restoration direction
d['restore_compadj']=-np.sign(d.bmi)*(d.comp_adj/(np.median(np.abs(d.comp_adj)) or 1.0))
d['y_dir_ca']=cellresid('restore_compadj')
def beta_tau(y):
    dd=d.dropna(subset=[y,'tau']); z=(dd.tau-dd.tau.mean())/dd.tau.std()
    r=sm.OLS(dd[y].values,sm.add_constant(pd.DataFrame({'tau':z.values}))).fit(cov_type='HC3')
    return r.params['tau'],r.pvalues['tau'],len(dd)
_rows=[]
for nm,y in [('reversal MAGNITUDE |Δ|','y_mag'),('toward-lean RESTORATION','y_dir'),('RESTORATION (composition-adj)','y_dir_ca')]:
    b,p,n=beta_tau(y); print(f'  tau -> {nm:32s}: beta={b:+.3f} p={p:.1e} (n={n})')
    _rows.append(dict(outcome=nm,col=y,tau_beta=b,p=p,n=n))
_op=S('results/orphan/tau_dissociation.tsv'); pd.DataFrame(_rows).to_csv(_op,sep='\t',index=False)
print(f'wrote {_op}')
print('\nINTERPRETATION: tau predicting MAGNITUDE but NOT direction => dynamic-range effect, NOT a restoration determinant.')
print('tau predicting RESTORATION (esp. composition-adjusted) => a real toward-health determinant.')

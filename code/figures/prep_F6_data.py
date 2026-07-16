#!/usr/bin/env python
"""prep_F6_data.py — source data for Figure 6 ("Reversibility is a property of context, not the molecule").
Reuses build_cross_context's 4 powered contexts. Writes:
  F6a_context_corr.csv     — 4x4 rev_beta Spearman matrix (tissue block structure)
  F6b_concordance.csv      — pairwise directional concordance, tagged same-/cross-tissue
  F6c_exemplars.csv        — OBJECTIVELY-selected context-specific genes x 4 contexts (the demonstration)
No cherry-picking: exemplars = top genes by context-specificity score (max|rev|-median|rev|), measured in >=3 contexts."""
import os, itertools, numpy as np, pandas as pd
from scipy import stats
HERE=os.path.dirname(os.path.abspath(__file__)); SIG=os.path.dirname(os.path.dirname(HERE)); PROJ=os.path.dirname(SIG)
OUT=os.path.join(HERE,'source_data'); os.makedirs(OUT,exist_ok=True)
PAN={'liver':'data/transcriptome/liver_reversal_table.tsv','adipose-CR':'data/transcriptome/adipose_reversal_table.tsv',
     'adipose-diet':'data/transcriptome/gse141221_adipose_reversal_table.tsv','blood':'data/transcriptome/gse273902_blood_reversal_table.tsv'}
TIS={'liver':'liver','adipose-CR':'adipose','adipose-diet':'adipose','blood':'blood'}
D={};SIG_={}
for k,f in PAN.items():
    d=pd.read_csv(os.path.join(PROJ,f),sep='\t')[['mark_id','rev_beta','rev_q']].dropna()
    D[k]=d.set_index('mark_id')['rev_beta']; SIG_[k]=set(d[d.rev_q<0.05]['mark_id'])
ctxs=list(PAN)

# a) rev_beta Spearman matrix. n_measured_both = the shared MEASURED universe each rho is computed on
# (the denominator; round-2 review asked for it explicitly alongside the shared SIGNIFICANT count in b).
rows=[]
for a in ctxs:
    for b in ctxs:
        sh=list(set(D[a].index)&set(D[b].index))
        rows.append({'ctx_a':a,'ctx_b':b,'rho':stats.spearmanr(D[a][sh],D[b][sh]).correlation,
                     'n_measured_both':len(sh)})
pd.DataFrame(rows).to_csv(os.path.join(OUT,'F6a_context_corr.csv'),index=False)

# b) directional concordance, same- vs cross-tissue. n_shared = shared SIGNIFICANT genes (the concordance
# denominator); n_measured_both = shared MEASURED genes (the rho denominator) — the two are different.
crows=[]
for a,b in itertools.combinations(ctxs,2):
    sh=list(SIG_[a]&SIG_[b])
    if sh:
        conc=float((np.sign(D[a][sh])==np.sign(D[b][sh])).mean())
        crows.append({'pair':f'{a} × {b}','n_shared':len(sh),'concordance':conc,
                      'kind':'same tissue' if TIS[a]==TIS[b] else 'cross tissue',
                      'n_measured_both':len(set(D[a].index)&set(D[b].index))})
pd.DataFrame(crows).to_csv(os.path.join(OUT,'F6b_concordance.csv'),index=False)

# c) exemplar genes: objective context-specificity selection (measured in >=3 contexts).
# EXCLUDE immunoglobulin + hemoglobin genes: in whole blood these are leukocyte/erythroid CELL-COMPOSITION
# markers, not tissue-intrinsic transcriptional reversal (project mandates composition control).
import re as _re
COMP=_re.compile(r'^IG[HKL]|^JCHAIN$|^HB[ABDGEZMQ]\d?$|^ALAS2$|^HBB$|^HBA')
allg=[g for g in set.union(*[set(D[k].index) for k in ctxs]) if not COMP.match(str(g))]
M=pd.DataFrame({k:D[k].reindex(allg) for k in ctxs})
n_meas=M.notna().sum(1)
M3=M[n_meas>=3].copy()
amax=M3.abs().max(1); amed=M3.abs().median(1)
M3['spec']=amax-amed; M3['home']=M3[ctxs].abs().idxmax(1).map(TIS)
M3['flips']=(M3[ctxs]>0).sum(1).between(1,n_meas[M3.index]-1) & (M3[ctxs].abs().max(1)>0.3)
# pick ~3 per home tissue by specificity + guarantee 2 sign-flippers
sel=[]
for h in ['liver','adipose','blood']:
    sel+=list(M3[M3.home==h].sort_values('spec',ascending=False).head(4).index)
flip=list(M3[M3.flips & (M3.spec>0.4)].sort_values('spec',ascending=False).head(3).index)
sel=list(dict.fromkeys(sel+flip))
ex=M3.loc[sel,ctxs+['home']].reset_index()
ex=ex.rename(columns={ex.columns[0]:'gene'})
ex.to_csv(os.path.join(OUT,'F6c_exemplars.csv'),index=False)

# d/e) discrete sig-set overlap (reversible & persistent), within-tissue vs cross-tissue (companions to the continuous matrix a)
from scipy.stats import hypergeom
bsd=pd.read_csv(os.path.join(SIG,'results','base_sd_transcriptome.tsv'),sep='\t')
bp={'liver':'liver_table','adipose-CR':'adipose_table','adipose-diet':'gse141221_adipose_table','blood':'gse273902_blood_table'}
TIS2={'liver':'liver','adipose-CR':'adipose','adipose-diet':'adipose','blood':'blood'}
PERS={}
for k,f in PAN.items():
    d=pd.read_csv(os.path.join(PROJ,f),sep='\t')[['mark_id','rev_beta']].dropna()
    b=bsd[bsd.panel==bp[k]].set_index('mark_id')['base_sd']
    d['bsd']=d.mark_id.map(b); d=d.dropna(subset=['bsd']); d['za']=d.rev_beta.abs()
    d['dec']=pd.qcut(d.bsd.rank(method='first'),20,labels=False,duplicates='drop')
    PERS[k]=set().union(*[set(g.loc[g.za<=g.za.quantile(1/3),'mark_id']) for _,g in d.groupby('dec')])
def overlap(setmap):
    rows=[]
    for a,b in itertools.combinations(ctxs,2):
        uni=set(D[a].index)&set(D[b].index); A=setmap[a]&uni; B=setmap[b]&uni
        ov=len(A&B); exp=len(A)*len(B)/len(uni) if len(uni) else np.nan
        rows.append({'pair':f'{a.replace("adipose","adi")}×{b.replace("adipose","adi")}','fold':ov/exp if exp else np.nan,
                     'p':float(hypergeom.sf(ov-1,len(uni),len(A),len(B))),
                     'kind':'within tissue' if TIS2[a]==TIS2[b] else 'cross tissue'})
    return pd.DataFrame(rows)
overlap(SIG_).to_csv(os.path.join(OUT,'F6d_reversible_overlap.csv'),index=False)
overlap(PERS).to_csv(os.path.join(OUT,'F6e_persistent_overlap.csv'),index=False)
print('F6 source data written.')
print(' matrix within-adipose rho=%.2f'%stats.spearmanr(D['adipose-CR'][list(set(D['adipose-CR'].index)&set(D['adipose-diet'].index))],D['adipose-diet'][list(set(D['adipose-CR'].index)&set(D['adipose-diet'].index))]).correlation)
print(' concordance:', [(r['pair'],round(r['concordance'],2),r['kind']) for r in crows])
print(' exemplars (%d):'%len(ex), list(ex.gene))

#!/usr/bin/env python
"""build_tissue_organization.py — the HARDENED tissue-as-governor test (Nature-Metabolism lever).
Now 4 tissues with MULTIPLE contexts each (adipose x8, muscle x3, liver x3, blood/immune x4), spanning
diet/CR/exercise/bariatric/drug/GLP-1RA. Tests: is weight-loss reversibility organized by TISSUE even
across very different interventions/cohorts within a tissue? Computes within-tissue vs cross-tissue mean
Spearman(rev_beta), per-tissue cohesion, and per-tissue consensus signatures. Excludes GSE117070
(batch-confounded signs, verifier-flagged). Output results/tissue_organization_summary.txt.
"""
import os, re, itertools, numpy as np, pandas as pd
from scipy import stats
SIG=os.path.dirname(os.path.dirname(os.path.abspath(__file__))); PROJ=os.path.dirname(SIG)
def P(*a): return os.path.join(PROJ,*a)
def S(*a): return os.path.join(SIG,*a)
NC=re.compile(r'-AS\d?|^LINC|orf\d|^RNU|^SNOR|^MIR\d|^LOC|^AC\d|^AL\d|^RP[0-9LS]|^IG[HKL]V')
# (label, tissue, intervention, file)  — GSE117070 EXCLUDED (batch-confounded directions)
PAN=[('liver-GSE83452','liver','diet+bariatric','liver_reversal_table'),
     ('liver-GSE48452','liver','bariatric','gse48452_liver_reversal_table'),
     ('liver-GSE106737','liver','bariatric','gse106737_liver_bariatric_reversal_table'),
     ('adi-CR','adipose','CR','adipose_reversal_table'),
     ('adi-LCD','adipose','diet-LCD','gse141221_adipose_reversal_table'),
     ('adi-dietEx','adipose','diet+exercise','gse43471_adipose_reversal_table'),
     ('adi-RYGB','adipose','bariatric','gse199063_adipose_reversal_table'),
     ('adi-diet2','adipose','diet-LCD/VLCD','gse77962_adipose_diet_reversal_table'),
     ('adi-dose','adipose','diet-graded','gse70529_adipose_dose_reversal_table'),
     ('adi-metformin','adipose','drug-metformin','gse107894_adipose_metformin_reversal_table'),
     ('adi-endur','adipose','exercise','gse224310_adipose_exercise_reversal_table'),
     ('mus-mtfRT','muscle','metformin+resistance','gse157585_muscle_reversal_table'),
     ('mus-exer','muscle','exercise','gse83352_muscle_exercise_reversal_table'),
     ('mus-bari','muscle','bariatric','gse161643_muscle_bariatric_reversal_table'),
     ('bld-bari','blood','bariatric','gse273902_blood_reversal_table'),
     ('bld-glp1','blood','GLP1RA','gse310742_blood_glp1ra_reversal_table'),
     ('bld-exer','blood','exercise','gse193771_blood_exercise_reversal_table'),
     ('bld-diet','blood','diet-Med','gse28358_pbmc_diet_reversal_table')]
D={}; TIS={}; SG={}
for lab,tis,iv,f in PAN:
    p=P('data/transcriptome',f+'.tsv')
    if not os.path.exists(p): print('  (skip',lab,'- missing)'); continue
    d=pd.read_csv(p,sep='\t').drop_duplicates('mark_id'); D[lab]=d.set_index('mark_id')['rev_beta']; TIS[lab]=tis
    s=d[pd.to_numeric(d.get('rev_q',1),errors='coerce')<0.05]; s=s[~s.mark_id.str.contains(NC,na=False)]
    SG[lab]=list(s.reindex(s.rev_beta.abs().sort_values(ascending=False).index).mark_id.head(20))
labs=list(D)
def rho(a,b):
    sh=list(set(D[a].index)&set(D[b].index)); return stats.spearmanr(D[a][sh],D[b][sh]).correlation if len(sh)>50 else np.nan
out=[]
def pr(s): out.append(s); print(s)
pr(f'Tissue-organization test: {len(labs)} transcriptome contexts across {len(set(TIS.values()))} tissues')
for t in ['liver','adipose','muscle','blood']:
    pr(f'  {t}: {sorted(l for l in labs if TIS[l]==t)}')
wt=[rho(a,b) for a,b in itertools.combinations(labs,2) if TIS[a]==TIS[b]]
xt=[rho(a,b) for a,b in itertools.combinations(labs,2) if TIS[a]!=TIS[b]]
pr(f'\n★ MEAN within-tissue r = {np.nanmean(wt):+.3f} (n={np.sum(np.isfinite(wt))} pairs)  vs  cross-tissue r = {np.nanmean(xt):+.3f} (n={np.sum(np.isfinite(xt))})')
pr(f'  Welch t(within vs cross): t={stats.ttest_ind(np.array(wt)[np.isfinite(wt)],np.array(xt)[np.isfinite(xt)],equal_var=False).statistic:.1f}  (>0 = tissue organizes reversibility)')
pr('\nper-tissue cohesion (mean within-tissue r):')
for t in ['liver','adipose','muscle','blood']:
    ls=[l for l in labs if TIS[l]==t]; rr=[rho(a,b) for a,b in itertools.combinations(ls,2)]
    pr(f'  {t:8s} ({len(ls)} contexts): mean within-r = {np.nanmean(rr):+.3f}')
pr('\nper-tissue CONSENSUS signature (genes recurrently top-reversed across a tissue contexts):')
for t in ['liver','adipose','muscle','blood']:
    from collections import Counter
    c=Counter(g for l in labs if TIS[l]==t for g in SG[l]); cons=[g for g,n in c.most_common(10) if n>=2]
    pr(f'  {t:8s}: '+', '.join(cons[:8]))
open(S('results/tissue_organization_summary.txt'),'w').write('\n'.join(out))
print('\nwrote results/tissue_organization_summary.txt')
# per-pair export for Figure 2a (distribution of within- vs cross-tissue pair correlations)
rows=[]
for a,b in itertools.combinations(labs,2):
    r=rho(a,b)
    if np.isfinite(r):
        rows.append({'a':a,'b':b,'tis_a':TIS[a],'tis_b':TIS[b],
                     'kind':'within' if TIS[a]==TIS[b] else 'cross',
                     'tissue':TIS[a] if TIS[a]==TIS[b] else 'cross','rho':r})
pd.DataFrame(rows).to_csv(S('results/tissue_pair_correlations.tsv'),sep='\t',index=False)
print('wrote results/tissue_pair_correlations.tsv (%d pairs)'%len(rows))

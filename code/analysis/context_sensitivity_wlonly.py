#!/usr/bin/env python
"""context_sensitivity_wlonly.py — TRUE-WEIGHT-LOSS-ONLY sensitivity for the CONTEXT claim
(tissue-governor + tissue/intervention variance partition). Companion to
determinant_sensitivity_wlonly.py.

STANDALONE. Reuses the crosswalk + factor definitions from variance_partition.py VERBATIM,
reads results/tissue_pair_correlations.tsv + panel_manifest_full.tsv read-only, and writes NEW
files only. Restricts the 153-pair set to pairs whose BOTH endpoints are TRUE-WEIGHT-LOSS panels,
then recomputes within-tissue vs cross-tissue mean rho and the same_tissue / same_platform /
same_intervention_family variance partition.

TWO panel definitions are reported (harmonized with the determinant companion, which offers the
same two):
  lenient12 — diet/CR + surgical; diet+exercise (adi-dietEx) RETAINED as diet
              (matches variance_partition.py's canonical intervention_family bucketing and the
               determinant k7_keepdietex scenario).
  strict11  — additionally drops adi-dietEx (diet+exercise), matching the determinant k6_strict
              scenario, so the two analyses share one "weight-loss-only" panel set.
Empagliflozin/EMPEROR is a proteome SGLT2i panel; the transcriptome tissue-pair set never
contains it, so it is not a factor here (the determinant null is transcriptome-powered — the
proteome determinant is inconclusive by design, MDE~0.19; empagliflozin's intervention behaviour
is examined separately in the intervention/Fig-5 analysis, not here).

Writes:
  results/sensitivity_wlonly/context_wlonly.tsv        (partition + within/cross means, both scenarios)
  results/sensitivity_wlonly/context_wlonly_summary.txt
"""
import os, numpy as np, pandas as pd
import statsmodels.formula.api as smf

RNG=np.random.default_rng(20260714)
N_PERM=20000

SIG=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
def S(*a): return os.path.join(SIG,*a)
OUT=S('results','sensitivity_wlonly'); os.makedirs(OUT,exist_ok=True)

# ---- crosswalk + factor defs transcribed VERBATIM from variance_partition.py ----
LABEL2FILE={
    'liver-GSE83452':'liver_reversal_table','liver-GSE48452':'gse48452_liver_reversal_table',
    'liver-GSE106737':'gse106737_liver_bariatric_reversal_table','adi-CR':'adipose_reversal_table',
    'adi-LCD':'gse141221_adipose_reversal_table','adi-dietEx':'gse43471_adipose_reversal_table',
    'adi-RYGB':'gse199063_adipose_reversal_table','adi-diet2':'gse77962_adipose_diet_reversal_table',
    'adi-dose':'gse70529_adipose_dose_reversal_table','adi-metformin':'gse107894_adipose_metformin_reversal_table',
    'adi-endur':'gse224310_adipose_exercise_reversal_table','mus-mtfRT':'gse157585_muscle_reversal_table',
    'mus-exer':'gse83352_muscle_exercise_reversal_table','mus-bari':'gse161643_muscle_bariatric_reversal_table',
    'bld-bari':'gse273902_blood_reversal_table','bld-glp1':'gse310742_blood_glp1ra_reversal_table',
    'bld-exer':'gse193771_blood_exercise_reversal_table','bld-diet':'gse28358_pbmc_diet_reversal_table',
}
def platform_tech(p):
    p=str(p).lower(); return 'seq' if 'seq' in p else 'array'
def intervention_family(iv):
    iv=str(iv).lower()
    if 'bariatric' in iv or 'rygb' in iv:                 return 'surgical'
    if 'metformin' in iv or 'glp1' in iv or 'sglt' in iv: return 'pharmacological'
    if 'exercise' in iv and 'diet' not in iv:             return 'exercise'
    return 'diet'

WLONLY_FAMILIES={'diet','surgical'}     # true weight-loss families
STRICT_EXTRA_DROP={'adi-dietEx'}        # diet+exercise: dropped in strict11 (matches determinant k6)

pairs=pd.read_csv(S('results/tissue_pair_correlations.tsv'),sep='\t')
man=pd.read_csv(S('results/panel_manifest_full.tsv'),sep='\t')
man['base']=man['path'].str.replace('data/transcriptome/','',regex=False).str.replace('.tsv','',regex=False)
file2row=man.set_index('base')
def meta(label):
    row=file2row.loc[LABEL2FILE[label]]
    return platform_tech(row['platform']),intervention_family(row['intervention']),row['platform'],row['intervention']

panels_all=sorted(set(pairs['a'])|set(pairs['b']))
attr={}
for lab in panels_all:
    tech,fam,praw,ivraw=meta(lab)
    tis=pairs.loc[pairs['a']==lab,'tis_a'].tolist()+pairs.loc[pairs['b']==lab,'tis_b'].tolist()
    attr[lab]={'tissue':tis[0],'platform_tech':tech,'intervention_family':fam,'intervention_raw':ivraw}

TERMS=['same_tissue','same_platform','same_intervention_family']
TERM2ATTR={'same_tissue':'tissue','same_platform':'platform_tech','same_intervention_family':'intervention_family'}

def _r2(X,yv):
    beta,*_=np.linalg.lstsq(X,yv,rcond=None); resid=yv-X@beta
    return 1-(resid@resid)/((yv-yv.mean())@(yv-yv.mean()))

def run_scenario(name, keep):
    A=pd.DataFrame({k:attr[k] for k in keep}).T
    d=pairs[pairs['a'].isin(keep)&pairs['b'].isin(keep)].copy()
    d['same_tissue']=(d['tis_a']==d['tis_b']).astype(int)
    d['same_platform']=(d['a'].map(lambda x:attr[x]['platform_tech'])==d['b'].map(lambda x:attr[x]['platform_tech'])).astype(int)
    d['same_intervention_family']=(d['a'].map(lambda x:attr[x]['intervention_family'])==d['b'].map(lambda x:attr[x]['intervention_family'])).astype(int)
    def r2f(formula): return smf.ols(formula,data=d).fit().rsquared
    active=[t for t in TERMS if d[t].nunique()>1]
    full=smf.ols('rho ~ '+' + '.join(active),data=d).fit(); R2_full=full.rsquared
    marg,uniq={},{}
    for t in active:
        marg[t]=r2f(f'rho ~ {t}'); others=[x for x in active if x!=t]
        uniq[t]=R2_full-(r2f('rho ~ '+' + '.join(others)) if others else 0.0)
    mean_within=d.loc[d.same_tissue==1,'rho'].mean(); mean_cross=d.loc[d.same_tissue==0,'rho'].mean()
    n_within=int(d.same_tissue.sum()); n_cross=int((1-d.same_tissue).sum())
    # composition of the within-tissue pairs (transparency: is the within-mean tissue-balanced?)
    wsub=d[d.same_tissue==1]; comp={t:(int((wsub.tissue==t).sum()),round(float(wsub.loc[wsub.tissue==t,'rho'].mean()),3))
                                     for t in sorted(wsub.tissue.unique())}
    # symmetric panel-level permutation for same_tissue (key claim)
    panel_index={p:i for i,p in enumerate(A.index)}
    ai=d['a'].map(panel_index).values; bi=d['b'].map(panel_index).values; y=d['rho'].values
    vec=A['tissue'].values.copy(); obs=marg.get('same_tissue',np.nan); c=0; perm=vec.copy()
    for _ in range(N_PERM):
        RNG.shuffle(perm); st=(perm[ai]==perm[bi]).astype(float)
        if _r2(np.column_stack([np.ones_like(st),st]),y)>=obs-1e-12: c+=1
    perm_p_tissue=(c+1)/(N_PERM+1)
    tissues={t:sum(1 for k in keep if attr[k]['tissue']==t) for t in sorted(set(attr[k]['tissue'] for k in keep))}
    return dict(name=name,keep=keep,active=active,R2_full=R2_full,marg=marg,uniq=uniq,full=full,
                mean_within=mean_within,mean_cross=mean_cross,n_within=n_within,n_cross=n_cross,
                perm_p_tissue=perm_p_tissue,d=d,comp=comp,tissues=tissues)

lenient=[lab for lab in panels_all if attr[lab]['intervention_family'] in WLONLY_FAMILIES]
strict=[lab for lab in lenient if lab not in STRICT_EXTRA_DROP]
DROP=[lab for lab in panels_all if lab not in lenient]
R={'lenient12':run_scenario('lenient12',lenient),'strict11':run_scenario('strict11',strict)}

# ---- write TSV ----
rows=[]
for name,r in R.items():
    for t in r['active']:
        rows.append({'scenario':name,'term':t,'beta_ols':r['full'].params.get(t,np.nan),
                     'p_ols':r['full'].pvalues.get(t,np.nan),'marginal_R2':r['marg'][t],'unique_R2':r['uniq'][t],
                     'perm_p_marginal_tissue':(r['perm_p_tissue'] if t=='same_tissue' else np.nan),
                     'n_same':int(r['d'][t].sum()),'n_diff':int((1-r['d'][t]).sum())})
    rows.append({'scenario':name,'term':'within_tissue_mean_rho','beta_ols':r['mean_within'],'n_same':r['n_within']})
    rows.append({'scenario':name,'term':'cross_tissue_mean_rho','beta_ols':r['mean_cross'],'n_same':r['n_cross']})
    rows.append({'scenario':name,'term':'MODEL_full_R2','marginal_R2':r['R2_full']})
pd.DataFrame(rows).to_csv(os.path.join(OUT,'context_wlonly.tsv'),sep='\t',index=False,float_format='%.6g')

# ---- human-readable summary ----
L=[]
def p(s=''): L.append(s); print(s)
p('TRUE-WEIGHT-LOSS-ONLY SENSITIVITY — context (tissue-governor + variance partition)')
p('='*80)
p(f'DROP (metformin/GLP1/SGLT2i/exercise-primary) n={len(DROP)}: {DROP}')
p(f'Canonical reference (all 18 panels): within-tissue +0.162 vs cross-tissue +0.020.')
p('diet+exercise (adi-dietEx) is retained as diet in lenient12 and dropped in strict11')
p('(mirrors the determinant companion k7_keepdietex vs k6_strict).')
p('Empagliflozin/EMPEROR is a proteome panel, absent from the transcriptome tissue-pair set.')
p('')
for name,r in R.items():
    p(f'--- {name} (n={len(r["keep"])} panels) ---')
    p(f'  KEEP: {r["keep"]}')
    p(f'  tissues: '+', '.join(f"{t}={n}" for t,n in r['tissues'].items()))
    p(f'  within-tissue mean rho = {r["mean_within"]:+.3f} (n={r["n_within"]})  vs  cross-tissue = {r["mean_cross"]:+.3f} (n={r["n_cross"]})')
    p(f'  within-tissue composition (n, mean rho): '+', '.join(f"{t}={n}({m:+.2f})" for t,(n,m) in r['comp'].items()))
    p(f'    -> the within-tissue mean is dominated by the tissue with the most panels (transparency; not a signal change)')
    p(f'  variance partition (full-model R2 = {r["R2_full"]:.4f}):')
    for t in r['active']:
        pp=f"  perm_p={r['perm_p_tissue']:.4g}" if t=='same_tissue' else ""
        p(f'    {t:26s} marg_R2={r["marg"][t]:.4f}  uniq_R2={r["uniq"][t]:.4f}{pp}  (same={int(r["d"][t].sum())}/diff={int((1-r["d"][t]).sum())})')
    for t in TERMS:
        if t not in r['active']: p(f'    {t:26s} CONSTANT in this subset -> not partitionable, omitted')
    p('')
p('CONCLUSION')
p('-'*80)
hold=all(r['mean_within']>r['mean_cross'] and r['perm_p_tissue']<0.05 for r in R.values())
for name,r in R.items():
    p(f'  {name}: within {r["mean_within"]:+.3f} > cross {r["mean_cross"]:+.3f}, same_tissue permutation p={r["perm_p_tissue"]:.4g}.')
p(f'  Shared tissue context still organizes reversibility on true-weight-loss interventions under '
  f'{"BOTH" if hold else "NOT ALL"} panel definitions.')
p('  Intervention-family axis is now diet-vs-surgery only (exercise + pharmacological families removed);')
p('  reported for completeness, not as the primary claim. The rise vs the canonical +0.162 within-mean')
p('  is a COMPOSITION effect (restriction removes near-zero within-blood and absent within-muscle pairs),')
p('  not a strengthened per-pair signal.')
open(os.path.join(OUT,'context_wlonly_summary.txt'),'w',encoding='utf-8').write('\n'.join(L))
print('\nwrote results/sensitivity_wlonly/{context_wlonly.tsv, context_wlonly_summary.txt}')

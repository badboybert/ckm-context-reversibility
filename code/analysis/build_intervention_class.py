#!/usr/bin/env python
"""build_intervention_class.py — tests the user's hypothesis: within ONE tissue (adipose), does the
reversibility signature differ by INTERVENTION CLASS (treatment[surgery/drug] vs lifestyle[diet/exercise])?
Tissue held fixed = SAT adipose. Groups all adipose panels by class, computes the pairwise rev_beta
Spearman matrix, contrasts WITHIN-class vs BETWEEN-class correlation, and lists per-class signatures.
Runs on whatever adipose panels are present (skips missing). CAVEAT: class is partly confounded with
timepoint/cohort/platform across panels — interpret as suggestive structure, not a clean experiment.
"""
import os, re, itertools, numpy as np, pandas as pd
from scipy import stats
SIG=os.path.dirname(os.path.dirname(os.path.abspath(__file__))); PROJ=os.path.dirname(SIG)
def P(*a): return os.path.join(PROJ,*a)
def S(*a): return os.path.join(SIG,*a)
NC=re.compile(r'-AS\d?|^LINC|orf\d|^RNU|^SNOR|^MIR\d|^LOC|^AC\d|^AL\d|^RP[0-9LS]')
# adipose panels: (label, class, file)
PANELS=[('CR','lifestyle-diet','adipose_reversal_table'),
        ('LCD','lifestyle-diet','gse141221_adipose_reversal_table'),
        ('diet+ex','lifestyle-mixed','gse43471_adipose_reversal_table'),
        ('endurance','lifestyle-exercise','gse224310_adipose_exercise_reversal_table'),
        ('metformin','treatment-drug','gse107894_adipose_metformin_reversal_table'),
        ('RYGB','treatment-surgery','gse199063_adipose_reversal_table')]
D={}; CLS={}; SIGS={}
for lab,cls,f in PANELS:
    path=P('data/transcriptome',f+'.tsv')
    if not os.path.exists(path): print(f'  (skip {lab}: not present)'); continue
    d=pd.read_csv(path,sep='\t').drop_duplicates('mark_id')
    D[lab]=d.set_index('mark_id')['rev_beta']; CLS[lab]=cls
    s=d[pd.to_numeric(d.get('rev_q',1),errors='coerce')<0.05]; s=s[~s.mark_id.str.contains(NC,na=False)]
    SIGS[lab]=list(s.reindex(s.rev_beta.abs().sort_values(ascending=False).index).mark_id.head(8))
labs=list(D)
print(f'\nadipose panels in test ({len(labs)}): '+', '.join(f'{l}[{CLS[l]}]' for l in labs))
print('\n=== pairwise rev_beta Spearman (tissue=adipose held fixed) ===')
print('             '+' '.join(f'{l[:9]:>9s}' for l in labs))
def rho(a,b): sh=list(set(D[a].index)&set(D[b].index)); return stats.spearmanr(D[a][sh],D[b][sh]).correlation if len(sh)>50 else np.nan
for a in labs: print(f'  {a:10s} '+' '.join(f'{rho(a,b):+9.2f}' for b in labs))
# within-class vs between-class
wc=[rho(a,b) for a,b in itertools.combinations(labs,2) if CLS[a]==CLS[b]]
bc=[rho(a,b) for a,b in itertools.combinations(labs,2) if CLS[a]!=CLS[b]]
# treatment vs lifestyle (broad)
def broad(c): return 'treatment' if c.startswith('treatment') else 'lifestyle'
wl=[rho(a,b) for a,b in itertools.combinations(labs,2) if broad(CLS[a])==broad(CLS[b])]
bl=[rho(a,b) for a,b in itertools.combinations(labs,2) if broad(CLS[a])!=broad(CLS[b])]
print(f'\n  MEAN within-CLASS r = {np.nanmean(wc):+.2f} (n={len(wc)})   vs between-CLASS r = {np.nanmean(bc):+.2f} (n={len(bc)})')
print(f'  MEAN within broad-group (treatment|lifestyle) r = {np.nanmean(wl):+.2f}   vs across-group r = {np.nanmean(bl):+.2f}')
print('  -> within > between supports intervention-class structuring adipose reversibility (beyond tissue).')
print('\n=== per-class top reversibility signatures ===')
for l in labs: print(f'  {l:10s} [{CLS[l]:18s}]: '+', '.join(SIGS[l]))
print('\nCAVEAT: class confounded with timepoint/cohort/platform; small n in drug/exercise arms -> suggestive.')

#!/usr/bin/env python
"""build_consensus_signatures.py — STRENGTHEN the reversible / persistent signatures at TISSUE and
INTERVENTION level using the expanded data (multiple contexts per tissue). A consensus signature =
genes recurrently top-reversed (or recurrently low-moving = persistent) across a tissue's independent
contexts (different interventions/cohorts/platforms) -> far more robust than single-panel. Reports
recurrence + direction-consistency. (Restorative/toward-lean is handled separately where a BMI axis
exists.) Excludes GSE117070 (batch-confounded). Output results/consensus_signatures.txt.
"""
import os, re, numpy as np, pandas as pd
from collections import Counter, defaultdict
SIG=os.path.dirname(os.path.dirname(os.path.abspath(__file__))); PROJ=os.path.dirname(SIG)
def P(*a): return os.path.join(PROJ,*a)
def S(*a): return os.path.join(SIG,*a)
NC=re.compile(r'-AS\d?|^LINC|orf\d|^RNU|^SNOR|^MIR\d|^LOC|^AC\d|^AL\d|^RP[0-9LS]|^IG[HKL][VJC]')
PAN=[('liver','diet+bariatric','liver_reversal_table'),('liver','bariatric','gse48452_liver_reversal_table'),
     ('liver','bariatric','gse106737_liver_bariatric_reversal_table'),
     ('adipose','CR','adipose_reversal_table'),('adipose','diet','gse141221_adipose_reversal_table'),
     ('adipose','diet+ex','gse43471_adipose_reversal_table'),('adipose','bariatric','gse199063_adipose_reversal_table'),
     ('adipose','diet','gse77962_adipose_diet_reversal_table'),('adipose','diet','gse70529_adipose_dose_reversal_table'),
     ('adipose','drug','gse107894_adipose_metformin_reversal_table'),('adipose','exercise','gse224310_adipose_exercise_reversal_table'),
     ('muscle','drug+ex','gse157585_muscle_reversal_table'),('muscle','exercise','gse83352_muscle_exercise_reversal_table'),
     ('muscle','bariatric','gse161643_muscle_bariatric_reversal_table'),
     ('blood','bariatric','gse273902_blood_reversal_table'),('blood','GLP1RA','gse310742_blood_glp1ra_reversal_table'),
     ('blood','exercise','gse193771_blood_exercise_reversal_table'),('blood','diet','gse28358_pbmc_diet_reversal_table')]
def load(f):
    p=P('data/transcriptome',f+'.tsv')
    if not os.path.exists(p): return None
    d=pd.read_csv(p,sep='\t').drop_duplicates('mark_id'); d=d[~d.mark_id.str.contains(NC,na=False)]
    d['q']=pd.to_numeric(d.get('rev_q',1),errors='coerce'); d['ab']=d.rev_beta.abs()
    d['toprev']=(d['q']<0.05)&(d['ab']>=d['ab'].quantile(0.95))      # strongly+sig reversed in this context
    d['lowmov']=d['ab']<=d['ab'].quantile(1/3)                        # low-mover (persistent candidate) in this context
    return d.set_index('mark_id')
out=[]
def pr(s): out.append(s); print(s)
# group panels by tissue and by intervention-class
byT=defaultdict(list); byI=defaultdict(list); LD={}
for tis,iv,f in PAN:
    d=load(f)
    if d is None: continue
    LD[f]=d; byT[tis].append(f); byI[iv.split('+')[0]].append((tis,f))
pr('===== PER-TISSUE CONSENSUS SIGNATURES (recurrence across a tissue contexts) =====')
for tis in ['liver','adipose','muscle','blood']:
    fs=byT[tis]; k=len(fs); need=max(2,(k+1)//2)
    rev=Counter(); per=Counter(); dirn=defaultdict(list)
    for f in fs:
        d=LD[f]
        for g in d.index[d['toprev']]: rev[g]+=1; dirn[g].append(np.sign(d.loc[g,'rev_beta']))
        for g in d.index[d['lowmov']]: per[g]+=1
    cons_rev=[g for g,n in rev.most_common() if n>=need]
    cons_per=[g for g,n in per.most_common(40) if n>=need]
    pr(f'\n[{tis}] {k} contexts (recurrence threshold >={need}):')
    pr(f'  REVERSIBLE consensus ({len(cons_rev)} genes, top): '+', '.join(cons_rev[:12]))
    pr(f'  PERSISTENT consensus ({len(cons_per)} genes, low-movers in >={need} contexts, top): '+', '.join(cons_per[:12]))
pr('\n===== PER-INTERVENTION CONSENSUS (genes reversed by an intervention ACROSS tissues) =====')
iv_rows=[]
for iv in ['bariatric','diet','exercise','drug','GLP1RA','CR']:
    items=byI.get(iv,[]); tissues=set(t for t,_ in items)
    if len(items)<2: continue
    rev=Counter()
    for t,f in items:
        for g in LD[f].index[LD[f]['toprev']]: rev[g]+=1
    need=max(2,(len(items)+1)//2)
    cons=[g for g,n in rev.most_common() if n>=need]
    iv_rows.append({'intervention':iv,'n_panels':len(items),'n_tissues':len(tissues),
                    'n_systemic_consensus':len(cons),'top_genes':','.join(cons[:12])})
    pr(f'  {iv:10s} ({len(items)} panels, tissues={sorted(tissues)}): '+(', '.join(cons[:10]) if cons else '(no gene reversed in >=half the panels — intervention reversal is tissue-specific)'))
pd.DataFrame(iv_rows).to_csv(S('results/intervention_systemic_consensus.tsv'),sep='\t',index=False)
open(S('results/consensus_signatures.txt'),'w').write('\n'.join(out))
print('\nwrote results/consensus_signatures.txt')

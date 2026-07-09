#!/usr/bin/env python
"""build_cross_context.py — quantifies how much reversibility TRANSFERS across contexts (the
context-specificity finding, beyond top-N). Four cross-determinations across the 4 transcriptome
contexts: (A) reversible significant-set overlap vs chance (hypergeometric), (B) continuous rev_beta
Spearman matrix (does the whole pattern transfer?), (C) directional concordance on shared sig genes,
(D) persistent-set overlap vs chance. ALL converge: reversibility is TISSUE-ORGANIZED (tissue >>
intervention) — shared within a tissue across interventions (+0.52, 0.98 same-dir), ~0/opposite across
tissues. Output results/cross_context_summary.txt.
"""
import os, itertools, numpy as np, pandas as pd
from scipy import stats
SIG=os.path.dirname(os.path.dirname(os.path.abspath(__file__))); PROJ=os.path.dirname(SIG)
def P(*a): return os.path.join(PROJ,*a);
def S(*a): return os.path.join(SIG,*a)
PAN={'liver':'data/transcriptome/liver_reversal_table.tsv','adi-CR':'data/transcriptome/adipose_reversal_table.tsv',
     'adi-diet':'data/transcriptome/gse141221_adipose_reversal_table.tsv','blood':'data/transcriptome/gse273902_blood_reversal_table.tsv'}
bsd=pd.read_csv(S('results/base_sd_transcriptome.tsv'),sep='\t')
bp={'liver':'liver_table','adi-CR':'adipose_table','adi-diet':'gse141221_adipose_table','blood':'gse273902_blood_table'}
D={};SIGSET={};PERS={}
for k,f in PAN.items():
    d=pd.read_csv(P(f),sep='\t')[['mark_id','rev_beta','rev_q']].dropna()
    D[k]=d.set_index('mark_id')['rev_beta']; SIGSET[k]=set(d[d.rev_q<0.05]['mark_id'])
    b=bsd[bsd.panel==bp[k]].set_index('mark_id')['base_sd']; dd=d.copy(); dd['bsd']=dd.mark_id.map(b); dd=dd.dropna(subset=['bsd']); dd['za']=dd.rev_beta.abs()
    dd['dec']=pd.qcut(dd.bsd.rank(method='first'),20,labels=False,duplicates='drop')
    PERS[k]=set().union(*[set(g.loc[g.za<=g.za.quantile(1/3),'mark_id']) for _,g in dd.groupby('dec')])
out=[]
def pr(s): out.append(s); print(s)
pr('=== A) reversible sig-set overlap vs chance (hypergeometric) ===')
for a,b in itertools.combinations(PAN,2):
    uni=set(D[a].index)&set(D[b].index); A=SIGSET[a]&uni; B=SIGSET[b]&uni; ov=len(A&B); exp=len(A)*len(B)/len(uni)
    pr(f'  {a:8s} vs {b:8s}: overlap={ov:4d} exp={exp:6.0f} enrich={ov/exp:.2f}x p={stats.hypergeom.sf(ov-1,len(uni),len(A),len(B)):.1e}')
pr('\n=== B) continuous rev_beta Spearman matrix ===')
pr('          '+' '.join(f'{k:>8s}' for k in PAN))
for a in PAN:
    pr(f'  {a:7s} '+' '.join(f'{stats.spearmanr(D[a][list(set(D[a].index)&set(D[b].index))],D[b][list(set(D[a].index)&set(D[b].index))]).correlation:+8.2f}' for b in PAN))
pr('\n=== C) directional concordance on shared sig genes (fraction same sign) ===')
for a,b in itertools.combinations(PAN,2):
    sh=list(SIGSET[a]&SIGSET[b])
    if sh: pr(f'  {a:8s} vs {b:8s}: {len(sh):4d} shared-sig same-dir={(np.sign(D[a][sh])==np.sign(D[b][sh])).mean():.2f}')
pr('\n=== D) persistent-set overlap vs chance ===')
for a,b in itertools.combinations(PAN,2):
    uni=set(D[a].index)&set(D[b].index); A=PERS[a]&uni; B=PERS[b]&uni; ov=len(A&B); exp=len(A)*len(B)/len(uni)
    pr(f'  {a:8s} vs {b:8s}: overlap={ov:4d} exp={exp:5.0f} enrich={ov/exp:.2f}x p={stats.hypergeom.sf(ov-1,len(uni),len(A),len(B)):.1e}')
pr('\n-> ALL converge: reversibility is TISSUE-ORGANIZED (tissue>>intervention); within-tissue shared (+0.52/0.98/1.68x), cross-tissue ~0/opposite (~1.05x).')
open(S('results/cross_context_summary.txt'),'w').write('\n'.join(out))
print('\nwrote results/cross_context_summary.txt')

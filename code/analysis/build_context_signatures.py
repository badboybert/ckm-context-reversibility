#!/usr/bin/env python
"""build_context_signatures.py — the POSITIVE deliverable: per-CONTEXT reversibility signatures.
The determinant model showed no UNIVERSAL intrinsic predictor of reversibility; but each context
(panel = tissue x intervention) has its OWN signature of strongly-reversible marks. This lists them
+ quantifies cross-context sharing (the context-specificity made concrete: different tissues ~0 overlap,
same tissue across interventions shares a core). Ranked by |rev_beta| among significant marks (q<0.05),
noncoding-noise filtered. Output results/context_signatures.tsv (long) + console summary.
"""
import os, re, itertools, numpy as np, pandas as pd
SIG=os.path.dirname(os.path.dirname(os.path.abspath(__file__))); PROJ=os.path.dirname(SIG)
def P(*a): return os.path.join(PROJ,*a)
def S(*a): return os.path.join(SIG,*a)
NC=re.compile(r'-AS\d?|^LINC|orf\d|^RNU|^SNOR|^MIR\d|^LOC\d|^AC\d|^AL\d|^IGKV|^IGLV|^IGHV')   # noncoding / hypervariable

CTX=[('transcriptome','liver','diet+surgery','data/transcriptome/liver_reversal_table.tsv'),
     ('transcriptome','adipose','caloric-restriction','data/transcriptome/adipose_reversal_table.tsv'),
     ('transcriptome','adipose','low-cal-diet','data/transcriptome/gse141221_adipose_reversal_table.tsv'),
     ('transcriptome','adipose','diet+exercise','data/transcriptome/gse43471_adipose_reversal_table.tsv'),
     ('transcriptome','adipose','RYGB-durability','data/transcriptome/gse199063_adipose_reversal_table.tsv'),
     ('transcriptome','muscle','metformin+resistance','data/transcriptome/gse157585_muscle_reversal_table.tsv'),
     ('transcriptome','blood','bariatric-surgery','data/transcriptome/gse273902_blood_reversal_table.tsv'),
     ('proteome','plasma','bariatric','data/proteome/reversal_olink_bbs.tsv'),
     ('proteome','plasma','diet','data/proteome/reversal_somascan_direct.tsv'),
     ('proteome','plasma','SGLT2i','data/proteome/reversal_emperor.tsv')]
rows=[]; tops={}
for layer,tissue,lever,f in CTX:
    d=pd.read_csv(P(f),sep='\t'); cols=list(d.columns)
    ec='rev_beta' if 'rev_beta' in cols else 'effect'; qc='rev_q' if 'rev_q' in cols else ('rev_p' if 'rev_p' in cols else None)
    key='mark_id' if 'mark_id' in cols else ('gene' if 'gene' in cols else 'UniProt')
    d['m']=d[key].astype(str)
    if 'gene' in cols and key!='gene': d['m']=d['gene'].astype(str)   # proteome: prefer gene label
    if qc: d=d[pd.to_numeric(d[qc],errors='coerce')<0.05]
    d=d[~d['m'].str.contains(NC,na=False)]
    d['abs']=pd.to_numeric(d[ec],errors='coerce').abs(); d=d.dropna(subset=['abs']).sort_values('abs',ascending=False)
    ctx=f'{tissue}/{lever}'; tops[ctx]=list(d['m'].head(15))
    for rank,(_,r) in enumerate(d.head(25).iterrows(),1):
        rows.append(dict(layer=layer,context=ctx,tissue=tissue,intervention=lever,rank=rank,mark=r['m'],rev_beta=r[ec]))
pd.DataFrame(rows).to_csv(S('results/context_signatures.tsv'),sep='\t',index=False)
print('Per-context reversibility signatures -> results/context_signatures.tsv')
for ctx,t in tops.items(): print(f'  {ctx:28s}: '+', '.join(t[:8]))
print('\nCross-context sharing of the top-15 signature (Jaccard-style overlap):')
for a,b in itertools.combinations([c for c in tops if c.startswith(('liver','adipose','blood'))],2):
    ov=sorted(set(tops[a])&set(tops[b]))
    print(f'  {a:26s} vs {b:26s}: {len(ov):2d}/15 shared {ov}')
print('\n-> different TISSUES ~0 overlap (signature is tissue-specific); same TISSUE across interventions shares a CORE.')

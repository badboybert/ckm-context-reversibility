#!/usr/bin/env python
"""add_base_sd_muscle.py — compute baseline variance for GSE83352 (skeletal muscle, exercise, n=42)
from its "Pre Training" samples, and APPEND to results/base_sd_transcriptome.tsv so muscle enters the
determinant meta (the new-tissue hardening lever). Probe->gene via GPL10558.annot, max-variance collapse.
"""
import os, gzip, re, numpy as np, pandas as pd
SIG=os.path.dirname(os.path.dirname(os.path.abspath(__file__))); PROJ=os.path.dirname(SIG)
def P(*a): return os.path.join(PROJ,*a)
def S(*a): return os.path.join(SIG,*a)
RAW=P('data/transcriptome/raw'); MX=os.path.join(RAW,'GSE83352_series_matrix.txt.gz'); ANN=os.path.join(RAW,'GPL10558.annot.gz')
# parse series matrix: source_name -> Pre/Post + subject ; keep Pre (baseline) sample columns
src=None; gsm=None; rows=[]; hdr=None; intab=False
with gzip.open(MX,'rt',errors='replace') as fh:
    for l in fh:
        if l.startswith('!Sample_geo_accession'): gsm=[x.strip().strip('"') for x in l.rstrip().split('\t')[1:]]
        elif l.startswith('!Sample_source_name_ch1') or (l.startswith('!Sample_title') and src is None):
            src=[x.strip().strip('"') for x in l.rstrip().split('\t')[1:]]
        elif l.startswith('!series_matrix_table_begin'): intab=True; continue
        elif l.startswith('!series_matrix_table_end'): break
        elif intab:
            parts=[x.strip().strip('"') for x in l.rstrip('\n').split('\t')]
            if hdr is None: hdr=parts; continue
            rows.append(parts)
pre_idx=[i for i,s in enumerate(src) if re.search(r'pre',str(s),re.I) and not re.search(r'post',str(s),re.I)]
pre_gsm=set(gsm[i] for i in pre_idx)
mat=pd.DataFrame(rows).set_index(0); mat.columns=hdr[1:]
pre_cols=[c for c in mat.columns if c in pre_gsm]
e=mat[pre_cols].apply(pd.to_numeric,errors='coerce')
psd=e.std(axis=1,ddof=1); pmean=e.mean(axis=1)
# probe -> symbol (GPL10558.annot: GEO SOFT/annot — find the real header row robustly)
import io
lines=gzip.open(ANN,'rt',errors='replace').read().split('\n')
hi=next((i for i,l in enumerate(lines) if (l.startswith('ID\t') or (l.split('\t')[0].strip()=='ID')) ), None)
if hi is None: hi=next(i for i,l in enumerate(lines) if 'platform_table_begin' in l)+1
end=next((i for i,l in enumerate(lines) if 'platform_table_end' in l), len(lines))
ann=pd.read_csv(io.StringIO('\n'.join(lines[hi:end])),sep='\t',dtype=str,low_memory=False,on_bad_lines='skip')
idc=ann.columns[0]; symc=([c for c in ann.columns if c.lower() in ('gene symbol','symbol','ilmn_gene','genesymbol')] or [c for c in ann.columns if 'symbol' in c.lower()])[0]
p2s=dict(zip(ann[idc],ann[symc]))
df=pd.DataFrame({'probe':e.index.astype(str),'base_sd':psd.values,'base_mean':pmean.values})
df['mark_id']=df.probe.map(p2s); df=df.dropna(subset=['mark_id','base_sd'])
df['var']=df.base_sd**2; df=df.sort_values('var',ascending=False).drop_duplicates('mark_id')
out=df[['mark_id','base_sd','base_mean']].copy(); out['panel']='gse83352_muscle_exercise_table'; out['n_base']=len(pre_cols)
print(f'GSE83352 muscle: {len(pre_cols)} Pre-Training baseline samples; base_sd for {len(out)} genes')
# append to base_sd_transcriptome (drop any prior rows for this panel)
bsd=pd.read_csv(S('results/base_sd_transcriptome.tsv'),sep='\t')
bsd=bsd[bsd.panel!='gse83352_muscle_exercise_table']
pd.concat([bsd,out],ignore_index=True).to_csv(S('results/base_sd_transcriptome.tsv'),sep='\t',index=False)
print('appended -> results/base_sd_transcriptome.tsv ; panels now:', sorted(pd.read_csv(S('results/base_sd_transcriptome.tsv'),sep='\t').panel.unique()))

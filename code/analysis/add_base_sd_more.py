#!/usr/bin/env python
"""add_base_sd_more.py — compute baseline variance for GSE106737 (liver) + GSE77962 (adipose) from
their baseline samples and APPEND to results/base_sd_transcriptome.tsv -> determinant meta k=5->7.
Reuses the muscle parse (quote-stripped series matrix); panel-specific baseline ID + probe->gene map.
"""
import os, gzip, re, json, io, numpy as np, pandas as pd
SIG=os.path.dirname(os.path.dirname(os.path.abspath(__file__))); PROJ=os.path.dirname(SIG)
def S(*a): return os.path.join(SIG,*a)
RAW=os.path.join(PROJ,'data/transcriptome/raw')

def parse_sm(path):
    meta={'title':None,'gsm':None,'char':[]}; rows=[]; hdr=None; intab=False
    with gzip.open(path,'rt',errors='replace') as fh:
        for l in fh:
            if l.startswith('!Sample_geo_accession'): meta['gsm']=[x.strip().strip('"') for x in l.rstrip().split('\t')[1:]]
            elif l.startswith('!Sample_title'): meta['title']=[x.strip().strip('"') for x in l.rstrip().split('\t')[1:]]
            elif l.startswith('!Sample_characteristics_ch1'): meta['char'].append([x.strip().strip('"') for x in l.rstrip().split('\t')[1:]])
            elif l.startswith('!series_matrix_table_begin'): intab=True; continue
            elif l.startswith('!series_matrix_table_end'): break
            elif intab:
                parts=[x.strip().strip('"') for x in l.rstrip('\n').split('\t')]
                if hdr is None: hdr=parts; continue
                rows.append(parts)
    mat=pd.DataFrame(rows).set_index(0); mat.columns=hdr[1:]
    return mat.apply(pd.to_numeric,errors='coerce'), meta

def collapse(mat, base_cols, p2s, panel):
    e=mat[base_cols]; psd=e.std(axis=1,ddof=1); pmean=e.mean(axis=1)
    d=pd.DataFrame({'probe':e.index.astype(str),'base_sd':psd.values,'base_mean':pmean.values})
    d['mark_id']=d.probe.map(p2s); d=d.dropna(subset=['mark_id','base_sd'])
    d['var']=d.base_sd**2; d=d.sort_values('var',ascending=False).drop_duplicates('mark_id')
    o=d[['mark_id','base_sd','base_mean']].copy(); o['panel']=panel; o['n_base']=len(base_cols); return o

new=[]
# --- GSE106737 liver: baseline = title contains 'baseline' ; probe->gene from json ---
mat,meta=parse_sm(os.path.join(RAW,'GSE106737_series_matrix.txt.gz'))
base=[meta['gsm'][i] for i,t in enumerate(meta['title']) if 'baseline' in str(t).lower()]
base=[c for c in base if c in mat.columns]
p2s=json.load(open(os.path.join(RAW,'gpl16686_probe2sym_reporter.json')))
o=collapse(mat,base,p2s,'gse106737_liver_bariatric_table'); new.append(o)
print(f'GSE106737 liver: {len(base)} baseline samples; base_sd for {len(o)} genes')
# --- GSE77962 adipose: baseline = char 'time point' contains 'baseline' ; probe->gene from GPL11532.annot ---
mat,meta=parse_sm(os.path.join(RAW,'GSE77962_series_matrix.txt.gz'))
tp=[r for r in meta['char'] if any('baseline' in str(v).lower() or 'time' in str(v).lower() for v in r)]
tpr=tp[0] if tp else None
base=[meta['gsm'][i] for i in range(len(meta['gsm'])) if tpr and 'baseline' in str(tpr[i]).lower()]
base=[c for c in base if c in mat.columns]
lines=gzip.open(os.path.join(RAW,'GPL11532.annot.gz'),'rt',errors='replace').read().split('\n')
hi=next((i for i,l in enumerate(lines) if l.split('\t')[0].strip()=='ID'),None)
end=next((i for i,l in enumerate(lines) if 'platform_table_end' in l),len(lines))
ann=pd.read_csv(io.StringIO('\n'.join(lines[hi:end])),sep='\t',dtype=str,low_memory=False,on_bad_lines='skip')
symc=([c for c in ann.columns if c.lower()=='gene symbol'] or [c for c in ann.columns if 'symbol' in c.lower()])[0]
p2s=dict(zip(ann[ann.columns[0]],ann[symc]))
o=collapse(mat,base,p2s,'gse77962_adipose_diet_table'); new.append(o)
print(f'GSE77962 adipose: {len(base)} baseline samples; base_sd for {len(o)} genes')

bsd=pd.read_csv(S('results/base_sd_transcriptome.tsv'),sep='\t')
bsd=bsd[~bsd.panel.isin(['gse106737_liver_bariatric_table','gse77962_adipose_diet_table'])]
pd.concat([bsd]+new,ignore_index=True).to_csv(S('results/base_sd_transcriptome.tsv'),sep='\t',index=False)
print('panels now:', sorted(pd.read_csv(S('results/base_sd_transcriptome.tsv'),sep='\t').panel.unique()))

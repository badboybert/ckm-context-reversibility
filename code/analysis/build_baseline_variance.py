#!/usr/bin/env python
"""build_baseline_variance.py — per-mark BASELINE (pre-intervention) variability for the
PERSISTENCE-eligible panels (SCORE_SPEC_FROZEN.md s3.2/s4). base_sd is used ONLY for DECILE
BINNING in the variance-conditional persistence call -> only its RANK matters, so each panel's
native scale is fine (no cross-panel harmonization needed here).

Persistence-eligible panels (manifest, n>=25 & measured/recompute base_sd):
  methylation : base_sd ALREADY on disk (directplus, drift2 reversal tables) -> read in classify step
  transcriptome: liver GSE83452 (array, baseline samples) + adipose GSE141221 (RNA-seq, CID1 baseline)
  metabolome  : Palau Biocrates (S8 T0 columns)
Proteome persistence = DEFERRED (no open individual baseline; manifest base_sd='proxy_*', not eligible).
Outputs results/base_sd_{transcriptome,metabolome}.tsv  (mark_id, panel, base_sd, n_base).
"""
import os, gzip, json, numpy as np, pandas as pd
SIG=os.path.dirname(os.path.dirname(os.path.abspath(__file__))); PROJ=os.path.dirname(SIG)
def P(*a): return os.path.join(PROJ,*a)
def S(*a): return os.path.join(SIG,*a)
RAW=P('data/transcriptome/raw')
def nm(s): return str(s).strip().rstrip('*').strip().lower()

# ============ TRANSCRIPTOME ============
tx=[]

# --- liver GSE83452 (HuGene array, log-RMA series matrix) ---
meta=pd.read_csv(P('data/transcriptome/GSE83452_sample_meta.tsv'),sep='\t')
base_gsm=set(meta.loc[meta['time'].astype(str).str.lower().str.contains('baseline'),'GSM'])
# parse series-matrix table block
rows=[]; hdr=None; inblk=False
with gzip.open(P('data/transcriptome/GSE83452_series_matrix.txt.gz'),'rt') as fh:
    for ln in fh:
        if ln.startswith('!series_matrix_table_begin'): inblk=True; continue
        if ln.startswith('!series_matrix_table_end'): break
        if inblk:
            parts=ln.rstrip('\n').split('\t')
            if hdr is None: hdr=[p.strip('"') for p in parts]; continue
            rows.append(parts)
mat=pd.DataFrame(rows,columns=hdr).set_index(hdr[0])
mat.index.name='probe'; mat.columns=[c.strip('"') for c in mat.columns]
base_cols=[c for c in mat.columns if c in base_gsm]
mat=mat[base_cols].apply(pd.to_numeric,errors='coerce')
probe_sd=mat.std(axis=1,ddof=1); probe_mean=mat.mean(axis=1)   # base_mean = abundance (detection-floor nuisance)
# join probe-SD/mean -> gene via the reversal table's probe column
lr=pd.read_csv(P('data/transcriptome/liver_reversal_table.tsv'),sep='\t')[['mark_id','probe']]
lr['probe']=lr['probe'].astype(str)
probe_sd.index=probe_sd.index.astype(str); probe_mean.index=probe_mean.index.astype(str)
lr['base_sd']=lr['probe'].map(probe_sd); lr['base_mean']=lr['probe'].map(probe_mean)
liver=lr.dropna(subset=['base_sd'])[['mark_id','base_sd','base_mean']].copy()
liver['panel']='liver_table'; liver['n_base']=len(base_cols); tx.append(liver)
print(f'liver GSE83452: {len(base_cols)} baseline samples; base_sd for {len(liver)}/{len(lr)} reversal genes')

# --- adipose GSE141221 (RNA-seq; replicate builder pairing, CID1=baseline) ---
gmeta=json.load(open(os.path.join(RAW,'gse141221_meta.json')))
def strip_b(k): return k[:-1] if k.endswith('b') else k
ic={(m['indiv'],m['cid']):strip_b(m['key']) for m in gmeta}
counts=pd.read_csv(gzip.open(os.path.join(RAW,'GSE141221_cleaned_count.csv.gz'),'rt'),index_col=0)
colset=set(counts.columns)
pre_cols=[ic.get((iv,'CID1')) for iv in sorted({m['indiv'] for m in gmeta})
          if ic.get((iv,'CID1')) in colset and ic.get((iv,'CID2')) in colset]
mat=counts[pre_cols].astype(float)                      # baseline samples only
lib=mat.sum(axis=0); cpm=mat.divide(lib,axis=1)*1e6
keep=(cpm>=1).sum(axis=1)>=0.2*mat.shape[1]
logcpm=np.log2(cpm[keep]+1.0)
gsd=logcpm.std(axis=1,ddof=1); gmean=logcpm.mean(axis=1)        # base_mean = abundance (detection-floor nuisance)
ens2sym=json.load(open(os.path.join(RAW,'gencode_v25_ens2sym.json')))
ad=pd.DataFrame({'ensg_versioned':logcpm.index,'base_sd':gsd.values,'base_mean':gmean.values})
ad['mark_id']=ad['ensg_versioned'].str.split('.').str[0].map(ens2sym)
ad=ad.dropna(subset=['mark_id'])
ad['var']=ad['base_sd']**2
ad=ad.sort_values('var',ascending=False).drop_duplicates('mark_id',keep='first')   # match builder collapse
ad=ad[['mark_id','base_sd','base_mean']]; ad['panel']='gse141221_adipose_table'; ad['n_base']=len(pre_cols); tx.append(ad)
print(f'adipose GSE141221: {len(pre_cols)} baseline samples; base_sd for {len(ad)} genes')

# --- adipose-CR GSE84046 (series matrix; 'before' baseline samples; probe->gene via reversal table) ---
import re
def _parse_sm(path):
    meta={'char':[]}; rows=[]; header=None; intab=False
    with gzip.open(path,'rt',encoding='utf-8',errors='replace') as f:
        for ln in f:
            if ln.startswith('!Sample_geo_accession'): meta['gsm']=[v.strip().strip('"') for v in ln.rstrip().split('\t')[1:]]
            elif ln.startswith('!Sample_characteristics_ch1'): meta['char'].append([v.strip().strip('"') for v in ln.rstrip().split('\t')[1:]])
            elif ln.startswith('!series_matrix_table_begin'): intab=True; continue
            elif ln.startswith('!series_matrix_table_end'): break
            elif intab:
                if header is None: header=[c.strip().strip('"') for c in ln.rstrip('\n').split('\t')]; continue
                rows.append(ln.rstrip('\n').split('\t'))
    expr=pd.DataFrame(rows).set_index(0); expr.columns=header[1:]; return expr.apply(pd.to_numeric,errors='coerce'),meta
expr,meta=_parse_sm(P('data/transcriptome/GSE84046_series_matrix.txt.gz'))
ba=[r for r in meta['char'] if any('before' in str(v).lower() for v in r)][0]
bef=[meta['gsm'][j] for j in range(len(meta['gsm'])) if 'before' in str(ba[j]).lower()]
e=expr[[c for c in bef if c in expr.columns]]; psd=e.std(axis=1,ddof=1); pmean=e.mean(axis=1)
ar=pd.read_csv(P('data/transcriptome/adipose_reversal_table.tsv'),sep='\t')[['mark_id','probe']]; ar['probe']=ar['probe'].astype(str)
psd.index=psd.index.astype(str); pmean.index=pmean.index.astype(str)
ar['base_sd']=ar['probe'].map(psd); ar['base_mean']=ar['probe'].map(pmean)
acr=ar.dropna(subset=['base_sd'])[['mark_id','base_sd','base_mean']].copy(); acr['panel']='adipose_table'; acr['n_base']=len(bef); tx.append(acr)
print(f'adipose-CR GSE84046: {len(bef)} before samples; base_sd for {len(acr)} genes')

# --- blood GSE273902 (counts; V1 baseline visit; ENSG->symbol) ---
si=pd.read_csv(gzip.open(P('data/transcriptome/raw/GSE273902_Sample_Info.tsv.gz'),'rt'),sep='\t',index_col=0)
cnt=pd.read_csv(gzip.open(P('data/transcriptome/raw/GSE273902_Counts_Table.tsv.gz'),'rt'),sep='\t',index_col=0)
v1=[s for s in si.index[si['Visit']=='V1'] if s in cnt.columns]
m=cnt[v1].astype(float); lib=m.sum(axis=0); cpm=m.divide(lib,axis=1)*1e6
keep=(cpm>=1).sum(axis=1)>=0.2*len(v1); lcpm=np.log2(cpm[keep]+1.0)
gi=pd.read_csv(gzip.open(P('data/transcriptome/raw/GSE273902_Gene_Info.tsv.gz'),'rt'),sep='\t').drop_duplicates('GENEID').set_index('GENEID')['SYMBOL']
bl=pd.DataFrame({'ensg':lcpm.index,'base_sd':lcpm.std(axis=1,ddof=1).values,'base_mean':lcpm.mean(axis=1).values})
bl['mark_id']=bl['ensg'].map(gi); bl=bl.dropna(subset=['mark_id']).sort_values('base_sd',ascending=False).drop_duplicates('mark_id')
bl=bl[['mark_id','base_sd','base_mean']]; bl['panel']='gse273902_blood_table'; bl['n_base']=len(v1); tx.append(bl)
print(f'blood GSE273902: {len(v1)} V1 baseline samples; base_sd for {len(bl)} genes')

tx=pd.concat(tx,ignore_index=True); tx.to_csv(S('results/base_sd_transcriptome.tsv'),sep='\t',index=False)
print(f'wrote results/base_sd_transcriptome.tsv ({len(tx)} rows)')

# ============ METABOLOME (Palau Biocrates, S8 T0 baseline) ============
NONMET={'MEDDM','MEDCOL','MEDINF','MEDHTA','GLU','INS','HOMA','HBA1C','HBA1C.mmol.mol','PESO','CC','CINT',
 'CAD','TAD','TAS','TG','COL','LDL','HDL','VLDL','PCR','LEP','ADIPO','GOT','GPT','GGT','URICO','CREAT','UREA',
 'HIERRO','TRANSF','FERR','bmi','Creatinine'}   # clinical chemistry, not metabolomics analytes
ps=pd.read_csv(P('data/metabolome/reversal_raw/palau_S8_raw.csv'))
t0=[c for c in ps.columns if c.endswith('_T0') and c[:-3] not in NONMET]
recs=[]
for c in t0:
    v=pd.to_numeric(ps[c],errors='coerce').dropna()
    if len(v)>=10: recs.append(dict(mark_key=nm(c[:-3]),base_sd=float(v.std(ddof=1)),n_base=len(v)))
mb=pd.DataFrame(recs).drop_duplicates('mark_key'); mb['panel']='palau'
mb.to_csv(S('results/base_sd_metabolome.tsv'),sep='\t',index=False)
# how many link to the metabolome score?
sc=pd.read_csv(S('results/reversibility_score_metabolome.tsv'),sep='\t')
link=mb['mark_key'].isin(set(sc['mark_key'])).sum()
print(f'metabolome Palau: base_sd for {len(mb)} metabolites; {link} link to the metabolome score (normalized-name match)')
print(f'wrote results/base_sd_metabolome.tsv')

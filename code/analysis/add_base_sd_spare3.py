#!/usr/bin/env python
"""add_base_sd_spare3.py — compute baseline variance for 3 spare on-disk transcriptome panels with
recoverable baseline matrices, and APPEND to results/base_sd_transcriptome.tsv to raise determinant k
above 6 (6 -> up to 9) and tighten the determinant-null MDE.

Panels (panel name = manifest panel that contributes_to_score, so build_determinants joins it):
  gse43471_adipose_table   (adipose, diet+exercise lifestyle, n=39 subjects) baseline = time point 'baseline'
  gse77962_adipose_diet_table (adipose, VLCD/LCD diet, n=48 subjects)        baseline = time point 'at study start'
  gse157585_muscle_table   (muscle, metPRT/plaPRT resistance training, n=47) baseline = biopsy '1' (P/M arms)

base_sd = per-gene SD across the BASELINE/pre samples ONLY, in the SAME value space and gene-collapse rule
(max-variance probe per symbol) used by each panel's existing reversal-table build, so the residualization
cells in build_determinants are consistent. Reuses the add_base_sd_more.py series-matrix parse (quote-strip).
"""
import os, gzip, re, io, numpy as np, pandas as pd
SIG=os.path.dirname(os.path.dirname(os.path.abspath(__file__))); PROJ=os.path.dirname(SIG)
def S(*a): return os.path.join(SIG,*a)
RAW=os.path.join(PROJ,'data/transcriptome/raw')
R=lambda *a: os.path.join(RAW,*a)

def parse_sm(path):
    """series-matrix parse (reused from add_base_sd_more.py): meta + numeric probe x sample matrix."""
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

def char_row(meta,key):
    """return the per-sample characteristic values for the row whose 'key:' prefix matches."""
    for r in meta['char']:
        v0=str(r[0])
        if v0.split(':')[0].strip().lower()==key.lower():
            return [str(x).split(':',1)[1].strip() if ':' in str(x) else str(x).strip() for x in r]
    return None

def annot_p2s(annot_path, sym_colnames=('gene symbol','symbol')):
    """GEO .annot probe->symbol; symbol = explicit named col, else 3rd column (idx 2) fallback."""
    lines=gzip.open(annot_path,'rt',errors='replace').read().split('\n')
    hi=next((i for i,l in enumerate(lines) if l.split('\t')[0].strip()=='ID'),None)
    if hi is None: hi=next(i for i,l in enumerate(lines) if 'platform_table_begin' in l)+1
    end=next((i for i,l in enumerate(lines) if 'platform_table_end' in l),len(lines))
    ann=pd.read_csv(io.StringIO('\n'.join(lines[hi:end])),sep='\t',dtype=str,low_memory=False,on_bad_lines='skip')
    idc=ann.columns[0]
    cands=[c for c in ann.columns if c.lower() in sym_colnames] or [c for c in ann.columns if 'symbol' in c.lower()]
    symc=cands[0] if cands else ann.columns[2]
    return dict(zip(ann[idc].astype(str),ann[symc].astype(str)))

def collapse_sd(e, p2s, panel, n_base):
    """e = probe x baseline-sample matrix (value space matching the panel's reversal build).
    per-gene base_sd = SD across baseline samples; max-variance probe per symbol; drop multi-symbol."""
    psd=e.std(axis=1,ddof=1); pmean=e.mean(axis=1)
    d=pd.DataFrame({'probe':e.index.astype(str),'base_sd':psd.values,'base_mean':pmean.values})
    d['mark_id']=d.probe.map(lambda p:p2s.get(p))
    d=d.dropna(subset=['mark_id','base_sd'])
    d=d[d.mark_id.astype(str).ne('')]
    d=d[~d.mark_id.astype(str).str.contains('///',regex=False)]
    d['var']=d.base_sd**2
    d=d.sort_values('var',ascending=False).drop_duplicates('mark_id')
    o=d[['mark_id','base_sd','base_mean']].copy(); o['panel']=panel; o['n_base']=n_base
    return o

new=[]

# ===== GSE43471 adipose (diet/exercise/diet+exercise lifestyle; control excluded) =====
mat,meta=parse_sm(R('GSE43471_series_matrix.txt.gz'))
tp=char_row(meta,'time point'); grp=char_row(meta,'sample group'); gsm=meta['gsm']
LIFE={'diet','exercise','diet+exercise'}
base_cols=[gsm[i] for i in range(len(gsm))
           if tp[i].lower()=='baseline' and grp[i].lower() in LIFE and gsm[i] in mat.columns]
e=np.log2(mat[base_cols].clip(lower=1))   # reversal build log2-stabilizes the non-log intensities
p2s=annot_p2s(R('GPL6947.annot.gz'))
o=collapse_sd(e,p2s,'gse43471_adipose_table',len(base_cols)); new.append(o)
print(f'GSE43471 adipose: {len(base_cols)} baseline (lifestyle) samples; base_sd for {len(o)} genes')

# ===== GSE77962 adipose (VLCD+LCD diet) baseline = 'at study start' =====
mat,meta=parse_sm(R('GSE77962_series_matrix.txt.gz'))
tp=char_row(meta,'time point'); gsm=meta['gsm']
base_cols=[gsm[i] for i in range(len(gsm)) if tp[i].lower()=='at study start' and gsm[i] in mat.columns]
e=mat[base_cols]   # already RMA log2-normalized intensities (reversal build uses raw values)
p2s=annot_p2s(R('GPL11532.annot.gz'))
o=collapse_sd(e,p2s,'gse77962_adipose_diet_table',len(base_cols)); new.append(o)
print(f'GSE77962 adipose: {len(base_cols)} baseline samples; base_sd for {len(o)} genes')

# ===== GSE157585 muscle (metPRT/plaPRT) baseline = Biopsy-1 of PAIRED P/M subjects =====
# raw counts; replicate reversal build: paired subjects only, CPM>=1 in >=20% of paired samples, log2(CPM+1)
counts=pd.read_csv(gzip.open(R('GSE157585_raw_counts.txt.gz'),'rt'),sep='\t',index_col=0)
counts.columns=[c.strip('"') for c in counts.columns]
subj={}
for c in counts.columns:
    m=re.match(r'([PMY])_(.+)_(\d+)$',c)
    if not m: raise ValueError('unparsed column '+c)
    subj.setdefault((m.group(1),m.group(2)),{})[m.group(3)]=c
pairs=[(a,s,d['1'],d['3']) for (a,s),d in subj.items() if a!='Y' and '1' in d and '3' in d]
pairs.sort()
pre_cols=[p[2] for p in pairs]; post_cols=[p[3] for p in pairs]
mat=counts[pre_cols+post_cols].astype(float)
lib=mat.sum(axis=0); cpm=mat.divide(lib,axis=1)*1e6
keep=(cpm>=1).sum(axis=1)>=(0.2*mat.shape[1])
logcpm=np.log2(cpm[keep]+1.0)
e=logcpm[pre_cols]   # BASELINE (Biopsy-1) only, same filtered gene set as the reversal build
e.index=[g.split('_',1)[1] if '_' in g else g for g in e.index]  # ENSG..._SYMBOL -> SYMBOL
e=e[[s!='' and s is not None for s in e.index]]
# symbol already embedded -> identity map; collapse_sd handles max-variance + multi-symbol drop
p2s={g:g for g in e.index}
o=collapse_sd(e,p2s,'gse157585_muscle_table',len(pre_cols)); new.append(o)
print(f'GSE157585 muscle: {len(pre_cols)} baseline (Biopsy-1, paired P/M) samples; base_sd for {len(o)} genes')

# ===== APPEND (drop any prior rows for these 3 panels) =====
TARGS={'gse43471_adipose_table','gse77962_adipose_diet_table','gse157585_muscle_table'}
bsd=pd.read_csv(S('results/base_sd_transcriptome.tsv'),sep='\t')
bsd=bsd[~bsd.panel.isin(TARGS)]
pd.concat([bsd]+new,ignore_index=True).to_csv(S('results/base_sd_transcriptome.tsv'),sep='\t',index=False)
final=pd.read_csv(S('results/base_sd_transcriptome.tsv'),sep='\t')
print('\nbase_sd panels now:')
print(final.panel.value_counts().to_string())

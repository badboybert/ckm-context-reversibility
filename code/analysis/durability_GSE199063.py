#!/usr/bin/env python
"""durability_GSE199063.py — DURABILITY/persistence-responder analysis (responder dataset #1).
GSE199063 = RYGB bariatric, SAT adipose, with Baseline / 2-year / 5-year timepoints. Of the reversals
seen at 2yr, which PERSIST to 5yr (durable = lasting benefit) vs REBOUND toward baseline (transient =
where molecular benefit is lost by 5yr = residual-risk candidates)? Uses subjects with ALL THREE
timepoints so 2yr and 5yr reversals are within the same people. Reuses the panel's probe->gene map.
"""
import os, gzip, numpy as np, pandas as pd
from scipy import stats
SIG=os.path.dirname(os.path.dirname(os.path.abspath(__file__))); PROJ=os.path.dirname(SIG)
def S(*a): return os.path.join(SIG,*a)
RAW=os.path.join(PROJ,'data/transcriptome/raw'); MX=os.path.join(RAW,'GSE199063_series_matrix.txt.gz'); PLAT=os.path.join(RAW,'GPL23126_platform_full.txt.gz')
def grab(l): return [p.strip().strip('"') for p in l.rstrip('\n').split('\t')]
gsm=subj=tp=None
with gzip.open(MX,'rt',errors='replace') as fh:
    for l in fh:
        if l.startswith('!Sample_geo_accession'): gsm=grab(l)[1:]
        elif l.startswith('!Sample_characteristics_ch1') and 'subject_id:' in l: subj=[v.replace('subject_id:','').strip() for v in grab(l)[1:]]
        elif l.startswith('!Sample_characteristics_ch1') and 'time_point:' in l: tp=[v.replace('time_point:','').strip() for v in grab(l)[1:]]
        elif l.startswith('!series_matrix_table_begin'): break
meta=pd.DataFrame({'gsm':gsm,'subject':subj,'tp':tp})
T={t:meta[meta.tp==t].set_index('subject')['gsm'] for t in ['Baseline','2_year','5_year']}
S3=sorted(set(T['Baseline'].index)&set(T['2_year'].index)&set(T['5_year'].index))
print(f'subjects with ALL THREE timepoints (Baseline/2yr/5yr): {len(S3)}')
mat=pd.read_csv(MX,sep='\t',comment='!',index_col=0,encoding_errors='replace')
mat.index=[str(i).strip('"') for i in mat.index]; mat.columns=[str(c).strip('"') for c in mat.columns]
mat=mat.loc[mat.index.to_series().str.match(r'^(TC|HTA)'),[c for c in mat.columns if c.startswith('GSM')]].apply(pd.to_numeric,errors='coerce')
base=mat[[T['Baseline'][s] for s in S3]].values; y2=mat[[T['2_year'][s] for s in S3]].values; y5=mat[[T['5_year'][s] for s in S3]].values
rev2=(y2-base).mean(1); rev5=(y5-base).mean(1); p2=stats.ttest_rel(y2,base,axis=1)[1]
var=mat[[T['Baseline'][s] for s in S3]].values.var(1)
def sym(ga):
    if not isinstance(ga,str) or ga in ('---',''): return None
    f=[x.strip() for x in ga.split('///')[0].split('//')]; return f[1] if len(f)>=2 and f[1] not in ('','---') else None
p2s=pd.read_csv(PLAT,sep='\t',skiprows=1,comment='!',dtype=str,encoding_errors='replace',usecols=['ID','gene_assignment'])
p2s['g']=p2s.gene_assignment.map(sym); m=dict(zip(p2s.ID,p2s.g))
r=pd.DataFrame({'probe':mat.index,'rev2':rev2,'rev5':rev5,'p2':np.where(np.isfinite(p2),p2,1),'var':var})
r['gene']=r.probe.map(m); r=r.dropna(subset=['gene']).sort_values('var',ascending=False).drop_duplicates('gene')

print(f'\nOverall durability: Spearman(2yr reversal, 5yr reversal) = {stats.spearmanr(r.rev2,r.rev5).correlation:+.2f}')
sig=r[(r.p2<0.05)&(r.rev2.abs()>0.2)].copy(); sig['ratio']=sig.rev5/sig.rev2   # 1=fully durable, 0=rebounded, <0=overshoot
print(f'among {len(sig)} genes reversed at 2yr: median 5yr/2yr retention ratio = {sig.ratio.median():.2f} (1=durable, 0=fully rebounded)')
print(f'  DURABLE (ratio>0.7, reversal persists to 5yr): {int((sig.ratio>0.7).sum())} ({(sig.ratio>0.7).mean()*100:.0f}%)')
print(f'  REBOUNDED (ratio<0.3, benefit lost by 5yr):    {int((sig.ratio<0.3).sum())} ({(sig.ratio<0.3).mean()*100:.0f}%)')
import re; NC=re.compile(r'-AS\d?|^LINC|^LOC|^AC\d|^AL\d|orf')
clean=sig[~sig.gene.astype(str).str.contains(NC,na=False)]
dur=clean[clean.ratio>0.7].reindex(clean[clean.ratio>0.7].rev2.abs().sort_values(ascending=False).index)
reb=clean[(clean.ratio<0.3)].reindex(clean[(clean.ratio<0.3)].rev2.abs().sort_values(ascending=False).index)
print('\n  top DURABLE reversers (persist to 5yr): '+', '.join(dur.gene.head(10)))
print('  top REBOUND reversers (lost by 5yr = residual-risk candidates): '+', '.join(reb.gene.head(10)))
r[['gene','rev2','rev5','p2']].assign(retention=r.rev5/r.rev2).to_csv(S('results/durability_GSE199063.tsv'),sep='\t',index=False)
print('\n  wrote results/durability_GSE199063.tsv  (CAVEAT: n=%d, microarray, RYGB single cohort -> hypothesis-generating)'%len(S3))

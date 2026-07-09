#!/usr/bin/env python
"""restoration_liver_blood.py
COMPLETE the restoration arm beyond adipose: compute per-context %toward-lean + net restoration for
LIVER and BLOOD transcriptome contexts, orienting each reversal panel by its tissue-matched BMI axis.

toward_lean = -sign(bmi_beta) * rev_beta   (rev_beta>0 = mark UP after intervention; convention).
A context "restores toward lean" when the intervention moves marks OPPOSITE to the BMI(obese) direction.

Axes:
  LIVER -> data/features/bmi_axis_transcriptome_liver.tsv  (GSE48452 in-vivo, VERIFIED 11/12 canonical)
  BLOOD -> data/features/bmi_axis_transcriptome_blood.tsv  (Homuth2015 KORA+SHIP, VERIFIED erythroid UP)
Reversal contexts (keyed by mark_id = gene symbol, rev_beta):
  liver:  liver_reversal_table (GSE83452), gse48452_liver_reversal_table, gse106737_liver_bariatric_reversal_table
  blood:  gse273902_blood_reversal_table
Compare to adipose 52-68% / plasma proteome 52-63%.
"""
import os, re, numpy as np, pandas as pd
PROJ=r"C:\Users\Bert\Downloads\CKM papers\weight-loss.omics"
SIG=os.path.join(PROJ,'signature-pivot')
NC=re.compile(r'-AS\d?|^LINC|orf\d|^RNU|^SNOR|^MIR\d|^LOC|^AC\d|^AL\d|^RP[0-9LS]|^IG[HKL][VJC]|rRNA|^7SK')

def load_axis(path):
    d=pd.read_csv(path,sep='\t')
    d['gene']=d['gene'].astype(str).str.upper()
    d['bmi_beta']=pd.to_numeric(d['bmi_beta'],errors='coerce')
    d=d.dropna(subset=['bmi_beta']).drop_duplicates('gene')
    return d.set_index('gene')['bmi_beta']

ax_liver=load_axis(os.path.join(SIG,'data','features','bmi_axis_transcriptome_liver.tsv'))
ax_blood=load_axis(os.path.join(SIG,'data','features','bmi_axis_transcriptome_blood.tsv'))

def ctx(path, axis, label):
    d=pd.read_csv(path,sep='\t')
    d=d.rename(columns={'mark_id':'gene'})
    d['gene']=d['gene'].astype(str).str.upper()
    d=d[~d['gene'].str.contains(NC,na=False)]
    d['rev_beta']=pd.to_numeric(d['rev_beta'],errors='coerce')
    d=d.dropna(subset=['rev_beta']).drop_duplicates('gene')
    j=d.join(axis,on='gene',how='inner').dropna(subset=['bmi_beta'])
    n=len(j)
    if n<20: return dict(label=label,n=n,note='insufficient overlap')
    j['toward_lean']=-np.sign(j['bmi_beta'])*j['rev_beta']
    pct=float((j['toward_lean']>0).mean())*100
    net=float(j['toward_lean'].mean())
    mag=float(j['rev_beta'].abs().median())
    # restrict to axis-CONFIDENT genes (bmi_p<0.05) for a cleaner orientation read
    return dict(label=label,n=n,pct_toward_lean=pct,net_toward_lean=net,resp_mag=mag)

T=os.path.join(PROJ,'data','transcriptome')
LIVER=[('liver GSE83452 (diet+BS)','liver_reversal_table'),
       ('liver GSE48452 (BS)','gse48452_liver_reversal_table'),
       ('liver GSE106737 (bariatric)','gse106737_liver_bariatric_reversal_table')]
BLOOD=[('blood GSE273902','gse273902_blood_reversal_table')]

out=[]
pr=lambda s:(out.append(s),print(s))
pr('================ RESTORATION (toward-lean) — LIVER + BLOOD transcriptome ================')
pr('orientation: toward_lean = -sign(bmi_beta)*rev_beta ; %>50 = majority of marks move toward lean')
pr(f'\n[LIVER]  axis = GSE48452 in-vivo BMI->expr (VERIFIED 11/12 canonical UP); n_axis={len(ax_liver)}')
pr(f'  {"context":30s} {"n_overlap":>9s} {"%toward_lean":>13s} {"net":>9s} {"resp_mag":>9s}')
for lab,f in LIVER:
    r=ctx(os.path.join(T,f+'.tsv'),ax_liver,lab)
    if 'pct_toward_lean' in r:
        pr(f'  {lab:30s} {r["n"]:9d} {r["pct_toward_lean"]:12.1f}% {r["net_toward_lean"]:+9.4f} {r["resp_mag"]:9.4f}')
    else:
        pr(f'  {lab:30s} {r["n"]:9d}   {r.get("note","")}')
pr(f'\n[BLOOD]  axis = Homuth2015 KORA+SHIP whole-blood BMI->expr (VERIFIED erythroid UP); n_axis={len(ax_blood)}')
pr(f'  {"context":30s} {"n_overlap":>9s} {"%toward_lean":>13s} {"net":>9s} {"resp_mag":>9s}')
for lab,f in BLOOD:
    r=ctx(os.path.join(T,f+'.tsv'),ax_blood,lab)
    if 'pct_toward_lean' in r:
        pr(f'  {lab:30s} {r["n"]:9d} {r["pct_toward_lean"]:12.1f}% {r["net_toward_lean"]:+9.4f} {r["resp_mag"]:9.4f}')
    else:
        pr(f'  {lab:30s} {r["n"]:9d}   {r.get("note","")}')
pr('\nBENCHMARK: adipose 52-68% toward-lean ; plasma proteome 52-63%.')
open(os.path.join(SIG,'results','restoration_liver_blood.txt'),'w',encoding='utf-8').write('\n'.join(out))
print('\nwrote results/restoration_liver_blood.txt')

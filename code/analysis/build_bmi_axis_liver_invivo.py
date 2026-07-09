#!/usr/bin/env python
"""build_bmi_axis_liver_invivo.py
Build an IN-VIVO liver BMI->expression axis from GSE48452 (Ahrens 2013, PMID 23803177):
human liver biopsies across control->obese->steatosis->NASH, n=73, with per-sample BMI.
Regress each gene's log2 expression on BMI (OLS), bmi_beta>0 = UP with BMI (matches convention).
This REPLACES the GTEx post-mortem liver axis (Signer2024) which was too weak/noise-dominated to
orient canonical NAFLD-up genes correctly. GSE48452 is in-vivo, obesity-spanning, and the SAME
platform (GPL11532) as the gse48452_liver reversal context.
Output: data/features/bmi_axis_transcriptome_liver.tsv (gene, bmi_beta, bmi_p) -- OVERWRITES GTEx version.
"""
import os, gzip, re, numpy as np
from scipy import stats
PROJ=r"C:\Users\Bert\Downloads\CKM papers\weight-loss.omics"
SIG=os.path.join(PROJ,'signature-pivot')
MAT=os.path.join(PROJ,'data','transcriptome','GSE48452_series_matrix.txt.gz')
P2S=os.path.join(SIG,'data','features','raw','gpl11532_probe2sym.tsv')
OUT=os.path.join(SIG,'data','features','bmi_axis_transcriptome_liver.tsv')

# ---- parse BMI per sample from series matrix header ----
bmi=None; samples=None; data_start=False; rows=[]
with gzip.open(MAT,'rt',encoding='utf-8',errors='replace') as f:
    for line in f:
        if line.startswith('!Sample_geo_accession'):
            samples=[x.strip().strip('"') for x in line.rstrip('\n').split('\t')[1:]]
        if line.startswith('!Sample_characteristics_ch1') and 'bmi:' in line.lower():
            vals=[x.strip().strip('"') for x in line.rstrip('\n').split('\t')[1:]]
            bmi=[]
            for v in vals:
                m=re.search(r'bmi:\s*([\d.]+)',v,re.I)
                bmi.append(float(m.group(1)) if m else np.nan)
        if line.startswith('!series_matrix_table_begin'):
            data_start=True; continue
        if line.startswith('!series_matrix_table_end'): break
        if data_start:
            parts=line.rstrip('\n').split('\t')
            if parts[0]=='"ID_REF"':
                hdr=[x.strip().strip('"') for x in parts[1:]]; continue
            rows.append(parts)
bmi=np.array(bmi,float)
print(f'samples={len(samples)} bmi non-nan={np.sum(~np.isnan(bmi))} bmi range {np.nanmin(bmi):.1f}-{np.nanmax(bmi):.1f}')

# align expression columns to bmi order (hdr == sample GSM order in table)
# samples from geo_accession should match hdr
idx=[samples.index(h) for h in hdr]
bmi_ord=bmi[idx]
keep=~np.isnan(bmi_ord)
bmi_use=bmi_ord[keep]
print(f'using {keep.sum()} samples with BMI')

# probe -> symbol
p2s={}
with open(P2S,encoding='utf-8') as f:
    for line in f:
        a=line.rstrip('\n').split('\t')
        if len(a)>=2 and a[1].strip():
            p2s[a[0].strip()]=a[1].split('///')[0].strip().upper()

# per-probe OLS expr~bmi, collapse to gene by min-p
best={}
n_used=0
for parts in rows:
    pid=parts[0].strip().strip('"')
    sym=p2s.get(pid)
    if not sym: continue
    try:
        expr=np.array([float(x) for x in parts[1:]],float)
    except: continue
    e=expr[keep]
    if np.sum(np.isfinite(e))<10 or np.std(e)==0: continue
    sl,ic,r,p,se=stats.linregress(bmi_use,e)
    n_used+=1
    if sym not in best or p<best[sym][1]:
        best[sym]=(sl,p)
print(f'probes regressed={n_used} -> genes={len(best)}')

# ---- self-check: canonical NAFLD/obesity liver UP genes should be bmi_beta>0 ----
UP=['LEP','FABP4','CD36','SPP1','LGALS3','CCL2','CD68','COL1A1','ANXA2','GPNMB','TREM2','CIDEC','PLIN1']
print('\n--- SELF-CHECK liver in-vivo (canonical obesity/NAFLD UP -> expect bmi_beta>0) ---')
ok=tot=0; flags=[]
for g in UP:
    if g in best:
        b,p=best[g]; tot+=1; m='OK' if b>0 else 'WRONG'
        if b>0: ok+=1
        else: flags.append(g)
        print(f'   {g:8s} bmi_beta={b:+.4f} p={p:.2e} {m}')
    else: print(f'   {g:8s} (not measured)')
bs=np.array([v[0] for v in best.values()])
print(f'   GLOBAL n={len(bs)} frac_pos={np.mean(bs>0):.2f} median|beta|={np.median(np.abs(bs)):.4f}')
print(f'   canonical-correct {ok}/{tot} flagged={flags}')

with open(OUT,'w',encoding='utf-8') as f:
    f.write('gene\tbmi_beta\tbmi_p\n')
    for g,(b,p) in sorted(best.items()):
        f.write(f'{g}\t{b}\t{p}\n')
print('wrote',OUT,len(best))
frac=ok/tot if tot else 0
print(f'\nLIVER IN-VIVO AXIS VERDICT: {ok}/{tot} ({frac:.0%}) -> {"OK" if frac>=0.6 else "STILL WEAK"}')

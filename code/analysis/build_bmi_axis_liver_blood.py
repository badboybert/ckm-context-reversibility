#!/usr/bin/env python
"""build_bmi_axis_liver_blood.py
Build BMI->EXPRESSION axes for LIVER and BLOOD to complete the restoration arm beyond adipose.

LIVER  source: Signer et al. 2024 medRxiv (PMID 39649609, DOI 10.1101/2024.11.26.24317923),
               GTEx BMI main-effect DE per tissue. Sheet 'S1 BMI DE Analyses', tissue=Liver.
               column Beta_BMImain = BMI->expression beta (positive = up with BMI). n=208 (Liver GTEx).
BLOOD  source: Homuth et al. 2015 BMC Med Genomics (PMID 26470795), KORA F4 + SHIP-TREND whole-blood
               mRNA, n=1977. Sheet 'Supp1', column Effect (col idx4) = meta BMI->expression beta;
               p_z(BH) (col idx3) = meta BH p. positive Effect = up with BMI.

Convention (matches METSIM adipose axis GSE70353_metsim_adipose_BMI.tsv): bmi_beta>0 => UP with BMI.
Output: data/features/bmi_axis_transcriptome_liver.tsv  (gene, bmi_beta, bmi_p)
        data/features/bmi_axis_transcriptome_blood.tsv  (gene, bmi_beta, bmi_p)
"""
import os, numpy as np, openpyxl
SIG=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW=os.path.join(SIG,'data','features','raw')
OUT=os.path.join(SIG,'data','features')

# ---- canonical BMI-UP / BMI-DOWN expected genes for self-check ----
# Blood whole-blood: inflammation/leptin-axis up with BMI is the textbook signature.
# (LEP is adipose-restricted; in blood use inflammation: SAA, CRP, S100A, leukocyte/ROS genes.)
EXPECT_UP_BLOOD   = ['SAA1','SAA2','S100A8','S100A9','TNFRSF1A','ALOX5AP']  # inflammation up w/ BMI in blood
EXPECT_UP_LIVER   = ['LEP','CRP','SAA1','SAA2','FABP4','CD68','COL1A1']     # leptin/inflammation/fibrosis up w/ BMI in liver
EXPECT_DOWN_LIVER = []  # do not hard-assume liver-down set

def collapse(rows):
    """rows: list of (gene,beta,p). collapse multi-probe by min-p per gene."""
    best={}
    for g,b,p in rows:
        if g is None: continue
        g=str(g).strip().upper()
        if g in ('','NONE','NA','---'): continue
        try: b=float(b); p=float(p)
        except: continue
        if not np.isfinite(b) or not np.isfinite(p): continue
        if g not in best or p<best[g][1]:
            best[g]=(b,p)
    return best

# ================= LIVER =================
wb=openpyxl.load_workbook(os.path.join(RAW,'media-4.xlsx'),read_only=True)
ws=wb['S1 BMI DE Analyses']
liver=[]
for i,row in enumerate(ws.iter_rows(values_only=True)):
    if i<2: continue
    if row[0]!='Liver': continue
    # gene_name=2, Beta_BMImain=7, pvalue_BMImain=9
    liver.append((row[2],row[7],row[9]))
wb.close()
lv=collapse(liver)
print(f'[LIVER] Signer2024 GTEx-Liver: {len(liver)} rows -> {len(lv)} unique genes')

# ================= BLOOD =================
wb=openpyxl.load_workbook(os.path.join(RAW,'homuth_S1.xlsx'),read_only=True)
ws=wb['Supp1']
blood=[]
for i,row in enumerate(ws.iter_rows(values_only=True)):
    if i<2: continue
    # gene=2, Effect(meta)=4, p_z(BH)=3
    blood.append((row[2],row[4],row[3]))
wb.close()
bl=collapse(blood)
print(f'[BLOOD] Homuth2015 KORA+SHIP whole-blood: {len(blood)} rows -> {len(bl)} unique genes')

def selfcheck(d,name,up,down):
    print(f'\n--- SELF-CHECK {name} (canonical BMI-UP should be bmi_beta>0) ---')
    ok=0; tot=0; flags=[]
    for g in up:
        if g in d:
            b,p=d[g]; tot+=1; sgn='UP' if b>0 else 'DOWN'
            mark='OK' if b>0 else 'WRONG'
            if b>0: ok+=1
            else: flags.append(g)
            print(f'   {g:10s} bmi_beta={b:+.4f} p={p:.2e}  [{sgn}] {mark}')
        else:
            print(f'   {g:10s} (not measured)')
    for g in down:
        if g in d:
            b,p=d[g]; tot+=1
            mark='OK' if b<0 else 'WRONG'
            if b<0: ok+=1
            else: flags.append(g)
            print(f'   {g:10s} bmi_beta={b:+.4f} p={p:.2e}  [expect DOWN] {mark}')
    # global sanity: median |beta| nonzero, fraction positive ~balanced or inflammation-skewed
    bs=np.array([v[0] for v in d.values()])
    print(f'   GLOBAL: n={len(bs)} frac_pos={np.mean(bs>0):.2f} median|beta|={np.median(np.abs(bs)):.4f}')
    print(f'   canonical-correct {ok}/{tot}; flagged={flags}')
    return ok,tot,flags

ok_l,tot_l,flag_l=selfcheck(lv,'LIVER',EXPECT_UP_LIVER,EXPECT_DOWN_LIVER)
ok_b,tot_b,flag_b=selfcheck(bl,'BLOOD',EXPECT_UP_BLOOD,[])

def write(d,path):
    with open(path,'w',encoding='utf-8') as f:
        f.write('gene\tbmi_beta\tbmi_p\n')
        for g,(b,p) in sorted(d.items()):
            f.write(f'{g}\t{b}\t{p}\n')
    print('wrote',path,len(d))

write(lv,os.path.join(OUT,'bmi_axis_transcriptome_liver.tsv'))
write(bl,os.path.join(OUT,'bmi_axis_transcriptome_blood.tsv'))

# decision: corrupt if <60% of MEASURED canonical-up genes are positive
def verdict(ok,tot,name):
    if tot==0: print(f'{name}: NO canonical genes measured -> cannot verify orientation'); return 'UNVERIFIED'
    frac=ok/tot
    v='OK' if frac>=0.6 else 'CORRUPT/WRONG-SIGNED'
    print(f'{name} AXIS VERDICT: {ok}/{tot} canonical correct ({frac:.0%}) -> {v}')
    return v
print()
verdict(ok_l,tot_l,'LIVER')
verdict(ok_b,tot_b,'BLOOD')

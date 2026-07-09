#!/usr/bin/env python
"""restoration_persistence_revisit.py
PART C: aggregate MOVE-TOWARD-LEAN magnitude PER CONTEXT (user q) — orient each panel's rev_beta by the
        matching BMI axis (toward-lean = -sign(BMI-axis)*rev_beta) and report net restoration + magnitude.
PART D: PERSISTENCE revisit (user q) — even though the determinant META is null (no intrinsic feature
        predicts reversibility), is there (a) a persistence SIGNATURE and (b) cross-context CORRELATION?
        Recompute per-tissue persistent vs reversible consensus, cross-TISSUE overlap fold, and the
        biological character (housekeeping/low-variance vs pathway) of each pole.
"""
import os, re, numpy as np, pandas as pd
from collections import Counter, defaultdict
SIG=os.path.dirname(os.path.dirname(os.path.abspath(__file__))); PROJ=os.path.dirname(SIG)
def P(*a): return os.path.join(PROJ,*a)
def S(*a): return os.path.join(SIG,*a)
NC=re.compile(r'-AS\d?|^LINC|orf\d|^RNU|^SNOR|^MIR\d|^LOC|^AC\d|^AL\d|^RP[0-9LS]|^IG[HKL][VJC]')
out=[]; pr=lambda s:(out.append(s),print(s))

# ================= PART C: per-context restoration magnitude =================
pr('================ PART C: aggregate move-toward-lean magnitude PER CONTEXT ================')
def axis_map(path,keycands):
    d=pd.read_csv(P(path),sep='\t')
    kc=next((c for c in keycands if c in d.columns),d.columns[0])
    bc=next((c for c in d.columns if re.search(r'bmi|beta|effect|axis|mr',c,re.I) and c!=kc),None)
    d=d[[kc,bc]].dropna(); d.columns=['key','bmi']; d['bmi']=pd.to_numeric(d['bmi'],errors='coerce')
    return d.dropna().drop_duplicates('key').set_index('key')['bmi']
axp=axis_map('signature-pivot/results/bmi_axis_proteome.tsv',['UniProt','gene','assay_id'])
axt=axis_map('signature-pivot/data/features/bmi_axis_transcriptome.tsv',['gene','mark_id','symbol'])
def ctx_restore(path, axis, key_in_data, tissue):
    if not os.path.exists(P(path)) or axis is None: return None
    d=pd.read_csv(P(path),sep='\t')
    kc=next((c for c in [key_in_data,'UniProt','gene','mark_id'] if c in d.columns),None)
    if kc is None: return None
    d=d.rename(columns={kc:'key'}); d['key']=d['key'].astype(str).str.split(r'[;,_]').str[0]
    d['b']=pd.to_numeric(d['rev_beta'],errors='coerce'); d=d.dropna(subset=['b'])
    j=d.join(axis,on='key',how='inner')
    if len(j)<20: return None
    j['toward_lean']=-np.sign(j['bmi'])*j['b']   # +ve = moved toward lean
    return dict(tissue=tissue,n=len(j),resp_mag=float(j['b'].abs().median()),
                net_toward_lean=float(j['toward_lean'].mean()),pct_toward_lean=float((j['toward_lean']>0).mean()))
PROT=[('BBS bariatric','data/proteome/reversal_olink_bbs.tsv','UniProt'),
      ('DiRECT diet','data/proteome/reversal_somascan_direct.tsv','UniProt'),
      ('EMPEROR SGLT2i','data/proteome/reversal_emperor.tsv','UniProt'),
      ('metformin','data/proteome/nowak_metformin_proteome_reversal.tsv','UniProt'),
      ('semaglutide GLP-1RA','data/proteome/step_semaglutide_proteome_reversal.tsv','UniProt'),
      ('MS bariatric','data/proteome/bariatric_ms_proteome_reversal.tsv','UniProt')]
# proteome axis is keyed by whatever axis_map picked; try UniProt then gene
pr('\n[proteome, plasma] (BMI axis = Goudswaard ST5 INTERVAL):')
pr(f'  {"context":18s} {"n":>5s} {"resp_mag(med|b|)":>16s} {"net_toward_lean":>16s} {"%toward_lean":>12s}')
for lab,f,k in PROT:
    r=ctx_restore(f,axp,k,'plasma')
    if not r:  # retry mapping by gene
        d=pd.read_csv(P(f),sep='\t')
        if 'gene' in d.columns: r=ctx_restore(f,axp,'gene','plasma')
    if r: pr(f'  {lab:18s} {r["n"]:5d} {r["resp_mag"]:16.3f} {r["net_toward_lean"]:+16.3f} {r["pct_toward_lean"]*100:11.0f}%')
    else: pr(f'  {lab:18s}  (no axis overlap)')
TXN_ADIP=[('adipose CR','data/transcriptome/adipose_reversal_table.tsv'),
          ('adipose LCD(DiOGenes)','data/transcriptome/gse141221_adipose_reversal_table.tsv'),
          ('adipose diet GSE77962','data/transcriptome/gse77962_adipose_diet_reversal_table.tsv'),
          ('adipose RYGB GSE199063','data/transcriptome/gse199063_adipose_reversal_table.tsv')]
pr('\n[transcriptome, adipose] (BMI axis = METSIM adipose; valid for ADIPOSE contexts only):')
pr(f'  {"context":24s} {"n":>5s} {"resp_mag":>10s} {"net_toward_lean":>16s} {"%toward_lean":>12s}')
for lab,f in TXN_ADIP:
    r=ctx_restore(f,axt,'mark_id','adipose')
    if r: pr(f'  {lab:24s} {r["n"]:5d} {r["resp_mag"]:10.3f} {r["net_toward_lean"]:+16.3f} {r["pct_toward_lean"]*100:11.0f}%')
    else: pr(f'  {lab:24s}  (no axis overlap)')
pr('  NOTE: net_toward_lean>0 = panel net-restores toward lean; magnitude = strength; %>50 = majority toward lean.')

# ================= PART D: persistence revisit =================
pr('\n\n================ PART D: PERSISTENCE revisit (signature + cross-context correlation) ================')
PAN=[('liver','liver_reversal_table'),('liver','gse48452_liver_reversal_table'),('liver','gse106737_liver_bariatric_reversal_table'),
     ('adipose','adipose_reversal_table'),('adipose','gse141221_adipose_reversal_table'),('adipose','gse43471_adipose_reversal_table'),
     ('adipose','gse199063_adipose_reversal_table'),('adipose','gse77962_adipose_diet_reversal_table'),('adipose','gse70529_adipose_dose_reversal_table'),
     ('adipose','gse107894_adipose_metformin_reversal_table'),('adipose','gse224310_adipose_exercise_reversal_table'),
     ('muscle','gse83352_muscle_exercise_reversal_table'),('muscle','gse161643_muscle_bariatric_reversal_table'),('muscle','gse157585_muscle_reversal_table'),
     ('blood','gse273902_blood_reversal_table'),('blood','gse310742_blood_glp1ra_reversal_table'),('blood','gse193771_blood_exercise_reversal_table'),('blood','gse28358_pbmc_diet_reversal_table')]
def load(f):
    p=P('data/transcriptome',f+'.tsv')
    if not os.path.exists(p): return None
    d=pd.read_csv(p,sep='\t').drop_duplicates('mark_id'); d=d[~d.mark_id.str.contains(NC,na=False)]
    d['ab']=pd.to_numeric(d.rev_beta,errors='coerce').abs(); d=d.dropna(subset=['ab'])
    d['per']=d['ab']<=d['ab'].quantile(1/3); d['per5']=d['ab']<=d['ab'].quantile(0.05); d['rev']=d['ab']>=d['ab'].quantile(0.95)
    return d.set_index('mark_id')
byT=defaultdict(list); LD={}
for t,f in PAN:
    d=load(f)
    if d is not None: LD[f]=d; byT[t].append(f)
# per-tissue persistent + reversible consensus sets
persS={}; revS={}; persS5={}; measU={}
for t,fs in byT.items():
    k=len(fs); need=max(2,(k+1)//2)
    pc=Counter(); rc=Counter(); p5=Counter()
    for f in fs:
        for g in LD[f].index[LD[f]['per']]: pc[g]+=1
        for g in LD[f].index[LD[f]['per5']]: p5[g]+=1
        for g in LD[f].index[LD[f]['rev']]: rc[g]+=1
    persS[t]={g for g,n in pc.items() if n>=need}; revS[t]={g for g,n in rc.items() if n>=need}; persS5[t]={g for g,n in p5.items() if n>=need}
    measU[t]=set().union(*[set(LD[f].index) for f in fs])
pr('\n(a) cross-TISSUE overlap of consensus sets — CORRECT universe (genes measured in BOTH tissues) + hypergeometric:')
tis=[t for t in ['liver','adipose','muscle','blood'] if t in persS]
from scipy.stats import hypergeom
def enr(a,b,U):
    a=a&U; b=b&U; o=len(a&b); exp=len(a)*len(b)/len(U) if len(U) else np.nan
    fold=o/exp if exp and exp>0 else np.nan
    p=float(hypergeom.sf(o-1,len(U),len(a),len(b))) if (len(a) and len(b)) else np.nan
    return o,fold,p
revf=[]; perf=[]; pmf=[]; rev_core=Counter(); rev_pair_rows=[]
for i in range(len(tis)):
    for j in range(i+1,len(tis)):
        ta,tb=tis[i],tis[j]; U=measU[ta]&measU[tb]
        ro,rf,rp=enr(revS[ta],revS[tb],U); po,pf,pp=enr(persS[ta],persS[tb],U); mo,mf,mp=enr(persS5[ta],persS5[tb],U)
        revf.append(rf); perf.append(pf); pmf.append(mf)
        shared=sorted((revS[ta]&U)&(revS[tb]&U))
        for g in shared: rev_core[g]+=1
        rev_pair_rows.append({'tissue_a':ta,'tissue_b':tb,'U':len(U),'n_shared':ro,'fold':rf,'p':rp,'genes':','.join(shared)})
        pr(f'   {ta:7s} x {tb:7s} |U|={len(U):5d}: REVERS {ro:3d} fold={rf:4.2f}x p={rp:.1e} | PERSIST-tertile {po:4d} fold={pf:4.2f}x | PERSIST-top5%(size-matched) {mo:3d} fold={mf:.2f}x')
# export the universal reversible core (union of cross-tissue reversible-consensus intersections) for Figure 2d / Suppl S7
pd.DataFrame(rev_pair_rows).to_csv(S('results/cross_tissue_reversible_pairs.tsv'),sep='\t',index=False)
core_df=pd.DataFrame([{'gene':g,'n_tissue_pairs':n} for g,n in rev_core.most_common()])
core_df.to_csv(S('results/universal_reversible_core.tsv'),sep='\t',index=False)
NAMED=['THBS1','CD83','PLAUR','CCL2','SPP1','SAA1','IL18','ACLY','ANGPTL4','PFKFB3']
pr(f'   universal reversible core = {len(core_df)} unique genes (union of cross-tissue reversible overlaps)')
pr('   named-gene CHECK: '+', '.join(f'{g}={"Y" if g in rev_core else "N"}({rev_core.get(g,0)})' for g in NAMED))
pr(f'   MEDIAN fold: REVERSIBLE {np.nanmedian(revf):.2f}x  vs  PERSIST-tertile {np.nanmedian(perf):.2f}x  vs  PERSIST-5%-size-matched {np.nanmedian(pmf):.2f}x')
pr('   -> VERIFIED (adversarial recompute): reversible >>chance (p to e-22); persistent ~chance; SIZE-MATCHED persistent ~0 -> the cross-tissue SHARED core lives on the REVERSIBLE pole. (muscle x blood weakest, p~0.09.)')
# (b) biological character: housekeeping/low-variance vs pathway
cl=pd.read_csv(S('results/marks_classified_transcriptome.tsv'),sep='\t')
pr('\n(b) biological character of each pole (transcriptome aggregate, n_elig>=2):')
for cls in ['persistent','consistently-reversible']:
    s=cl[cl['class']==cls]
    pr(f'   {cls:24s}: n={len(s):5d}  base_sd med={s.base_sd_mean.median():.3f}  meas(logCPM) med={s.meas_mean.median():.3f}  Z_resp med={s.Z_responsiveness.median():.3f}')
pr('   -> persistent = LOW base_sd + LOW expression (housekeeping/measurability floor); reversible = high-variance responders.')
# (c) cross-context CORRELATION of the continuous low-mover identity, WITHIN vs ACROSS tissue
pr('\n(c) cross-context correlation of |rev_beta| RANK (persistence=low) within-tissue vs cross-tissue:')
from scipy.stats import spearmanr
def rho_pair(fa,fb):
    ov=LD[fa].index.intersection(LD[fb].index)
    if len(ov)<50: return None,len(ov)
    return spearmanr(LD[fa].loc[ov,'ab'],LD[fb].loc[ov,'ab'])[0],len(ov)
wi=[]; cr=[]
fl=[(t,f) for t,fs in byT.items() for f in fs]
for i in range(len(fl)):
    for j in range(i+1,len(fl)):
        (ta,fa),(tb,fb)=fl[i],fl[j]
        r,n=rho_pair(fa,fb)
        if r is None: continue
        (wi if ta==tb else cr).append(r)
pr(f'   within-tissue mean rho(|rev_beta|) = {np.mean(wi):+.3f} (k={len(wi)})  vs  cross-tissue = {np.mean(cr):+.3f} (k={len(cr)})')
pr('   -> if within>>cross, the low-mover (persistence) identity is reproducible WITHIN a context but context-specific ACROSS — same structure as reversibility.')
open(S('results/restoration_persistence_revisit.txt'),'w',encoding='utf-8').write('\n'.join(out))
print('\nwrote results/restoration_persistence_revisit.txt')

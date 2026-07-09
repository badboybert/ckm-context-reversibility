#!/usr/bin/env python
"""build_local_features.py — mark-intrinsic features available LOCALLY (genetic architecture +
CKM-causal tier) for the gene-keyed determinant layers (proteome, transcriptome). External
gene features (LOEUF, druggability, disease-assoc, tau, functional) come from the acquisition
workflow -> features/*.tsv. Causal status enters as a GRADED FEATURE (leave-eGFR-out), never a filter.
Outputs signature-pivot/data/features/{proteome,transcriptome}_local.tsv.
"""
import os, glob, numpy as np, pandas as pd
SIG=os.path.dirname(os.path.dirname(os.path.abspath(__file__))); PROJ=os.path.dirname(SIG)
def P(*a): return os.path.join(PROJ,*a)
def S(*a): return os.path.join(SIG,*a)
NONEGFR=['T2D','NAFLD','CAD','AF','stroke','HF','CKD']   # leave-eGFR-out causal tier

# ---------- PROTEOME (key UniProt; gene = cis_gene) ----------
pm=pd.read_csv(P('data/proteome/trackA_pqtl_causal_master.tsv'),sep='\t')
for o in NONEGFR+['eGFR']:
    pm[f'{o}_causal']=pm[f'{o}_causal'].astype(str).str.lower().isin(['true','1','1.0'])
pm['causal_nonEGFR']=pm[[f'{o}_causal' for o in NONEGFR]].sum(axis=1)
nsig=pd.read_csv(P('data/proteome/ukbppp_cis_instruments.tsv'),sep='\t').groupby('UniProt').size().rename('n_cis_signals')
prot=pm.groupby('UniProt').agg(gene=('cis_gene','first'),cis_F=('F','max'),
     pav=('pav_flag','max'),causal_nonEGFR=('causal_nonEGFR','max')).reset_index()
prot['pav']=prot['pav'].astype(str).str.lower().isin(['true','1','1.0']).astype(int)
prot['causal_any']=(prot.causal_nonEGFR>0).astype(int)
prot=prot.merge(nsig,left_on='UniProt',right_index=True,how='left')
prot['n_cis_signals']=prot['n_cis_signals'].fillna(1)
prot.to_csv(S('data/features/proteome_local.tsv'),sep='\t',index=False)
print(f'proteome_local: {len(prot)} proteins | causal_any={int(prot.causal_any.sum())} pav={int(prot.pav.sum())} '
      f'medF={prot.cis_F.median():.0f}')

# ---------- TRANSCRIPTOME (key gene symbol) ----------
GCST2OUT={'GCST006414':'AF','GCST90103634':'eGFR','GCST90104540':'stroke','GCST90162626':'HF',
          'GCST90476128':'CKD','ebi-a-GCST006867':'T2D','NAFLD':'NAFLD','ebi-a-GCST003116':'CAD'}
caus={}; Fmax={}
for f in glob.glob(P('data/transcriptome/liver_causal_*_final.tsv')):
    key=os.path.basename(f).replace('liver_causal_','').replace('_final.tsv',''); out=GCST2OUT.get(key,key)
    d=pd.read_csv(f,sep='\t')
    isc=d['is_causal'].astype(str).str.lower().isin(['true','1','1.0']) if 'is_causal' in d else (pd.to_numeric(d.get('q_fdr',1),errors='coerce')<0.10)
    caus[out]=set(d.loc[isc,'mark_id'])
    for _,r in d[['mark_id','F']].dropna().iterrows(): Fmax[r['mark_id']]=max(Fmax.get(r['mark_id'],0),r['F'])
allg=set().union(*caus.values())
tcaus=pd.DataFrame([(g,sum(1 for o,s in caus.items() if o!='eGFR' and g in s)) for g in allg],
                   columns=['gene','causal_nonEGFR'])
# architecture from GTEx-liver eQTL credible sets
cs=pd.read_csv(P('data/eqtl/QTD000266.credible_sets.tsv.gz'),sep='\t')
arch=cs.groupby('gene_id').agg(eqtl_n_cs=('cs_id','nunique'),eqtl_max_pip=('pip','max')).reset_index()
arch['ensg_base']=arch['gene_id'].str.split('.').str[0]
br=pd.read_csv(P('data/transcriptome/liver_causal_table.tsv'),sep='\t')[['mark_id','ensg']].drop_duplicates()
br['ensg_base']=br['ensg'].astype(str).str.split('.').str[0]
arch=arch.merge(br,on='ensg_base',how='inner').groupby('mark_id').agg(
     eqtl_n_cs=('eqtl_n_cs','max'),eqtl_max_pip=('eqtl_max_pip','max')).reset_index().rename(columns={'mark_id':'gene'})
tx=arch.merge(tcaus,on='gene',how='outer')
tx['causal_nonEGFR']=tx['causal_nonEGFR'].fillna(0); tx['causal_any']=(tx.causal_nonEGFR>0).astype(int)
tx['has_eqtl']=tx.eqtl_n_cs.notna().astype(int)
tx['eqtl_F']=tx['gene'].map(Fmax)
tx.to_csv(S('data/features/transcriptome_local.tsv'),sep='\t',index=False)
cset=set(tx.loc[tx.causal_any==1,'gene']); pc=[g for g in ['SORT1','LPA','TSPAN8','CELSR2','PSRC1'] if g in cset]
print(f'transcriptome_local: {len(tx)} genes | has_eqtl={int(tx.has_eqtl.sum())} causal_any={int(tx.causal_any.sum())} | positive-control causal genes present: {pc}')

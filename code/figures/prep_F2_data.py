#!/usr/bin/env python
"""prep_F2_data.py — tidy per-panel source-data CSVs for Figure 2 (tissue-governor + reversible-pole asymmetry).
Reads verified engine outputs in results/; writes source_data/F2{a,b,c,d}_*.csv. No re-analysis; pure reshape."""
import os, pandas as pd
HERE=os.path.dirname(os.path.abspath(__file__)); SIG=os.path.dirname(os.path.dirname(HERE))
RES=os.path.join(SIG,'results'); OUT=os.path.join(HERE,'source_data'); os.makedirs(OUT,exist_ok=True)
R=lambda p: pd.read_csv(os.path.join(RES,p),sep='\t')

# a) within- vs cross-tissue pair correlations (distribution)
pair=R('tissue_pair_correlations.tsv')
pair[['kind','tissue','rho']].to_csv(os.path.join(OUT,'S9_pair_correlations.csv'),index=False)

# b) per-tissue cohesion (verified values from tissue_organization_summary.txt)
coh=pd.DataFrame({'tissue':['liver','adipose','muscle','blood'],
                  'mean_within_r':[0.407,0.190,0.056,-0.035],'n_contexts':[3,8,3,4]})
coh.to_csv(os.path.join(OUT,'F2c_cohesion.csv'),index=False)

# c) reversible-pole asymmetry: fold-enrichment of cross-tissue overlap, reversible vs persistent, both layers
rev=R('cross_tissue_reversible_pairs.tsv')          # transcriptome reversible folds (6 pairs)
meth=R('cross_tissue_methylation_persistence_edgecond.tsv')  # methylome 3 pairs (edge-conditioned)
rows=[]
for _,x in rev.iterrows():
    rows.append({'layer':'transcriptome','pair':f"{x.tissue_a}×{x.tissue_b}",'pole':'reversible','fold':x.fold})
    rows.append({'layer':'transcriptome','pair':f"{x.tissue_a}×{x.tissue_b}",'pole':'persistent','fold':0.0})  # size-matched top-5%, all 0
for _,x in meth.iterrows():
    rows.append({'layer':'methylome','pair':x['pair'],'pole':'reversible','fold':x.rev_ec_fold})
    rows.append({'layer':'methylome','pair':x['pair'],'pole':'persistent','fold':x.pers_ec_fold})
pd.DataFrame(rows).to_csv(os.path.join(OUT,'F2d_pole_asymmetry.csv'),index=False)

# d) universal reversible core (100 genes); flag the cross-shared (>=2 pairs) + named inflammatory members
core=R('universal_reversible_core.tsv')
named=['THBS1','CD83','PLAUR','CCL2','SPP1','SAA1','IL18','ACLY','ANGPTL4','PFKFB3']
core['named']=core['gene'].isin(named)
core.to_csv(os.path.join(OUT,'F2e_universal_core.csv'),index=False)
print('F2 source data written:', len(pair),'pairs |', len(rows),'asym rows |', len(core),'core genes')

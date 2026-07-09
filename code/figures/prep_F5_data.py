#!/usr/bin/env python
"""prep_F5_data.py — source-data for Figure 5 (intervention acts through physiology; drug-class; dissociation).
Reads verified engine outputs + raw proteome reversal tables (direction-only). Writes source_data/F5*."""
import os, re, numpy as np, pandas as pd
HERE=os.path.dirname(os.path.abspath(__file__)); SIG=os.path.dirname(os.path.dirname(HERE)); PROJ=os.path.dirname(SIG)
RES=os.path.join(SIG,'results'); OUT=os.path.join(HERE,'source_data'); os.makedirs(OUT,exist_ok=True)
R=lambda p: pd.read_csv(os.path.join(RES,p),sep='\t')

# a) intervention systemic vs tissue-specific (cross-tissue consensus gene count)
iv=R('intervention_systemic_consensus.tsv')
iv['kind']=np.where(iv.n_systemic_consensus>0,'systemic','tissue-specific')
iv[['intervention','n_panels','n_tissues','n_systemic_consensus','kind']].to_csv(
    os.path.join(OUT,'F5a_intervention_systemic.csv'),index=False)

# b) semaglutide tracks weight loss, not drug class
dc=R('drug_class_corr.tsv')
def get(a,b):
    r=dc[((dc.a==a)&(dc.b==b))|((dc.a==b)&(dc.b==a))]
    return (float(r.iloc[0].rho),float(r.iloc[0].p),int(r.iloc[0].n)) if len(r) and pd.notna(r.iloc[0].rho) else (np.nan,np.nan,0)
comp=[('BBS','surgery','cross'),('DiRECT','diet','same'),('MS_bariatric','surgery','cross'),
      ('metformin','drug (biguanide)','same'),('empagliflozin','drug (SGLT2i)','cross')]
rows=[]
for c,cls,plat in comp:
    rho,p,n=get('semaglutide',c)
    rows.append({'comparator':{'MS_bariatric':'MS bariatric','empagliflozin':'SGLT2i','metformin':'metformin'}.get(c,c),
                 'class':cls,'platform':plat,'rho':rho,'p':p,'n':n,'sig':p<0.05})
pd.DataFrame(rows).to_csv(os.path.join(OUT,'F5b_semaglutide_corr.csv'),index=False)
# platform control (from drug_class_robustness): same +0.128 (k=8) vs cross +0.146 (k=11)
pd.DataFrame([{'platform':'same','mean_rho':0.128,'k':8},{'platform':'cross','mean_rho':0.146,'k':11}]
            ).to_csv(os.path.join(OUT,'F5b_platform_control.csv'),index=False)

# c) layered plasma core: marker x intervention direction grid (rev_beta sign; direction-only harmonization)
panels={'surgery\n(BBS)':'reversal_olink_bbs.tsv','diet\n(DiRECT)':'reversal_somascan_direct.tsv',
        'SGLT2i\n(EMPEROR)':'reversal_emperor.tsv','metformin':'nowak_metformin_proteome_reversal.tsv',
        'sema-\nglutide':'step_semaglutide_proteome_reversal.tsv'}
markers={'P08833':'IGFBP1','P18065':'IGFBP2','P05231':'IL6','P41159':'LEP','P02786':'TFRC'}  # GDF15 omitted (guard: diet timepoint-confound)
grid=[]
for lab,f in panels.items():
    d=pd.read_csv(os.path.join(PROJ,'data/proteome',f),sep='\t')
    d['k']=d['UniProt'].astype(str).str.split(r'[;,_]').str[0]
    for uid,g in markers.items():
        r=d[d['k']==uid]
        b=float(r.iloc[0]['rev_beta']) if len(r) else np.nan
        grid.append({'intervention':lab,'marker':g,'rev_beta':b})
pd.DataFrame(grid).to_csv(os.path.join(OUT,'F5c_layered_core.csv'),index=False)

# e) intervention-pair reversal-correlation matrix (surgery/diet agree; drugs disagree)
IVN=['BBS','MS_bariatric','DiRECT','semaglutide','metformin','empagliflozin']  # powered set (n>=30 pairs)
LAB={'BBS':'BBS (surg)','MS_bariatric':'MS (surg)','DiRECT':'DiRECT (diet)','semaglutide':'sema',
     'metformin':'metformin','empagliflozin':'SGLT2i'}
def getrho(a,b):
    if a==b: return 1.0
    r=dc[((dc.a==a)&(dc.b==b))|((dc.a==b)&(dc.b==a))]
    return float(r.iloc[0].rho) if len(r) and pd.notna(r.iloc[0].rho) else np.nan
mrows=[]
for a in IVN:
    for b in IVN:
        mrows.append({'a':LAB[a],'b':LAB[b],'rho':getrho(a,b)})
pd.DataFrame(mrows).to_csv(os.path.join(OUT,'F5e_intervention_matrix.csv'),index=False)

# d) clinical dissociation (caveated, hypothesis-generating)
hd=R('responder_dissociation_hardened.tsv'); bo=R('responder_dissociation_bootstrap.tsv').set_index('metric')['value']
dd=hd[hd.stratum.isin(['remission','no_remission'])][['stratum','n_pairs','median_abs_sig05']].copy()
dd['label']={'remission':'remitter','no_remission':'non-remitter'}[dd.stratum.iloc[0]] if False else dd.stratum.map({'remission':'remitter','no_remission':'non-remitter'})
dd.to_csv(os.path.join(OUT,'F5d_dissociation_magnitude.csv'),index=False)
pd.DataFrame([{'resembles_surgery_no_remission':float(bo['spearman_no_remission_vs_no_diabetes']),
              'resembles_surgery_remission':float(bo['spearman_remission_vs_no_diabetes']),
              'boot_ci_lo':float(bo['boot_diff_ci_lo']),'boot_ci_hi':float(bo['boot_diff_ci_hi']),
              'boot_diff_mean':float(bo['boot_diff_mean']),'point_diff':float(bo['diff_nor_minus_rem_point'])}]
            ).to_csv(os.path.join(OUT,'F5d_dissociation_stats.csv'),index=False)
print('F5 source data written.')
print(' intervention systemic:', dict(zip(iv.intervention, iv.n_systemic_consensus)))
print(' sema corr:', [(r['comparator'],round(r['rho'],2)) for r in rows])

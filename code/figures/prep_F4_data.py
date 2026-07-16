#!/usr/bin/env python
"""prep_F4_data.py — tidy source-data for Figure 4 (restoration toward lean, scaled by intensity; durability).
Reads verified engine outputs; writes source_data/F4{a,b,c,d}_*.csv. Reshape + verified-constant values only."""
import os, re, numpy as np, pandas as pd
from scipy import stats
HERE=os.path.dirname(os.path.abspath(__file__)); SIG=os.path.dirname(os.path.dirname(HERE)); PROJ=os.path.dirname(SIG)
RES=os.path.join(SIG,'results'); OUT=os.path.join(HERE,'source_data'); os.makedirs(OUT,exist_ok=True)
R=lambda p: pd.read_csv(os.path.join(RES,p),sep='\t')

# a) restoration distribution: R_restoration (toward-lean) for proteome + metabolome + canonical marker points
pro=R('restoration_score_proteome.tsv'); met=R('restoration_score_metabolome.tsv')
dist=pd.concat([pro[['R_restoration']].assign(layer='proteome'),
                met[['R_restoration']].assign(layer='metabolome')])
# clip extreme tails for plotting only (medians/stats computed on full set, recorded below)
dist.to_csv(os.path.join(OUT,'F4a_restoration_dist.csv'),index=False)
canon={'P08833':'IGFBP1','P41159':'LEP','P15090':'FABP4','P02741':'CRP','P18065':'IGFBP2'}
rows=[]
for uid,g in canon.items():
    r=pro[pro.mark_key==uid]
    if len(r): x=r.iloc[0]; rows.append({'gene':g,'R':x.R_restoration,'ci_lo':x.R_ci_lo,'ci_hi':x.R_ci_hi,'k':x.k_panels,'ci_excl0':bool(x.restorative)})
pd.DataFrame(rows).to_csv(os.path.join(OUT,'F4a_canonical_markers.csv'),index=False)
summ=pd.DataFrame({'layer':['proteome','metabolome'],
                   'median_R':[pro.R_restoration.median(),met.R_restoration.median()],
                   'pct_restorative':[100*pro.restorative.mean(),100*met.restorative.mean()]})
summ.to_csv(os.path.join(OUT,'F4a_summary.csv'),index=False)

# b) %toward-lean per context. Values transcribed from the canonical engines:
#   plasma + adipose -> results/restoration_persistence_revisit.txt PART C (printed to 0 dp)
#   liver + blood    -> results/restoration_liver_blood.txt (printed to 1 dp: 53.2 / 54.4 / 51.5)
# NB blood is 51.5 (canonical 51.47%), NOT 52 — it had been double-rounded up (51.5 -> 52), which
# nudged the one near-chance context the round-2 review asked about away from the 50% line.
ctx=[ # context, tissue, intervention, pct_toward_lean, weight-loss tier (for ordering/shading)
 ('DiRECT','plasma','diet',63,'high'),('BBS','plasma','surgery',59,'high'),
 ('semaglutide','plasma','GLP-1RA',57,'mod'),('MS bariatric','plasma','surgery',55,'high'),
 ('metformin','plasma','drug',54,'low'),('EMPEROR','plasma','SGLT2i',52,'low'),
 ('RYGB','adipose','surgery',68,'high'),('CR','adipose','diet',57,'high'),
 ('diet','adipose','diet',54,'mod'),('LCD','adipose','diet',52,'mod'),
 ('GSE83452','liver','diet+surgery',53.2,'high'),('GSE106737','liver','surgery',54.4,'high'),
 ('GSE273902','blood','surgery',51.5,'high')]
_b=pd.DataFrame(ctx,columns=['context','tissue','intervention','pct_toward_lean','wl_tier'])
# Round-2 review (Figure 4 row): "mark near-chance contexts". ns_vs_chance flags the contexts whose
# restoration is NOT distinguishable from the 50% chance line on an exact binomial test — all four are
# the small-n plasma/blood panels. Source: results/restoration_uncertainty/restoration_uncertainty.tsv
# (binom_p_vs_50): MS bariatric 0.124, metformin 0.225, empagliflozin/EMPEROR 0.314, GSE273902 0.101;
# every other context has P < 1e-3. Read from the engine so the flag cannot drift from its evidence.
_u=pd.read_csv(os.path.join(SIG,'results','restoration_uncertainty','restoration_uncertainty.tsv'),sep='\t')
_KEY={'DiRECT':'DiRECT (plasma)','BBS':'By-Band-Sleeve (plasma)','semaglutide':'semaglutide (plasma)',
      'MS bariatric':'MS bariatric (plasma)','metformin':'metformin (plasma)','EMPEROR':'empagliflozin (plasma)',
      'RYGB':'RYGB (adipose)','CR':'CR (adipose)','diet':'diet GSE77962 (adipose)','LCD':'LCD DiOGenes (adipose)',
      'GSE83452':'GSE83452 (liver)','GSE106737':'GSE106737 (liver)','GSE273902':'GSE273902 (blood)'}
_p=_u.set_index('context')['binom_p_vs_50'].to_dict()
_b['binom_p_vs_50']=_b['context'].map(lambda c:_p[_KEY[c]])
_b['ns_vs_chance']=_b['binom_p_vs_50']>=0.05
assert set(_b.loc[_b.ns_vs_chance,'context'])=={'MS bariatric','metformin','EMPEROR','GSE273902'}, \
    f"ns set changed: {sorted(_b.loc[_b.ns_vs_chance,'context'])}"
_b.to_csv(os.path.join(OUT,'F4b_pct_toward_lean.csv'),index=False)

# c) durability: adipose 2yr vs 5yr (RYGB), CANONICAL genes only.
# GSE199063 = Affymetrix Clariom-D (GPL23126); ~32% of probes carry AceView novel-transcript names
# (lowercase, "Transcript Identified by AceView") = novel/non-coding, not interpretable genes. We restrict
# durability to annotated HGNC symbols (drops AceView novel + lncRNA/pseudogene noise classes).
NC=re.compile(r'-AS\d?|^LINC|orf\d|^RNU|^SNOR|^MIR\d|^LOC|^AC\d|^AL\d|^RP[0-9LS]|^IG[HKL][VJC]')
d=R('durability_GSE199063.tsv').dropna(subset=['rev2','rev5'])
g=d['gene'].astype(str)
d=d[g.str.match(r'^[A-Z][A-Z0-9-]*$') & ~g.str.contains(NC,na=False)].copy()   # canonical only
sp_all=stats.spearmanr(d.rev2,d.rev5).correlation
sig=d[d.p2<0.05].copy()                     # '96%' = sign-retained among p<0.05 canonical genes
sign_ret=float((np.sign(sig.rev5)==np.sign(sig.rev2)).mean())
samp=sig.sample(min(4000,len(sig)),random_state=1)[['gene','rev2','rev5']]
samp.to_csv(os.path.join(OUT,'F4c_durability_adipose_points.csv'),index=False)
pd.DataFrame([{'spearman_all':sp_all,'n_all':len(d),'sign_retained_sig':sign_ret*100,'n_sig':len(sig),
              'prot_spearman_2_12yr':0.6706,'prot_sign_retained':83.74,'prot_mag_ratio':1.0389,'prot_n_sig':910}]
            ).to_csv(os.path.join(OUT,'F4c_durability_summary.csv'),index=False)

# e) proteome 2yr->12yr durability (Yousri RYGB SomaLogic) — robust-z scatter (12-year horizon)
def _loadY(f):
    d=pd.read_csv(os.path.join(PROJ,'data/proteome',f),sep='\t')
    d['UniProt']=d['UniProt'].astype(str).str.split(r'[;,_]').str[0]
    return d.sort_values('rev_p',na_position='last').drop_duplicates('UniProt').set_index('UniProt')
y2=_loadY('yousri_rygb_2yr_reversal.tsv'); y12=_loadY('yousri_rygb_12yr.tsv')
ov=y2.index.intersection(y12.index)
def _rz(s): return s/(1.4826*(s-s.median()).abs().median())
yj=pd.DataFrame({'gene':y2.loc[ov,'gene'].astype(str),
                 'z2':_rz(y2.loc[ov,'rev_beta']).values,'z12':_rz(y12.loc[ov,'rev_beta']).values}).dropna()
yj.sample(min(1200,len(yj)),random_state=2).to_csv(os.path.join(OUT,'F4e_yousri_points.csv'),index=False)

# d) methylation restoration = age-drift-confounded caution (axis-correlation rises with sampling interval)
ad=R('orphan/methylation_age_drift.tsv')
ad['interval_mo']=ad['panel'].map({'drift2':3,'directplus':18})
ad['cohort']=ad['panel'].map({'drift2':'DRIFT2','directplus':'DIRECT-PLUS'})
ad[['cohort','interval_mo','spearman_delta_bmi']].to_csv(os.path.join(OUT,'F4d_methyl_age_drift.csv'),index=False)
print('F4 source data written. medR pro %.3f met %.3f | adipose dur %.3f sign %.1f%% | prot dur 0.671/83.7%%/1.04x'%(
    pro.R_restoration.median(),met.R_restoration.median(),sp_all,sign_ret*100))

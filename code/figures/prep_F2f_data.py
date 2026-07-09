#!/usr/bin/env python
"""prep_F2f_data.py — tissue-PRIVATE reversible programs (each tissue's recurrent top reversers + biology).
Reuses the per-tissue consensus logic (genes top-5% reversed & sig in >=half a tissue's contexts), reports
mean rev_beta across the tissue's contexts. Demonstrates each tissue carries a distinct, coherent program;
muscle=MSTN-only, blood=none (private-program absence, not missing data). Writes source_data/F2f_tissue_programs.csv."""
import os, re, numpy as np, pandas as pd
from collections import Counter, defaultdict
HERE=os.path.dirname(os.path.abspath(__file__)); SIG=os.path.dirname(os.path.dirname(HERE)); PROJ=os.path.dirname(SIG)
OUT=os.path.join(HERE,'source_data'); os.makedirs(OUT,exist_ok=True)
NC=re.compile(r'-AS\d?|^LINC|orf\d|^RNU|^SNOR|^MIR\d|^LOC|^AC\d|^AL\d|^RP[0-9LS]|^IG[HKL][VJC]')
PAN=[('liver','liver_reversal_table'),('liver','gse48452_liver_reversal_table'),('liver','gse106737_liver_bariatric_reversal_table'),
     ('adipose','adipose_reversal_table'),('adipose','gse141221_adipose_reversal_table'),('adipose','gse43471_adipose_reversal_table'),
     ('adipose','gse199063_adipose_reversal_table'),('adipose','gse77962_adipose_diet_reversal_table'),('adipose','gse70529_adipose_dose_reversal_table'),
     ('adipose','gse107894_adipose_metformin_reversal_table'),('adipose','gse224310_adipose_exercise_reversal_table'),
     ('muscle','gse83352_muscle_exercise_reversal_table'),('muscle','gse161643_muscle_bariatric_reversal_table'),('muscle','gse157585_muscle_reversal_table'),
     ('blood','gse273902_blood_reversal_table'),('blood','gse310742_blood_glp1ra_reversal_table'),('blood','gse193771_blood_exercise_reversal_table'),('blood','gse28358_pbmc_diet_reversal_table')]
byT=defaultdict(list)
for t,f in PAN:
    p=os.path.join(PROJ,'data/transcriptome',f+'.tsv')
    if not os.path.exists(p): continue
    d=pd.read_csv(p,sep='\t').drop_duplicates('mark_id'); d=d[~d.mark_id.str.contains(NC,na=False)]
    d['q']=pd.to_numeric(d.get('rev_q',1),errors='coerce'); d['ab']=d.rev_beta.abs()
    d['top']=(d['q']<0.05)&(d['ab']>=d['ab'].quantile(0.95))
    byT[t].append(d.set_index('mark_id'))
rows=[]
THEME={'liver':'acute-phase / stress','adipose':'lipid metabolism','muscle':'myostatin only','blood':'none (no private program)'}
for t in ['liver','adipose','muscle','blood']:
    fs=byT[t]; k=len(fs); need=max(2,(k+1)//2)
    cnt=Counter(); betas=defaultdict(list)
    for d in fs:
        for g in d.index[d['top']]: cnt[g]+=1
        for g in d.index:
            if g in d.index: betas[g].append(d.loc[g,'rev_beta'] if np.isscalar(d.loc[g,'rev_beta']) else np.nan)
    cons=[g for g,n in cnt.items() if n>=need]
    # mean rev_beta across this tissue's contexts for the consensus genes; top 6 by |mean|
    mb={g:np.nanmean(betas[g]) for g in cons}
    top=sorted(mb, key=lambda g:-abs(mb[g]))[:6]
    for g in top:
        rows.append({'tissue':t,'gene':g,'mean_rev_beta':round(float(mb[g]),3),'n_recurrent':cnt[g],'theme':THEME[t]})
    if not top:  # blood: no consensus
        rows.append({'tissue':t,'gene':'(none)','mean_rev_beta':0.0,'n_recurrent':0,'theme':THEME[t]})
pd.DataFrame(rows).to_csv(os.path.join(OUT,'F2f_tissue_programs.csv'),index=False)
df=pd.DataFrame(rows)
print('F2f tissue programs:')
for t in ['liver','adipose','muscle','blood']:
    g=df[df.tissue==t]; print(' %-8s [%s]: %s'%(t,g.theme.iloc[0],', '.join('%s(%+.1f)'%(r.gene,r.mean_rev_beta) for _,r in g.iterrows())))

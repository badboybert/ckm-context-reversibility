#!/usr/bin/env python
"""prep_ED5_data.py — Extended Data 5: cross-layer (meth-vs-expr) node-concordance NO-GO.
Re-emits the PER-GENE meth-Delta vs expr-Delta vectors (the committed cross_layer_concordance.tsv has only
2 summary rows). Reuses src/cross_layer_concordance.py logic verbatim. Recomputes Spearman to confirm
muscle -0.19/-0.35 and blood -0.01 ns. Writes source_data/ED5_{muscle,blood}_scatter.csv + ED5_summary.csv."""
import os, csv, collections, numpy as np, pandas as pd
from scipy import stats
HERE=os.path.dirname(os.path.abspath(__file__)); OUT=os.path.join(HERE,'source_data'); os.makedirs(OUT,exist_ok=True)
BASE=os.path.dirname(os.path.dirname(os.path.dirname(HERE)))  # weight-loss.omics
CPG_GENE_MAP=BASE+"/data/methylation/cpg_gene_map.tsv"
PAIRS={"muscle":dict(meth=BASE+"/data/methylation/gse60655_muscle_exercise_reversal.tsv",
                     expr=BASE+"/data/transcriptome/gse60655_muscle_expression_reversal_table.tsv",
                     meth_acc="GSE60655(450K)",expr_acc="GSE60590(RNA-seq)"),
       "blood":dict(meth=BASE+"/data/methylation/gse193730_blood_exercise_reversal.tsv",
                    expr=BASE+"/data/transcriptome/gse193771_blood_exercise_reversal_table.tsv",
                    meth_acc="GSE193730(EPIC)",expr_acc="GSE193771(RNA-seq)")}   # GEO GSE193771 = Expression profiling by high throughput sequencing (GPL11154); corrected from "array" 2026-07-17
nk=lambda s:s.strip().upper()
def load_map(p):
    m={}
    with open(p) as fh:
        for r in csv.DictReader(fh,delimiter="\t"):
            cpg=r["cpg"].strip(); gf=(r.get("genes") or "").strip()
            if cpg and gf: m[cpg]={nk(g) for g in gf.split(";") if g.strip()}
    return m
def load_meth(p):
    d={}
    with open(p) as fh:
        for r in csv.DictReader(fh,delimiter="\t"):
            try: d[r["cpg"].strip()]=float(r["delta_beta"])
            except: pass
    return d
def load_expr(p):
    d={}
    with open(p) as fh:
        for r in csv.DictReader(fh,delimiter="\t"):
            try: b=float(r["rev_beta"])
            except: continue
            try: pv=float(r["rev_p"])
            except: pv=float("nan")
            d[nk(r["mark_id"])]=(b,pv)
    return d
cpg2g=load_map(CPG_GENE_MAP)
summ=[]
for name,cfg in PAIRS.items():
    meth=load_meth(cfg["meth"]); expr=load_expr(cfg["expr"])
    acc=collections.defaultdict(list)
    for cpg,db in meth.items():
        for g in cpg2g.get(cpg,()): acc[g].append(db)
    gmeth={g:float(np.mean(v)) for g,v in acc.items()}
    genes=sorted(set(gmeth)&set(expr))
    df=pd.DataFrame({'gene':genes,'meth_delta':[gmeth[g] for g in genes],
                     'expr_delta':[expr[g][0] for g in genes],'expr_p':[expr[g][1] for g in genes]})
    df['expr_sig']=df.expr_p<0.05
    df.to_csv(os.path.join(OUT,f'ED5_{name}_scatter.csv'),index=False)
    rho_all=stats.spearmanr(df.meth_delta,df.expr_delta).correlation
    s=df[df.expr_sig]; rho_sig=stats.spearmanr(s.meth_delta,s.expr_delta).correlation if len(s)>3 else float('nan')
    summ.append({'pair':name,'meth':cfg['meth_acc'],'expr':cfg['expr_acc'],'n_all':len(df),
                 'rho_all':round(rho_all,4),'n_sig':int(df.expr_sig.sum()),'rho_sig':round(rho_sig,4)})
pd.DataFrame(summ).to_csv(os.path.join(OUT,'ED5_summary.csv'),index=False)
for s in summ: print('%-7s n=%d rho_all=%+.3f | sig n=%d rho_sig=%+.3f'%(s['pair'],s['n_all'],s['rho_all'],s['n_sig'],s['rho_sig']))

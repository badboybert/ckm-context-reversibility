import csv, glob, os
from collections import defaultdict
ROOT=os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))  # weight-loss.omics
DP=os.path.join(ROOT,"data"); SIG=os.path.join(ROOT,"signature-pivot"); PAPERA=os.path.join(ROOT,"4-layers null")
OUTDIR=os.path.join(SIG,"results"); os.makedirs(OUTDIR,exist_ok=True)
GCST2OUT={'GCST006414':'AF','GCST90103634':'eGFR','GCST90104540':'stroke','GCST90162626':'HF',
          'GCST90476128':'CKD','ebi-a-GCST006867':'T2D','NAFLD':'NAFLD','ebi-a-GCST003116':'CAD'}
NONEGFR=['T2D','NAFLD','CAD','AF','stroke','HF','CKD']
def rd(p):
    with open(p,encoding="utf-8") as f: return list(csv.DictReader(f,delimiter="\t"))

# ---- Paper B substrate ----
txloc=rd(os.path.join(SIG,"data/features/transcriptome_local.tsv"))
tx_causal_val={r["gene"]:float(r["causal_nonEGFR"]) for r in txloc if float(r["causal_nonEGFR"])>=1}
prloc=rd(os.path.join(SIG,"data/features/proteome_local.tsv"))
pr_causal_val={r["UniProt"]:(float(r["causal_nonEGFR"]),r["gene"],r["n_cis_signals"],r["pav"]) for r in prloc if float(r["causal_nonEGFR"])>=1}
pr107={u for u,(v,_,_,_) in pr_causal_val.items() if v==1}
pr16={u for u,(v,_,_,_) in pr_causal_val.items() if v>=2}
print("tx causal (>=1):",len(tx_causal_val)," prot causal_any(>=1):",len(pr_causal_val),
      "| ==1:",len(pr107)," >=2:",len(pr16))

# ---- transcriptome is_causal detail ----
tx_causal={}
for f in glob.glob(os.path.join(DP,"transcriptome","liver_causal_*_final.tsv")):
    key=os.path.basename(f).replace("liver_causal_","").replace("_final.tsv","")
    out=GCST2OUT.get(key,key)
    for r in rd(f):
        if str(r.get("is_causal")).lower() in ("true","1","1.0"):
            tx_causal.setdefault(r["mark_id"],{})[out]=r
txcoloc={}
for r in rd(os.path.join(DP,"transcriptome","coloc_results.tsv")):
    txcoloc[(r["symbol"],r["outcome"])]=r; txcoloc[(r["ensg"],r["outcome"])]=r

# ---- proteome trackA ----
trackA=rd(os.path.join(DP,"proteome","trackA_pqtl_causal_master.tsv"))
prot_by_uni=defaultdict(list)
for r in trackA: prot_by_uni[r["UniProt"]].append(r)
prcoloc={}
for r in rd(os.path.join(DP,"proteome","coloc_results_primary.tsv")):
    prcoloc[(r["UniProt"],r["outcome"])]=r

# ---- panel membership (determinant panels) ----
def keyset(path,key):
    s=set()
    with open(os.path.join(ROOT,path),encoding="utf-8") as f:
        rd2=csv.DictReader(f,delimiter="\t")
        for r in rd2:
            s.add(str(r.get(key)))
    return s
TXPANELS={"adipose_table":"data/transcriptome/adipose_reversal_table.tsv",
 "liver_table":"data/transcriptome/liver_reversal_table.tsv",
 "gse141221_adipose_table":"data/transcriptome/gse141221_adipose_reversal_table.tsv",
 "gse273902_blood_table":"data/transcriptome/gse273902_blood_reversal_table.tsv",
 "gse157585_muscle_table":"data/transcriptome/gse157585_muscle_reversal_table.tsv",
 "gse43471_adipose_table":"data/transcriptome/gse43471_adipose_reversal_table.tsv",
 "gse83352_muscle_exercise_table":"data/transcriptome/gse83352_muscle_exercise_reversal_table.tsv",
 "gse106737_liver_bariatric_table":"data/transcriptome/gse106737_liver_bariatric_reversal_table.tsv",
 "gse77962_adipose_diet_table":"data/transcriptome/gse77962_adipose_diet_reversal_table.tsv"}
PRPANELS={"reversal_olink_bbs":"data/proteome/reversal_olink_bbs.tsv",
 "reversal_somascan_direct":"data/proteome/reversal_somascan_direct.tsv",
 "reversal_emperor":"data/proteome/reversal_emperor.tsv"}
txpanel_keys={p:keyset(pa,"mark_id") for p,pa in TXPANELS.items()}
prpanel_keys={p:keyset(pa,"UniProt") for p,pa in PRPANELS.items()}
def tx_npanel(g): return sum(1 for s in txpanel_keys.values() if g in s)
def pr_npanel(u): return sum(1 for s in prpanel_keys.values() if u in s)

# ---- S2 sets ----
S2=rd(os.path.join(PAPERA,"supp/S2_nominated_features.tsv"))
s2_prot={r["feature"] for r in S2 if r["layer"]=="Proteome"}
s2_prot_uni={r["stable_id"] for r in S2 if r["layer"]=="Proteome"}
s2_tx_coloc={r["feature"] for r in S2 if r["layer"]=="Transcriptome" and r["tier"].startswith("coloc")}
s2_tx_mr={r["feature"] for r in S2 if r["layer"]=="Transcriptome" and r["tier"].startswith("MR")}
s2_tx_all=s2_tx_coloc|s2_tx_mr

def h2r(pp4,pp3):
    try:
        pp4=float(pp4);pp3=float(pp3); return "%.3g"%(pp4/pp3) if pp3>0 else "Inf"
    except: return "NA"

rows=[]
# transcriptome
for g in sorted(tx_causal_val):
    npan=tx_npanel(g)
    for o,r in sorted(tx_causal[g].items()):
        if o=="eGFR": continue
        c=txcoloc.get((g,o)) or txcoloc.get((r["ensg"],o))
        pp4=c["PP4"] if c else ""; pp3=c["PP3"] if c else ""
        s2flag="Y (coloc-tier)" if g in s2_tx_coloc else ("Y (MR-tier)" if g in s2_tx_mr else "N")
        rows.append([g,"Transcriptome",o,r["rsid"],"1 (finemapped lead)","%.1f"%float(r["F"]),
            "%+.4f"%float(r["mr_beta"]),"%.2e"%float(r["mr_p"]),"%.2e"%float(r["q_fdr"]),
            ("%.3f"%float(pp4)) if pp4 else "NA (not coloc-staged)",
            h2r(pp4,pp3) if pp4 else "NA","%.3f"%float(r["pip"]),
            "cis-MR Bonferroni-p & Steiger-fwd, F>=10 (GTEx-liver eQTL); is_causal",
            int(tx_causal_val[g]),"Y","%d/9"%npan,s2flag,r["mhc"]])
# proteome (all 123 causal_any)
for uni,(val,gene,ncs,pav) in sorted(pr_causal_val.items(), key=lambda x:x[1][1]):
    npan=pr_npanel(uni)
    assays=prot_by_uni.get(uni,[])
    outs=set()
    for a in assays:
        for o in NONEGFR:
            if a.get(f"{o}_causal")=="1": outs.add(o)
    for o in sorted(outs):
        best=None
        for a in assays:
            if a.get(f"{o}_causal")=="1":
                qv=float(a[f"{o}_q"])
                if best is None or qv<best[0]: best=(qv,a)
        qv,a=best
        c=prcoloc.get((uni,o)); pp4=c["PP4"] if c else ""; pp3=c["PP3"] if c else ""
        s2flag="Y" if (gene in s2_prot or uni in s2_prot_uni) else "N"
        rows.append([gene,"Proteome",o,a["rsid"],str(int(float(ncs))),"%.1f"%float(a["F"]),
            "%+.4f"%float(a[f"{o}_mrb"]),"%.2e"%float(a[f"{o}_p"]),"%.2e"%float(a[f"{o}_q"]),
            ("%.3f"%float(pp4)) if pp4 else "NA (not coloc-staged)",
            h2r(pp4,pp3) if pp4 else "NA","NA",
            "cis-MR BH-FDR q<0.10 & non-MHC (UKB-PPP/Olink pQTL); Wald-ratio lead cis variant",
            int(val),"Y","%d/3"%npan,s2flag,a["mhc"]])

cols=["gene","layer","outcome","instrument_rsid","n_instr","F","MR_beta","MR_p","q",
      "coloc_PP_H4","PP_H4_over_H3","finemap_pip","tier","n_nonEGFR_causal_outcomes",
      "included_in_determinant","determinant_panels_measured","in_PaperA_S2","mhc_region"]
outpath=os.path.join(OUTDIR,"causal_status_supp_table.tsv")
with open(outpath,"w",newline="",encoding="utf-8") as f:
    w=csv.writer(f,delimiter="\t"); w.writerow(cols); w.writerows(rows)
print("WROTE",outpath,"rows",len(rows))
ntx=sum(1 for r in rows if r[1]=="Transcriptome"); npr=sum(1 for r in rows if r[1]=="Proteome")
print("  transcriptome rows",ntx," proteome rows",npr)
print("  unique proteome genes",len(set(r[0] for r in rows if r[1]=="Proteome")))

# reconciliation prints
print("\n== PROTEOME reconciliation ==")
prg={gene for uni,(val,gene,_,_) in pr_causal_val.items()}
print("causal_any UniProt 123; in S2(deCODE) by uni/gene:",
      len({u for u,(v,g,_,_) in pr_causal_val.items() if (u in s2_prot_uni or g in s2_prot)}))

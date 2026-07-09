#!/usr/bin/env python
"""
Cross-layer concordance: methylation-Delta vs expression-Delta, gene-level.

QUESTION
--------
Within a tissue, do the genes whose PROMOTER/BODY CpGs change methylation with an
exercise intervention move OPPOSITE in expression? Classic coordinated regulation
predicts a NEGATIVE Spearman (promoter hypomethylation -> up-regulation, and vice
versa), strongest among expression-significant genes. A null (~0) means the two
layers are uncoordinated at the gene level.

PAIRS (each a methylation arm x a SAME-tissue, same-intervention expression arm)
--------------------------------------------------------------------------------
MUSCLE (endurance training, vastus lateralis; Lindholm SuperSeries GSE60591,
        PMID 25484259):
  meth = data/methylation/gse60655_muscle_exercise_reversal.tsv         (450K Delta-beta)
  expr = data/transcriptome/gse60655_muscle_expression_reversal_table.tsv
         ^ this file is the GSE60590 RNA-seq EXPRESSION arm of the SAME study
           (its build script build_reversal_60590.py downloads GSE60590 FPKM and
            writes this table). GSE60655(meth) x GSE60590(expr) = the same-subject
            endurance pair, per the project memory.

BLOOD (exercise; pediatric exercise cohorts):
  meth = data/methylation/gse193730_blood_exercise_reversal.tsv
  expr = data/transcriptome/gse193771_blood_exercise_reversal_table.tsv

CAVEATS (persisted)
-------------------
1. NOT strictly same-subject. The muscle meth (GSE60655) and expr (GSE60590) arms
   come from the same SuperSeries/cohort but the gene-level Delta vectors are
   aggregated independently, not paired per individual. The blood meth (GSE193730)
   and expr (GSE193771) are DIFFERENT GEO accessions from related but not identical
   sample sets. So this is a cohort/tissue-level concordance, not a within-person test.
2. NOT promoter-resolved. The CpG->gene map (data/methylation/cpg_gene_map.tsv, the
   EPIC/450K UCSC_RefGene annotation) is gene-level (any CpG annotated to the gene),
   not restricted to promoter/TSS200 windows. CpG Delta-beta is aggregated to the
   gene by MEAN over all mapped CpGs. A true promoter-biased aggregation would sharpen
   the negative coupling; this mean-over-all-CpGs version is conservative.

METHOD
------
- methylation Delta = delta_beta (POST-PRE; >0 = hyper after intervention).
- expression  Delta = rev_beta  (POST-PRE log2; >0 = up after intervention).
- CpG -> gene via cpg_gene_map.tsv; ';'-delimited multi-gene rows are exploded so a
  CpG mapping to k genes contributes to all k.
- gene methylation Delta = MEAN(delta_beta over its CpGs).
- merge on gene symbol (case-normalized).
- Spearman(meth-Delta, expr-Delta): (a) OVERALL all merged genes; (b) among
  EXPRESSION-SIGNIFICANT genes (expr rev_p < 0.05).
- persist n_genes, both Spearman rho + p, per pair, plus the caveats.

OUTPUT
------
results/cross_layer_concordance.tsv
"""
import csv
import collections
import numpy as np
from scipy import stats

BASE = r"C:/Users/Bert/Downloads/CKM papers/weight-loss.omics"
CPG_GENE_MAP = BASE + "/data/methylation/cpg_gene_map.tsv"
OUT = BASE + "/signature-pivot/results/cross_layer_concordance.tsv"

PAIRS = {
    "MUSCLE": dict(
        meth=BASE + "/data/methylation/gse60655_muscle_exercise_reversal.tsv",
        expr=BASE + "/data/transcriptome/gse60655_muscle_expression_reversal_table.tsv",
        meth_acc="GSE60655(450K)", expr_acc="GSE60590(RNA-seq)",
    ),
    "BLOOD": dict(
        meth=BASE + "/data/methylation/gse193730_blood_exercise_reversal.tsv",
        expr=BASE + "/data/transcriptome/gse193771_blood_exercise_reversal_table.tsv",
        meth_acc="GSE193730(EPIC)", expr_acc="GSE193771(array)",
    ),
}

EXPR_SIG_P = 0.05  # expression-significant = nominal rev_p < 0.05


def nk(s):
    return s.strip().upper()


def load_cpg_gene_map(path):
    """cpg -> list of gene symbols (explode ';'-delimited, dedup)."""
    m = {}
    with open(path) as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        for r in rd:
            cpg = r["cpg"].strip()
            genes_field = (r.get("genes") or "").strip()
            if not cpg or not genes_field:
                continue
            genes = {nk(g) for g in genes_field.split(";") if g.strip()}
            if genes:
                m[cpg] = genes
    return m


def load_meth(path):
    """cpg -> delta_beta."""
    d = {}
    with open(path) as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        for r in rd:
            try:
                d[r["cpg"].strip()] = float(r["delta_beta"])
            except (KeyError, ValueError, TypeError):
                pass
    return d


def load_expr(path):
    """gene -> (rev_beta, rev_p)."""
    d = {}
    with open(path) as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        for r in rd:
            try:
                beta = float(r["rev_beta"])
            except (KeyError, ValueError, TypeError):
                continue
            try:
                p = float(r["rev_p"])
            except (KeyError, ValueError, TypeError):
                p = float("nan")
            d[nk(r["mark_id"])] = (beta, p)
    return d


def gene_meth_delta(meth, cpg2genes):
    """gene -> mean(delta_beta over its mapped CpGs)."""
    acc = collections.defaultdict(list)
    for cpg, db in meth.items():
        genes = cpg2genes.get(cpg)
        if not genes:
            continue
        for g in genes:
            acc[g].append(db)
    return {g: float(np.mean(v)) for g, v in acc.items()}


def spearman(x, y):
    if len(x) < 3:
        return float("nan"), float("nan")
    rho, p = stats.spearmanr(x, y)
    return float(rho), float(p)


def main():
    print("loading CpG->gene map ...")
    cpg2genes = load_cpg_gene_map(CPG_GENE_MAP)
    print("  mapped CpGs:", len(cpg2genes))

    rows = []
    for name, cfg in PAIRS.items():
        print(f"\n=== {name} ===")
        meth = load_meth(cfg["meth"])
        expr = load_expr(cfg["expr"])
        print(f"  meth CpGs: {len(meth)} | expr genes: {len(expr)}")

        gmeth = gene_meth_delta(meth, cpg2genes)
        print(f"  genes with meth Delta: {len(gmeth)}")

        # merge on gene
        genes = sorted(set(gmeth) & set(expr))
        mvals = np.array([gmeth[g] for g in genes])
        evals = np.array([expr[g][0] for g in genes])
        epvals = np.array([expr[g][1] for g in genes])
        n_genes = len(genes)
        print(f"  merged genes (meth & expr): {n_genes}")

        rho_all, p_all = spearman(mvals, evals)

        sig = np.isfinite(epvals) & (epvals < EXPR_SIG_P)
        n_sig = int(sig.sum())
        rho_sig, p_sig = spearman(mvals[sig], evals[sig])
        print(f"  Spearman OVERALL    rho={rho_all:+.4f} p={p_all:.3g} (n={n_genes})")
        print(f"  Spearman EXPR-SIG   rho={rho_sig:+.4f} p={p_sig:.3g} (n={n_sig}, expr rev_p<{EXPR_SIG_P})")

        rows.append(dict(
            pair=name,
            meth_dataset=cfg["meth_acc"],
            expr_dataset=cfg["expr_acc"],
            n_genes_overall=n_genes,
            spearman_overall=round(rho_all, 4),
            spearman_overall_p=p_all,
            n_genes_expr_sig=n_sig,
            spearman_expr_sig=round(rho_sig, 4),
            spearman_expr_sig_p=p_sig,
            expr_sig_threshold=f"rev_p<{EXPR_SIG_P}",
            caveat_same_subject="cohort/tissue-level, NOT within-person paired",
            caveat_promoter="gene-level CpG->gene (UCSC_RefGene), NOT promoter/TSS-resolved; CpG Delta aggregated by mean",
        ))

    fields = list(rows[0].keys())
    with open(OUT, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow(r)
    print("\nWROTE", OUT)


if __name__ == "__main__":
    main()

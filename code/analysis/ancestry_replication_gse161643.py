#!/usr/bin/env python3
"""
Ancestry-replication check (small-n, NOT a powered ancestry contrast).

GSE161643 = muscle bariatric reversal, BLACK FEMALES, n=12 -- the only
non-European transcriptome panel on disk. Compare its muscle reversal
program to the EUROPEAN muscle panels:
  GSE157585 (MASTERS, bariatric/aging-master-athletes muscle, n=47)
  GSE83352  (STRRIDE, exercise muscle, n=42)

Outputs:
  - Spearman(rev_beta) on shared genes (BLACK vs each EUR panel, and vs EUR mean)
  - sig-gene overlap fold-enrichment (rev_q<0.1) vs hypergeometric expectation
  - directional concordance on shared sig genes
  - canonical muscle-consensus replication: MSTN suppression; contractile/ECM
    (MYH1, MYH2, COL1A1, COL3A1, PPARGC1A) directions across cohorts

Convention: rev_beta>0 = mark UP after intervention.
Honest framing: n=12 is too weak to *conclude* an ancestry difference; this is a
replication-direction check only.
"""
import sys
import numpy as np
import pandas as pd
from scipy import stats

import os
_WLO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))  # weight-loss.omics
DATA = os.path.join(_WLO, "data", "transcriptome")
OUT = os.path.join(_WLO, "signature-pivot", "results", "ancestry_replication_gse161643.tsv")

PANELS = {
    "GSE161643_BLACK_muscle_bariatric": f"{DATA}/gse161643_muscle_bariatric_reversal_table.tsv",
    "GSE157585_EUR_muscle_masters":     f"{DATA}/gse157585_muscle_reversal_table.tsv",
    "GSE83352_EUR_muscle_strride":      f"{DATA}/gse83352_muscle_exercise_reversal_table.tsv",
}
BLACK = "GSE161643_BLACK_muscle_bariatric"
EUR = ["GSE157585_EUR_muscle_masters", "GSE83352_EUR_muscle_strride"]

QSIG = 0.10  # significance threshold on rev_q (matches project convention elsewhere)
CANON = ["MSTN", "MYH1", "MYH2", "COL1A1", "COL3A1", "PPARGC1A"]


def load(path):
    df = pd.read_csv(path, sep="\t")
    df = df[["mark_id", "rev_beta", "rev_q", "n_pairs"]].copy()
    df = df.dropna(subset=["rev_beta"])
    df = df.drop_duplicates(subset=["mark_id"])
    return df.set_index("mark_id")


def spearman_pair(a, b):
    shared = a.index.intersection(b.index)
    if len(shared) < 10:
        return len(shared), np.nan, np.nan
    rho, p = stats.spearmanr(a.loc[shared, "rev_beta"], b.loc[shared, "rev_beta"])
    return len(shared), rho, p


def sig_overlap(a, b):
    """Hypergeometric fold-enrichment of co-significant genes + directional concordance."""
    shared = a.index.intersection(b.index)
    N = len(shared)
    sa = set(shared[a.loc[shared, "rev_q"] < QSIG])
    sb = set(shared[b.loc[shared, "rev_q"] < QSIG])
    k = len(sa & sb)
    expected = len(sa) * len(sb) / N if N > 0 else np.nan
    fold = k / expected if expected and expected > 0 else np.nan
    # hypergeometric p (P[X>=k])
    if N > 0 and len(sa) > 0 and len(sb) > 0:
        hp = stats.hypergeom.sf(k - 1, N, len(sa), len(sb))
    else:
        hp = np.nan
    # directional concordance among co-sig genes
    co = list(sa & sb)
    if co:
        conc = np.mean(
            np.sign(a.loc[co, "rev_beta"].values) == np.sign(b.loc[co, "rev_beta"].values)
        )
    else:
        conc = np.nan
    return N, len(sa), len(sb), k, expected, fold, hp, conc


def main():
    P = {name: load(path) for name, path in PANELS.items()}
    rows = []

    # --- 1. pairwise Spearman + sig overlap: BLACK vs each EUR ---
    for eur in EUR:
        n_sh, rho, sp = spearman_pair(P[BLACK], P[eur])
        N, sa, sb, k, exp, fold, hp, conc = sig_overlap(P[BLACK], P[eur])
        rows.append({
            "comparison": f"{BLACK} vs {eur}",
            "metric": "spearman_rev_beta",
            "n_shared": n_sh, "rho": round(rho, 4), "p": sp,
            "black_n_sig": sa, "eur_n_sig": sb, "co_sig": k,
            "exp_co_sig": round(exp, 3) if exp == exp else np.nan,
            "fold_enrich": round(fold, 3) if fold == fold else np.nan,
            "hypergeom_p": hp, "dir_concord_cosig": round(conc, 4) if conc == conc else np.nan,
            "value": round(rho, 4), "note": "BLACK vs EUR panel",
        })

    # --- 2. BLACK vs EUR consensus (mean rev_beta over EUR panels on shared genes) ---
    eur_shared = P[EUR[0]].index.intersection(P[EUR[1]].index)
    eur_mean = pd.DataFrame({
        "rev_beta": (P[EUR[0]].loc[eur_shared, "rev_beta"] + P[EUR[1]].loc[eur_shared, "rev_beta"]) / 2,
        # consensus-sig = significant in BOTH EUR panels
        "rev_q": np.where(
            (P[EUR[0]].loc[eur_shared, "rev_q"] < QSIG) & (P[EUR[1]].loc[eur_shared, "rev_q"] < QSIG),
            0.0, 1.0),
    }, index=eur_shared)
    n_sh, rho, sp = spearman_pair(P[BLACK], eur_mean)
    N, sa, sb, k, exp, fold, hp, conc = sig_overlap(P[BLACK], eur_mean)
    rows.append({
        "comparison": f"{BLACK} vs EUR_muscle_consensus(mean)",
        "metric": "spearman_rev_beta",
        "n_shared": n_sh, "rho": round(rho, 4), "p": sp,
        "black_n_sig": sa, "eur_n_sig": sb, "co_sig": k,
        "exp_co_sig": round(exp, 3) if exp == exp else np.nan,
        "fold_enrich": round(fold, 3) if fold == fold else np.nan,
        "hypergeom_p": hp, "dir_concord_cosig": round(conc, 4) if conc == conc else np.nan,
        "value": round(rho, 4), "note": "BLACK vs EUR consensus(mean over 2 EUR muscle panels)",
    })

    # --- 3. EUR vs EUR (internal benchmark: how reproducible are EUR panels themselves?) ---
    n_sh, rho, sp = spearman_pair(P[EUR[0]], P[EUR[1]])
    N, sa, sb, k, exp, fold, hp, conc = sig_overlap(P[EUR[0]], P[EUR[1]])
    rows.append({
        "comparison": f"{EUR[0]} vs {EUR[1]}",
        "metric": "spearman_rev_beta",
        "n_shared": n_sh, "rho": round(rho, 4), "p": sp,
        "black_n_sig": sa, "eur_n_sig": sb, "co_sig": k,
        "exp_co_sig": round(exp, 3) if exp == exp else np.nan,
        "fold_enrich": round(fold, 3) if fold == fold else np.nan,
        "hypergeom_p": hp, "dir_concord_cosig": round(conc, 4) if conc == conc else np.nan,
        "value": round(rho, 4), "note": "EUR-vs-EUR internal benchmark (ceiling for replication)",
    })

    # --- 4. canonical muscle-consensus replication table ---
    for g in CANON:
        rec = {"comparison": f"CANON::{g}", "metric": "rev_beta_by_cohort"}
        for name in PANELS:
            if g in P[name].index:
                rec[name + "__beta"] = round(float(P[name].loc[g, "rev_beta"]), 4)
                rec[name + "__q"] = float(P[name].loc[g, "rev_q"])
            else:
                rec[name + "__beta"] = np.nan
                rec[name + "__q"] = np.nan
        # replication verdict: does BLACK match EUR sign?
        bb = rec.get(BLACK + "__beta", np.nan)
        eb = [rec.get(e + "__beta", np.nan) for e in EUR]
        eb = [x for x in eb if x == x]
        if bb == bb and eb:
            eur_sign = np.sign(np.mean(eb))
            rec["black_matches_eur_sign"] = bool(np.sign(bb) == eur_sign)
        else:
            rec["black_matches_eur_sign"] = np.nan
        rec["note"] = "MSTN expect suppression(<0); MYH1/contractile+ECM directionality"
        rows.append(rec)

    out = pd.DataFrame(rows)
    out.to_csv(OUT, sep="\t", index=False)

    # --- console summary ---
    print("=== ANCESTRY-REPLICATION CHECK: GSE161643 (BLACK FEMALES, n=12) ===")
    print("Small-n replication-direction check; NOT a powered ancestry contrast.\n")
    for r in rows:
        if r["metric"] == "spearman_rev_beta":
            print(f"{r['comparison']}")
            print(f"   Spearman rho={r['rho']} (p={r['p']:.2e}), n_shared={r['n_shared']}")
            print(f"   sig overlap: black={r['black_n_sig']} eur={r['eur_n_sig']} co={r['co_sig']} "
                  f"fold={r['fold_enrich']} hyperP={r['hypergeom_p']:.2e} dirConcord={r['dir_concord_cosig']}")
    print("\n--- canonical muscle consensus (rev_beta; *=q<0.1) ---")
    hdr = f"{'gene':8s} " + " ".join([f"{n.split('_')[0]:>10s}" for n in PANELS])
    print(hdr + "   BLACK_matches_EUR_sign")
    for r in rows:
        if r["comparison"].startswith("CANON::"):
            g = r["comparison"].split("::")[1]
            cells = []
            for n in PANELS:
                b = r.get(n + "__beta", np.nan)
                q = r.get(n + "__q", np.nan)
                star = "*" if (q == q and q < QSIG) else " "
                cells.append(f"{b:>9.3f}{star}" if b == b else f"{'NA':>10s}")
            print(f"{g:8s} " + " ".join(cells) + f"   {r['black_matches_eur_sign']}")
    print(f"\nWritten: {OUT}")


if __name__ == "__main__":
    main()

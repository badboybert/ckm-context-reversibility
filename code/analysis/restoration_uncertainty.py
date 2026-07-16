#!/usr/bin/env python
"""restoration_uncertainty.py — ROUND-2 item: what uncertainty can honestly be put on the
"% of marks moving toward lean" statistic for the near-chance (52-54%) restoration contexts?

The round-2 review offered three options for the 52-54% contexts:
  (i) confidence intervals, (ii) binomial/permutation P values, or (iii) "near-random" wording.
This engine tests which instrument is admissible instead of asserting it, and reports the answer
even where it is unflattering.

Orientation is identical to the canonical engines (src/restoration_liver_blood.py,
src/restoration_persistence_revisit.py): toward_lean = -sign(bmi_beta) * rev_beta.

WHAT IS COMPUTED, per context
  1. pct_toward_lean = mean(toward_lean > 0).
     ANCHOR: reproduces every canonical published percentage to <0.5 pp. It is NOT bit-identical:
     this engine de-duplicates mark keys on both sides, the canonical proteome path de-duplicates
     only the axis, so n differs in some contexts (reported as n_delta_vs_canonical). No conclusion
     depends on it; the difference is disclosed rather than hidden.

  2. THE BINOMIAL INSTRUMENT (the review's option i/ii), vs the 50% line: Wilson 95% CI + exact P.

  3. A SECOND, MARGINAL-CONDITIONED BASELINE. Panels differ enormously in the marginal sign balance
     of their reversal vector (EMPEROR: 90% of plasma proteins rise; metformin: 84% fall;
     semaglutide: centred). Crossed with a sign-imbalanced BMI axis, a RANDOM pairing does not give
     50%. Permuting the BMI-axis sign labels (both margins fixed) has expectation available in
     closed form, p0 = P(b>0)P(r<0) + P(b<0)P(r>0) = 0.5*(1 - A*R), with A and R the sign
     asymmetries of the two vectors. Confirmed against a permutation null to <=0.02 pp.

     ★ p0 IS A SECOND ESTIMAND, NOT A CORRECTION TO 50%. P(rev>0) is a POST-TREATMENT quantity and
     is not ancillary: against a one-sided BMI axis a genuine global restorative shift IS a shift in
     the reversal marginal, so conditioning on it conditions away part of the signal of interest.
     The 50% line answers "did marks move toward lean more often than not?"; p0 answers "did they do
     so beyond what each vector's own marginal forces?". Both are reported. Only conclusions that
     hold under BOTH baselines are treated as robust (agrees_across_baselines).

  4. CENTRING SENSITIVITY. A global shift may be technical (a concentration/normalisation effect)
     rather than restorative. Every statistic is recomputed after removing each panel's median
     rev_beta. Conclusions that flip under centring are marked and are NOT relied on.

  5. WHY NO CORRELATION-RESPECTING INTERVAL IS AVAILABLE. The draft response letter offered, as a
     fallback, "permutation-based intervals that preserve the inter-mark correlation structure".
     That instrument does not exist here, and NOT for the reason of a simulation result: a gene-label
     permutation randomises gene identity, so the inter-mark correlation of the reversal vector
     cannot enter its null BY CONSTRUCTION. Its SD is deterministic given the margins:
         perm_SD / binom_SD = 4*sqrt(a(1-a)b(1-b)),  a = P(bmi>0), b = P(rev<0)
     (verified here to <0.03 across all contexts) — it equals ~1 whenever both margins are balanced,
     so the observed "ratio ~= 1" is an arithmetic identity, not evidence about independence.
     A correlation-respecting interval would need the replication unit the question is actually
     about — the SUBJECT — and these panels are published per-mark summary effects, so subject-level
     resampling is not possible for them. Hence: report the binomial interval for what it is, and
     describe the near-chance contexts as near-random (the review's option iii).

STANDALONE + CANONICAL READ-ONLY: reads only committed data + canonical axes; writes only to
results/restoration_uncertainty/. Touches no canonical output.
"""
import os, re, json
import numpy as np
import pandas as pd
from scipy import stats

SIG = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PROJ = os.path.dirname(SIG)
def P(*a): return os.path.join(PROJ, *a)
def S(*a): return os.path.join(SIG, *a)

OUT = S('results', 'restoration_uncertainty')
os.makedirs(OUT, exist_ok=True)
RNG = np.random.default_rng(11)
B = 10000

NC = re.compile(r'-AS\d?|^LINC|orf\d|^RNU|^SNOR|^MIR\d|^LOC|^AC\d|^AL\d|^RP[0-9LS]|^IG[HKL][VJC]|rRNA|^7SK')

# canonical published percentages (results/restoration_persistence_revisit.txt PART C;
# results/restoration_liver_blood.txt) — the anchor this engine must reproduce.
CANON = {'DiRECT (plasma)': 63, 'By-Band-Sleeve (plasma)': 59, 'semaglutide (plasma)': 57,
         'MS bariatric (plasma)': 55, 'metformin (plasma)': 54, 'empagliflozin (plasma)': 52,
         'RYGB (adipose)': 68, 'CR (adipose)': 57, 'diet GSE77962 (adipose)': 54,
         'LCD DiOGenes (adipose)': 52, 'GSE83452 (liver)': 53.2, 'GSE106737 (liver)': 54.4,
         'GSE273902 (blood)': 51.5}

def axis_gene(path):
    d = pd.read_csv(path, sep='\t')
    d['gene'] = d['gene'].astype(str).str.upper()
    d['bmi_beta'] = pd.to_numeric(d['bmi_beta'], errors='coerce')
    return d.dropna(subset=['bmi_beta']).drop_duplicates('gene').set_index('gene')['bmi_beta']

def axis_map(path, keycands):
    d = pd.read_csv(P(path), sep='\t')
    kc = next((c for c in keycands if c in d.columns), d.columns[0])
    bc = next((c for c in d.columns if re.search(r'bmi|beta|effect|axis|mr', c, re.I) and c != kc), None)
    d = d[[kc, bc]].dropna(); d.columns = ['key', 'bmi']
    d['bmi'] = pd.to_numeric(d['bmi'], errors='coerce')
    return d.dropna().drop_duplicates('key').set_index('key')['bmi']

ax_liver = axis_gene(S('data', 'features', 'bmi_axis_transcriptome_liver.tsv'))
ax_blood = axis_gene(S('data', 'features', 'bmi_axis_transcriptome_blood.tsv'))
ax_adip  = axis_map('signature-pivot/data/features/bmi_axis_transcriptome.tsv', ['gene', 'mark_id', 'symbol'])
ax_prot  = axis_map('signature-pivot/results/bmi_axis_proteome.tsv', ['UniProt', 'gene', 'assay_id'])

def vec_transcriptome(path, axis, drop_nc):
    d = pd.read_csv(P(path), sep='\t').rename(columns={'mark_id': 'gene'})
    d['gene'] = d['gene'].astype(str).str.upper()
    if drop_nc:
        d = d[~d['gene'].str.contains(NC, na=False)]
    d['rev_beta'] = pd.to_numeric(d['rev_beta'], errors='coerce')
    d = d.dropna(subset=['rev_beta']).drop_duplicates('gene')
    j = d.join(axis.rename('bmi'), on='gene', how='inner').dropna(subset=['bmi'])
    return j['rev_beta'].to_numpy(float), j['bmi'].to_numpy(float)

def vec_proteome(path, axis):
    d = pd.read_csv(P(path), sep='\t')
    kc = next((c for c in ['UniProt', 'gene', 'mark_id'] if c in d.columns), None)
    d = d.rename(columns={kc: 'key'})
    d['key'] = d['key'].astype(str).str.split(r'[;,_]').str[0]
    d['rev_beta'] = pd.to_numeric(d['rev_beta'], errors='coerce')
    d = d.dropna(subset=['rev_beta']).drop_duplicates('key')
    j = d.join(axis.rename('bmi'), on='key', how='inner').dropna(subset=['bmi'])
    return j['rev_beta'].to_numpy(float), j['bmi'].to_numpy(float)

CTX = [
    ('DiRECT (plasma)',         'plasma',  'diet',      lambda: vec_proteome('data/proteome/reversal_somascan_direct.tsv', ax_prot)),
    ('By-Band-Sleeve (plasma)', 'plasma',  'surgery',   lambda: vec_proteome('data/proteome/reversal_olink_bbs.tsv', ax_prot)),
    ('semaglutide (plasma)',    'plasma',  'GLP-1RA',   lambda: vec_proteome('data/proteome/step_semaglutide_proteome_reversal.tsv', ax_prot)),
    ('MS bariatric (plasma)',   'plasma',  'surgery',   lambda: vec_proteome('data/proteome/bariatric_ms_proteome_reversal.tsv', ax_prot)),
    ('metformin (plasma)',      'plasma',  'drug',      lambda: vec_proteome('data/proteome/nowak_metformin_proteome_reversal.tsv', ax_prot)),
    ('empagliflozin (plasma)',  'plasma',  'drug',      lambda: vec_proteome('data/proteome/reversal_emperor.tsv', ax_prot)),
    ('RYGB (adipose)',          'adipose', 'surgery',   lambda: vec_transcriptome('data/transcriptome/gse199063_adipose_reversal_table.tsv', ax_adip, False)),
    ('CR (adipose)',            'adipose', 'diet',      lambda: vec_transcriptome('data/transcriptome/adipose_reversal_table.tsv', ax_adip, False)),
    ('diet GSE77962 (adipose)', 'adipose', 'diet',      lambda: vec_transcriptome('data/transcriptome/gse77962_adipose_diet_reversal_table.tsv', ax_adip, False)),
    ('LCD DiOGenes (adipose)',  'adipose', 'diet',      lambda: vec_transcriptome('data/transcriptome/gse141221_adipose_reversal_table.tsv', ax_adip, False)),
    ('GSE83452 (liver)',        'liver',   'diet+surg', lambda: vec_transcriptome('data/transcriptome/liver_reversal_table.tsv', ax_liver, True)),
    ('GSE106737 (liver)',       'liver',   'surgery',   lambda: vec_transcriptome('data/transcriptome/gse106737_liver_bariatric_reversal_table.tsv', ax_liver, True)),
    ('GSE273902 (blood)',       'blood',   'surgery',   lambda: vec_transcriptome('data/transcriptome/gse273902_blood_reversal_table.tsv', ax_blood, True)),
]

def analytic_p0(sb, sr):
    """Closed-form permutation baseline: p0 = P(b>0)P(r<0) + P(b<0)P(r>0) = 0.5*(1 - A*R)."""
    return 100.0 * ((sb > 0).mean() * (sr < 0).mean() + (sb < 0).mean() * (sr > 0).mean())

def stat_block(rev, bmi):
    sb, sr = np.sign(bmi), np.sign(rev)
    n = len(rev)
    k = int(((-sb * rev) > 0).sum())
    pct = 100.0 * k / n
    p0 = analytic_p0(sb, sr)
    bt = stats.binomtest(k, n, 0.5)
    lo, hi = bt.proportion_ci(confidence_level=0.95, method='wilson')
    null = np.array([100.0 * ((-RNG.permutation(sb) * rev) > 0).mean() for _ in range(B)])
    nm = float(null.mean())
    perm_p = float((np.abs(null - nm) >= abs(pct - nm)).mean())
    a, b = (sb > 0).mean(), (sr < 0).mean()
    ratio_closed = 4.0 * np.sqrt(a * (1 - a) * b * (1 - b))
    return dict(n=n, pct=pct, p0=p0, perm_mean=nm, perm_sd=float(null.std(ddof=1)),
                perm_p=perm_p, binom_p=bt.pvalue, ci_lo=100 * lo, ci_hi=100 * hi,
                binom_sd=100 * np.sqrt(0.25 / n), ratio_closed=float(ratio_closed),
                pct_rev_up=100 * (sr > 0).mean(), pct_bmi_up=100 * (sb > 0).mean())

rows = []
print(f"{'context':26s} {'n':>6s} {'pct':>6s} {'binomP50':>9s} {'p0':>6s} {'permP':>7s} {'agree':>6s} {'centred pct/P':>16s}")
print('-' * 104)
for label, tissue, fam, load in CTX:
    rev, bmi = load()
    if len(rev) < 20:
        print(f'  {label:26s} SKIP (n={len(rev)})'); continue
    s = stat_block(rev, bmi)
    c = stat_block(rev - np.median(rev), bmi)          # centring sensitivity
    sig50, sigp0 = s['binom_p'] < 0.05, s['perm_p'] < 0.05
    agree = (sig50 == sigp0)
    canon = CANON.get(label)
    rows.append(dict(
        context=label, tissue=tissue, family=fam, n_marks=s['n'],
        pct_toward_lean=round(s['pct'], 2),
        canonical_published_pct=canon, anchor_delta_pp=round(s['pct'] - canon, 2) if canon else None,
        pct_rev_up=round(s['pct_rev_up'], 1), pct_bmi_up=round(s['pct_bmi_up'], 1),
        # baseline 1: the 50% line (the manuscript's)
        wilson_ci_lo=round(s['ci_lo'], 2), wilson_ci_hi=round(s['ci_hi'], 2),
        binom_p_vs_50=s['binom_p'], sig_vs_50=bool(sig50),
        # baseline 2: marginal-conditioned (a DIFFERENT estimand, reported for robustness only)
        chance_baseline_p0=round(s['p0'], 2), perm_null_mean=round(s['perm_mean'], 2),
        excess_over_p0_pp=round(s['pct'] - s['p0'], 2), perm_p_vs_p0=s['perm_p'], sig_vs_p0=bool(sigp0),
        agrees_across_baselines=bool(agree),
        # centring sensitivity
        centred_pct=round(c['pct'], 2), centred_binom_p_vs_50=c['binom_p'],
        centred_sig_vs_50=bool(c['binom_p'] < 0.05),
        flips_on_centring=bool((c['binom_p'] < 0.05) != sig50),
        # the perm/binom SD ratio is an arithmetic identity, kept only to document that
        perm_over_binom_sd=round(s['perm_sd'] / s['binom_sd'], 3),
        perm_over_binom_sd_closed_form=round(s['ratio_closed'], 3),
    ))
    print(f"  {label:24s} {s['n']:6d} {s['pct']:6.2f} {s['binom_p']:9.2e} {s['p0']:6.2f} {s['perm_p']:7.4f} "
          f"{str(agree):>6s} {c['pct']:8.2f}/{c['binom_p']:.3f}")

df = pd.DataFrame(rows)
df.to_csv(os.path.join(OUT, 'restoration_uncertainty.tsv'), sep='\t', index=False)

# closed-form identity check for the SD ratio
ratio_err = float((df.perm_over_binom_sd - df.perm_over_binom_sd_closed_form).abs().max())
anchor_err = float(df.anchor_delta_pp.abs().max())
nearc = df[df.pct_toward_lean < 55]
robust_ns = df[(~df.sig_vs_50) & (~df.sig_vs_p0)]

summary = {
    'B_permutations': B,
    'n_contexts': int(len(df)),
    'anchor_max_abs_delta_vs_canonical_pp': round(anchor_err, 2),
    'perm_binom_SD_ratio_closed_form_max_err': round(ratio_err, 3),
    'contexts_lt_55pct': sorted(nearc.context.tolist()),
    'NOT_distinguishable_from_50pct_line': {r.context: {'pct': r.pct_toward_lean,
        'wilson_ci': [r.wilson_ci_lo, r.wilson_ci_hi], 'binom_p': round(r.binom_p_vs_50, 3)}
        for r in df[~df.sig_vs_50].itertuples()},
    'ROBUST_ns_under_BOTH_baselines': sorted(robust_ns.context.tolist()),
    'baseline_disagreement': sorted(df[~df.agrees_across_baselines].context.tolist()),
    'flips_on_centring': sorted(df[df.flips_on_centring].context.tolist()),
    'FINDING_1': ('The review\'s option (i)/(ii) is admissible but WEAK, and it works against the '
                  'manuscript where it matters: on the exact binomial against the 50% line, '
                  f'{len(df[~df.sig_vs_50])} of {len(df)} contexts are NOT distinguishable from chance '
                  f'({", ".join(f"{r.context} {r.pct_toward_lean}%, P={r.binom_p_vs_50:.2f}" for r in df[~df.sig_vs_50].itertuples())}). '
                  'Every one of them is a small-n plasma or blood panel.'),
    'FINDING_2': ('A second baseline that additionally conditions on each panel\'s marginal sign '
                  f'balance (p0, range {df.chance_baseline_p0.min():.1f}-{df.chance_baseline_p0.max():.1f}%, closed form 0.5*(1-A*R), '
                  'confirmed against permutation) is a DIFFERENT estimand, not a correction: P(rev>0) is '
                  'post-treatment and not ancillary, so conditioning on it conditions away part of a '
                  'genuine global restorative shift. Only conclusions agreeing under BOTH baselines are '
                  f'relied on. Robustly not-distinguishable-from-chance under both: {sorted(robust_ns.context.tolist())}. '
                  f'Baseline-dependent (do NOT claim): {sorted(df[~df.agrees_across_baselines].context.tolist())}.'),
    'FINDING_3': ('Empagliflozin is the only context whose observed percentage falls below its own p0 '
                  '(90% of EMPEROR plasma proteins rise). That comparison is NOT relied on: it is '
                  'sensitive to whether the global rise is restorative biology or a concentration effect '
                  f'— under median-centring the contexts that change significance are {sorted(df[df.flips_on_centring].context.tolist())}. '
                  'The baseline-independent statement is the weaker one: empagliflozin is not '
                  'distinguishable from chance (P = '
                  f'{df.loc[df.context=="empagliflozin (plasma)","binom_p_vs_50"].iloc[0]:.2f} vs the 50% line; P = '
                  f'{df.loc[df.context=="empagliflozin (plasma)","perm_p_vs_p0"].iloc[0]:.2f} vs p0).'),
    'FINDING_4': ('No correlation-respecting interval is available for these contexts, and the '
                  '"correlation-preserving permutation intervals" offered as a fallback in the draft '
                  'response letter cannot be delivered as described. A gene-label permutation randomises '
                  'gene identity, so inter-mark correlation cannot enter its null BY CONSTRUCTION; its SD '
                  'is deterministic given the margins, perm_SD/binom_SD = 4*sqrt(a(1-a)b(1-b)) '
                  f'(closed form matches simulation to {ratio_err:.3f} here), equal to ~1 only when both '
                  'margins are balanced. The observed ratio is therefore an arithmetic identity, not '
                  'evidence. A correlation-respecting interval needs subject-level resampling, and these '
                  'panels are published per-mark summary effects.'),
}
with open(os.path.join(OUT, 'restoration_uncertainty_summary.json'), 'w', encoding='utf-8') as f:
    json.dump(summary, f, indent=2)

print(f'\nanchor: max |delta vs canonical| = {anchor_err:.2f} pp')
print(f'SD-ratio closed form vs simulation: max err {ratio_err:.3f}')
print('\n--- summary ---')
print(json.dumps(summary, indent=2))
print(f'\nwrote {OUT}')

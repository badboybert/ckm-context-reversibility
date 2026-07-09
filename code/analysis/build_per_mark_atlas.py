#!/usr/bin/env python
"""build_per_mark_atlas.py  (TASK B9) -- reviewer-accessible per-mark reproducibility atlas.

STANDALONE. Does NOT modify build_determinants.py nor overwrite any canonical result.
Reads existing summary-level outputs and emits ONE tidy per-mark table:
    results/per_mark_atlas.tsv

One row = one MARK per LAYER (the marks_classified aggregation granularity), joining:
  * results/marks_classified_{transcriptome,methylation,metabolome}.tsv
        -> reversible/persistent CLASS, n_elig/n_persistent/n_reversible, base_sd_mean, Z_responsiveness
        (proteome has NO marks_classified file: proteome panels are persistence_eligible=False in
         panel_manifest, so the frozen persistence spec assigns NO persistent/reversible class. Proteome
         rows are therefore sourced from results/reversibility_score_proteome.tsv and carry class=NA.)
  * results/reversibility_score_{layer}.tsv
        -> signed pooled effect (pooled_signed), score (Z_responsiveness), Z CI, I2, sign_concord,
           contributing panels, k_panels.  se_or_q = SE of pooled Z derived from the committed CI:
           SE = (Z_ci_hi - Z_ci_lo) / (2 * 1.959964).
  * intrinsic gene/UniProt features (GRADED, never a filter):
        data/features/loeuf.tsv (loeuf), data/features/tissue_tau.tsv (tau),
        data/features/{transcriptome,proteome}_local.tsv (architecture arch_nsig/arch_str + causal_nonEGFR).
        gene-keyed features only populate for transcriptome (mark_key = gene symbol) and proteome
        (UniProt -> gene via proteome_local); methylation (CpG) and metabolome have NA -- expected.
  * results/context_signatures.tsv
        -> in_context_signature flag + comma list of contexts where the mark is a top-ranked signature
        (matched on gene symbol; available for transcriptome + proteome only).
  * results/panel_manifest_full.tsv is the panel-metadata companion for the `panel` lists (not row-joined;
     see results/panel_manifest_full.tsv for platform/intervention/tissue/n_pairs per panel).

DUA: every input is already SUMMARY-LEVEL (per-mark statistics: pooled Z, betas, base_sd, counts).
No individual-level / per-subject values are read or written.

Baseline anchor: rows with source=='marks_classified' == exact sum of the three marks_classified
files, per layer (no silent drop/duplication). Proteome rows are ADDITIVE and separately counted.
"""
import os, numpy as np, pandas as pd

SIG = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
def S(*a): return os.path.join(SIG, *a)
Z = 1.959964

# ---------------- intrinsic feature tables (gene / UniProt keyed) ----------------
loeuf = pd.read_csv(S('data/features/loeuf.tsv'), sep='\t')[['gene', 'loeuf']]
loeuf = loeuf.dropna(subset=['gene']).drop_duplicates('gene').set_index('gene')['loeuf']
tau = pd.read_csv(S('data/features/tissue_tau.tsv'), sep='\t')[['gene', 'tau']]
tau = tau.dropna(subset=['gene']).groupby('gene')['tau'].max()

tx_loc = pd.read_csv(S('data/features/transcriptome_local.tsv'), sep='\t')
tx_loc = tx_loc.dropna(subset=['gene']).drop_duplicates('gene').set_index('gene')
prot_loc = pd.read_csv(S('data/features/proteome_local.tsv'), sep='\t')
prot_loc = prot_loc.drop_duplicates('UniProt').set_index('UniProt')

# ---------------- context signatures (gene-symbol keyed, per layer) --------------
ctx = pd.read_csv(S('results/context_signatures.tsv'), sep='\t')
ctx_map = (ctx.groupby(['layer', 'mark'])['context']
           .apply(lambda s: ','.join(sorted(set(s)))).to_dict())

REV_COLS = ['mark_key', 'k_panels', 'n_total_pairs', 'pooled_signed', 'Z_responsiveness',
            'Z_ci_lo', 'Z_ci_hi', 'I2', 'panels', 'sign_concord']

OUT_COLS = ['mark_id', 'gene_or_uniprot_or_cpg_or_metabolite', 'gene', 'layer',
            'panel', 'k_panels', 'n_elig', 'effect', 'se_or_q', 'score',
            'z_ci_lo', 'z_ci_hi', 'i2', 'sign_concord', 'reversible_persistent_label',
            'n_persistent', 'n_reversible', 'base_sd_or_abundance', 'meas_mean',
            'eligibility_flags', 'loeuf', 'tau', 'arch_nsig', 'arch_str', 'causal_status',
            'in_context_signature', 'context_signatures', 'source']

REAL_CLASSES = {'persistent', 'consistently-reversible', 'intermediate'}


def rev_score(layer):
    df = pd.read_csv(S(f'results/reversibility_score_{layer}.tsv'), sep='\t', usecols=REV_COLS)
    df['se_or_q'] = (df['Z_ci_hi'] - df['Z_ci_lo']) / (2 * Z)
    return df.rename(columns={'k_panels': 'k_panels', 'pooled_signed': 'effect',
                              'Z_responsiveness': 'score', 'Z_ci_lo': 'z_ci_lo',
                              'Z_ci_hi': 'z_ci_hi', 'I2': 'i2', 'panels': 'panel'})


def attach_gene_features(d, gene_col):
    d['loeuf'] = d[gene_col].map(loeuf)
    d['tau'] = d[gene_col].map(tau)
    return d


def build_classified(layer, gene_is_markkey=True):
    """transcriptome / methylation / metabolome: rows defined by marks_classified (anchor)."""
    mc = pd.read_csv(S(f'results/marks_classified_{layer}.tsv'), sep='\t')
    rv = rev_score(layer)
    d = mc.merge(rv, on='mark_key', how='left')
    d['mark_id'] = d['mark_key']
    d['gene_or_uniprot_or_cpg_or_metabolite'] = d['mark_key']
    d['reversible_persistent_label'] = d['class']
    d['base_sd_or_abundance'] = d['base_sd_mean']
    d['n_elig'] = d['n_elig']
    d['source'] = 'marks_classified'
    # gene resolution + gene-keyed intrinsic features
    if layer == 'transcriptome':
        d['gene'] = d['mark_key']
        d = attach_gene_features(d, 'gene')
        d['arch_nsig'] = d['gene'].map(tx_loc['eqtl_n_cs'])
        d['arch_str'] = d['gene'].map(tx_loc['eqtl_max_pip'])
        d['causal_status'] = d['gene'].map(tx_loc['causal_nonEGFR'])
    else:  # methylation (CpG) / metabolome -> no gene mapping in these inputs
        d['gene'] = np.nan
        for c in ['loeuf', 'tau', 'arch_nsig', 'arch_str', 'causal_status']:
            d[c] = np.nan
    # context signatures (gene-symbol keyed)
    d['context_signatures'] = [ctx_map.get((layer, g)) for g in d['mark_key']]
    d['in_context_signature'] = d['context_signatures'].notna()
    # eligibility flags
    multi = d['n_elig'] >= 2
    classified = d['class'].isin(REAL_CLASSES)
    d['eligibility_flags'] = np.where(multi, 'multi_panel', 'single_panel_excluded')
    d['eligibility_flags'] = d['eligibility_flags'] + np.where(classified, ';persistence_classified', '')
    return d.reindex(columns=OUT_COLS)


def build_proteome():
    """No marks_classified (persistence_eligible=False). Rows from reversibility_score_proteome."""
    d = rev_score('proteome')
    d['layer'] = 'proteome'
    d['mark_id'] = d['mark_key']
    d['gene_or_uniprot_or_cpg_or_metabolite'] = d['mark_key']  # UniProt
    d['gene'] = d['mark_key'].map(prot_loc['gene'])
    d = attach_gene_features(d, 'gene')
    d['arch_nsig'] = d['mark_key'].map(prot_loc['n_cis_signals'])
    d['arch_str'] = np.log10(d['mark_key'].map(prot_loc['cis_F']).clip(lower=1))
    d['causal_status'] = d['mark_key'].map(prot_loc['causal_nonEGFR'])
    d['reversible_persistent_label'] = np.nan   # not persistence-eligible
    d['n_elig'] = np.nan
    d['n_persistent'] = np.nan
    d['n_reversible'] = np.nan
    d['base_sd_or_abundance'] = np.nan           # proteome has no base_sd
    d['meas_mean'] = np.nan
    d['source'] = 'reversibility_score_proteome'
    # context signatures keyed on gene symbol for proteome
    d['context_signatures'] = [ctx_map.get(('proteome', g)) if pd.notna(g) else None
                               for g in d['gene']]
    d['in_context_signature'] = d['context_signatures'].notna()
    d['eligibility_flags'] = np.where(d['k_panels'] >= 2, 'multi_panel', 'single_panel') + ';proteome_no_persistence_class'
    return d.reindex(columns=OUT_COLS)


parts = {
    'transcriptome': build_classified('transcriptome'),
    'methylation':   build_classified('methylation'),
    'metabolome':    build_classified('metabolome'),
    'proteome':      build_proteome(),
}
atlas = pd.concat(parts.values(), ignore_index=True)
atlas.to_csv(S('results/per_mark_atlas.tsv'), sep='\t', index=False)
# compressed companion (atlas is large: ~0.9M rows, methylation-dominated)
try:
    atlas.to_parquet(S('results/per_mark_atlas.parquet'), index=False, compression='zstd')
except Exception as e:  # engine/codec missing -> TSV remains the primary deliverable
    print(f'(parquet companion skipped: {e})')

# ------------------------------- report -----------------------------------------
CLASSIFIED_LAYERS = ['transcriptome', 'methylation', 'metabolome']
mc_counts = {L: len(pd.read_csv(S(f'results/marks_classified_{L}.tsv'), sep='\t', usecols=['mark_key']))
             for L in CLASSIFIED_LAYERS}
mc_sum = sum(mc_counts.values())
atlas_mc = int((atlas['source'] == 'marks_classified').sum())

print('=== per_mark_atlas.tsv written ===')
print(f'total rows           : {len(atlas)}')
print(f'columns              : {len(atlas.columns)}')
print('rows by source       :', dict(atlas['source'].value_counts()))
print('\n=== BASELINE ANCHOR (marks_classified: no silent drop/duplication) ===')
for L in CLASSIFIED_LAYERS:
    n = int(((atlas.layer == L) & (atlas.source == 'marks_classified')).sum())
    print(f'  {L:14s}: atlas={n:>7d}  marks_classified={mc_counts[L]:>7d}  match={n==mc_counts[L]}')
print(f'  SUM classified atlas rows = {atlas_mc}  ==  sum(marks_classified) = {mc_sum}  -> {atlas_mc==mc_sum}')
n_prot = int((atlas.source == 'reversibility_score_proteome').sum())
print(f'  proteome (additive, no marks_classified) = {n_prot}')

print('\n=== unique marks per layer ===')
for L in ['transcriptome', 'proteome', 'methylation', 'metabolome']:
    sub = atlas[atlas.layer == L]
    print(f'  {L:14s}: rows={len(sub):>7d}  unique mark_id={sub.mark_id.nunique():>7d}  dup={int(sub.mark_id.duplicated().sum())}')

print('\n=== column coverage (non-NA / total, NA count) ===')
N = len(atlas)
for c in atlas.columns:
    nna = int(atlas[c].isna().sum())
    print(f'  {c:38s} non-NA={N-nna:>7d}  NA={nna:>7d}  ({100*(N-nna)/N:5.1f}%)')

# integrity: no full-row duplicates on (layer, mark_id)
dups = int(atlas.duplicated(subset=['layer', 'mark_id']).sum())
print(f'\n(layer,mark_id) duplicate rows: {dups}')

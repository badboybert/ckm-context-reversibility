# A multi-omic atlas of weight-loss molecular reversibility: context, not molecular identity (CKM)

Analysis code and frozen source data for:

> **A multi-omic atlas defines context, not molecular identity, as the determinant of weight-loss reversibility.**

The study asks whether the reversibility of an obesity-associated molecular mark with weight loss is an intrinsic property of the mark or a property of its biological context. Across a 40-panel human multi-omic intervention atlas (plasma/serum proteome, tissue transcriptome, blood/tissue methylome, serum metabolome; four tissues and plasma; six intervention types), each mark is summarized into a per-mark reversibility score and regressed on mark-intrinsic features (genetic constraint, druggability, disease-association burden, cis-regulatory architecture, tissue-specificity and genetically nominated causal status) and against contextual determinants (tissue, intervention, molecular layer). The powered transcriptome determinant analysis finds no continuous mark-intrinsic feature predicts a practically meaningful amount of reversibility; genetically nominated causal status is reported on the interpretable raw scale as an underpowered contrast. Reversibility is instead organized by shared biological context, is reproducible within and across cohorts, restores toward the lean state in proportion to weight loss, and is durable where it occurs.

## Repository layout

```
code/
  analysis/  Python analysis engines: reversibility scoring; mark classification; the
             determinant regression + robustness / sensitivity / equivalence (TOST) / mode-sweep /
             per-panel / complete-case variants; causal-status contrast; predictive null;
             variance partition; tissue organization; restoration and durability
             (Yousri, GSE199063); intervention- and drug-class analyses and reconciliation;
             cross-layer concordance; split-half reproducibility; ancestry replication;
             context / consensus / cross-context signatures; the per-mark atlas and panel manifest.
  figures/   R scripts that render Figures 1-6, Extended Data ED1-ED6 and Supplementary
             Figures S7-S9, the Python prep_*.py source-data builders (incl. render_F3.py),
             and the locked theme (theme_publication_ckm.R; theme_ckm.py).
  tables/    build_tables.py (Table 1 + Table 2 + Supplementary Tables workbook),
             build_paperB_docx.py (manuscript docx), assemble_manuscript_GM.py (Genome
             Medicine assembly), build_gm_package.py (submission-package assembler).
source_data/          frozen, machine-readable inputs for every figure panel (+ manifest).
supplementary_tables/ results tables backing Supplementary Tables S1-S14 + the per-mark
                      reproducibility atlas (per_mark_atlas.parquet).
```

## Reproducing the display items

Figures render in R from `source_data/` (frozen per-panel CSVs); the `code/figures/prep_*.py`
builders regenerate those CSVs from the tables in `supplementary_tables/` (Figure 3's exporter is
`render_F3.py`, and `render_F3.R` reads the determinant/predictive-null tables directly).
Tables build in Python from `supplementary_tables/`. **No upstream genetic or individual-level
data is needed to reproduce the display items.**

| Display item | Render script | Source-data builder | Key source files |
|---|---|---|---|
| Figure 1 (score is reproducible & power-independent) | `code/figures/render_F1.R` | `prep_F1ef_data.py` | `source_data/F1e_split_half.csv`, `F1f_ancestry.csv` (panels a-d programmatic / from results) |
| Figure 2 (reversibility organized by shared context) | `code/figures/render_F2.R` | `prep_F2_data.py`, `prep_F2f_data.py` | `F2a_pair_correlations.csv`, `F2b_cohesion.csv`, `F2c_pole_asymmetry.csv`, `F2d_universal_core.csv`, `F2f_tissue_programs.csv` |
| Figure 3 (no intrinsic feature predicts reversibility) | `code/figures/render_F3.R` | `render_F3.py` | `F3a_determinant_transcriptome.csv`, `F3c_determinant_proteome.csv`, `F3d_predictive_null.csv`, `F3e_ablation.csv`, `F3f_tau.csv` (renderer also reads `supplementary_tables/{determinant_meta,predictive_null,predictive_null_ablation,causal_status_effect}.tsv`) |
| Figure 4 (restoration toward lean + durability) | `code/figures/render_F4.R` | `prep_F4_data.py` | `F4a_*.csv`, `F4b_pct_toward_lean.csv`, `F4c_durability_*.csv`, `F4e_yousri_points.csv` (final panel **d**), `F4d_methyl_age_drift.csv` (final panel **e**) |
| Figure 5 (systemic response tracks weight-loss type, not drug label) | `code/figures/render_F5.R` | `prep_F5_data.py` | `F5a_intervention_systemic.csv`, `F5b_platform_control.csv`, `F5b_semaglutide_corr.csv`, `F5c_layered_core.csv`, `F5d_dissociation_*.csv`, `F5e_intervention_matrix.csv` |
| Figure 6 (context-not-molecule synthesis) | `code/figures/render_F6.R` | `prep_F6_data.py` | `F6a_context_corr.csv`, `F6b_concordance.csv`, `F6c_exemplars.csv`, `F6d_reversible_overlap.csv`, `F6e_persistent_overlap.csv` |
| Extended Data ED1-ED6 | `render_ED1.R` … `render_ED6.R` | `prep_ED1_data.py`, `prep_ED2_ED6_data.py`, `prep_ED3_data.py`, `prep_ED4_data.py`, `prep_ED5_data.py` | `source_data/ED1*.csv` … `ED6*.csv` |
| Suppl. Figs S7-S9 | `render_S7.R`, `render_S8.R`, `render_S9.R` | (self-contained) | `S7_drug_class_reconciliation.csv`; `S8_determinant_tost_forest.csv`; `S9a_variance_partition.csv`, `S9b_cell_means.csv` |
| Table 1 (determinant), Table 2 (cohort inventory), Supp. Tables S1-S14 | `code/tables/build_tables.py` | — | `supplementary_tables/*.tsv` (e.g. `determinant_meta.tsv`, `panel_manifest_full.tsv`, `causal_status_*.tsv`, `determinant_{tost,modes,per_panel,completecase}.tsv`, `variance_partition.tsv`, `downsampling_sensitivity.tsv`, `drug_class_*.tsv`, `predictive_null*.tsv`, `tissue_pair_correlations.tsv`, `context_signatures.tsv`, `intervention_systemic_consensus.tsv`, `universal_reversible_core.tsv`, `causal_nominatable_universe.tsv`, `intrinsic_missingness.tsv`) |
| Per-mark reproducibility atlas (Additional file 12) | — | `code/analysis/build_per_mark_atlas.py` | `supplementary_tables/per_mark_atlas.parquet` (903,025 mark x layer records; the 207 MB `.tsv` is not bundled and regenerates from this script) |

`source_data/source_data_manifest.tsv` maps every panel to its file and both scripts.

## Setup

**R** (>= 4.4): `ggplot2`, `cowplot`, `scales`, `patchwork`, `dplyr` (see `R_packages.txt`). The figure
scripts `source()` `code/figures/theme_publication_ckm.R` (the lab-locked publication theme, bundled here).

**Python** (>= 3.10): see `requirements.txt` (`pandas`, `numpy`, `scipy`, `statsmodels`, `scikit-learn`,
`openpyxl`, `python-docx`; `pyarrow` to read the atlas parquet).

> **Paths.** The scripts were written for the original environment and set an absolute `ROOT`/`RES`
> path near the top; the render scripts additionally `source()` the theme via an absolute path. After
> cloning, edit those paths to your local clone (the theme is bundled at
> `code/figures/theme_publication_ckm.R`), or run from the repo root.

## Data provenance

`source_data/` and `supplementary_tables/` contain **processed, summary-level** results only; **no
individual-level data are included.** All inputs are summary-level. The atlas was assembled from public
repositories and published per-feature effect estimates (accessions as stated in the manuscript's
*Availability of data and materials*):

- **Transcriptome (GEO):** GSE84046, GSE141221, GSE77962, GSE43471, GSE199063, GSE70529, GSE107894,
  GSE224310 (subcutaneous adipose); GSE83452, GSE106737, GSE48452 (liver); GSE83352, GSE157585,
  GSE161643, GSE60590 (skeletal muscle); GSE273902, GSE28358, GSE193771, GSE310742 (whole blood/PBMC).
  GSE117070 was excluded from signed analyses (batch-inverted directionality) and GSE114763 was retained
  only as a labelled epigenetic-memory (detraining) comparator.
- **DNA methylation (GEO):** GSE171140, GSE193730, GSE60655 (skeletal-muscle methylation series,
  distinct from the GSE60590 expression series).
- **Methylome (ArrayExpress):** CENTRAL cohort under E-MTAB-8956.
- **Metabolome (MetaboLights):** bariatric metabolome panel MTBLS218.
- **Proteome (published summary statistics, not re-derived from individual-level data):** DiRECT (SomaScan)
  and By-Band-Sleeve (Olink) from Goudswaard et al.; semaglutide (STEP, SomaScan) from Maretty et al.;
  empagliflozin (EMPEROR, Olink) from Zannad et al.; long-term RYGB (SomaLogic, 2- and 12-year) from
  Yousri et al.; metformin (Olink) from Connolly et al. The mass-spectrometry bariatric proteome panel
  was re-derived from the publicly deposited MaxQuant output of Wewer Albrechtsen et al. at
  **ProteomeXchange PXD009348 (MassIVE MSV000084344)**.
- **Genetic nomination inputs (upstream, not bundled):** cis-QTL atlases (UKB-PPP / deCODE pQTL;
  GTEx / eQTLGen / eQTL Catalogue eQTL; GoDMC mQTL; Chen et al. CLSA metaboQTL) and CKM outcome GWAS.

The **upstream cis-Mendelian-randomization / colocalization pipeline** (raw GWAS/QTL → causal
nomination) is documented in the `code/analysis/` scripts but is **not bundled**, because it depends on
these large external QTL/GWAS datasets; the frozen `supplementary_tables/` and `source_data/` let a
reviewer reproduce every reported number without it. Some `code/analysis/` scripts reference external
GEO/summary inputs by path and are provided to document the analysis logic; they require the corresponding
public dataset to re-run from scratch, but no such data are copied into this repository.

## License & citation

Code released under the MIT License (see `LICENSE`). If you use this code or data, please cite the
manuscript (Genome Medicine, 2026; citation to be added on publication) and this archived software
release (Zenodo DOI to be added on release — see `.zenodo.json`, `CITATION.cff`, and
`../SETUP_GITHUB_ZENODO.md`).

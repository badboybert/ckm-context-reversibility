#!/usr/bin/env python
"""Assemble the Genome Medicine submission package into 'submission package v3/' (snapshot COPIES;
canonical files in manuscript/ untouched). GM layout: manuscript + cover/declarations +
main figures + main tables (Table 1 + Table 2) + Additional files (1-6 ex-ED figs S1-S6, 7 supp tables,
8 STROBE-MR, 9-11 new supp figs S7-S9, 12 per-mark atlas) + source data (per-panel CSVs + the
causal-status effect-size table backing Supplementary Table 14) + author-todo.
Re-run to regenerate."""
import os, shutil, glob
SIG = os.path.dirname(os.path.abspath(__file__))
MAN = os.path.join(SIG, "manuscript"); FIG = os.path.join(MAN, "figures"); SD = os.path.join(FIG, "source_data")
RES = os.path.join(SIG, "results")
PKG = os.path.join(SIG, "submission package v3")
PKG_PREV = os.path.join(SIG, "submission package v2")  # carry _author_todo forward (v2 == v1, verbatim)
if os.path.isdir(PKG): shutil.rmtree(PKG)
def cp(src, reldst):
    dst = os.path.join(PKG, reldst); os.makedirs(os.path.dirname(dst), exist_ok=True)
    if os.path.exists(src): shutil.copy2(src, dst); return True
    print("  MISSING:", src); return False
def wr(rel, text):
    p = os.path.join(PKG, rel); os.makedirs(os.path.dirname(p), exist_ok=True)
    open(p, "w", encoding="utf-8").write(text)

README_V3 = """# Submission package v3 - Genome Medicine (post-review revision)

**Manuscript:** *A multi-omic atlas defines context, not molecular identity, as the determinant of weight-loss reversibility*
**Target journal:** Genome Medicine (BMC / Springer Nature), Research article.
**Companion:** Paper A - *in preparation, to be submitted to Communications Medicine* (cited as unpublished; provided on request).

> Snapshot copies. Canonical working files stay in `signature-pivot/manuscript/` and `signature-pivot/results/` - nothing moved. Regenerate via `signature-pivot/build_gm_package.py`. Author-metadata placeholders remain intended for the review round (see `_author_todo/`).

**v3 changelog:** v3: post-review — genetic causal-status result reframed from the standardization-compressed −0.004 SD to the interpretable raw effect (−0.18 SD, 95% CI −0.37..+0.01, n=9, underpowered); powered practical-equivalence null retained for the continuous features; +Supplementary Table 14; Figure 2/3b re-rendered; ED Table 1→Table 2; abstract ≤350w.

**v2 changelog:** v2: reviewer-driven revision - determinant-null reframed as practical-equivalence (SESOI/TOST), tissue+intervention context reframe, Fig 5 drug-class reconciliation, Table 1 heterogeneity, causal-status Supplementary Table 8, 6 new supp tables + 3 new supp figures + per-mark atlas.

**Display items:** 6 main figures + 2 tables (Table 1 + Table 2) + 12 additional files + 14 supplementary-table sheets.

## Contents
| Folder | Files |
|---|---|
| `01_manuscript/` | `MANUSCRIPT.docx` (integrated Genome Medicine format: Abstract -> Keywords -> Background -> **Methods** -> Results -> Discussion -> Conclusions -> List of abbreviations -> Declarations -> References -> Figures/Tables/Additional-file legends; 6 main figures embedded), `MANUSCRIPT.md` |
| `02_cover_declarations/` | `COVER_LETTER_AND_TITLE_PAGE` (Genome-Medicine-addressed), `DECLARATIONS` (8-item GM order) |
| `03_main_figures/` | `Figure_1`-`Figure_6` (PDF + PNG) |
| `04_main_tables/` | `Paper_B_Tables.xlsx` - **Table 1** (determinant meta-analysis) + **Table 2** (cohort inventory, promoted from the former Extended Data Table 1) + **14 supplementary-table sheets** (S1-S14, including S8 causal-status, S9 equivalence/TOST, and S14 causal-status effect sizes) |
| `05_additional_files/` | **Additional files 1-6** = Figures S1-S6 (former Extended Data figures; PDF + PNG); **Additional file 7** = Supplementary Tables (in the workbook under `04_main_tables/`); **Additional file 8** = STROBE-MR reporting checklist; **Additional files 9-11** = new Figures S7-S9 (PDF + PNG); **Additional file 12** = per-mark reversibility atlas (`per_mark_atlas.parquet`; the full 207 MB TSV lives in the code repository) |
| `06_source_data/` | per-panel machine-readable CSVs (one per plotted panel) + `causal_status_effect.tsv` (source for Supplementary Table 14) + manifest + provenance |
| `_author_todo/` | `AUTHOR_ACTION_FLAGS`, `PRESUBMISSION_TODO`, `CITATIONS_RESOLVED` |

## Genome Medicine format compliance
- **Structured abstract** (Background/Methods/Results/Conclusions) + keywords.
- GM section order with **Methods in the middle**; standalone **Conclusions**; **List of abbreviations**; 8-item **Declarations**.
- **No Extended Data tier** - former ED figures are Additional files 1-6; ED Table 1 promoted to main **Table 2**; new revision figures are Additional files 9-11 (Figures S7-S9).
- **References: numbered Vancouver by order of appearance**, contiguous, all cited.
- STROBE-MR checklist supplied (Additional file 8). No GM word cap on main text.

## Status
Author-metadata fields (co-author list / ORCIDs / CRediT roles; the code Zenodo DOI + GitHub URL; per-cohort data-use-agreement contacts) remain author-pending - see `_author_todo/`. Corresponding author Bertrand Chin-Ming Tan (ORCID 0000-0002-2218-7115); funding grants filled (NSTC 112-2320-B-182-011-MY3; 114-2320-B-182-013-MY3; CGMH CMRPD1P0221; BMRP960); data accessions incl. PXD009348.
"""

for d in ["01_manuscript","02_cover_declarations","03_main_figures","04_main_tables",
          "05_additional_files","06_source_data","_author_todo"]:
    os.makedirs(os.path.join(PKG, d), exist_ok=True)

# 01 manuscript
cp(os.path.join(MAN,"MANUSCRIPT_GM.docx"), "01_manuscript/MANUSCRIPT.docx")
cp(os.path.join(MAN,"MANUSCRIPT_GM.md"),   "01_manuscript/MANUSCRIPT.md")
# 02 cover + declarations
for f,dst in [("COVER_LETTER_GM","COVER_LETTER_AND_TITLE_PAGE"),("DECLARATIONS_GM","DECLARATIONS")]:
    cp(os.path.join(MAN,f+".docx"), f"02_cover_declarations/{dst}.docx")
    cp(os.path.join(MAN,f+".md"),   f"02_cover_declarations/{dst}.md")
# 03 main figures
for n in range(1,7):
    for e in ("png","pdf"): cp(os.path.join(FIG,f"Figure_{n}.{e}"), f"03_main_figures/Figure_{n}.{e}")
# 04 main tables (Table 1 determinant + Table 2 cohort [promoted] + supp = the workbook)
cp(os.path.join(MAN,"tables","Paper_B_Tables.xlsx"), "04_main_tables/Paper_B_Tables.xlsx")
# 05 additional files: 1-6 = former ED figures (S1-S6); 7 = supp tables (in workbook, noted); 8 = STROBE-MR
for n in range(1,7):
    for e in ("png","pdf"): cp(os.path.join(FIG,f"Figure_ED{n}.{e}"), f"05_additional_files/Additional_file_{n}_Figure_S{n}.{e}")
cp(os.path.join(MAN,"STROBE_MR_CHECKLIST_GM.docx"), "05_additional_files/Additional_file_8_STROBE-MR_checklist.docx")
cp(os.path.join(MAN,"STROBE_MR_CHECKLIST_GM.md"),   "05_additional_files/Additional_file_8_STROBE-MR_checklist.md")
# 9-11 = new v2 supplementary figures S7-S9; 12 = per-mark reversibility atlas (Supplementary Data)
for supp, af in [("S7",9),("S8",10),("S9",11)]:
    for e in ("png","pdf"): cp(os.path.join(FIG,f"Figure_{supp}.{e}"), f"05_additional_files/Additional_file_{af}_Figure_{supp}.{e}")
# atlas: ship the compact parquet only; the 207 MB .tsv stays in the code repo (noted in PROVENANCE)
cp(os.path.join(RES,"per_mark_atlas.parquet"), "05_additional_files/Additional_file_12_per_mark_atlas.parquet")
# 06 source data
sd = sorted(glob.glob(os.path.join(SD,"*.csv")))
for f in sd: cp(f, "06_source_data/"+os.path.basename(f))
# causal-status effect-size table = source for Supplementary Table 14
cp(os.path.join(RES,"causal_status_effect.tsv"), "06_source_data/causal_status_effect.tsv")
wr("06_source_data/source_data_manifest.tsv",
   "file\tbacks_panel\n"
   + "\n".join(f"{os.path.basename(f)}\t{os.path.basename(f).split('_')[0]}" for f in sd)
   + "\ncausal_status_effect.tsv\tSupplementary Table 14 (genetic causal-status effect sizes)"
   + "\nper_mark_atlas.parquet\tAdditional file 12 (per-mark reversibility atlas; full 207 MB TSV in the code repository, not shipped)\n")
wr("06_source_data/PROVENANCE.md",
   "# Source data provenance\n\n"
   "Per-panel machine-readable data for every main figure (F1-6) and every Additional-file figure "
   "(Additional files 1-11 = Figures S1-S9), one CSV per plotted panel. "
   "`causal_status_effect.tsv` is the source for Supplementary Table 14 (genetic causal-status effect sizes: "
   "raw group SMD, multivariable-adjusted, nominatable-universe, and standardized-anchor estimates on the nine causal transcripts). "
   "Additional file 12 is the per-mark reversibility atlas (../05_additional_files/Additional_file_12_per_mark_atlas.parquet); "
   "the full per-mark table is also available as a 207 MB TSV in the code repository (not shipped here to keep the package lightweight). "
   "GSE199063 AceView novel transcripts excluded from gene-level analyses; adipose 2->5yr durability on canonical protein-coding genes. "
   "Table values are in ../04_main_tables/Paper_B_Tables.xlsx. Accessions in ../02_cover_declarations/DECLARATIONS.md.\n")

# carry over author-todo (verbatim) from the previous package + README (updated for v3)
for f in sorted(glob.glob(os.path.join(PKG_PREV, "_author_todo", "*"))):
    cp(f, "_author_todo/"+os.path.basename(f))
wr("README_SUBMISSION.md", README_V3)

print("GM package v3 assembled:", PKG)
print("source_data CSVs:", len(sd))

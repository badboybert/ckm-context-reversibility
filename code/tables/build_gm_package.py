#!/usr/bin/env python
"""Assemble the Genome Medicine submission package into 'submission package v4/' (snapshot COPIES;
canonical files in manuscript/ untouched). GM layout: manuscript + cover/declarations +
main figures + main tables (Table 1 + Table 2) + Additional files (1-6 ex-ED figs S1-S6, 7 supp tables,
8 STROBE-MR, 9-11 new supp figs S7-S9; the per-mark atlas is deposited, not shipped) + source data (per-panel CSVs + the
causal-status effect-size table backing Supplementary Table 14) + author-todo.
Re-run to regenerate."""
import os, re, shutil, glob, subprocess, sys
SIG = os.path.dirname(os.path.abspath(__file__))
MAN = os.path.join(SIG, "manuscript"); FIG = os.path.join(MAN, "figures"); SD = os.path.join(FIG, "source_data")
RES = os.path.join(SIG, "results")
PKG = os.path.join(SIG, "submission package v11")
# (no PKG_PREV: _author_todo is generated from its canonical source in manuscript/_author_todo — see below)
if os.path.isdir(PKG): shutil.rmtree(PKG)
def cp(src, reldst):
    dst = os.path.join(PKG, reldst); os.makedirs(os.path.dirname(dst), exist_ok=True)
    if os.path.exists(src): shutil.copy2(src, dst); return True
    print("  MISSING:", src); return False
def wr(rel, text):
    p = os.path.join(PKG, rel); os.makedirs(os.path.dirname(p), exist_ok=True)
    open(p, "w", encoding="utf-8").write(text)

README_V5 = """# Submission package v11 - Genome Medicine

**Manuscript:** *A multi-omic atlas reveals context-dependent molecular reversibility after weight loss and cardiometabolic intervention*
**Target journal:** Genome Medicine (BMC / Springer Nature), Research article.
**Companion:** the companion study is a public medRxiv preprint (doi:10.64898/2026.07.27.26358659), cited as reference [17].

> Snapshot copies. Canonical working files stay in `signature-pivot/manuscript/` and `signature-pivot/results/` - nothing moved. Regenerate via `signature-pivot/build_gm_package.py`.
> Author-metadata placeholders are listed in `_author_todo/`.

**Display items:** 6 main figures + 2 main tables (Word tables in the manuscript) + 11 additional files + 16 supplementary-table sheets.

## Contents
| Folder | Files |
|---|---|
| `01_manuscript/` | `MANUSCRIPT.docx` (integrated Genome Medicine format: Abstract -> Keywords -> Background -> **Methods** -> Results -> Discussion -> Conclusions -> List of abbreviations -> Declarations -> References -> Figures/Tables/Additional-file legends; 6 main figures embedded), `MANUSCRIPT.md` |
| `02_cover_declarations/` | `COVER_LETTER_ONLY.docx` — **upload this one into the portal's cover-letter field**; `COVER_LETTER_AND_TITLE_PAGE` (same letter plus a title page, for filling the portal's author/affiliation fields); `DECLARATIONS` (8-item GM order) |
| `03_main_figures/` | `Figure_1`-`Figure_6` (PDF + PNG) |
| `05_additional_files/` | **Additional files 1-6** = Figures S1-S6 (former Extended Data figures; PDF + PNG); **Additional file 7** = the supplementary-table workbook (`Table 2 - full registry` + S1-S16; main Tables 1-2 are Word tables in the manuscript); **Additional file 8** = STROBE-MR reporting checklist; **Additional files 9-11** = new Figures S7-S9 (PDF + PNG). The per-mark atlas is NOT an additional file: at 35.81 MiB it is over the 20 MB limit and is deposited with the code under the Zenodo DOI |
| `06_source_data/` | per-panel machine-readable CSVs (one per plotted panel) + `causal_status_effect.tsv` (source for Supplementary Table 14) + manifest + provenance |
| `_author_todo/` | `AUTHOR_ACTION_FLAGS`, `PRESUBMISSION_TODO`, `CITATIONS_RESOLVED` |

## Genome Medicine format compliance
- **Structured abstract** (Background/Methods/Results/Conclusions) + keywords.
- GM section order with **Methods in the middle**; standalone **Conclusions**; **List of abbreviations**; 8-item **Declarations**.
- **No Extended Data tier** - former ED figures are Additional files 1-6; ED Table 1 promoted to main **Table 2**; new revision figures are Additional files 9-11 (Figures S7-S9).
- **References: numbered Vancouver by order of appearance**, contiguous, all cited.
- STROBE-MR checklist supplied (Additional file 8). No GM word cap on main text.

## Status
Corresponding author Bertrand Chin-Ming Tan (ORCID 0000-0002-2218-7115); funding grants filled (NSTC 112-2320-B-182-011-MY3; 114-2320-B-182-013-MY3; CGMH CMRPD1P0221; BMRP960); data accessions incl. PXD009348.
Author-metadata fields (co-author list / ORCIDs / CRediT roles; the code Zenodo DOI) remain author-pending - see `_author_todo/`.
"""

# ---------------------------------------------------------------------------------------------------
# STALENESS GATE. The package is a snapshot of generated files, so it can ship outputs that are older
# than the code that makes them: this happened — build_tables.py was corrected and the package, built
# earlier, kept copying the PREVIOUS workbook, so a fixed verdict rule shipped unfixed. Nothing in the
# copy step can notice that. Fail loudly instead of snapshotting stale bytes.
def _stale(out, srcs):
    if not os.path.exists(out): return f"MISSING OUTPUT {out}"
    o = os.path.getmtime(out)
    old = [s for s in srcs if os.path.exists(s) and os.path.getmtime(s) > o]
    return f"{os.path.basename(out)} is OLDER than {[os.path.basename(x) for x in old]}" if old else None
_checks = [
    (os.path.join(MAN, "tables", "Paper_B_Tables.xlsx"), [os.path.join(MAN, "tables", "build_tables.py")]),
    (os.path.join(MAN, "MANUSCRIPT_GM.docx"), [os.path.join(MAN, "MANUSCRIPT_GM.md")]),
    (os.path.join(MAN, "DECLARATIONS_GM.docx"), [os.path.join(MAN, "DECLARATIONS_GM.md")]),
    (os.path.join(MAN, "COVER_LETTER_GM.docx"), [os.path.join(MAN, "COVER_LETTER_GM.md")]),
    (os.path.join(MAN, "COVER_LETTER_ONLY.md"), [os.path.join(MAN, "COVER_LETTER_GM.md")]),
    (os.path.join(MAN, "COVER_LETTER_ONLY.docx"), [os.path.join(MAN, "COVER_LETTER_ONLY.md")]),
    (os.path.join(MAN, "STROBE_MR_CHECKLIST_GM.docx"), [os.path.join(MAN, "STROBE_MR_CHECKLIST_GM.md")]),
] + [(os.path.join(FIG, f"Figure_{n}.png"), [os.path.join(FIG, f"render_F{n}.R")]) for n in range(1, 7)] \
  + [(os.path.join(FIG, f"Figure_ED{n}.png"), [os.path.join(FIG, f"render_ED{n}.R")]) for n in range(1, 7)] \
  + [(os.path.join(FIG, f"Figure_S{n}.png"), [os.path.join(FIG, f"render_S{n}.R")]) for n in (7, 8, 9)]
_stale_msgs = [m for m in (_stale(o, s) for o, s in _checks) if m]
if _stale_msgs:
    raise SystemExit("REFUSING TO BUILD — regenerate these first, the package would ship stale bytes:\n  "
                     + "\n  ".join(_stale_msgs))
# ---------------------------------------------------------------------------------------------------
# CONTENT GATE. The staleness gate catches OUT-OF-DATE bytes; this catches WRONG bytes. Round 7 item 3:
# the Methods sentence shipped "19 analysed panels" over an enumeration of 13, contiguous and internally
# consistent but untrue. Refuse to package a Methods transcriptome enumeration that has drifted from the
# ANALYZED transcriptome set in results/panel_manifest_full.tsv (the check self-tests with --self-test).
_pc = subprocess.run([sys.executable, os.path.join(SIG, "qa", "check_methods_panels.py")],
                     capture_output=True, text=True)
if _pc.returncode != 0:
    raise SystemExit("REFUSING TO BUILD — Methods panel enumeration drifted from the manifest:\n"
                     + _pc.stdout + _pc.stderr)
# ---------------------------------------------------------------------------------------------------

for d in ["01_manuscript","02_cover_declarations","03_main_figures",
          "05_additional_files","06_source_data","_author_todo"]:
    os.makedirs(os.path.join(PKG, d), exist_ok=True)

# 01 manuscript
cp(os.path.join(MAN,"MANUSCRIPT_GM.docx"), "01_manuscript/MANUSCRIPT.docx")
cp(os.path.join(MAN,"MANUSCRIPT_GM.md"),   "01_manuscript/MANUSCRIPT.md")
# 02 cover + declarations
for f,dst in [("COVER_LETTER_GM","COVER_LETTER_AND_TITLE_PAGE"),("DECLARATIONS_GM","DECLARATIONS")]:
    cp(os.path.join(MAN,f+".docx"), f"02_cover_declarations/{dst}.docx")
    cp(os.path.join(MAN,f+".md"),   f"02_cover_declarations/{dst}.md")
# Letter-only version for the portal's cover-letter field: the combined file appends a title page the
# manuscript already carries in full, so uploading it makes the editor read a duplicate before the
# letter. Derived from the same .md by manuscript/split_cover_letter.py, so the two cannot diverge.
cp(os.path.join(MAN,"COVER_LETTER_ONLY.docx"), "02_cover_declarations/COVER_LETTER_ONLY.docx")
# 03 main figures
for n in range(1,7):
    for e in ("png","pdf"): cp(os.path.join(FIG,f"Figure_{n}.{e}"), f"03_main_figures/Figure_{n}.{e}")
# 04/05 tables. The workbook ships ONCE, as Additional file 7.
#
# It used to be copied to BOTH 04_main_tables/ and 05_additional_files/, and round 6 caught the two
# files being byte-identical — two names for one object is two sources of truth for an editor to
# reconcile. Main Tables 1 and 2 are now Word table objects inside MANUSCRIPT.docx (generated from
# this same workbook by manuscript/sync_main_tables.py), so 04_main_tables/ has nothing left to hold
# and the workbook is purely supplementary.
cp(os.path.join(MAN,"tables","Paper_B_Tables.xlsx"), "05_additional_files/Additional_file_7_Supplementary_Tables.xlsx")
# 05 additional files: 1-6 = former ED figures (S1-S6); 7 = supp tables (in workbook, noted); 8 = STROBE-MR
for n in range(1,7):
    for e in ("png","pdf"): cp(os.path.join(FIG,f"Figure_ED{n}.{e}"), f"05_additional_files/Additional_file_{n}_Figure_S{n}.{e}")
cp(os.path.join(MAN,"STROBE_MR_CHECKLIST_GM.docx"), "05_additional_files/Additional_file_8_STROBE-MR_checklist.docx")
cp(os.path.join(MAN,"STROBE_MR_CHECKLIST_GM.md"),   "05_additional_files/Additional_file_8_STROBE-MR_checklist.md")
# 9-11 = new v2 supplementary figures S7-S9; 12 = per-mark reversibility atlas (Supplementary Data)
for supp, af in [("S7",9),("S8",10),("S9",11)]:
    for e in ("png","pdf"): cp(os.path.join(FIG,f"Figure_{supp}.{e}"), f"05_additional_files/Additional_file_{af}_Figure_{supp}.{e}")
# The per-mark atlas is NOT shipped. At 35.81 MiB it exceeds the 20 MB per-additional-file limit;
# lossless recompression reaches only 32.3 MiB and a per-layer split still leaves the methylation
# layer at 29.6 MiB, so it is deposited with the analysis code under the Zenodo DOI and cited in
# "Availability of data and materials" instead. There is no Additional file 12.
# 06 source data
sd = sorted(glob.glob(os.path.join(SD,"*.csv")))
for f in sd: cp(f, "06_source_data/"+os.path.basename(f))
# causal-status effect-size table = source for Supplementary Table 14
cp(os.path.join(RES,"causal_status_effect.tsv"), "06_source_data/causal_status_effect.tsv")
# true-weight-loss-only sensitivity = source for Supplementary Table 15
for f in ("determinant_meta_wlonly.tsv","context_wlonly.tsv"):
    cp(os.path.join(RES,"sensitivity_wlonly",f), "06_source_data/"+f)
# restoration uncertainty (binomial CI / permutation calibration) = source for Supplementary Table 16
# and for the open-circle "not distinguishable from chance" flags in Figure 4b
for f in ("restoration_uncertainty.tsv","between_context_spread.tsv","restoration_uncertainty_summary.json"):
    cp(os.path.join(RES,"restoration_uncertainty",f), "06_source_data/"+f)
# --- source-data manifest, naming the GENERATING SCRIPT per file ------------------------------------
# The Availability statement promises "a manifest naming the script that generates each" panel file. The
# render_*.{R,py} scripts READ these CSVs (rd()/read.csv) — they do NOT write them; the WRITERS are the
# prep_*.py panels plus a few render scripts that also emit their own data, and the src/ engines for the
# supplementary-table TSVs. We scan those writers' write idioms so the mapping cannot silently drift, then
# ASSERT every plotted-panel file has a named generator and refuse to build otherwise.
_PY_W = re.compile(r"""(?:to_csv|to_json)\(\s*os\.path\.join\([^,]+,\s*f?['"]([\w./{}-]+\.(?:csv|tsv|json))['"]""")
_PY_J = re.compile(r"""open\(\s*os\.path\.join\([^,]+,\s*f?['"]([\w./{}-]+\.json)['"]""")
_R_W = re.compile(r"""(?:write\.csv|fwrite|write_tsv)\([\s\S]*?file\.path\(SD,\s*f?['"]([\w./{}-]+\.(?:csv|tsv))['"]""")
def _gen_matchers():
    out = []
    writers = (glob.glob(os.path.join(FIG, "prep_*.py"))
               + [os.path.join(FIG, s) for s in ("render_F2.R", "render_F3.py", "render_S7.R", "render_S8.R")]
               + [os.path.join(SIG, "src", s) for s in ("restoration_uncertainty.py",
                  "determinant_sensitivity_wlonly.py", "context_sensitivity_wlonly.py")])
    for w in writers:
        if not os.path.exists(w):
            continue
        txt = open(w, encoding="utf-8").read(); rel = os.path.relpath(w, SIG).replace("\\", "/")
        for tmpl in set(_PY_W.findall(txt)) | set(_PY_J.findall(txt)) | set(_R_W.findall(txt)):
            rx = "^" + re.sub(r"\\\{[^}]*\\\}", "[^/]+", re.escape(os.path.basename(tmpl))) + "$"
            out.append((re.compile(rx), rel))
    out.append((re.compile(r"^F3.+\.csv$"), "manuscript/figures/render_F3.py"))  # F3 written via nm+'.csv'
    return out
_GEN_M = _gen_matchers()
_GEN_EXPLICIT = {  # verified writers the os.path.join idiom scan does not reach, or files with no in-tree generator
    "causal_status_effect.tsv": "src/causal_status_effect.py",           # writes via a S('results/...') path helper
    "per_mark_atlas.parquet": "src/build_per_mark_atlas.py",
    "between_context_spread.tsv": "(precomputed restoration intermediate; no generating script in this release)",
}
def _gen_for(fn):
    hits = sorted({rel for rx, rel in _GEN_M if rx.match(fn)})
    if len(hits) > 1:
        raise SystemExit(f"REFUSING TO BUILD — source-data file {fn!r} has >1 candidate generator: {hits}")
    return hits[0] if hits else _GEN_EXPLICIT.get(fn)
_BACKS = {
    "causal_status_effect.tsv": "Supplementary Table 14 (genetic causal-status effect sizes)",
    "restoration_uncertainty.tsv": "Supplementary Table 16 + the Figure 4b not-distinguishable-from-chance flags",
    "restoration_uncertainty_summary.json": "Supplementary Table 16 (restoration-uncertainty summary)",
    "between_context_spread.tsv": "restoration between-context spread (supporting Supplementary Table 16 / Fig. 4b)",
    "determinant_meta_wlonly.tsv": "Supplementary Table 15 (weight-loss-only determinant sensitivity)",
    "context_wlonly.tsv": "Supplementary Table 15 (weight-loss-only tissue/variance-partition sensitivity)",
}
_shipped_sd = sorted(os.path.basename(f) for f in glob.glob(os.path.join(PKG, "06_source_data", "*"))
                     if os.path.basename(f) not in ("source_data_manifest.tsv", "PROVENANCE.md"))
_panel = re.compile(r"(F\d|ED\d|S\d)")
_rows, _miss = [], []
for _fn in _shipped_sd:
    _g = _gen_for(_fn)
    if _g is None and _panel.match(_fn):
        _miss.append(_fn)
    _rows.append((_fn, _BACKS.get(_fn, _fn.split("_")[0]), _g or "n/a"))
if _miss:
    raise SystemExit("REFUSING TO BUILD — plotted-panel source-data file(s) without a named generating "
                     "script (the Availability statement promises one): " + ", ".join(_miss))
# per_mark_atlas is deposited (not shipped in 06_source_data) but named in the manifest for provenance
_rows.append(("per_mark_atlas.parquet", "deposited with the analysis code under the Zenodo DOI (over the "
              "20 MB additional-file limit); full 207 MB TSV also in the code repository",
              _GEN_EXPLICIT["per_mark_atlas.parquet"]))
wr("06_source_data/source_data_manifest.tsv",
   "file\tbacks_panel\tgenerating_script\n" + "\n".join(f"{f}\t{b}\t{g}" for f, b, g in _rows) + "\n")
wr("06_source_data/PROVENANCE.md",
   "# Source data provenance\n\n"
   "Per-panel machine-readable data for every main figure (F1-6) and every Additional-file figure "
   "(Additional files 1-11 = Figures S1-S9), one CSV per plotted panel. "
   "`causal_status_effect.tsv` is the source for Supplementary Table 14 (genetic causal-status effect sizes: "
   "raw group SMD, multivariable-adjusted, nominatable-universe, and standardized-anchor estimates on the nine causal transcripts). "
   "`determinant_meta_wlonly.tsv` and `context_wlonly.tsv` are the source for Supplementary Table 15 (true-weight-loss-only sensitivity: determinant meta and tissue/variance partition recomputed on diet/CR + bariatric panels only). "
   "The per-mark reversibility atlas (per_mark_atlas.parquet) is deposited with the analysis code under the Zenodo DOI, "
   "not as an Additional file, because it exceeds the 20 MB per-file limit; the full 207 MB TSV is in the code repository. "
   "GSE199063 AceView novel transcripts excluded from gene-level analyses; adipose 2->5yr durability on canonical protein-coding genes. "
   "Table values are in ../05_additional_files/Additional_file_7_Supplementary_Tables.xlsx. Accessions in ../02_cover_declarations/DECLARATIONS.md.\n")

# carry over author-todo (verbatim) from the previous package + README (updated for v3)
# _author_todo now has a CANONICAL source (manuscript/_author_todo) rather than being carried forward
# from the previous package: that indirection silently shipped an EMPTY directory the moment the previous
# package was archived, while README kept claiming three files. Assert, don't hope.
_todo = sorted(glob.glob(os.path.join(MAN, "_author_todo", "*")))
for f in _todo:
    cp(f, "_author_todo/"+os.path.basename(f))
_want = {"AUTHOR_ACTION_FLAGS.md", "PRESUBMISSION_TODO.md", "CITATIONS_RESOLVED.md"}
_got = set(os.listdir(os.path.join(PKG, "_author_todo")))
assert _got == _want, f"_author_todo mismatch: missing {_want-_got}, unexpected {_got-_want}"
wr("README_SUBMISSION.md", README_V5)

# ---------------------------------------------------------------------------
# CLEAN UPLOAD VARIANT. Excludes the internal _author_todo notes from the journal upload. The clean
# README is DERIVED from the same README_V5 string, never maintained as a second copy, and the
# derivation is asserted, so the two READMEs cannot drift apart.
# ---------------------------------------------------------------------------
UPLOAD = PKG + " (upload)"
if os.path.isdir(UPLOAD): shutil.rmtree(UPLOAD)
shutil.copytree(PKG, UPLOAD, ignore=shutil.ignore_patterns("_author_todo"))

_clean = "\n".join(l for l in README_V5.split("\n")
                   if "_author_todo" not in l and "AUTHOR_ACTION_FLAGS" not in l)
assert "_author_todo" not in _clean and "AUTHOR_ACTION_FLAGS" not in _clean, \
    "clean README still references the internal to-do files"
assert len(_clean.split("\n")) < len(README_V5.split("\n")), "clean README derivation removed nothing"

# The line-strip above once emptied a whole section: the entire '## Status' body was a single line
# mentioning _author_todo, so the upload README ended on a bare heading. Assert that no heading in
# EITHER README is left without a body, and that neither carries the retired post-review framing.
def _check_readme(text, which):
    lines = text.split("\n")
    for i, l in enumerate(lines):
        if not l.startswith("#"):
            continue
        body = [b for b in lines[i+1:] if b.strip()]
        assert body and not body[0].startswith("#"), f"{which} README: section {l!r} has an empty body"
    for bad in ("post-review", "in response to peer review", "response-to-reviewers"):
        assert bad not in text, f"{which} README carries retired framing {bad!r} (the paper has never been reviewed)"
_check_readme(README_V5, "full")
_check_readme(_clean, "upload")
open(os.path.join(UPLOAD, "README_SUBMISSION.md"), "w", encoding="utf-8").write(_clean)

assert not os.path.isdir(os.path.join(UPLOAD, "_author_todo")), "_author_todo leaked into the upload variant"
_pkg_n = sum(len(f) for _, _, f in os.walk(PKG))
_upl_n = sum(len(f) for _, _, f in os.walk(UPLOAD))
assert _upl_n == _pkg_n - 3, f"upload variant should drop exactly the 3 author-todo files: {_pkg_n} -> {_upl_n}"

print("GM package assembled:", PKG, f"({_pkg_n} files)")
print("clean upload variant:", UPLOAD, f"({_upl_n} files; _author_todo excluded)")
print("source_data CSVs:", len(sd))


#!/usr/bin/env python
"""Assemble the GENOME MEDICINE manuscript (MANUSCRIPT_GM.md) from the GM-format pieces + reordered body.
Deterministic transforms only: section reorder (Background -> Methods -> Results -> Discussion ->
Conclusions -> List of abbreviations -> Declarations -> References), Extended Data -> Additional files,
ED Table 1 -> main Table 2. In-text author-date/name/PMID citations are LEFT AS-IS here; the
author-date -> numbered-Vancouver conversion is a separate verified pass (convert_citations step).
Canonical Nature files are untouched; GM output = *_GM files.
"""
import os, re
HERE = os.path.dirname(os.path.abspath(__file__))
def R(fn):
    p = os.path.join(HERE, fn); return open(p, encoding="utf-8").read() if os.path.exists(p) else ""

TITLE = "A multi-omic atlas defines context, not molecular identity, as the determinant of weight-loss reversibility"

def strip_lead(body):
    """Drop a leading H1 title + leading italic editorial meta lines."""
    out, started = [], False
    for ln in body.splitlines():
        if not started:
            if ln.startswith("# "): continue
            if ln.strip().startswith("*") and ln.strip().endswith("*") and len(ln.strip()) > 2: continue
            if ln.strip() == "": continue
            started = True
        out.append(ln)
    return "\n".join(out).strip()

def relabel_ed(text):
    """Extended Data -> Additional files; ED Table 1 -> Table 2."""
    text = re.sub(r"Extended Data Table\s*1", "Table 2", text)
    text = re.sub(r"Extended Data Fig(?:ure)?\.?\s*(\d+)", r"Additional file \1", text)
    text = re.sub(r"Extended Data Fig(?:ure)?s?\.?\s*(\d+)\s*[–-]\s*(\d+)", r"Additional files \1-\2", text)
    return text

def demote(text):
    """Demote in-section ## / ### by one level so they sit under GM H2 section heads."""
    text = re.sub(r"^### ", "#### ", text, flags=re.M)
    text = re.sub(r"^## ", "### ", text, flags=re.M)
    return text

def section_body(md_file):
    return demote(relabel_ed(strip_lead(R(md_file))))

def piece_body(gm_file, drop_first_h2=None):
    """Return a GM piece (already ## structured), relabel ED, keep its ## as section content."""
    t = relabel_ed(strip_lead(R(gm_file)))
    return t

parts = []
parts.append(f"# {TITLE}\n")
parts.append("**Authors:** [AUTHOR — full co-author list, order, degrees and ORCIDs]; Bertrand Chin-Ming Tan (ORCID 0000-0002-2218-7115), corresponding author.\n")
parts.append("**Affiliations:** 1. Department of Biomedical Sciences, College of Medicine, Chang Gung University, Taoyuan City 33302, Taiwan. 2. Graduate Institute of Biomedical Sciences, College of Medicine, Chang Gung University, Taoyuan City 33302, Taiwan. 3. Center for Bioenergetics and Metabolic Translational Research, Chang Gung University, Taoyuan City 33302, Taiwan. 4. Department of Neurosurgery, Lin-Kou Medical Center, Chang Gung Memorial Hospital, Taoyuan 333, Taiwan.\n")
parts.append("**Correspondence:** Bertrand Chin-Ming Tan, btan@mail.cgu.edu.tw.\n")
parts.append("\n---\n")

# Abstract + Keywords (already GM-structured piece; keep its ## Abstract / ## Keywords as H2)
parts.append(R("ABSTRACT_GM.md").strip() + "\n")

# Body in GM order
parts.append("\n## Background\n")
parts.append(section_body("INTRODUCTION.md") + "\n")
parts.append("\n## Methods\n")
parts.append(section_body("METHODS.md") + "\n")
parts.append("\n## Results\n")
parts.append(section_body("RESULTS.md") + "\n")
parts.append("\n## Discussion\n")
parts.append(section_body("DISCUSSION.md") + "\n")

# Conclusions (GM piece already has "## Conclusions")
parts.append("\n" + R("CONCLUSIONS_GM.md").strip() + "\n")

# List of abbreviations (piece has "## List of abbreviations")
parts.append("\n" + R("ABBREVIATIONS_GM.md").strip() + "\n")

# Declarations (piece has "## Declarations" + 8 subsections)
parts.append("\n" + relabel_ed(R("DECLARATIONS_GM.md").strip()) + "\n")

# References (numbered list produced by convert step; fallback to placeholder)
parts.append("\n## References\n")
refs_gm = R("REFERENCES_GM.md").strip()
parts.append((refs_gm if refs_gm else "[REFERENCES — numbered-Vancouver conversion pending (convert_citations step)]") + "\n")

# Figure legends (main F1-6) + Table legends (Table 1 + promoted Table 2) + Additional-file legends
parts.append("\n## Figures\n")
leg = relabel_ed(strip_lead(R("LEGENDS.md")))
# split main (Figure 1-6) vs additional-file (former ED) legends
main_fig = []
addl = []
for block in re.split(r"\n(?=\*\*(?:Figure|Additional file) )", leg):
    b = block.strip()
    if not b: continue
    (addl if b.startswith("**Additional file") else main_fig).append(b)
parts.append("\n\n".join(main_fig) + "\n")
parts.append("\n## Tables\n")
parts.append(relabel_ed(strip_lead(R("TABLE_LEGENDS.md"))) + "\n")
parts.append("\n## Additional files\n")
parts.append("Additional files 1-6 (Figures S1-S6) — legends below. Additional file 7 — Supplementary Tables S1-S7 (workbook). Additional file 8 — STROBE-MR reporting checklist.\n")
parts.append("\n".join(addl) + "\n")

out = "\n".join(parts).strip() + "\n"
open(os.path.join(HERE, "MANUSCRIPT_GM.md"), "w", encoding="utf-8").write(out)

def wc(s): return len(re.findall(r"\S+", re.sub(r"[#*`|>_-]", " ", s)))
main = section_body("INTRODUCTION.md")+section_body("RESULTS.md")+section_body("DISCUSSION.md")+R("CONCLUSIONS_GM.md")
print("wrote MANUSCRIPT_GM.md")
print("Background+Results+Discussion+Conclusions ~words:", wc(main))
print("Additional-file (former ED) legend blocks:", len(addl), "| main figure legend blocks:", len(main_fig))
print("residual 'Extended Data' mentions:", out.count("Extended Data"))

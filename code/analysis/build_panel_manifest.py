#!/usr/bin/env python
"""Build data/signature/panel_manifest.tsv — the authoritative per-panel registry
for the reversibility-signature pivot. Reads each reversal table on disk, derives
n_marks / sig-count / has_se from the actual file, and stamps the frozen metadata
(intervention, tissue, BMI-axis source, baseline-SD source) + the n>=30 power-floor
gate (SCORE_SPEC_FROZEN.md sections 4 & 7). No scoring here; manifest only.
"""
import os, pandas as pd, numpy as np

SIG = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))   # signature-pivot/
PROJ = os.path.dirname(SIG)                                          # weight-loss.omics/ (shared root)
def P(*a): return os.path.join(PROJ, *a)   # SHARED inputs (referenced, never moved)
def S(*a): return os.path.join(SIG, *a)    # signature-pivot outputs

# Frozen per-panel metadata. n_pairs: None => read from a table column; else hardcoded from MAP.
# bmi_axis / base_sd sources encode the data-availability matrix (spec section 7).
PANELS = [
 # proteome (baseline individual data NOT open -> base_sd = population-variance proxy)
 dict(file='data/proteome/reversal_olink_bbs.tsv',       layer='proteome', platform='Olink',    intervention='bariatric', tissue='plasma', n_pairs=118,  bmi_axis='goudswaard_ST5_INTERVAL', base_sd='proxy_ukbppp'),
 dict(file='data/proteome/reversal_somascan_direct.tsv', layer='proteome', platform='SomaScan', intervention='diet',      tissue='plasma', n_pairs=292,  bmi_axis='goudswaard_ST5_INTERVAL', base_sd='proxy_decode'),
 dict(file='data/proteome/reversal_emperor.tsv',         layer='proteome', platform='Olink',    intervention='SGLT2i_empagliflozin', tissue='plasma', n_pairs=1134, bmi_axis='goudswaard_ST5_INTERVAL', base_sd='proxy_ukbppp', note='treatment-effect not paired'),
 # transcriptome (baseline recomputable from GSE series matrices)
 dict(file='data/transcriptome/liver_reversal_table.tsv',   layer='transcriptome', platform='HuGene', intervention='diet+bariatric',      tissue='liver',        n_pairs=None, bmi_axis='acquire', base_sd='recompute_GSE83452'),
 dict(file='data/transcriptome/adipose_reversal_table.tsv', layer='transcriptome', platform='microarray', intervention='caloric_restriction', tissue='SAT_adipose', n_pairs=None, bmi_axis='acquire', base_sd='recompute_GSE84046'),
 # GEO expansion. GSE95640 dropped: build failed + redundant with GSE141221 (SAME DiOGenes cohort).
 dict(file='data/transcriptome/gse141221_adipose_reversal_table.tsv', layer='transcriptome', platform='RNAseq', intervention='diet_LCD',  tissue='SAT_adipose',  n_pairs=None, bmi_axis='acquire', base_sd='recompute_GSE141221', note='DiOGenes n=220; LEP -0.68 q1e-26; use INSTEAD of GSE95640'),
 dict(file='data/transcriptome/gse273902_blood_reversal_table.tsv',   layer='transcriptome', platform='RNAseq', intervention='bariatric', tissue='whole_blood',  n_pairs=None, bmi_axis='acquire', base_sd='recompute_GSE273902', note='PMMO bariatric n=24'),
 dict(file='data/transcriptome/gse48452_liver_reversal_table.tsv',    layer='transcriptome', platform='HuGene', intervention='bariatric', tissue='liver',        n_pairs=None, bmi_axis='acquire', base_sd='recompute_GSE48452',  note='Ahrens2013 n=13 -> corroboration-only (below n>=20)'),
 # MAP-expansion panels. base_sd='pending' => score/signature/cross-context only (NOT persistence/determinant until base_sd recomputed).
 dict(file='data/transcriptome/gse157585_muscle_reversal_table.tsv',  layer='transcriptome', platform='RNAseq', intervention='metformin+resistance', tissue='skeletal_muscle', n_pairs=None, bmi_axis='acquire', base_sd='pending', note='MASTERS n=47 — NEW TISSUE muscle (treatment+lifestyle)'),
 dict(file='data/transcriptome/gse43471_adipose_reversal_table.tsv',  layer='transcriptome', platform='microarray', intervention='diet+exercise', tissue='SAT_adipose', n_pairs=None, bmi_axis='acquire', base_sd='pending', note='n=39 lifestyle (diet/exercise)'),
 dict(file='data/transcriptome/gse199063_adipose_reversal_table.tsv', layer='transcriptome', platform='microarray', intervention='bariatric_RYGB', tissue='SAT_adipose', n_pairs=None, bmi_axis='acquire', base_sd='pending', note='n=49 RYGB treatment-surgery; 0/2yr/5yr durability'),
 # MAP-exhaustion wave. GSE117070 EXCLUDED (batch-confounded signs). base_sd computed only for muscle GSE83352.
 dict(file='data/transcriptome/gse83352_muscle_exercise_reversal_table.tsv', layer='transcriptome', platform='array', intervention='exercise', tissue='skeletal_muscle', n_pairs=None, bmi_axis='acquire', base_sd='recompute_GSE83352', note='STRRIDE n=42 — MUSCLE in determinant meta (5th tissue)'),
 dict(file='data/transcriptome/gse106737_liver_bariatric_reversal_table.tsv', layer='transcriptome', platform='array', intervention='bariatric_RYGB', tissue='liver', n_pairs=None, bmi_axis='acquire', base_sd='pending', note='n=41 — 2nd liver context (CRP/SAA/SERPINE1 down)'),
 dict(file='data/transcriptome/gse77962_adipose_diet_reversal_table.tsv', layer='transcriptome', platform='array', intervention='diet-LCD/VLCD', tissue='SAT_adipose', n_pairs=None, bmi_axis='acquire', base_sd='pending', note='n=48 — diet + maintenance'),
 dict(file='data/transcriptome/gse28358_pbmc_diet_reversal_table.tsv', layer='transcriptome', platform='array', intervention='diet-Mediterranean', tissue='PBMC', n_pairs=None, bmi_axis='acquire', base_sd='pending', note='n=22 PREDIMED PBMC'),
 # methylation (base_sd already on disk)
 dict(file='data/methylation/directplus_reversal.tsv', layer='methylation', platform='EPIC', intervention='diet',       tissue='blood',   n_pairs=None, bmi_axis='acquire', base_sd='measured'),
 dict(file='data/methylation/drift2_reversal.tsv',     layer='methylation', platform='EPIC', intervention='behavioral', tissue='blood',   n_pairs=None, bmi_axis='acquire', base_sd='measured'),
 dict(file='data/methylation/predimed_reversal.tsv',   layer='methylation', platform='450K', intervention='diet',       tissue='blood',   n_pairs=None, bmi_axis='acquire', base_sd='measured'),
 dict(file='data/methylation/leipzig_reversal.tsv',    layer='methylation', platform='EPIC', intervention='bariatric',  tissue='blood',   n_pairs=None, bmi_axis='acquire', base_sd='measured'),
 dict(file='data/methylation/benton_sat_reversal.tsv', layer='methylation', platform='450K', intervention='RYGB',       tissue='SAT_adipose', n_pairs=None, bmi_axis='acquire', base_sd='measured'),
 # metabolome
 dict(file='data/metabolome/bbs_reversal_metabolon.tsv', layer='metabolome', platform='Metabolon', intervention='bariatric', tissue='serum', n_pairs=264, bmi_axis='ondisk_bmi_mr', base_sd='missing_summary'),
 dict(file='data/metabolome/reversal_raw/bagheri_reversal.tsv', layer='metabolome', platform='Metabolon', intervention='bariatric', tissue='plasma', n_pairs=104, bmi_axis='borrow_bbs', base_sd='missing_summary'),
 dict(file='data/metabolome/reversal_raw/palau_reversal.tsv',   layer='metabolome', platform='Biocrates', intervention='bariatric', tissue='serum', n_pairs=None, bmi_axis='borrow_bbs', base_sd='recompute_palau_S8'),
 dict(file='data/metabolome/reversal_raw/mtbls218_reversal.tsv',layer='metabolome', platform='LC-HRMS',  intervention='bariatric', tissue='serum', n_pairs=44,  bmi_axis='borrow_bbs', base_sd='recompute_mtbls218_maf'),

 # ===== FULL-SET EXTENSION (analyzed datasets beyond the 24 core scoring panels) =====
 # --- proteome drug panels (published affinity summary stats; treatment-effect betas) ---
 dict(file='data/proteome/step_semaglutide_proteome_reversal.tsv', layer='proteome', platform='SomaScan', intervention='GLP1RA_semaglutide', tissue='plasma', n_pairs=1133, bmi_axis='goudswaard_ST5_INTERVAL', base_sd='proxy_decode',  role='drug-panel', status='ANALYZED', note='STEP semaglutide Nat Med 2025 (PMID 39753963); treatment-effect'),
 dict(file='data/proteome/nowak_metformin_proteome_reversal.tsv',  layer='proteome', platform='Olink',    intervention='metformin',            tissue='plasma', n_pairs=98,   bmi_axis='goudswaard_ST5_INTERVAL', base_sd='proxy_ukbppp', role='drug-panel', status='ANALYZED', note='Nowak metformin; treatment-effect'),
 dict(file='data/proteome/bariatric_ms_proteome_reversal.tsv',     layer='proteome', platform='MS',       intervention='bariatric',            tissue='plasma', n_pairs=44,   bmi_axis='goudswaard_ST5_INTERVAL', base_sd='proxy_ukbppp', role='drug-panel', status='ANALYZED', note='open MS proteomics fallback; per-protein n_pairs<=44'),
 dict(file='data/proteome/liraglutide_proteome_reversal.tsv',      layer='proteome', platform='Olink',    intervention='GLP1RA_liraglutide',   tissue='plasma', n_pairs=20,   bmi_axis='goudswaard_ST5_INTERVAL', base_sd='proxy_ukbppp', role='drug-panel', status='corroboration-only', note='n=20 corroboration-only'),
 dict(file='data/proteome/dapagliflozin_proteome_reversal.tsv',    layer='proteome', platform='Olink',    intervention='SGLT2i_dapagliflozin', tissue='plasma', n_pairs=-1,   bmi_axis='goudswaard_ST5_INTERVAL', base_sd='proxy_ukbppp', role='drug-panel', status='corroboration-only', note='corroboration-only; small panel'),
 # --- transcriptome context / durability panels (small-n; analyzed for cross-context, not scoring) ---
 dict(file='data/transcriptome/gse199063_adipose_reversal_table.tsv', layer='transcriptome', platform='microarray', intervention='bariatric_RYGB', tissue='SAT_adipose', n_pairs=None, bmi_axis='acquire', base_sd='pending', role='context-durability', status='ANALYZED', note='RYGB 0/2yr/5yr durability series (already a scoring panel; durability role)') if False else None,
 dict(file='data/transcriptome/gse107894_adipose_metformin_reversal_table.tsv', layer='transcriptome', platform='RNAseq', intervention='metformin', tissue='SAT_adipose', n_pairs=None, bmi_axis='acquire', base_sd='pending', role='context-durability', status='ANALYZED', note='adipose metformin context; platform corrected array->RNAseq 2026-07-17 (GEO GSE107894 !Series_type = Expression profiling by high throughput sequencing; GPL16791 Illumina HiSeq 2500)'),
 dict(file='data/transcriptome/gse161643_muscle_bariatric_reversal_table.tsv',  layer='transcriptome', platform='array', intervention='bariatric', tissue='skeletal_muscle', n_pairs=None, bmi_axis='acquire', base_sd='pending', role='context-durability', status='ANALYZED', note='muscle bariatric context; platform corrected RNAseq->array 2026-07-17 (GEO GSE161643 !Series_type = Expression profiling by array; GPL570 Affymetrix HG-U133_Plus_2)'),
 dict(file='data/transcriptome/gse224310_adipose_exercise_reversal_table.tsv',  layer='transcriptome', platform='array', intervention='exercise', tissue='SAT_adipose', n_pairs=None, bmi_axis='acquire', base_sd='pending', role='context-durability', status='ANALYZED', note='adipose exercise context; platform corrected RNAseq->array 2026-07-17 (GEO GSE224310 !Series_type = Expression profiling by array; GPL23159 Affymetrix Clariom S)'),
 dict(file='data/transcriptome/gse310742_blood_glp1ra_reversal_table.tsv',      layer='transcriptome', platform='RNAseq', intervention='GLP1RA', tissue='blood_neutrophil', n_pairs=None, bmi_axis='acquire', base_sd='pending', role='context-durability', status='ANALYZED', note='blood GLP1RA context; tissue corrected whole_blood->blood_neutrophil 2026-07-17 (GEO GSE310742 !Series_overall_design = "Neutrophils (isolated by PolymorphPrep) from patients (n=11) with obesity and cardiovascular disease ... before and after semaglutide treatment for 6 months" -- an isolated polymorphonuclear fraction, i.e. NOT whole blood and NOT PBMC). Grouped with the blood-family panels in the coarse 4-tissue analysis axis (as PBMC already is); see build_tissue_organization.py PAN.'),
 dict(file='data/transcriptome/gse70529_adipose_dose_reversal_table.tsv',       layer='transcriptome', platform='array', intervention='caloric_restriction_dose', tissue='SAT_adipose', n_pairs=None, bmi_axis='acquire', base_sd='pending', role='context-durability', status='ANALYZED', note='adipose CR dose context'),
 # --- transcriptome batch-confounded EXCLUSION ---
 dict(file='data/transcriptome/gse117070_muscle_exercise_reversal_table.tsv',   layer='transcriptome', platform='array', intervention='exercise', tissue='skeletal_muscle', n_pairs=None, bmi_axis='acquire', base_sd='pending', role='context-durability', status='EXCLUDED', note='EXCLUDED — batch-inverted directionality'),
 # --- cross-layer panels (same GSE measured across layers; analyzed for cross-layer concordance) ---
 dict(file='data/transcriptome/gse60655_muscle_expression_reversal_table.tsv',  layer='transcriptome', platform='RNAseq', intervention='exercise', tissue='skeletal_muscle', n_pairs=None, bmi_axis='acquire', base_sd='pending', role='cross-layer', status='ANALYZED', note='GSE60590 expression arm; same endurance-training study as the GSE60655 methylation panel (cross-layer). NOTE the filename says gse60655 for historical reasons but the DATA are GSE60590; platform corrected array->RNAseq 2026-07-17 (GEO GSE60590 !Series_type = Expression profiling by high throughput sequencing; GPL11154 Illumina HiSeq 2000). Not a member of the 18-panel tissue-pair set (see variance_partition.py LABEL2FILE), so this correction does not affect the variance partition.'),
 dict(file='data/transcriptome/gse193771_blood_exercise_reversal_table.tsv',    layer='transcriptome', platform='RNAseq', intervention='exercise', tissue='whole_blood', n_pairs=None, bmi_axis='acquire', base_sd='pending', role='cross-layer', status='ANALYZED', note='GSE193771 expression arm (cross-layer)'),
 # --- methylation cross-layer / comparator panels ---
 dict(file='data/methylation/gse60655_muscle_exercise_reversal.tsv',   layer='methylation', platform='450K', intervention='exercise', tissue='skeletal_muscle', n_pairs=None, bmi_axis='acquire', base_sd='measured', role='cross-layer', status='ANALYZED', note='GSE60655 methylation arm (cross-layer w/ the GSE60590 expression panel)'),
 dict(file='data/methylation/gse193730_blood_exercise_reversal.tsv',   layer='methylation', platform='EPIC', intervention='exercise', tissue='blood', n_pairs=None, bmi_axis='acquire', base_sd='measured', role='cross-layer', status='ANALYZED', note='GSE193730 blood exercise (cross-layer)'),
 dict(file='data/methylation/gse171140_muscle_exercise_reversal.tsv',  layer='methylation', platform='EPIC', intervention='exercise', tissue='skeletal_muscle', n_pairs=None, bmi_axis='acquire', base_sd='measured', role='context-durability', status='ANALYZED', note='muscle exercise methylation context'),
 dict(file='data/methylation/gse114763_muscle_detraining_reversal.tsv',layer='methylation', platform='EPIC', intervention='detraining', tissue='skeletal_muscle', n_pairs=None, bmi_axis='acquire', base_sd='measured', role='comparator', status='RETAINED-memory-comparator', note='RETAINED-memory-comparator (detraining = epigenetic memory)'),
 # --- proteome durability anchor (Yousri 2022, PMID 34796696; RYGB SomaLogic 2yr + 12yr arms) ---
 dict(file='data/proteome/yousri_rygb_2yr_reversal.tsv', layer='proteome', platform='SomaScan', intervention='bariatric_RYGB', tissue='plasma', n_pairs=200, bmi_axis='goudswaard_ST5_INTERVAL', base_sd='proxy_decode', role='durability', status='ANALYZED', note='Yousri 2022 (PMID 34796696) RYGB SomaLogic 2yr arm; durability anchor (2yr->12yr rho=0.67)'),
 dict(file='data/proteome/yousri_rygb_12yr.tsv',         layer='proteome', platform='SomaScan', intervention='bariatric_RYGB', tissue='plasma', n_pairs=200, bmi_axis='goudswaard_ST5_INTERVAL', base_sd='proxy_decode', role='durability', status='ANALYZED', note='Yousri 2022 (PMID 34796696) RYGB SomaLogic 12yr arm; 84% same-dir vs 2yr, |12yr|/|2yr|=1.04'),
 # --- methylation cross-cohort replication (E-MTAB-8956 blood diet n=120) ---
 dict(file='data/methylation/central_emtab8956_reversal.tsv', layer='methylation', platform='EPIC', intervention='diet', tissue='blood', n_pairs=120, bmi_axis='acquire', base_sd='measured', role='cross-cohort-replication', status='ANALYZED', note='E-MTAB-8956 blood diet n=120; cross-cohort methylation replication panel'),
]
PANELS = [m for m in PANELS if m is not None]

SIG_Q = ['q','rev_q','fdr','rev_q_adj','fdr_adj']         # FDR cols -> threshold 0.10
SIG_P = ['rev_p','p','rev_p_adj','p_adj','pval']           # p-only fallback -> 0.05
SE    = ['rev_se','se','rev_se_adj','se_adj']
EFF   = ['rev_beta','effect','delta_beta','rev_beta_adj','effect_adj']

def first(cols, cands):
    for c in cands:
        if c in cols: return c
    return None

rows=[]
for m in PANELS:
    f=P(m['file'])
    if not os.path.exists(f):
        rows.append({**m, 'status':'MISSING'}); continue
    df=pd.read_csv(f, sep=None, engine='python')
    cols=list(df.columns)
    n_marks=len(df)
    se_col=first(cols,SE); has_se=se_col is not None
    qcol=first(cols,SIG_Q); pcol=first(cols,SIG_P)
    if qcol: sig=int((pd.to_numeric(df[qcol],errors='coerce')<0.10).sum()); sig_basis=f'{qcol}<0.10'
    elif pcol: sig=int((pd.to_numeric(df[pcol],errors='coerce')<0.05).sum()); sig_basis=f'{pcol}<0.05'
    else: sig=-1; sig_basis='none'
    # n_pairs: read from table if metadata None
    npairs=m['n_pairs']
    if npairs is None:
        npc=first(cols,['n_pairs','n'])
        npairs=int(pd.to_numeric(df[npc],errors='coerce').dropna().iloc[0]) if npc else -1
    eff_col=first(cols,EFF)
    # SPLIT power floor (SCORE_SPEC_FROZEN.md s4): score n>=20 & >=1 sig ; persistence n>=25 & measured/recompute base_sd
    contributes = (npairs>=20) and (sig>=1)
    base_measured = str(m['base_sd']).startswith(('measured','recompute_'))
    persistence_eligible = contributes and (npairs>=25) and base_measured
    # role/status: default the 24 core scoring panels to 'primary'/'ANALYZED'; new panels carry explicit labels.
    role   = m.get('role','primary')
    status = m.get('status','ANALYZED')
    if status=='EXCLUDED':            # excluded panels never count toward score/persistence
        contributes=False; persistence_eligible=False
    rows.append(dict(panel=os.path.basename(m['file']).replace('_reversal','').replace('.tsv',''),
        layer=m['layer'], platform=m['platform'], intervention=m['intervention'], tissue=m['tissue'],
        n_pairs=npairs, n_marks=n_marks, has_se=has_se, eff_col=eff_col, sig_marks=sig, sig_basis=sig_basis,
        bmi_axis=m['bmi_axis'], base_sd=m['base_sd'], role=role, status=status,
        contributes_to_score=contributes,
        persistence_eligible=persistence_eligible, note=m.get('note',''), path=m['file']))

out=pd.DataFrame(rows)
op=S('results/panel_manifest_full.tsv')
out.to_csv(op, sep='\t', index=False)
# console summary
cols_show=['panel','layer','intervention','n_pairs','n_marks','sig_marks','role','status','base_sd','contributes_to_score','persistence_eligible']
with pd.option_context('display.width',240,'display.max_columns',50,'display.max_colwidth',26):
    print(out[cols_show].to_string(index=False))
print(f'\nwrote {op}')

# ---- FULL analyzed-dataset accounting (for Table 1 / abstract) ----
analyzed = out[out.status!='MISSING'].copy()
included = analyzed[analyzed.status!='EXCLUDED']   # EXCLUDED panels are reported but not counted as analyzed-in
print('\n=== PER-LAYER analyzed-dataset counts (status != EXCLUDED) ===')
print(included.groupby('layer').size().to_string())
print(f'\nTOTAL analyzed datasets (incl. all roles, excl. EXCLUDED): {len(included)}')
print(f'TOTAL incl. EXCLUDED (all rows on disk): {len(analyzed)}')
print('\n=== by role (analyzed, status != EXCLUDED) ===')
print(included.groupby('role').size().to_string())
print('\n=== EXCLUDED ===', ', '.join(out.loc[out.status=="EXCLUDED","panel"].tolist()) or 'none')
print('=== RETAINED-memory-comparator ===', ', '.join(out.loc[out.status=="RETAINED-memory-comparator","panel"].tolist()) or 'none')
print('=== MISSING ===', ', '.join(out.loc[out.status=="MISSING","panel"].tolist()) if (out.status=='MISSING').any() else 'none')
print('\nSCORE-contributing panels per layer (n>=20 & >=1 sig):')
print(out[out.contributes_to_score].groupby('layer').size().to_string())
print('\nPERSISTENCE-eligible panels per layer:')
pe=out[out.persistence_eligible]
print(pe.groupby('layer').size().to_string() if len(pe) else '  (none)')

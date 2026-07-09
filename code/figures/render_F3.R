# render_F3.R — Figure 3 (HEADLINE): determinant null + constructive predictive null.
# R/ggplot following publication-figures-tables skill + LOCKED theme (theme_publication_ckm.R).
# v2: equal 2x3 grid (column alignment), tag/title separation, abbreviated labels, decluttered panels.
suppressMessages({library(ggplot2); library(patchwork); library(scales)})
source("C:/Users/Bert/Downloads/CKM papers/paper 3.eas/analysis/tools/theme_publication_ckm.R")
set_publication_defaults_ckm()
RES <- "C:/Users/Bert/Downloads/CKM papers/weight-loss.omics/signature-pivot/results"
rd  <- function(p) read.delim(file.path(RES, p), stringsAsFactors = FALSE)
TXN<-"#009E73"; PRO<-"#0072B2"; NEU<-"#333333"; UP<-"#D55E00"; DN<-"#56B4E9"; BAND<-"#E4E4E4"

dm <- rd("determinant_meta.tsv"); dt <- subset(dm, layer=="transcriptome"); dp <- subset(dm, layer=="proteome")
pn <- rd("predictive_null.tsv");  ab <- rd("predictive_null_ablation.tsv"); tau <- rd("orphan/tau_dissociation.tsv")
NAME <- c(loeuf="LOEUF", n_drug_log="druggability", n_gwas_log="GWAS burden",
 tau="tissue-spec.", arch_nsig="cis signals", arch_str="cis F",
 has_arch="has cis-QTL", causal_nonEGFR="causal status", is_enzyme="enzyme", is_membrane="membrane",
 is_secreted="secreted", loeuf_miss="LOEUF NA", tau_miss="tissue-spec. NA")
mde <- median(dt$mde)
ebh <- function(...) geom_errorbar(..., orientation="y", width=0)   # horizontal error bar (ggplot2 4.0)

# a) determinant forest (uniform — no feature is a determinant)
dt2 <- dt[order(dt$pooled_beta),]; dt2$lab <- factor(NAME[dt2$feature], levels=NAME[dt2$feature])
pa <- ggplot(dt2, aes(pooled_beta, lab)) +
  annotate("rect", xmin=-mde, xmax=mde, ymin=-Inf, ymax=Inf, fill=BAND) +
  geom_vline(xintercept=0, linewidth=0.4) +
  ebh(aes(xmin=pooled_beta-ci_halfwidth, xmax=pooled_beta+ci_halfwidth), color=NEU, linewidth=0.5) +
  geom_point(color=TXN, size=1.6) +
  scale_x_continuous(limits=c(-0.08,0.08), breaks=c(-0.05,0,0.05)) +
  annotate("text", x=0.078, y=1.4, label="all |β| < 0.029", hjust=1, size=2.3, fontface="italic") +
  labs(x="pooled standardized β", y=NULL, title="Determinant model (k=9)")

# b) causal status on the RAW / interpretable scale — the standardized coef in (a) is
#    compressed ~50x by this feature's rare-binary SD (9 of ~18,600 transcripts), so it is
#    not an effect size; show the raw group difference (Supplementary Table 14). Underpowered.
cse <- rd("causal_status_effect.tsv")
csr <- cse[match(c("raw_group_smd_unadjusted","multivariable_raw_coef","nominatable_universe_smd"), cse$estimate_name),]
csr$lab <- factor(c("unadjusted","adjusted","nominatable\nuniverse"),
                  levels=c("nominatable\nuniverse","adjusted","unadjusted"))
pb <- ggplot(csr, aes(pooled, lab)) +
  geom_vline(xintercept=0, linewidth=0.4) +
  ebh(aes(xmin=ci_lo, xmax=ci_hi), color=NEU, linewidth=0.6) +
  geom_point(color=TXN, size=2.4) +
  geom_text(aes(label=sprintf("%+.2f", pooled)), vjust=-1.05, size=2.2) +
  scale_x_continuous(limits=c(-0.46,0.20), breaks=c(-0.4,-0.2,0)) +
  annotate("text", x=-0.45, y=0.58, hjust=0, vjust=1, size=1.75, fontface="italic", color="grey35", lineheight=0.95,
           label="(a) std. coef is compressed;\n(b) raw effect, see Table S14") +
  labs(x="raw group difference (SD)", y=NULL, title="Genetic causal status:\nunderpowered (n = 9)")

# c) layer power
lp <- rbind(data.frame(layer="transcriptome\n(k=9)", ci=dt$ci_halfwidth),
            data.frame(layer="proteome\n(k=3)", ci=dp$ci_halfwidth))
lp$layer <- factor(lp$layer, levels=c("transcriptome\n(k=9)","proteome\n(k=3)"))
pc <- ggplot(lp, aes(layer, ci, color=layer)) +
  geom_jitter(width=0.16, height=0, size=1.4, alpha=0.9) +
  scale_color_manual(values=c("transcriptome\n(k=9)"=TXN, "proteome\n(k=3)"=PRO), guide="none") +
  scale_y_log10() +
  labs(x=NULL, y="95% CI half-width (SD)", title="Layer power", subtitle="powered vs inconclusive")

# d) predictive null
pnl <- data.frame(ctx=factor(pn$context, levels=pn$context), within=pn$within_spearman, loco=pn$loco_spearman)
pdn <- ggplot(pnl) +
  geom_segment(aes(x=ctx, xend=ctx, y=within, yend=loco), color="grey70", linewidth=0.5) +
  geom_hline(yintercept=0, linewidth=0.4) +
  geom_point(aes(ctx, loco, fill="leave-one-out (LOCO)"), shape=21, color=TXN, size=2, stroke=0.8) +
  geom_point(aes(ctx, within, fill="within-context"), shape=21, color=TXN, size=2.1, stroke=0.5) +
  scale_fill_manual(values=c("within-context"=TXN, "leave-one-out (LOCO)"="white"), name=NULL,
                    breaks=c("within-context","leave-one-out (LOCO)")) +
  scale_y_continuous(limits=c(-0.06,0.31), breaks=c(0,0.1,0.2,0.3)) +
  annotate("text", x=3, y=-0.048, label="LOCO mean +0.002", size=2.2, fontface="italic") +
  labs(x=NULL, y="Spearman (pred vs actual)", title="Predictive null") +
  theme(axis.text.x=element_text(angle=35, hjust=1),
        legend.position="inside", legend.position.inside=c(0.30,1.0), legend.justification=c(0,1),
        legend.key.size=unit(7,"pt"), legend.text=element_text(size=6.3),
        legend.background=element_rect(fill=alpha("white",0.8), color=NA))

# e) ablation
abm <- rbind(data.frame(ctx=ab$context, set="all", auc=ab$loco_auc_full),
             data.frame(ctx=ab$context, set="intrinsic", auc=ab$loco_auc_intrinsic_only),
             data.frame(ctx=ab$context, set="measurability", auc=ab$loco_auc_measurability_only))
abm$ctx <- factor(abm$ctx, levels=ab$context); abm$set <- factor(abm$set, levels=c("all","intrinsic","measurability"))
pe <- ggplot(abm, aes(ctx, auc, fill=set)) +
  geom_col(position=position_dodge(0.8), width=0.74, color="black", linewidth=0.2) +
  geom_hline(yintercept=0.5, linetype="dashed", linewidth=0.4) +
  scale_fill_manual(values=c("all"=NEU, "intrinsic"=TXN, "measurability"=PRO), name=NULL) +
  coord_cartesian(ylim=c(0.45,0.85)) +
  labs(x=NULL, y="cross-context AUC", title="Detectability ablation") +
  theme(axis.text.x=element_text(angle=35, hjust=1),
        legend.position="inside", legend.position.inside=c(0.40,0.99), legend.justification=c(0,1),
        legend.key.size=unit(7,"pt"), legend.text=element_text(size=6.3),
        legend.background=element_rect(fill=alpha("white",0.7), color=NA))

# f) tau dissociation
tf <- data.frame(x=factor(c("magnitude","restoration","restoration (adj.)"),
                          levels=c("magnitude","restoration","restoration (adj.)")),
                 v=c(tau$tau_beta[tau$col=="y_mag"], tau$tau_beta[tau$col=="y_dir"], tau$tau_beta[tau$col=="y_dir_ca"]))
tf$lbl <- sprintf("%+.3f", tf$v); tf$vj <- ifelse(tf$v>=0, -0.5, 1.3)
pf <- ggplot(tf, aes(x, v, fill=x)) +
  geom_col(width=0.62, color="black", linewidth=0.3) +
  geom_hline(yintercept=0, linewidth=0.4) +
  geom_text(aes(label=lbl), vjust=tf$vj, size=2.3) +
  scale_fill_manual(values=setNames(c(UP,DN,"#9ecae1"), levels(tf$x)), guide="none") +
  coord_cartesian(ylim=c(-0.08,0.13)) +
  labs(x=NULL, y="τ → outcome (β)", title="τ: magnitude, not direction") +
  theme(axis.text.x=element_text(angle=35, hjust=1))

fig <- wrap_plots(pa, pb, pc, free(pdn, type="space", side="l"), pe, pf, ncol=3) +
  plot_annotation(tag_levels="a") &
  theme(plot.tag = element_text(face="bold", size=11, hjust=0, vjust=1),
        plot.tag.position = c(0.01, 1.02),
        plot.title = element_text(size=8.5, hjust=0, margin=margin(t=2, b=4)),
        plot.subtitle = element_text(size=7, color="grey25", margin=margin(b=3)),
        plot.margin = margin(11, 6, 10, 6, "pt"))
save_figure_ckm(fig, "Figure_3", width_mm=183, height_mm=170, output_dir=".")
cat("done\n")

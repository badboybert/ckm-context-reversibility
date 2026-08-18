# render_F3.R — Figure 3 (HEADLINE): determinant null + constructive predictive null.
# R/ggplot, LOCKED publication theme (theme_publication_ckm.R).
# v2: equal 2x3 grid (column alignment), tag/title separation, abbreviated labels, decluttered panels.
suppressMessages({library(ggplot2); library(patchwork); library(scales)})
source("theme_publication_ckm.R")  # bundled alongside this script
set_publication_defaults_ckm()
RES <- "../../results"
rd  <- function(p) read.delim(file.path(RES, p), stringsAsFactors = FALSE)
TXN<-"#009E73"; PRO<-"#0072B2"; NEU<-"#333333"; UP<-"#D55E00"; DN<-"#56B4E9"; BAND<-"#E4E4E4"

dm <- rd("determinant_meta.tsv"); dt <- subset(dm, layer=="transcriptome"); dp <- subset(dm, layer=="proteome")
pn <- rd("predictive_null.tsv");  ab <- rd("predictive_null_ablation.tsv"); tau <- rd("orphan/tau_dissociation.tsv")
NAME <- c(loeuf="LOEUF", n_drug_log="druggability", n_gwas_log="GWAS burden",
 tau="tissue specificity", arch_nsig="cis-eQTL signals", arch_str="cis-instrument F",
 has_arch="has cis-QTL", causal_nonEGFR="causal status", is_enzyme="enzyme", is_membrane="membrane",
 is_secreted="secreted", loeuf_miss="LOEUF missing", tau_miss="tissue specificity missing")
mde <- median(dt$mde)
ebh <- function(...) geom_errorbar(..., orientation="y", width=0)   # horizontal error bar (ggplot2 4.0)

# a) determinant forest. The RARE-BINARY genetic causal status
#    be visually separated from the continuous features: its standardized coefficient is compressed
#    ~50x by SD_x ~ 0.020 (9 of ~18,600 transcripts), so it is not on the same footing as the others
#    and must not read as "another equivalent-to-zero feature". It is therefore sorted to the bottom,
#    drawn as an OPEN point, separated by a rule, and labelled with its n. Panel b gives its raw scale.
dt$is_rare_binary <- dt$feature == "causal_nonEGFR"
dt2 <- dt[order(dt$is_rare_binary, dt$pooled_beta, decreasing=c(TRUE, FALSE), method="radix"),]
dt2$lab <- ifelse(dt2$is_rare_binary, "causal status\n(rare binary, n=9)", NAME[dt2$feature])
dt2$lab <- factor(dt2$lab, levels=dt2$lab)
y_sep <- sum(dt2$is_rare_binary) + 0.5          # rule between the binary and the continuous block
pa <- ggplot(dt2, aes(pooled_beta, lab)) +
  # nested bands: outer light = ±0.05 SD SESOI (the equivalence bound, the headline test);
  # inner darker = ±0.013 SD median MDE (detectability). Dashed guides at the SESOI edge.
  annotate("rect", xmin=-0.05, xmax=0.05, ymin=-Inf, ymax=Inf, fill="#F2F2F2") +
  annotate("rect", xmin=-mde, xmax=mde, ymin=-Inf, ymax=Inf, fill=BAND) +
  geom_vline(xintercept=c(-0.05, 0.05), linetype="dashed", linewidth=0.35, color="grey60") +
  geom_vline(xintercept=0, linewidth=0.4) +
  geom_hline(yintercept=y_sep, linetype="dotted", linewidth=0.4, color="grey45") +
  ebh(aes(xmin=pooled_beta-ci_halfwidth, xmax=pooled_beta+ci_halfwidth), color=NEU, linewidth=0.5) +
  geom_point(aes(shape=is_rare_binary), color=TXN, fill="white", size=1.6, stroke=0.6) +
  scale_shape_manual(values=c(`FALSE`=16, `TRUE`=21), guide="none") +
  scale_x_continuous(limits=c(-0.085,0.085), breaks=c(-0.05,0,0.05)) +
  # No in-panel band key. Panel a spans only 0.17 SD, so any key either clips at the axis or lands
  # on a data row -- it did both. The outer (+/-0.05 SD SESOI) and inner (+/-0.013 SD MDE) bands are
  # defined in the figure legend, which is itself "outside the data region".
  annotate("text", x=0.083, y=0.62, label="rare binary — see b", hjust=1, size=2.5,
           fontface="italic", color="grey35") +
  labs(x="pooled standardized β", y=NULL)

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
  geom_text(aes(label=sprintf("%+.2f", pooled)), vjust=-1.05, size=2.8) +
  scale_x_continuous(limits=c(-0.46,0.20), breaks=c(-0.4,-0.2,0)) +
  # (in-panel footnote REMOVED: the legend already states that the standardized coefficient is
  #  compressed, is not an effect size, and is reported on the raw scale in Supplementary Table 14.
  #  Any text long enough to say it ran across the x=0 reference line, which struck through it;
  #  a white label box would instead occlude that line. The legend carries it.)
  labs(x="raw group difference (SD)", y=NULL)

# c) layer power
lp <- rbind(data.frame(layer="transcriptome\n(k=9)", ci=dt$ci_halfwidth),
            data.frame(layer="proteome\n(k=3)", ci=dp$ci_halfwidth))
lp$layer <- factor(lp$layer, levels=c("transcriptome\n(k=9)","proteome\n(k=3)"))
set.seed(11)   # geom_jitter below is stochastic; seed so the figure is reproducible
pc <- ggplot(lp, aes(layer, ci, color=layer)) +
  geom_jitter(width=0.16, height=0, size=1.4, alpha=0.9) +
  scale_color_manual(values=c("transcriptome\n(k=9)"=TXN, "proteome\n(k=3)"=PRO), guide="none") +
  scale_y_log10() +
  labs(x=NULL, y="95% CI half-width (SD)") +
  theme(axis.text.x=element_text(angle=35, hjust=1, lineheight=0.9))  # 35° = same tilt as d/e/f

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
  annotate("text", x=3, y=-0.048, label="LOCO mean +0.002", size=2.8, fontface="italic") +
  labs(x=NULL, y="Spearman (pred vs actual)") +
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
  labs(x=NULL, y="cross-context AUC") +
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
  geom_text(aes(label=lbl), vjust=tf$vj, size=2.8) +
  scale_fill_manual(values=setNames(c(UP,DN,"#9ecae1"), levels(tf$x)), guide="none") +
  coord_cartesian(ylim=c(-0.08,0.13)) +
  labs(x=NULL, y="τ → outcome (β)") +
  theme(axis.text.x=element_text(angle=35, hjust=1))

# free(pc, side="b"): panel c's tilted two-line x-labels are taller than a/b's, and patchwork
# aligns the shared bottom-space row, which pushed a/b's x-axis titles far below their axes.
fig <- wrap_plots(pa, pb, free(pc, type="space", side="b"), free(pdn, type="space", side="l"), pe, pf, ncol=3) +
  plot_annotation(tag_levels="a") &
  theme(plot.tag = element_text(face="bold", size=11, hjust=0, vjust=1),
        plot.tag.position = c(0.01, 1.02),
        plot.margin = margin(11, 6, 10, 6, "pt"))
save_figure_ckm(fig, "Figure_3", width_mm=183, height_mm=170, output_dir=".")
cat("done\n")

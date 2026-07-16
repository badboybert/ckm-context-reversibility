# render_F4.R — Figure 4: Weight loss restores the molecular landscape toward lean, scaled by intensity;
# where reversal occurs it is durable. R/ggplot, LOCKED theme. Source: prep_F4_data.py -> source_data/F4*.
suppressMessages({library(ggplot2); library(patchwork); library(scales)})
source("C:/Users/Bert/Downloads/CKM papers/paper 3.eas/analysis/tools/theme_publication_ckm.R")
set_publication_defaults_ckm()
SD <- "source_data"; rd <- function(p) read.csv(file.path(SD, p), stringsAsFactors = FALSE)
PRO<-"#0072B2"; METAB<-"#CC79A7"; METH<-"#D55E00"; NEU<-"#333333"
TIS <- c(plasma="#4D4D4D", liver="#8C6D31", adipose="#E1A100", muscle="#B2182B", blood="#2166AC")

# ---- a) Restoration toward lean: distribution of R (proteome + metabolome) ----
dist <- rd("F4a_restoration_dist.csv"); summ <- rd("F4a_summary.csv")
dist$layer <- factor(dist$layer, levels=c("metabolome","proteome"))
summ$layer <- factor(summ$layer, levels=c("metabolome","proteome"))
summ$lab <- sprintf("median %+.2f · %.0f%% CI>0", summ$median_R, summ$pct_restorative)
pa <- ggplot(dist, aes(R_restoration, layer, fill=layer)) +
  geom_vline(xintercept=0, linewidth=0.4, color="grey55") +
  geom_violin(scale="width", color=NA, alpha=0.55, width=0.8) +
  geom_point(data=summ, aes(median_R, layer), shape=23, size=2.4, fill="white", color="black", stroke=0.6) +
  geom_text(data=summ, aes(x=median_R, y=layer, label=lab), vjust=-1.6, size=2.4, color=NEU) +
  scale_fill_manual(values=c(proteome=PRO, metabolome=METAB), guide="none") +
  scale_x_continuous(limits=c(-5,6), breaks=c(-4,-2,0,2,4,6)) +
  annotate("text", x=-4.6, y=2.5, label="toward lean →", hjust=0, size=2.3, fontface="italic", color="grey35") +
  labs(x="restoration toward lean (R)", y=NULL)

# ---- b) % toward lean, scaled by weight-loss intensity ----
pb_d <- rd("F4b_pct_toward_lean.csv")
pb_d$tissue <- factor(pb_d$tissue, levels=names(TIS))
pb_d$lab <- sprintf("%s (%s)", pb_d$context, pb_d$tissue)
pb_d <- pb_d[order(pb_d$pct_toward_lean),]; pb_d$lab <- factor(pb_d$lab, levels=pb_d$lab)
# open circles = NOT distinguishable from the 50% chance line (exact binomial P >= 0.05) — the
# round-2 review asked for the near-chance contexts to be marked.
pb_d$ns <- as.logical(pb_d$ns_vs_chance)
pb <- ggplot(pb_d, aes(pct_toward_lean, lab, color=tissue)) +
  geom_vline(xintercept=50, linetype="dashed", linewidth=0.4, color="grey40") +
  geom_segment(aes(x=50, xend=pct_toward_lean, yend=lab), linewidth=0.5) +
  geom_point(aes(shape=ns, fill=ns), size=2, stroke=0.7) +
  scale_color_manual(values=TIS, name=NULL) +
  scale_shape_manual(values=c(`FALSE`=16, `TRUE`=21), guide="none") +
  scale_fill_manual(values=c(`FALSE`=NA, `TRUE`="white"), guide="none") +
  scale_x_continuous(limits=c(48,70), breaks=c(50,55,60,65,70)) +
  annotate("text", x=50, y=0.4, label="chance", size=2.1, color="grey35", fontface="italic", vjust=0) +
  # placed mid-right: rows 5-8 all end below 55%, and the tissue legend occupies the bottom-right
  # (inside c(0.74,0.04)) — anchoring here keeps the note clear of both the data and the legend.
  annotate("text", x=69.6, y=6.9, label="open = not distinguishable\nfrom chance (P ≥ 0.05)",
           size=1.95, color="grey35", fontface="italic", hjust=1, vjust=0.5, lineheight=0.95) +
  labs(x="% of marks moving toward lean", y=NULL) +
  theme(axis.text.y=element_text(size=9),
        legend.position="inside", legend.position.inside=c(0.74,0.04), legend.justification=c(0,0),
        legend.key.size=unit(7,"pt"), legend.text=element_text(size=6.3),
        legend.background=element_rect(fill=alpha("white",0.7), color=NA))

# ---- c) Durability: adipose 2yr vs 5yr (RYGB) ----
pc_d <- rd("F4c_durability_adipose_points.csv"); ds <- rd("F4c_durability_summary.csv")
lim <- 2.2
pcc <- pmin(pmax(pc_d[,c("rev2","rev5")], -lim), lim)
pc_d$rev2 <- pcc$rev2; pc_d$rev5 <- pcc$rev5
pc <- ggplot(pc_d, aes(rev2, rev5)) +
  annotate("rect", xmin=0, xmax=lim, ymin=0, ymax=lim, fill="#E8F0E8") +
  annotate("rect", xmin=-lim, xmax=0, ymin=-lim, ymax=0, fill="#E8F0E8") +
  geom_hline(yintercept=0, linewidth=0.3, color="grey70") + geom_vline(xintercept=0, linewidth=0.3, color="grey70") +
  geom_abline(slope=1, intercept=0, linetype="dashed", linewidth=0.4, color="grey45") +
  geom_point(size=0.35, alpha=0.25, color="#1B7837") +
  scale_x_continuous(limits=c(-lim,lim), breaks=c(-2,-1,0,1,2)) +
  scale_y_continuous(limits=c(-lim,lim), breaks=c(-2,-1,0,1,2)) +
  annotate("label", x=-2.1, y=2.05, hjust=0, vjust=1, size=2.3, fill="white", linewidth=0, lineheight=0.95,
           label=sprintf("ρ = %+.2f\n%.0f%% sign-retained", ds$spearman_all, ds$sign_retained_sig)) +
  labs(x="2-year reversal", y="5-year reversal")

# ---- e) proteome 2yr->12yr durability (Yousri RYGB SomaLogic) ----
pe_d <- rd("F4e_yousri_points.csv"); lim2 <- 4   # heavy-tailed robust-z; out-of-range points dropped (not clamped) -> no boundary pile-up
pe <- ggplot(pe_d, aes(z2, z12)) +
  annotate("rect", xmin=0, xmax=lim2, ymin=0, ymax=lim2, fill="#E8EEF6") +
  annotate("rect", xmin=-lim2, xmax=0, ymin=-lim2, ymax=0, fill="#E8EEF6") +
  geom_hline(yintercept=0, linewidth=0.3, color="grey70") + geom_vline(xintercept=0, linewidth=0.3, color="grey70") +
  geom_abline(slope=1, intercept=0, linetype="dashed", linewidth=0.4, color="grey45") +
  geom_point(size=0.35, alpha=0.28, color="#0072B2") +
  scale_x_continuous(limits=c(-lim2,lim2), breaks=c(-3,0,3)) +
  scale_y_continuous(limits=c(-lim2,lim2), breaks=c(-3,0,3)) +
  annotate("label", x=-3.9, y=3.8, hjust=0, vjust=1, size=2.3, fill="white", linewidth=0, lineheight=0.95,
           label=sprintf("ρ = %+.2f\n%.0f%% same direction\nmagnitude 1.04×", ds$prot_spearman_2_12yr, ds$prot_sign_retained)) +
  labs(x="2-year reversal (z)", y="12-year reversal (z)")

# ---- d) Methylation restoration = age-drift confounded (caution) ----
pd_d <- rd("F4d_methyl_age_drift.csv")
pdd <- ggplot(pd_d, aes(interval_mo, spearman_delta_bmi)) +
  geom_line(linewidth=0.5, color=METH) +
  geom_point(size=2.4, color=METH) +
  geom_text(aes(label=cohort), vjust=-1.1, size=2.3) +
  scale_x_continuous(limits=c(0,21), breaks=c(3,18)) +
  scale_y_continuous(limits=c(0,0.7), breaks=c(0,0.2,0.4,0.6)) +
  annotate("text", x=10.5, y=0.10, label="epigenetic-ageing drift\ntracks the BMI axis", size=2.2, fontface="italic", color="grey30", lineheight=0.9) +
  labs(x="sampling interval (months)", y="Δ vs BMI-axis (Spearman)")

# layout: top row = restoration (a wide dist, b wide gradient); bottom row = durability pair (c,d) + methylome (e).
# tag order a..e = restoration-dist / %toward-lean / adipose-durability / proteome-durability / methylome-caution.
fig <- (pa + pb) / (pc + pe + pdd) +
  plot_layout(heights=c(1,1)) +
  plot_annotation(tag_levels="a") &
  theme(plot.tag = element_text(face="bold", size=11, hjust=0, vjust=1),
        plot.tag.position = c(0.01, 1.02),
        plot.margin = margin(11, 6, 8, 6, "pt"))
save_figure_ckm(fig, "Figure_4", width_mm=183, height_mm=150, output_dir=".")
cat("done\n")

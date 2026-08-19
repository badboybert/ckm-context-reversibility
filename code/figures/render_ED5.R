# render_ED5.R — Extended Data 5: methylation-expression reversal shows limited concordance in two
# cohort-level tissue comparisons.
# (a) muscle meth-Δ vs expr-Δ (weak negative); (b) blood (no coupling); (c) what the comparison can
# and cannot support. Panels a/b now carry their tissue in the panel subtitle: without it a reader
# cannot tell which scatter is which.
suppressMessages({library(ggplot2); library(patchwork); library(scales)})
source("theme_publication_ckm.R")  # bundled alongside this script
set_publication_defaults_ckm()
SD <- "source_data"; rd <- function(p) read.csv(file.path(SD, p), stringsAsFactors=FALSE)
SIGc<-"#B2182B"; NSc<-"grey65"; NEU<-"#333333"
sm <- rd("ED5_summary.csv")

mkscatter <- function(file, pair) {
  d <- rd(file); d$expr_sig <- d$expr_sig %in% c("True", TRUE, "TRUE")
  qx <- quantile(d$meth_delta, c(.01,.99), na.rm=TRUE); qy <- quantile(d$expr_delta, c(.01,.99), na.rm=TRUE)
  s <- sm[sm$pair==pair,]
  ggplot(d, aes(meth_delta, expr_delta)) +
    geom_hline(yintercept=0, linewidth=0.3, color="grey80") + geom_vline(xintercept=0, linewidth=0.3, color="grey80") +
    geom_point(data=subset(d,!expr_sig), size=0.4, alpha=0.18, color=NSc) +
    geom_point(data=subset(d, expr_sig), size=0.5, alpha=0.5, color=SIGc) +
    geom_smooth(method="lm", se=FALSE, color="black", linewidth=0.5, formula=y~x) +
    coord_cartesian(xlim=qx, ylim=qy) +
    annotate("label", x=qx[1], y=qy[2], hjust=0, vjust=1, size=2.8, fill=NA, linewidth=0, lineheight=0.95,
             label=sprintf("ρ(all) = %+.2f  (n=%d)\nρ(expr-sig) = %+.2f  (n=%d)", s$rho_all, s$n_all, s$rho_sig, s$n_sig)) +
    labs(x="methylation Δβ (post−pre)", y="expression Δ (post−pre)",
         subtitle=sprintf("%s (cohort-level)", pair)) +
    theme(plot.subtitle=element_text(size=6.8, color="grey30"))
}
pa <- mkscatter("ED5_muscle_scatter.csv","muscle")
pb <- mkscatter("ED5_blood_scatter.csv","blood")

# c) scope of the cross-layer comparison — neutral result summary (round 7: replaced a rebuttal-style slide)
wrapc <- function(s, w=30) paste(strwrap(s, width=w), collapse="\n")
pc <- ggplot() + xlim(0,10) + ylim(0,10) +
  annotate("text", x=0.3, y=9.8, hjust=0, vjust=1, size=2.55, fontface="bold", label="Cross-layer concordance") +
  annotate("text", x=0.3, y=9.0, hjust=0, vjust=1, size=2.15, lineheight=1.2, label=paste0(
    wrapc("Gene-level methylation and expression reversal are weakly to non-concordant: muscle ρ = −0.19 / −0.35, blood ρ = −0.01 (ns)."), "\n\n",
    wrapc("Correlations are cohort-level and gene-level (UCSC RefGene), not within-person or promoter-resolved."), "\n\n",
    wrapc("A within-person, promoter-resolved test needs same-subject, multi-layer cohorts."))) +
  theme_void(base_family=BASE_FONT) +
  theme(plot.background=element_rect(fill="white", color=NA),
        panel.background=element_rect(fill="white", color=NA))

fig <- wrap_plots(pa, pb, pc, ncol=3) +
  plot_annotation(tag_levels="a",
    theme=theme(plot.background=element_rect(fill="white", color=NA))) &
  theme(plot.tag=element_text(face="bold", size=11, hjust=0, vjust=1), plot.tag.position=c(0.01,1.04),
        plot.margin=margin(13,6,8,6,"pt"))
save_figure_ckm(fig, "Figure_ED5", width_mm=183, height_mm=82, output_dir=".")
cat("done\n")

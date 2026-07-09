# render_ED5.R — Extended Data 5: cross-layer node concordance is a NO-GO (disclosed honest negative).
# (a) muscle meth-Δ vs expr-Δ (weak negative); (b) blood (no coupling); (c) caveat/disclosure of the abandoned capstone.
suppressMessages({library(ggplot2); library(patchwork); library(scales)})
source("C:/Users/Bert/Downloads/CKM papers/paper 3.eas/analysis/tools/theme_publication_ckm.R")
set_publication_defaults_ckm()
SD <- "source_data"; rd <- function(p) read.csv(file.path(SD, p), stringsAsFactors=FALSE)
SIGc<-"#B2182B"; NSc<-"grey65"; NEU<-"#333333"
sm <- rd("ED5_summary.csv")

mkscatter <- function(file, pair, title, sub) {
  d <- rd(file); d$expr_sig <- d$expr_sig %in% c("True", TRUE, "TRUE")
  qx <- quantile(d$meth_delta, c(.01,.99), na.rm=TRUE); qy <- quantile(d$expr_delta, c(.01,.99), na.rm=TRUE)
  s <- sm[sm$pair==pair,]
  ggplot(d, aes(meth_delta, expr_delta)) +
    geom_hline(yintercept=0, linewidth=0.3, color="grey80") + geom_vline(xintercept=0, linewidth=0.3, color="grey80") +
    geom_point(data=subset(d,!expr_sig), size=0.4, alpha=0.18, color=NSc) +
    geom_point(data=subset(d, expr_sig), size=0.5, alpha=0.5, color=SIGc) +
    geom_smooth(method="lm", se=FALSE, color="black", linewidth=0.5, formula=y~x) +
    coord_cartesian(xlim=qx, ylim=qy) +
    annotate("label", x=qx[1], y=qy[2], hjust=0, vjust=1, size=2.1, fill="white", label.size=0, lineheight=0.95,
             label=sprintf("ρ(all) = %+.2f  (n=%d)\nρ(expr-sig) = %+.2f  (n=%d)", s$rho_all, s$n_all, s$rho_sig, s$n_sig)) +
    labs(x="methylation Δβ (post−pre)", y="expression Δ (post−pre)", title=title, subtitle=sub)
}
pa <- mkscatter("ED5_muscle_scatter.csv","muscle","Muscle: weak negative coupling","GSE60655 meth × GSE60590 expr")
pb <- mkscatter("ED5_blood_scatter.csv","blood","Blood: no coupling","GSE193730 meth × GSE193771 expr")

# c) the NO-GO disclosure (honest negative; the abandoned capstone)
pc <- ggplot() + xlim(0,10) + ylim(0,10) +
  annotate("rect", xmin=0.2, xmax=9.8, ymin=6.4, ymax=9.6, fill="#FBEEE4", color="#D55E00", linewidth=0.4) +
  annotate("text", x=5, y=8.7, size=2.5, fontface="bold", label="Cross-layer node concordance:") +
  annotate("text", x=5, y=7.9, size=2.5, fontface="bold", color="#B2182B", label="NO-GO") +
  annotate("text", x=5, y=7.0, size=2.0, lineheight=1.0, label="the original capstone is NOT a positive claim") +
  annotate("text", x=0.5, y=5.2, hjust=0, size=2.05, lineheight=1.1, label=
    "• Same gene's CpG vs transcript:\n   weak (muscle −0.19 / −0.35) to\n   absent (blood −0.01, ns).") +
  annotate("text", x=0.5, y=2.9, hjust=0, size=2.05, lineheight=1.1, label=
    "• Cohort-level, not within-person.\n• Gene-level (UCSC_RefGene),\n   not promoter-resolved.") +
  annotate("text", x=0.5, y=0.7, hjust=0, size=2.05, fontface="italic", color="grey30", lineheight=1.1, label=
    "→ a limitation; 4-layer node atlas\n   needs same-subject (OCEAN/CKB).") +
  theme_void(base_family=BASE_FONT) +
  theme(plot.title=element_text(size=8.5,hjust=0,margin=margin(b=4)),
        plot.background=element_rect(fill="white", color=NA),
        panel.background=element_rect(fill="white", color=NA)) +
  labs(title="Why the capstone was not pursued")

fig <- wrap_plots(pa, pb, pc, ncol=3) +
  plot_annotation(tag_levels="a",
    title="Additional file 5 | Cross-layer node concordance is a NO-GO, disclosed as an honest negative",
    theme=theme(plot.background=element_rect(fill="white", color=NA))) &
  theme(plot.tag=element_text(face="bold", size=11, hjust=0, vjust=1), plot.tag.position=c(0.01,1.04),
        plot.title=element_text(size=8.5, hjust=0, margin=margin(t=2,b=4)),
        plot.subtitle=element_text(size=7, color="grey25", margin=margin(b=3)),
        plot.margin=margin(13,6,8,6,"pt"))
save_figure_ckm(fig, "Figure_ED5", width_mm=183, height_mm=82, output_dir=".")
cat("done\n")

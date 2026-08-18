# render_ED2.R — Extended Data 2: cross-ancestry transfer is directional-only, no genome-wide correlation.
suppressMessages({library(ggplot2); library(patchwork); library(scales)})
source("theme_publication_ckm.R")  # bundled alongside this script
set_publication_defaults_ckm()
SD <- "source_data"; rd <- function(p) read.csv(file.path(SD, p), stringsAsFactors=FALSE)
TXN<-"#009E73"; UP<-"#B2182B"; DN<-"#2166AC"; NEU<-"#333333"; GREY<-"grey62"

# a) genome-wide rho: Black-cohort-vs-European null vs European-vs-European ceiling
ga <- rd("ED2a_genomewide.csv"); ga <- ga[order(ga$rho),]
ga$label <- sprintf("%s
(n = %s shared genes)", ga$comparison, format(ga$n_shared, big.mark=",", trim=TRUE))
ga$comparison <- factor(ga$label, levels=ga$label)
pa <- ggplot(ga, aes(rho, comparison, fill=kind)) +
  geom_vline(xintercept=0, linewidth=0.4, color="grey55") +
  geom_col(width=0.62, color="black", linewidth=0.2) +
  geom_text(aes(label=sprintf("%+.2f", rho), hjust=ifelse(rho>=0,-0.25,1.2)), size=2.8) +
  scale_fill_manual(values=c("Black cohort"=GREY, "EUR ceiling"=TXN), guide="none") +
  scale_x_continuous(limits=c(-0.18,0.40), breaks=c(-0.1,0,0.1,0.2,0.3)) +
  labs(x="genome-wide Spearman ρ", y=NULL) +
  theme(axis.text.y=element_text(size=6.9, lineheight=0.95))

# b) canonical genes: biology-anchored directional concordance
cb <- rd("ED2b_canonical.csv"); sm <- rd("ED2_summary.csv")
cb$gene <- factor(cb$gene, levels=cb$gene[order(cb$black_beta)])
pb <- ggplot(cb, aes(black_beta, gene)) +
  geom_vline(xintercept=0, linewidth=0.4, color="grey55") +
  geom_segment(aes(x=0, xend=black_beta, yend=gene), linewidth=0.5, color=NEU) +
  geom_point(aes(color=black_beta<0), size=2.4) +
  scale_color_manual(values=c(`TRUE`=DN, `FALSE`=UP), guide="none") +
  scale_x_continuous(limits=c(-0.8,0.7), breaks=c(-0.5,0,0.5)) +
  annotate("label", x=0.68, y=1.05, hjust=1, vjust=0, size=2.55, fill="white", linewidth=0, lineheight=1.0,
           label=sprintf("n = 12 participants\ndirectional concordance %.0f%%\n(p = 0.009, %d genes)", sm$concordance[1]*100, sm$n[1])) +
  labs(x="reversal β (Black women's muscle)", y=NULL) +
  theme(axis.text.y=element_text(face="italic", size=8))

fig <- wrap_plots(pa, pb, ncol=2) +
  plot_annotation(tag_levels="a",
    theme=theme(plot.background=element_rect(fill="white", color=NA))) &
  theme(plot.tag=element_text(face="bold", size=11, hjust=0, vjust=1), plot.tag.position=c(0.01,1.04),
        plot.margin=margin(13,6,8,6,"pt"))
save_figure_ckm(fig, "Figure_ED2", width_mm=170, height_mm=80, output_dir=".")
cat("done\n")

# render_ED1.R — Extended Data 1: reversibility-score validity.
# (a) per-layer Z distributions (folded magnitude); (b) biology face-validity (signed Z); (c) not a measurability proxy.
suppressMessages({library(ggplot2); library(patchwork); library(scales)})
source("C:/Users/Bert/Downloads/CKM papers/paper 3.eas/analysis/tools/theme_publication_ckm.R")
set_publication_defaults_ckm()
SD <- "source_data"; rd <- function(p) read.csv(file.path(SD, p), stringsAsFactors=FALSE)
PRO<-"#0072B2"; TXN<-"#009E73"; METAB<-"#CC79A7"; METH<-"#D55E00"; UP<-"#B2182B"; DN<-"#2166AC"; NEU<-"#333333"
LC <- c(proteome=PRO, transcriptome=TXN, metabolome=METAB, methylation=METH)

# a) per-layer Z_responsiveness distributions
za <- rd("ED1a_z_dist.csv"); za$layer <- factor(za$layer, levels=c("methylation","metabolome","transcriptome","proteome"))
pa <- ggplot(za, aes(Z, layer, fill=layer)) +
  geom_violin(scale="width", color=NA, alpha=0.7, width=0.85) +
  scale_fill_manual(values=LC, guide="none") +
  scale_x_continuous(limits=c(0,4), breaks=c(0,1,2,3,4)) +
  geom_vline(xintercept=1.96, linetype="dashed", linewidth=0.35, color="grey45") +
  annotate("text", x=1.96, y=4.6, label="Z=1.96", size=2.0, color="grey40", vjust=0) +
  labs(x="responsiveness Z (folded magnitude)", y=NULL, title="Per-layer score distributions",
       subtitle="folded: piles at 0, heavy responder tail")

# b) biology face-validity (signed Z) — score recovers known weight-loss biology
bi <- rd("ED1b_biology.csv"); bi <- bi[bi$gene!="IGFBP1",]
bi$gene <- factor(bi$gene, levels=bi$gene[order(bi$signed_Z)])
pb <- ggplot(bi, aes(signed_Z, gene, color=signed_Z>0)) +
  geom_vline(xintercept=0, linewidth=0.4, color="grey55") +
  geom_segment(aes(x=0, xend=signed_Z, yend=gene), linewidth=0.6) +
  geom_point(size=2.4) +
  scale_color_manual(values=c(`TRUE`=UP, `FALSE`=DN), labels=c("down (toward lean)","up"), name=NULL,
                     breaks=c(FALSE,TRUE)) +
  scale_x_continuous(limits=c(-3,1.2), breaks=c(-2,-1,0,1)) +
  labs(x="signed responsiveness Z", y=NULL, title="Score recovers known biology",
       subtitle="SAA1/SAA2/CRP/LEP ↓, GDF15 ↑") +
  theme(axis.text.y=element_text(face="italic"),
        legend.position="inside", legend.position.inside=c(0.55,0.04), legend.justification=c(0,0),
        legend.key.size=unit(7,"pt"), legend.text=element_text(size=6.3),
        legend.background=element_rect(fill=alpha("white",0.7), color=NA))

# c) score is not a measurability proxy
mc <- rd("ED1c_measurability.csv"); cc <- rd("ED1c_corr.csv")
pc <- ggplot(mc, aes(meas_mean, Z_responsiveness)) +
  geom_point(size=0.4, alpha=0.18, color=TXN) +
  geom_smooth(method="lm", se=FALSE, color="black", linewidth=0.5, formula=y~x) +
  scale_y_continuous(limits=c(0,4)) +
  annotate("label", x=Inf, y=4, hjust=1, vjust=1, size=2.3, fill="white", label.size=0,
           label=sprintf("ρ = %+.2f  (n=%s)", cc$corr[1], format(cc$n[1],big.mark=","))) +
  labs(x="measurability (mean logCPM)", y="responsiveness Z", title="Not a measurability proxy",
       subtitle="independent of detectability")

fig <- wrap_plots(pa, pb, pc, ncol=3) +
  plot_annotation(tag_levels="a",
    title="Additional file 1 | The reversibility score is valid: folded-magnitude distributions, biological face-validity, measurability-independent",
    theme=theme(plot.background=element_rect(fill="white", color=NA))) &
  theme(plot.tag=element_text(face="bold", size=11, hjust=0, vjust=1), plot.tag.position=c(0.01,1.04),
        plot.title=element_text(size=8.3, hjust=0, margin=margin(t=2,b=4)),
        plot.subtitle=element_text(size=7, color="grey25", margin=margin(b=3)),
        plot.margin=margin(13,6,8,6,"pt"))
save_figure_ckm(fig, "Figure_ED1", width_mm=183, height_mm=82, output_dir=".")
cat("done\n")

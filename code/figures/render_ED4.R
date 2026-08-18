# render_ED4.R — Extended Data 4: persistence is a variance-matched measurability/noise floor, NOT memory.
# (a) base_sd vs Z by class; (b) variance-conditional class counts; (c) cross-tissue persistence COLLAPSES on
# base-mean conditioning while the reversible pole holds (the falsification). LOCKED theme; backs F2c/F3/F6e.
suppressMessages({library(ggplot2); library(patchwork); library(scales)})
source("theme_publication_ckm.R")  # bundled alongside this script
set_publication_defaults_ckm()
SD <- "source_data"; rd <- function(p) read.csv(file.path(SD, p), stringsAsFactors=FALSE)
PERS<-"#B2182B"; REV<-"#1B7837"; TXN<-"#009E73"; METH<-"#D55E00"; NEU<-"#333333"

# a) base_sd vs Z_responsiveness by class (transcriptome) — persistent = low-variance, low-responsiveness floor
sa <- rd("ED4a_class_scatter.csv"); sa$class <- factor(sa$class, levels=c("reversible","persistent"))
med <- aggregate(cbind(base_sd_mean,Z_responsiveness)~class, sa, median)
pa <- ggplot(sa, aes(base_sd_mean, Z_responsiveness, color=class)) +
  geom_point(size=0.5, alpha=0.35) +
  geom_point(data=med, size=3, shape=23, aes(fill=class), color="black", stroke=0.5) +
  scale_color_manual(values=c(reversible=REV, persistent=PERS), name=NULL) +
  scale_fill_manual(values=c(reversible=REV, persistent=PERS), guide="none") +
  scale_x_continuous(limits=c(0,0.8)) + scale_y_continuous(limits=c(0,1.2)) +
  annotate("label", x=0.78, y=1.15, hjust=1, vjust=1, size=2.35, fill=NA, linewidth=0, lineheight=0.95,
           label="persistent: base_sd 0.24, Z 0.15\nreversible: base_sd 0.27, Z 0.41\nbase_sd ratio 0.90 (median)") +
  labs(x="baseline SD", y="responsiveness Z") +
  theme(legend.position="inside", legend.position.inside=c(0.98,0.03), legend.justification=c(1,0),
        legend.key.size=unit(7,"pt"), legend.text=element_text(size=7),
        legend.background=element_rect(fill=alpha("white",0.7), color=NA))

# b) variance-conditional class counts (transcriptome + methylation), log scale
cb <- rd("ED4b_class_counts.csv")
cb$class <- factor(cb$class, levels=c("persistent","reversible","intermediate"))
cb$layer <- factor(cb$layer, levels=c("transcriptome","methylation"))
pb <- ggplot(cb, aes(class, n, fill=layer)) +
  geom_col(position=position_dodge(0.75), width=0.7, color="black", linewidth=0.2) +
  geom_text(aes(label=scales::comma(n)), position=position_dodge(0.75), vjust=-0.4, size=2.8) +
  scale_fill_manual(values=c(transcriptome=TXN, methylation=METH), name=NULL) +
  scale_y_log10(limits=c(1,2e6), breaks=c(1e2,1e3,1e4,1e5,1e6), labels=scales::comma) +
  labs(x=NULL, y="marks (log)") +
  theme(axis.text.x=element_text(size=8, angle=20, hjust=1),
        legend.position="inside", legend.position.inside=c(0.02,0.99), legend.justification=c(0,1),
        legend.key.size=unit(7,"pt"), legend.text=element_text(size=7),
        legend.background=element_rect(fill=alpha("white",0.7), color=NA))

# c) THE FALSIFICATION: cross-tissue persistence collapses on conditioning; reversible holds
fc <- rd("ED4c_falsification.csv")
fc$cond <- factor(fc$cond, levels=c("raw","edge-conditioned"))
fc$grp <- paste(fc$pole, fc$pair)
pc <- ggplot(fc, aes(cond, fold, group=grp, color=pole)) +
  geom_hline(yintercept=1, linetype="dashed", linewidth=0.4, color="grey45") +
  geom_line(linewidth=0.5, alpha=0.9) + geom_point(size=1.6) +
  scale_color_manual(values=c(persistent=PERS, reversible=REV), name=NULL) +
  scale_y_continuous(limits=c(0.7,5.3), breaks=c(1,2,3,4,5)) +
  annotate("text", x=2.0, y=0.85, label="chance", size=2.8, color="grey35", fontface="italic", hjust=1) +
  labs(x=NULL, y="cross-tissue overlap (fold)") +
  theme(axis.text.x=element_text(size=8),
        legend.position="inside", legend.position.inside=c(0.40,0.99), legend.justification=c(0,1),
        legend.key.size=unit(7,"pt"), legend.text=element_text(size=7),
        legend.background=element_rect(fill=alpha("white",0.7), color=NA))

fig <- wrap_plots(pa, pb, pc, ncol=3) +
  plot_annotation(tag_levels="a",
    theme=theme(plot.background=element_rect(fill="white", color=NA))) &
  theme(plot.tag=element_text(face="bold", size=11, hjust=0, vjust=1), plot.tag.position=c(0.01,1.04),
        plot.margin=margin(13,6,8,6,"pt"))
save_figure_ckm(fig, "Figure_ED4", width_mm=183, height_mm=82, output_dir=".")
cat("done\n")

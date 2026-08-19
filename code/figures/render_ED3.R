# render_ED3.R — Extended Data 3: the determinant null is POWERED, not under-powered.
# (a) per-panel coefficient forest (k=9) + pooled meta + MDE band; (b) τ sign-flips by panel (dynamic-range);
# (c) causal status NS in every panel. LOCKED theme; backs F3a/b/f.
suppressMessages({library(ggplot2); library(patchwork); library(scales)})
source("theme_publication_ckm.R")  # bundled alongside this script
set_publication_defaults_ckm()
SD <- "source_data"; rd <- function(p) read.csv(file.path(SD, p), stringsAsFactors=FALSE)
TXN<-"#009E73"; UP<-"#D55E00"; DN<-"#2166AC"; NEU<-"#333333"; BAND<-"#E6E6E6"
PMAP <- c(liver_table="liver", adipose_table="adipose (CR)", gse141221_adipose_table="adipose (LCD)",
          gse273902_blood_table="blood", gse157585_muscle_table="muscle (metformin)", gse43471_adipose_table="adipose (diet+exercise)",
          gse83352_muscle_exercise_table="muscle (exercise)", gse106737_liver_bariatric_table="liver (bariatric)", gse77962_adipose_diet_table="adipose (diet)")

# a) per-panel determinant forest
a <- rd("ED3a_perpanel.csv"); mt <- rd("ED3a_meta.csv")
# round 7 9e: open the abbreviated noun label to match "tissue specificity" used everywhere else
a$feat <- gsub("tissue-spec.","tissue spec.",a$feat,fixed=TRUE); mt$feat <- gsub("tissue-spec.","tissue spec.",mt$feat,fixed=TRUE)
lev <- mt$feat[order(abs(mt$pooled_beta))]
a$feat <- factor(a$feat, levels=lev); mt$feat <- factor(mt$feat, levels=lev)
mde <- median(mt$mde)
pa <- ggplot() +
  annotate("rect", xmin=-mde, xmax=mde, ymin=0.4, ymax=length(lev)+0.6, fill=BAND) +
  geom_vline(xintercept=0, linewidth=0.4) +
  geom_jitter(data=a, aes(beta, feat), width=0, height=0.16, size=0.7, alpha=0.5, color="grey45") +
  geom_point(data=mt, aes(pooled_beta, feat), shape=23, size=2.1, fill=TXN, color="black", stroke=0.5) +
  scale_x_continuous(limits=c(-0.12,0.12), breaks=c(-0.1,0,0.1)) +
  labs(x="standardized β (grey = 9 panels, ◆ = pooled)", y=NULL) +
  theme(axis.text.y=element_text(size=8))

# b) tau sign-flip by panel (dynamic-range, not a determinant)
tb <- rd("ED3b_tau.csv"); tb$panel <- factor(tb$panel, levels=tb$panel[order(tb$beta)])
pb <- ggplot(tb, aes(beta, panel, color=beta>0)) +
  geom_vline(xintercept=0, linewidth=0.4) +
  geom_point(size=2) +
  scale_color_manual(values=c(`TRUE`=UP, `FALSE`=DN), guide="none") +
  scale_y_discrete(labels=PMAP) +
  scale_x_continuous(limits=c(-0.12,0.12), breaks=c(-0.1,0,0.1)) +
  labs(x="τ β per panel", y=NULL) +
  theme(axis.text.y=element_text(size=7.5))

# c) causal status NS in every panel
cc <- rd("ED3c_causal.csv"); cc$panel <- factor(cc$panel, levels=cc$panel[order(cc$beta)])
pc <- ggplot(cc, aes(beta, panel)) +
  annotate("rect", xmin=-mde, xmax=mde, ymin=0.4, ymax=nrow(cc)+0.6, fill=BAND) +
  geom_vline(xintercept=0, linewidth=0.4) +
  geom_errorbar(aes(xmin=beta-1.96*se, xmax=beta+1.96*se), orientation="y", width=0, color="grey45", linewidth=0.4) +
  geom_point(size=1.8, color=NEU) +
  scale_y_discrete(labels=PMAP) +
  scale_x_continuous(limits=c(-0.06,0.06), breaks=c(-0.05,0,0.05)) +
  labs(x="causal-status β", y=NULL) +
  theme(axis.text.y=element_text(size=7.5))

fig <- wrap_plots(pa, pb, pc, ncol=3) +
  plot_annotation(tag_levels="a",
    theme=theme(plot.background=element_rect(fill="white", color=NA))) &
  theme(plot.tag=element_text(face="bold", size=11, hjust=0, vjust=1), plot.tag.position=c(0.01,1.04),
        plot.margin=margin(13,6,8,6,"pt"))
save_figure_ckm(fig, "Figure_ED3", width_mm=183, height_mm=88, output_dir=".")
cat("done\n")

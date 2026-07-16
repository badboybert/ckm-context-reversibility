# render_F2.R — Figure 2: Reversibility is organized by shared context; cross-tissue sharing lives on the reversible pole.
# R/ggplot, LOCKED theme. Sources: prep_F2_data.py -> source_data/F2{a,b,c,d}_*.csv.
# Legend note: tissue-organization = KNOWN backdrop (MoTrPAC PMID 38693863; Hinte PMID 39558077);
# the NOVEL contribution is the reversible-pole cross-tissue asymmetry (panels c,d).
suppressMessages({library(ggplot2); library(patchwork); library(scales)})
source("C:/Users/Bert/Downloads/CKM papers/paper 3.eas/analysis/tools/theme_publication_ckm.R")
set_publication_defaults_ckm()
SD <- "source_data"; rd <- function(p) read.csv(file.path(SD, p), stringsAsFactors = FALSE)
TXN<-"#009E73"; METH<-"#D55E00"; NEU<-"#333333"; BAND<-"#ECECEC"
TIS <- c(liver="#8C6D31", adipose="#E1A100", muscle="#B2182B", blood="#2166AC", cross="grey70")

# ---- a) Organized by shared context: within- vs cross-tissue pair correlations ----
pa_d <- rd("F2a_pair_correlations.csv")
pa_d$kind <- factor(pa_d$kind, levels=c("within","cross"),
                    labels=c("within\ntissue","cross\ntissue"))
pa_d$tissue <- factor(pa_d$tissue, levels=names(TIS))
mu <- aggregate(rho~kind, pa_d, mean)
set.seed(7)
pa <- ggplot(pa_d, aes(kind, rho)) +
  geom_hline(yintercept=0, linewidth=0.4, color="grey55") +
  geom_jitter(aes(fill=tissue), width=0.16, height=0, shape=21, color="white",
              stroke=0.25, size=1.7, alpha=0.9) +
  geom_crossbar(data=mu, aes(y=rho, ymin=rho, ymax=rho), width=0.45, linewidth=0.5, color=NEU) +
  scale_fill_manual(values=TIS, breaks=c("liver","adipose","muscle","blood"), name=NULL) +
  scale_y_continuous(limits=c(-0.6,0.85), breaks=c(-0.5,0,0.5)) +
  # annotate("text"), NOT annotate("label", fill="white"): a white label box can erase data —
  # the old "mean +0.020" box, sitting at y=0.40, painted over the cross-tissue pair at rho=0.398
  # (one of the 113 the legend counts). Text has no fill, so it cannot erase a point; the
  # cross label is raised to y=0.58, clear of that column's maximum (rho=0.448).
  annotate("text", x=1, y=0.74, label="mean +0.162", size=2.4) +
  annotate("text", x=2, y=0.58, label="mean +0.020", size=2.4) +
  annotate("text", x=1.5, y=0.84, label="Welch t = 4.4 (descriptive)", size=2.5, fontface="italic") +
  labs(x=NULL, y="Spearman r (panel pair)") +
  theme(legend.position="inside", legend.position.inside=c(0.62,0.04), legend.justification=c(0,0),
        legend.key.size=unit(7,"pt"), legend.text=element_text(size=6.3),
        legend.background=element_rect(fill=alpha("white",0.7), color=NA))

# ---- b) Per-tissue cohesion ----
pb_d <- rd("F2b_cohesion.csv")
pb_d$tissue <- factor(pb_d$tissue, levels=c("liver","adipose","muscle","blood"))
# single-line label: when tilted, a two-line "tissue\n(k=N)" pushes its second line into the
# neighbouring label. One line tilts cleanly.
pb_d$lab <- sprintf("%s (k=%d)", pb_d$tissue, pb_d$n_contexts)
pb_d$lab <- factor(pb_d$lab, levels=pb_d$lab[order(-pb_d$mean_within_r)])
pb <- ggplot(pb_d, aes(lab, mean_within_r, fill=tissue)) +
  geom_col(width=0.66, color="black", linewidth=0.25) +
  geom_hline(yintercept=0, linewidth=0.4) +
  geom_text(aes(label=sprintf("%+.2f", mean_within_r), vjust=ifelse(mean_within_r>=0,-0.5,1.3)), size=2.4) +
  scale_fill_manual(values=TIS, guide="none") +
  scale_y_continuous(limits=c(-0.10,0.50), breaks=c(0,0.2,0.4)) +
  labs(x=NULL, y="mean within-tissue r") +
  # tilt (locked recipe for long x-labels): at 0° "adipose" and "muscle" abutted with zero gap
  # and read as "adiposemuscle". Tilting keeps the uniform 9 pt rather than shrinking the text.
  theme(axis.text.x=element_text(angle=30, hjust=1, lineheight=0.9))

# ---- c) Reversible-pole asymmetry (the contribution) ----
pc_d <- rd("F2c_pole_asymmetry.csv")
pc_d$pole <- factor(pc_d$pole, levels=c("reversible","persistent"))
pc_d$layer <- factor(pc_d$layer, levels=c("transcriptome","methylome"))
med <- aggregate(fold~pole+layer, pc_d, median)
set.seed(3)
pc <- ggplot(pc_d, aes(pole, fold, color=layer)) +
  geom_hline(yintercept=1, linetype="dashed", linewidth=0.4, color="grey40") +
  geom_point(position=position_jitterdodge(jitter.width=0.18, dodge.width=0.55),
             size=1.8, alpha=0.9) +
  geom_crossbar(data=med, aes(ymin=fold, ymax=fold, group=layer),
                position=position_dodge(0.55), width=0.42, linewidth=0.45) +
  scale_color_manual(values=c(transcriptome=TXN, methylome=METH), name=NULL) +
  scale_y_continuous(limits=c(-0.2,5.2), breaks=c(0,1,2,3,4,5)) +
  annotate("text", x=2, y=1.45, label="chance", size=2.2, color="grey35", fontface="italic") +
  annotate("text", x=0.78, y=5.05, label="shared", size=2.3, fontface="italic") +
  labs(x=NULL, y="cross-tissue overlap (fold)") +
  theme(legend.position="inside", legend.position.inside=c(0.42,0.78), legend.justification=c(0,0),
        legend.key.size=unit(7,"pt"), legend.text=element_text(size=6.3),
        legend.background=element_rect(fill=alpha("white",0.7), color=NA))

# ---- d) The universal reversible core (n=100; inflammation/innate-immune) ----
core <- rd("F2d_universal_core.csv")
show <- c("CD83","PLAUR","THBS1","PFKFB3","CCL2","SPP1","SAA1","IL18","ACLY","ANGPTL4")
THEME <- c(CD83="inflammation / innate-immune", PLAUR="inflammation / innate-immune",
           THBS1="inflammation / innate-immune", CCL2="inflammation / innate-immune",
           SPP1="inflammation / innate-immune", SAA1="inflammation / innate-immune",
           IL18="inflammation / innate-immune", PFKFB3="metabolic remodelling",
           ACLY="metabolic remodelling", ANGPTL4="metabolic remodelling")
pd_d <- core[match(show, core$gene),]
pd_d$theme <- THEME[pd_d$gene]
pd_d <- pd_d[order(pd_d$theme, pd_d$n_tissue_pairs, pd_d$gene),]
pd_d$gene <- factor(pd_d$gene, levels=pd_d$gene)
pdd <- ggplot(pd_d, aes(n_tissue_pairs, gene, color=theme)) +
  geom_segment(aes(x=0, xend=n_tissue_pairs, yend=gene), linewidth=0.4) +
  geom_point(size=2.2) +
  scale_color_manual(values=c("inflammation / innate-immune"="#762A83",
                              "metabolic remodelling"="#1B7837"), name=NULL) +
  scale_x_continuous(limits=c(0,3.4), breaks=c(1,3)) +
  labs(x="tissue pairs shared", y=NULL) +
  theme(legend.position="inside", legend.position.inside=c(0.46,0.10), legend.justification=c(0,0),
        legend.key.size=unit(7,"pt"), legend.text=element_text(size=6.3),
        legend.background=element_rect(fill=alpha("white",0.7), color=NA),
        axis.text.y=element_text(face="italic"))

# ---- e) tissue-PRIVATE reversible programs (each tissue's recurrent reversers + biology) ----
pe_d <- rd("F2f_tissue_programs.csv"); pe_d <- pe_d[pe_d$gene != "(none)", ]
pe_d$tissue <- factor(pe_d$tissue, levels=c("liver","adipose","muscle"))
pe_d <- pe_d[order(pe_d$tissue, pe_d$mean_rev_beta), ]
pe_d$gene <- factor(pe_d$gene, levels=pe_d$gene)
pe <- ggplot(pe_d, aes(mean_rev_beta, gene, color=tissue)) +
  geom_vline(xintercept=0, linewidth=0.4, color="grey55") +
  geom_segment(aes(x=0, xend=mean_rev_beta, yend=gene), linewidth=0.5) +
  geom_point(size=2) +
  scale_color_manual(values=TIS, breaks=c("liver","adipose","muscle"), name=NULL) +
  scale_x_continuous(limits=c(-1.85,1.05), breaks=c(-1.5,-1,-0.5,0,0.5,1)) +
  annotate("text", x=1.0, y=3, label="muscle: MSTN only\nblood: none", size=2.0, hjust=1, vjust=0.5, lineheight=0.9, fontface="italic", color="grey35") +
  labs(x="mean reversal β (within tissue)", y=NULL) +
  theme(axis.text.y=element_text(size=9, face="italic"),
        legend.position="inside", legend.position.inside=c(0.03,0.55), legend.justification=c(0,0.5),
        legend.key.size=unit(7,"pt"), legend.text=element_text(size=6.3),
        legend.background=element_rect(fill=alpha("white",0.7), color=NA))

fig <- (pa + pb + pc) / (pdd + pe) + plot_layout(heights=c(1,1.15)) +
  plot_annotation(tag_levels="a") &
  theme(plot.tag = element_text(face="bold", size=11, hjust=0, vjust=1),
        plot.tag.position = c(0.01, 1.02),
        plot.margin = margin(11, 6, 8, 6, "pt"))
save_figure_ckm(fig, "Figure_2", width_mm=183, height_mm=155, output_dir=".")
cat("done\n")

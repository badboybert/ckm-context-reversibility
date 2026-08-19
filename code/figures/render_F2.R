# render_F2.R — Figure 2 (ADOPTED 2026-07-21): variance partition promoted to panels a/b; context panels c–f.
# The permutation-based variance partition (2x2 tissue×intervention cell means + unique-R² bars) leads the
# figure; the descriptive within-vs-cross Welch-t jitter moved to Figure S9. Source: prep_F2 + variance_partition.tsv.

suppressMessages({library(ggplot2); library(patchwork); library(scales)})
source("theme_publication_ckm.R")
set_publication_defaults_ckm()
SD <- "source_data"; rd <- function(p) read.csv(file.path(SD, p), stringsAsFactors = FALSE)
RES <- "../../results"
TXN<-"#009E73"; METH<-"#D55E00"; NEU<-"#333333"
TIS <- c(liver="#8C6D31", adipose="#E1A100", muscle="#B2182B", blood="#2166AC", cross="grey70")

# ============ PROMOTED variance partition (from S9) ============
v <- read.delim(file.path(RES, "variance_partition.tsv"), stringsAsFactors = FALSE)
.sup <- c("0"="⁰","1"="¹","2"="²","3"="³","4"="⁴","5"="⁵","6"="⁶","7"="⁷","8"="⁸","9"="⁹")
fmtP <- function(p){ if(is.na(p)) return("NA"); if(p>=0.01) return(sprintf("%.2f (NS)",p))
  e<-floor(log10(p)); m<-p/10^e; ms<-if(abs(m-round(m))<0.05) sprintf("%.0f",round(m)) else sprintf("%.1f",m)
  ds<-paste0(.sup[strsplit(as.character(abs(e)),"")[[1]]],collapse=""); paste0(ms,"×10⁻",ds) }
fmtR2 <- function(r) if(r<0.001) sprintf("%.4f",r) else sprintf("%.3f",r)

# --- a1: 2x2 tissue × intervention-family cell means ---
cval <- function(t) v$beta_ols[v$term==t]; cn <- function(t) v$n_same[v$term==t]
cm <- data.frame(
  tissue = factor(c("different","different","same","same"), levels=c("different","same")),
  family = factor(c("different","same","different","same"), levels=c("different","same")),
  rho = c(cval("cell_diffT_diffF_mean"),cval("cell_diffT_sameF_mean"),cval("cell_sameT_diffF_mean"),cval("cell_sameT_sameF_mean")),
  n   = c(cn("cell_diffT_diffF_mean"),cn("cell_diffT_sameF_mean"),cn("cell_sameT_diffF_mean"),cn("cell_sameT_sameF_mean")))
write.csv(cm, file.path(SD,"F2a_cell_means.csv"), row.names=FALSE)
pA <- ggplot(cm, aes(family, tissue, fill=rho)) +
  geom_tile(color="white", linewidth=1.4) +
  geom_text(aes(label=sprintf("ρ %+.2f\nn %d", rho, n), color=rho>0.20), size=2.5, lineheight=0.95, show.legend=FALSE) +
  scale_color_manual(values=c(`TRUE`="white",`FALSE`="grey10")) +
  scale_fill_gradient(low="white", high="#1B7837", limits=c(0,0.36), breaks=c(0,0.1,0.2,0.3), name="mean ρ") +
  labs(x="intervention family", y="tissue", subtitle="Between-context agreement (2×2)") +
  theme(legend.position="right", legend.key.width=unit(6,"pt"), legend.key.height=unit(11,"pt"),
        legend.title=element_text(size=7.5), legend.text=element_text(size=7),
        plot.subtitle=element_text(size=8.2, face="bold"))

# --- a2: unique variance explained (jackknife + permutation P; platform null) ---
f <- v[v$term %in% c("same_tissue","same_intervention_family","same_platform"),]
disp <- c(same_tissue="shared tissue", same_intervention_family="shared intervention family", same_platform="shared assay platform")
f$disp <- factor(disp[f$term], levels=c("shared assay platform","shared tissue","shared intervention family"))
f$grp  <- ifelse(f$term=="same_platform","assay platform (control)","biological context")
f$lab  <- sprintf("R² = %s · P = %s", vapply(f$unique_R2, fmtR2, ""), vapply(f$perm_p_unique, fmtP, ""))
write.csv(f[,c("term","unique_R2","perm_p_unique","unique_R2_jack_min","unique_R2_jack_max")], file.path(SD,"F2b_variance_partition.csv"), row.names=FALSE)
pB <- ggplot(f, aes(unique_R2, disp, fill=grp)) +
  geom_col(width=0.58, color="black", linewidth=0.3) +
  geom_errorbar(aes(xmin=unique_R2_jack_min, xmax=unique_R2_jack_max), orientation="y", width=0.20, linewidth=0.4, color="grey25") +
  geom_text(aes(label=lab), x=0.004, hjust=0, position=position_nudge(y=0.42), size=2.3, color="grey15") +
  annotate("text", x=0.02, y=1, label="no detectable unique assay-platform contribution", hjust=0, size=2.3, fontface="italic", color="grey35") +
  scale_fill_manual(values=c("biological context"="#1B7837","assay platform (control)"="grey72"), guide="none") +
  scale_x_continuous(limits=c(0,0.46), breaks=c(0,0.1,0.2,0.3), expand=expansion(mult=c(0,0.02))) +
  labs(x="unique variance explained (unique R², jackknife range)", y=NULL,
       subtitle="Tissue and intervention co-govern") +
  theme(axis.text.y=element_text(size=8.6), plot.subtitle=element_text(size=8.2, face="bold"))

# ============ retained context panels (current Fig 2 b–e) ============
pb_d <- rd("F2c_cohesion.csv"); pb_d$tissue <- factor(pb_d$tissue, levels=c("liver","adipose","muscle","blood"))
pb_d$lab <- sprintf("%s (k=%d)", pb_d$tissue, pb_d$n_contexts); pb_d$lab <- factor(pb_d$lab, levels=pb_d$lab[order(-pb_d$mean_within_r)])
pC <- ggplot(pb_d, aes(lab, mean_within_r, fill=tissue)) +
  geom_col(width=0.66, color="black", linewidth=0.25) + geom_hline(yintercept=0, linewidth=0.4) +
  geom_text(aes(label=sprintf("%+.2f", mean_within_r), vjust=ifelse(mean_within_r>=0,-0.5,1.3)), size=2.8) +
  scale_fill_manual(values=TIS, guide="none") + scale_y_continuous(limits=c(-0.10,0.50), breaks=c(0,0.2,0.4)) +
  labs(x=NULL, y="mean within-tissue ρ") + theme(axis.text.x=element_text(angle=30, hjust=1, lineheight=0.9))

pc_d <- rd("F2d_pole_asymmetry.csv"); pc_d$pole <- factor(pc_d$pole, levels=c("reversible","persistent")); pc_d$layer <- factor(pc_d$layer, levels=c("transcriptome","methylome"))
med <- aggregate(fold~pole+layer, pc_d, median); set.seed(3)
pD <- ggplot(pc_d, aes(pole, fold, color=layer)) +
  geom_hline(yintercept=1, linetype="dashed", linewidth=0.4, color="grey40") +
  geom_point(position=position_jitterdodge(jitter.width=0.18, dodge.width=0.55), size=1.8, alpha=0.9) +
  geom_crossbar(data=med, aes(ymin=fold, ymax=fold, group=layer), position=position_dodge(0.55), width=0.42, linewidth=0.45) +
  scale_color_manual(values=c(transcriptome=TXN, methylome=METH), name=NULL) + scale_y_continuous(limits=c(-0.2,5.2), breaks=c(0,1,2,3,4,5)) +
  annotate("text", x=2, y=1.45, label="chance", size=2.8, color="grey35", fontface="italic") +
  annotate("text", x=0.78, y=5.05, label="shared", size=2.8, fontface="italic") +
  labs(x=NULL, y="cross-tissue overlap (fold)") +
  theme(legend.position="inside", legend.position.inside=c(0.42,0.78), legend.justification=c(0,0),
        legend.key.size=unit(7,"pt"), legend.text=element_text(size=6.3), legend.background=element_rect(fill=alpha("white",0.7), color=NA))

core <- rd("F2e_universal_core.csv")
show <- c("CD83","PLAUR","THBS1","PFKFB3","CCL2","SPP1","SAA1","IL18","ACLY","ANGPTL4")
THEME <- c(CD83="inflammation / innate-immune",PLAUR="inflammation / innate-immune",THBS1="inflammation / innate-immune",CCL2="inflammation / innate-immune",SPP1="inflammation / innate-immune",SAA1="inflammation / innate-immune",IL18="inflammation / innate-immune",PFKFB3="metabolic remodelling",ACLY="metabolic remodelling",ANGPTL4="metabolic remodelling")
pd_d <- core[match(show, core$gene),]; pd_d$theme <- THEME[pd_d$gene]; pd_d <- pd_d[order(pd_d$theme, pd_d$n_tissue_pairs, pd_d$gene),]; pd_d$gene <- factor(pd_d$gene, levels=pd_d$gene)
pE <- ggplot(pd_d, aes(n_tissue_pairs, gene, color=theme)) +
  geom_segment(aes(x=0, xend=n_tissue_pairs, yend=gene), linewidth=0.4) + geom_point(size=2.2) +
  scale_color_manual(values=c("inflammation / innate-immune"="#762A83","metabolic remodelling"="#1B7837"), name=NULL) +
  scale_x_continuous(limits=c(0,3.4), breaks=c(1,3)) + labs(x="tissue pairs shared", y=NULL) +
  theme(legend.position="inside", legend.position.inside=c(0.40,0.10), legend.justification=c(0,0),
        legend.key.size=unit(7,"pt"), legend.text=element_text(size=6.3), legend.background=element_rect(fill=alpha("white",0.7), color=NA), axis.text.y=element_text(face="italic"))

pe_d <- rd("F2f_tissue_programs.csv"); pe_d <- pe_d[pe_d$gene != "(none)", ]; pe_d$tissue <- factor(pe_d$tissue, levels=c("liver","adipose","muscle")); pe_d <- pe_d[order(pe_d$tissue, pe_d$mean_rev_beta), ]; pe_d$gene <- factor(pe_d$gene, levels=pe_d$gene)
pF <- ggplot(pe_d, aes(mean_rev_beta, gene, color=tissue)) +
  geom_vline(xintercept=0, linewidth=0.4, color="grey55") + geom_segment(aes(x=0, xend=mean_rev_beta, yend=gene), linewidth=0.5) + geom_point(size=2) +
  scale_color_manual(values=TIS, breaks=c("liver","adipose","muscle"), name=NULL) + scale_x_continuous(limits=c(-1.85,1.05), breaks=c(-1.5,-1,-0.5,0,0.5,1)) +
  labs(x="mean reversal β (within tissue)", y=NULL) +
  theme(axis.text.y=element_text(size=9, face="italic"), legend.position="inside", legend.position.inside=c(0.03,0.55), legend.justification=c(0,0.5),
        legend.key.size=unit(7,"pt"), legend.text=element_text(size=6.3), legend.background=element_rect(fill=alpha("white",0.7), color=NA))

# ============ layout: full-width promoted partition on top, context below ============
design <- "AAABBBBB\nCCCCDDDD\nEEEEFFFF"
fig <- pA + pB + pC + pD + pE + pF +
  plot_layout(design=design, heights=c(1.0,1.0,1.12)) +
  plot_annotation(tag_levels="a") &
  theme(plot.tag=element_text(face="bold", size=11, hjust=0, vjust=1), plot.tag.position=c(0.01,1.03),
        plot.margin=margin(10,6,7,6,"pt"))
save_figure_ckm(fig, "Figure_2", width_mm=183, height_mm=200, output_dir=".")
cat("done F2\n")

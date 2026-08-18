# render_F5.R — Figure 5: Reversibility tracks the physiological change, not the drug; uncoupled from clinical benefit.
# R/ggplot, LOCKED theme. Source: prep_F5_data.py -> source_data/F5*.
suppressMessages({library(ggplot2); library(patchwork); library(scales)})
source("theme_publication_ckm.R")  # bundled alongside this script
set_publication_defaults_ckm()
SD <- "source_data"; rd <- function(p) read.csv(file.path(SD, p), stringsAsFactors = FALSE)
NEU<-"#333333"; UP<-"#B2182B"; DN<-"#2166AC"
CLS <- c("surgery"="#8C6D31","diet"="#E1A100","drug (biguanide)"="#6A51A3","drug (SGLT2i)"="#2166AC")

# ---- a) Intervention acts through physiology: systemic vs tissue-specific ----
pa_d <- rd("F5a_intervention_systemic.csv")
pa_d <- pa_d[order(pa_d$n_systemic_consensus),]
pa_d$lab <- sprintf("%s (%dt)", pa_d$intervention, pa_d$n_tissues)
pa_d$lab <- factor(pa_d$lab, levels=pa_d$lab)
pa_d$fill <- ifelse(pa_d$n_systemic_consensus>0,"systemic","tissue-specific (0)")
pa <- ggplot(pa_d, aes(lab, n_systemic_consensus, fill=fill)) +
  geom_col(width=0.66, color="black", linewidth=0.25) +
  geom_text(aes(label=n_systemic_consensus), vjust=-0.5, size=2.8) +
  scale_fill_manual(values=c("systemic"="#1B7837","tissue-specific (0)"="grey75"), name=NULL) +
  scale_y_continuous(limits=c(0,72), breaks=c(0,20,40,60)) +
  labs(x=NULL, y="cross-tissue consensus genes") +
  theme(axis.text.x=element_text(size=9, angle=30, hjust=1),
        legend.position="inside", legend.position.inside=c(0.03,0.98), legend.justification=c(0,1),
        legend.key.size=unit(7,"pt"), legend.text=element_text(size=6.3),
        legend.background=element_rect(fill=alpha("white",0.7), color=NA))

# ---- b) Semaglutide tracks weight loss, not drug class ----
pb_d <- rd("F5b_semaglutide_corr.csv"); plc <- rd("F5b_platform_control.csv")
CLAB <- c("BBS"="BBS (surgery)","DiRECT"="DiRECT (diet)","metformin"="metformin (drug)",
          "MS bariatric"="MS (surgery)","empagliflozin"="empagliflozin")
pb_d$ylab <- CLAB[pb_d$comparator]
pb_d$ylab <- factor(pb_d$ylab, levels=pb_d$ylab[order(pb_d$rho)])
pb_d$class <- factor(pb_d$class, levels=names(CLS))
pb <- ggplot(pb_d, aes(rho, ylab)) +
  geom_vline(xintercept=0, linewidth=0.4, color="grey55") +
  geom_segment(aes(x=0, xend=rho, yend=ylab, color=class), linewidth=0.6) +
  geom_point(aes(color=class), size=2.4) +
  # left-align the "ns" and clear the marker: at nudge_x=0.024 with the default hjust=0.5 the label's
  # leading glyph sat on the point and was clipped by it (the same occlusion class as Fig 2a).
  geom_text(aes(label=ifelse(sig=="True"|sig==TRUE,"","ns")), nudge_x=0.030, hjust=0, size=2.8, color="grey45") +
  scale_color_manual(values=CLS, guide="none") +
  scale_x_continuous(limits=c(-0.02,0.42), breaks=c(0,0.1,0.2,0.3,0.4)) +
  # No in-panel argument: the previous 3-line label was anchored inside the data region and
  # covered the MS (surgery) marker outright, and it stated "not a platform artifact" as an
  # absolute. The cross-platform control (+0.146, k = 11 vs +0.128, k = 8) is in the legend.
  labs(x="Spearman vs semaglutide", y=NULL) +
  theme(axis.text.y=element_text(size=9))

# ---- c) Layered plasma core: shared IGF axis + mechanism-specific overlays ----
pc_d <- rd("F5c_layered_core.csv")
pc_d$marker <- factor(pc_d$marker, levels=rev(c("IGFBP1","IGFBP2","IL6","LEP","TFRC")))
ivl <- c("surgery\n(BBS)","diet\n(DiRECT)","sema-\nglutide","metformin","empagliflozin\n(EMPEROR)")
xlab <- c("surgery\n(BBS)"="surgery","diet\n(DiRECT)"="diet","sema-\nglutide"="semaglutide",
          "metformin"="metformin","empagliflozin\n(EMPEROR)"="empagliflozin")
pc_d$intervention <- factor(pc_d$intervention, levels=ivl)
pc <- ggplot(pc_d, aes(intervention, marker, fill=rev_beta)) +
  geom_tile(color="white", linewidth=0.6) +
  geom_text(aes(label=ifelse(is.na(rev_beta),"·",sprintf("%+.1f",rev_beta))), size=2.5,
            color=ifelse(!is.na(pc_d$rev_beta) & abs(pc_d$rev_beta)>0.55,"white","black")) +
  scale_fill_gradient2(low=DN, mid="white", high=UP, midpoint=0, limits=c(-1,1), oob=scales::squish,
                       na.value="grey88", name=NULL, breaks=c(-1,0,1), labels=c("down","0","up")) +
  scale_x_discrete(labels=xlab) +
  labs(x=NULL, y=NULL) +
  theme(axis.text.x=element_text(size=9, angle=35, hjust=1), axis.text.y=element_text(face="italic"),
        legend.position="right", legend.key.width=unit(6,"pt"), legend.key.height=unit(10,"pt"),
        legend.title=element_text(size=6.3), legend.text=element_text(size=6.0))

# ---- d) Molecular reversal != clinical benefit: difference estimate with uncertainty (not bars) ----
ds <- rd("F5d_dissociation_stats.csv"); lo <- rd("F5d_loso.csv")
pdd <- ggplot() +
  geom_vline(xintercept=0, linetype="dashed", linewidth=0.4, color="grey55") +
  # leave-one-subject-out difference estimates (each drops one patient) — leverage spread, all positive
  geom_jitter(data=lo, aes(diff_nor_minus_rem, 1), width=0, height=0.11, size=1.3, alpha=0.45, color="#5AAE61") +
  # bootstrap 95% CI (crosses 0) + the point difference
  geom_errorbarh(aes(y=1, xmin=ds$boot_ci_lo, xmax=ds$boot_ci_hi), height=0.14, linewidth=0.7, color="grey30") +
  geom_point(aes(ds$point_diff, 1), size=3.1, shape=18, color="grey15") +
  scale_x_continuous(limits=c(-0.62,0.82), breaks=c(-0.4,0,0.4)) +
  scale_y_continuous(limits=c(0.5,1.7)) +
  annotate("text", x=0.1, y=1.46, size=2.7, label=sprintf("Δ %+.2f  [%.2f, %+.2f]", ds$point_diff, ds$boot_ci_lo, ds$boot_ci_hi)) +
  annotate("text", x=-0.58, y=1.28, hjust=0, size=2.5, fontface="italic", color="grey35",
           label="CI crosses 0 → ns") +
  annotate("text", x=-0.58, y=0.62, hjust=0, size=2.5, fontface="italic", color="grey35",
           label="green: leave-one-subject-out") +
  labs(x="Δ median |reversal| (non-remitter − remitter)", y=NULL) +
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.line.y=element_blank())
# ---- e) intervention-pair reversal-correlation matrix (no drug-class signature) ----
pe_d <- rd("F5e_intervention_matrix.csv")
ivord <- c("BBS (surg)","DiRECT (diet)","MS (surg)","sema","metformin","empagliflozin")
pe_d$ra <- match(pe_d$a, ivord); pe_d$rb <- match(pe_d$b, ivord)
# lower triangle only (one cell per unordered pair) — frees the other half for larger labels
pe_d$show <- ifelse(pe_d$ra > pe_d$rb, pe_d$rho, NA)
pe_d$a <- factor(pe_d$a, levels=ivord); pe_d$b <- factor(pe_d$b, levels=rev(ivord))
DLAB <- c("BBS (surg)"="By-Band-Sleeve","DiRECT (diet)"="DiRECT","MS (surg)"="MS bariatric",
          "sema"="semaglutide","metformin"="metformin","empagliflozin"="empagliflozin")
pe <- ggplot(pe_d, aes(a, b, fill=show)) +
  geom_tile(aes(alpha=is.na(show)), color="white", linewidth=0.6) +
  geom_text(aes(label=ifelse(is.na(show),"",sprintf("%+.2f",show))), size=2.75,
            color=ifelse(!is.na(pe_d$show) & abs(pe_d$show)>0.22,"white","black")) +
  scale_fill_gradient2(low=DN, mid="white", high=UP, midpoint=0, limits=c(-0.35,0.35), oob=scales::squish,
                       na.value="grey96", name="ρ", breaks=c(-0.3,0,0.3)) +
  scale_alpha_manual(values=c(`FALSE`=1,`TRUE`=0), guide="none") +
  scale_x_discrete(labels=DLAB) + scale_y_discrete(labels=DLAB) +
  annotate("text", x=1.35, y=2.5, hjust=0, size=2.5, fontface="italic", color="grey35", lineheight=0.95,
           label="the four weight-loss\npanels cohere;\nempagliflozin\nstands apart") +
  labs(x=NULL, y=NULL) +
  theme(axis.text.x=element_text(size=9, angle=35, hjust=1), axis.text.y=element_text(size=9),
        legend.key.width=unit(5,"pt"), legend.key.height=unit(9,"pt"),
        legend.title=element_text(size=9), legend.text=element_text(size=7.5))

# 2x3 grid (a,b,c / d,e-wide): a..e reading order = RESULTS callouts; panel e spans 2 cols so the 6x6 matrix is legible.
fig <- wrap_plots(pa, pb, pc, pdd, pe, design="ABC\nDEE") +
  plot_annotation(tag_levels="a") &
  theme(plot.tag = element_text(face="bold", size=11, hjust=0, vjust=1),
        plot.tag.position = c(0.01, 1.02),
        plot.margin = margin(11, 6, 8, 6, "pt"))
save_figure_ckm(fig, "Figure_5", width_mm=183, height_mm=150, output_dir=".")
cat("done\n")

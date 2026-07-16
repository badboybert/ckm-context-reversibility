# render_F5.R — Figure 5: Reversibility tracks the physiological change, not the drug; uncoupled from clinical benefit.
# R/ggplot, LOCKED theme. Source: prep_F5_data.py -> source_data/F5*.
suppressMessages({library(ggplot2); library(patchwork); library(scales)})
source("C:/Users/Bert/Downloads/CKM papers/paper 3.eas/analysis/tools/theme_publication_ckm.R")
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
  geom_text(aes(label=n_systemic_consensus), vjust=-0.5, size=2.6) +
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
          "MS bariatric"="MS (surgery)","SGLT2i"="empagliflozin")
pb_d$ylab <- CLAB[pb_d$comparator]
pb_d$ylab <- factor(pb_d$ylab, levels=pb_d$ylab[order(pb_d$rho)])
pb_d$class <- factor(pb_d$class, levels=names(CLS))
pb <- ggplot(pb_d, aes(rho, ylab)) +
  geom_vline(xintercept=0, linewidth=0.4, color="grey55") +
  geom_segment(aes(x=0, xend=rho, yend=ylab, color=class), linewidth=0.6) +
  geom_point(aes(color=class), size=2.4) +
  # left-align the "ns" and clear the marker: at nudge_x=0.024 with the default hjust=0.5 the label's
  # leading glyph sat on the point and was clipped by it (the same occlusion class as Fig 2a).
  geom_text(aes(label=ifelse(sig=="True"|sig==TRUE,"","ns")), nudge_x=0.030, hjust=0, size=2.0, color="grey45") +
  scale_color_manual(values=CLS, guide="none") +
  scale_x_continuous(limits=c(-0.02,0.42), breaks=c(0,0.1,0.2,0.3,0.4)) +
  annotate("label", x=0.42, y=1.55, hjust=1, vjust=0.5, size=2.05, fill="white", linewidth=0, lineheight=0.95, color="grey25",
           label="cross-platform +0.15 ≥\nsame +0.13\n(not a platform artifact)") +
  labs(x="Spearman vs semaglutide", y=NULL) +
  theme(axis.text.y=element_text(size=9))

# ---- c) Layered plasma core: shared IGF axis + mechanism-specific overlays ----
pc_d <- rd("F5c_layered_core.csv")
pc_d$marker <- factor(pc_d$marker, levels=rev(c("IGFBP1","IGFBP2","IL6","LEP","TFRC")))
ivl <- c("surgery\n(BBS)","diet\n(DiRECT)","sema-\nglutide","metformin","SGLT2i\n(EMPEROR)")
xlab <- c("surgery\n(BBS)"="surgery","diet\n(DiRECT)"="diet","sema-\nglutide"="sema",
          "metformin"="metformin","SGLT2i\n(EMPEROR)"="empa")
pc_d$intervention <- factor(pc_d$intervention, levels=ivl)
pc <- ggplot(pc_d, aes(intervention, marker, fill=rev_beta)) +
  geom_tile(color="white", linewidth=0.6) +
  geom_text(aes(label=ifelse(is.na(rev_beta),"·",sprintf("%+.1f",rev_beta))), size=2.1,
            color=ifelse(!is.na(pc_d$rev_beta) & abs(pc_d$rev_beta)>0.55,"white","black")) +
  scale_fill_gradient2(low=DN, mid="white", high=UP, midpoint=0, limits=c(-1,1), oob=scales::squish,
                       na.value="grey88", name=NULL, breaks=c(-1,0,1), labels=c("down","0","up")) +
  scale_x_discrete(labels=xlab) +
  labs(x=NULL, y=NULL) +
  theme(axis.text.x=element_text(size=9, angle=35, hjust=1), axis.text.y=element_text(face="italic"),
        legend.position="right", legend.key.width=unit(6,"pt"), legend.key.height=unit(10,"pt"),
        legend.title=element_text(size=6.3), legend.text=element_text(size=6.0))

# ---- d) Molecular reversal != clinical benefit (caveated) ----
md <- rd("F5d_dissociation_magnitude.csv"); ds <- rd("F5d_dissociation_stats.csv")
md$label <- factor(md$label, levels=c("remitter","non-remitter"))
pdd <- ggplot(md, aes(label, median_abs_sig05, fill=label)) +
  geom_col(width=0.6, color="black", linewidth=0.25) +
  geom_text(aes(label=sprintf("%.2f", median_abs_sig05)), vjust=-0.5, size=2.5) +
  geom_text(aes(label=sprintf("n=%d", n_pairs)), y=0.04, size=2.1, color="white") +
  scale_fill_manual(values=c("remitter"="#9970AB","non-remitter"="#5AAE61"), guide="none") +
  scale_x_discrete(labels=c("remitter"="remitter","non-remitter"="non-\nremitter")) +
  scale_y_continuous(limits=c(0,1.02), breaks=c(0,0.25,0.5,0.75,1.0)) +
  annotate("label", x=1.5, y=0.97, size=2.05, fill="white", linewidth=0, lineheight=0.95, color="grey20",
           label=sprintf("resembles generic surgery:\nnon-rem %+.2f vs rem %+.2f\nΔ CI [%+.2f,%+.2f] crosses 0",
                         ds$resembles_surgery_no_remission, ds$resembles_surgery_remission, ds$boot_ci_lo, ds$boot_ci_hi)) +
  labs(x=NULL, y="median |reversal| (sig genes)")
# ---- e) intervention-pair reversal-correlation matrix (no drug-class signature) ----
pe_d <- rd("F5e_intervention_matrix.csv")
ivord <- c("BBS (surg)","DiRECT (diet)","MS (surg)","sema","metformin","SGLT2i")
pe_d$a <- factor(pe_d$a, levels=ivord); pe_d$b <- factor(pe_d$b, levels=rev(ivord))
pe_d$show <- ifelse(as.character(pe_d$a)==as.character(pe_d$b), NA, pe_d$rho)
pe <- ggplot(pe_d, aes(a, b, fill=show)) +
  geom_tile(color="white", linewidth=0.6) +
  geom_text(aes(label=ifelse(is.na(show),"",sprintf("%+.2f",show))), size=1.95,
            color=ifelse(!is.na(pe_d$show) & abs(pe_d$show)>0.22,"white","black")) +
  scale_fill_gradient2(low=DN, mid="white", high=UP, midpoint=0, limits=c(-0.35,0.35), oob=scales::squish,
                       na.value="grey92", name="ρ", breaks=c(-0.3,0,0.3)) +
  scale_x_discrete(labels=c("BBS (surg)"="BBS","DiRECT (diet)"="DiRECT","MS (surg)"="MS","sema"="sema","metformin"="metf","SGLT2i"="empa")) +
  scale_y_discrete(labels=c("BBS (surg)"="BBS","DiRECT (diet)"="DiRECT","MS (surg)"="MS","sema"="sema","metformin"="metf","SGLT2i"="empa")) +
  annotate("rect", xmin=0.5, xmax=4.5, ymin=2.5, ymax=6.5, fill=NA, color="black", linewidth=0.5) +
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

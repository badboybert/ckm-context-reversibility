# render_S8.R — Supplementary Fig. S8 (provisional): SESOI / TOST equivalence forest (workstream B5; review §3.3/§9-123).
# Per-feature pooled determinant coefficient with 90% CI vs the SESOI bands (+/-0.05 primary, +/-0.03 secondary).
# Equivalent = 90% CI inside the band; inconclusive = CI spills over. Source: results/determinant_tost.tsv.
suppressMessages({library(ggplot2); library(scales)})
source("C:/Users/Bert/Downloads/CKM papers/paper 3.eas/analysis/tools/theme_publication_ckm.R")
set_publication_defaults_ckm()
FIGDIR <- "C:/Users/Bert/Downloads/CKM papers/weight-loss.omics/signature-pivot/manuscript/figures"
RES    <- "C:/Users/Bert/Downloads/CKM papers/weight-loss.omics/signature-pivot/results"
SD     <- file.path(FIGDIR, "source_data")

t <- read.delim(file.path(RES, "determinant_tost.tsv"), stringsAsFactors = FALSE)
fmap <- c(causal_nonEGFR="genetic causal status", arch_str="cis-QTL architecture (strength)",
          has_arch="has cis-QTL architecture", n_drug_log="druggability (log drug count)",
          is_enzyme="enzyme", loeuf="LOEUF constraint", tau_miss="τ missing",
          n_gwas_log="GWAS burden (log)", is_secreted="secreted protein", is_membrane="membrane protein",
          arch_nsig="cis-QTL architecture (n signals)", loeuf_miss="LOEUF missing", tau="tissue-specificity τ")
t$flabel <- fmap[t$feature]
t$equiv05 <- t$equiv_05 %in% c("True","TRUE",TRUE)
t$equiv03 <- t$equiv_03 %in% c("True","TRUE",TRUE)
t$excl0   <- (t$ci90_lo > 0) | (t$ci90_hi < 0)
t$status  <- ifelse(!t$equiv05, "inconclusive",
              ifelse(t$equiv03, "equivalent (±0.03 SD)", "equivalent (±0.05 SD only)"))
t$status  <- factor(t$status, levels=c("equivalent (±0.03 SD)","equivalent (±0.05 SD only)","inconclusive"))
t <- t[order(t$pooled_beta),]
t$flabel <- factor(t$flabel, levels=t$flabel)     # ascending beta, bottom -> top
statcol <- c("equivalent (±0.03 SD)"="#009E73","equivalent (±0.05 SD only)"="#7FBF9B","inconclusive"="#D55E00")

# ---- source CSV ----
write.csv(t[,c("feature","flabel","pooled_beta","se","ci90_lo","ci90_hi","tost_p_05","equiv05",
               "tost_p_03","equiv03","excl0","status","I2")],
          file.path(SD,"S8_determinant_tost_forest.csv"), row.names=FALSE)

yN <- nlevels(t$flabel)
p <- ggplot(t, aes(pooled_beta, flabel, color=status)) +
  annotate("rect", xmin=-0.05, xmax=0.05, ymin=0.4, ymax=yN+0.6, fill="grey93") +
  annotate("rect", xmin=-0.03, xmax=0.03, ymin=0.4, ymax=yN+0.6, fill="grey84") +
  geom_vline(xintercept=0, linewidth=0.45, color="grey35") +
  geom_vline(xintercept=c(-0.05,0.05), linewidth=0.35, linetype="22", color="grey50") +
  geom_vline(xintercept=c(-0.03,0.03), linewidth=0.35, linetype="12", color="grey50") +
  geom_errorbar(aes(xmin=ci90_lo, xmax=ci90_hi), orientation="y", width=0, linewidth=0.55) +
  geom_point(size=2.3) +
  geom_text(data=subset(t, excl0), aes(x=0.062, label="CI ≠ 0"), hjust=1, size=1.9, color="grey25") +
  scale_color_manual(values=statcol, name="90% CI vs SESOI") +
  scale_x_continuous(limits=c(-0.066,0.066), breaks=c(-0.05,-0.03,0,0.03,0.05),
                     labels=c("−0.05","−0.03","0","+0.03","+0.05")) +
  labs(x="pooled determinant coefficient  (standardized β, 90% CI)", y=NULL,
       title="Determinant coefficients are statistically equivalent to zero within the SESOI (practical-equivalence test)",
       subtitle="SESOI bands: ±0.05 SD (light) and ±0.03 SD (dark). 12/13 features equivalent at ±0.05, 10/13 at ±0.03.") +
  theme(axis.text.y=element_text(size=8.2),
        plot.title=element_text(size=8.3, margin=margin(b=2)),
        plot.subtitle=element_text(size=7, color="grey30", lineheight=1.02, margin=margin(b=5)),
        legend.position="right",
        legend.title=element_text(size=7.5), legend.text=element_text(size=7),
        legend.key.size=unit(9,"pt"))

# callouts for the two carrier/covariate cases (positions tuned for the beta-ordered y)
yid <- setNames(seq_len(yN), levels(t$flabel))
yc  <- yid[["genetic causal status"]]; yt <- yid[["tissue-specificity τ"]]
p <- p +
  geom_curve(data=data.frame(1), aes(x=-0.058,y=yc+2.4,xend=-0.0036,yend=yc), inherit.aes=FALSE,
             curvature=0.2, linewidth=0.3, color="grey45", arrow=arrow(length=unit(3.5,"pt"),type="closed")) +
  geom_label(data=data.frame(1), aes(x=-0.063,y=yc+2.7), inherit.aes=FALSE, hjust=0, vjust=0.5, size=2.0,
             label="genetic causal status:\n90% CI [−0.006,−0.001] excludes 0\nyet inside SESOI = practical equivalence",
             fill="white", label.size=0, lineheight=0.92, color="grey15") +
  geom_curve(data=data.frame(1), aes(x=-0.004,y=yt-0.9,xend=0.0246,yend=yt-0.12), inherit.aes=FALSE,
             curvature=-0.25, linewidth=0.3, color="grey45", arrow=arrow(length=unit(3.5,"pt"),type="closed")) +
  geom_label(data=data.frame(1), aes(x=-0.006,y=yt-0.9), inherit.aes=FALSE, hjust=1, vjust=0.5, size=2.0,
             label="τ inconclusive — CI crosses the SESOI\n(dynamic-range covariate, not a determinant)",
             fill="white", label.size=0, lineheight=0.92, color="#8C3A00")

save_figure_ckm(p, "Figure_S8", width_mm=183, height_mm=105, output_dir=FIGDIR)
cat("done S8\n")

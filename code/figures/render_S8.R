# render_S8.R — Supplementary Fig. S8 (provisional): SESOI / TOST equivalence forest.
# Per-feature pooled determinant coefficient with 90% CI vs the SESOI bands (+/-0.05 primary, +/-0.03 secondary).
# Equivalent = 90% CI inside the band; inconclusive = CI spills over. Source: results/determinant_tost.tsv.
suppressMessages({library(ggplot2); library(scales)})
source("theme_publication_ckm.R")  # bundled alongside this script
set_publication_defaults_ckm()
FIGDIR <- "."
RES    <- "../../results"
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
# Genetic causal status is a rare 0/1 feature (9 transcripts): its standardized beta is
# compressed ~50x and is NOT an effect size, so it is shown as a separate underpowered
# rare-binary exception rather than an "equivalent" dot (raw contrast -0.18 SD; see Table S14).
t$status[t$feature=="causal_nonEGFR"] <- "rare binary (underpowered)"
t$status  <- factor(t$status, levels=c("equivalent (±0.03 SD)","equivalent (±0.05 SD only)","inconclusive","rare binary (underpowered)"))
t <- t[order(t$pooled_beta),]
t$flabel <- factor(t$flabel, levels=t$flabel)     # ascending beta, bottom -> top
statcol <- c("equivalent (±0.03 SD)"="#009E73","equivalent (±0.05 SD only)"="#7FBF9B","inconclusive"="#D55E00","rare binary (underpowered)"="#762A83")

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
  geom_text(data=subset(t, excl0), aes(x=0.062, label="CI ≠ 0"), hjust=1, size=2.8, color="grey25") +
  scale_color_manual(values=statcol, name="90% CI vs SESOI") +
  scale_x_continuous(limits=c(-0.066,0.066), breaks=c(-0.05,-0.03,0,0.03,0.05),
                     labels=c("−0.05","−0.03","0","+0.03","+0.05")) +
  labs(x="pooled determinant coefficient  (standardized β, 90% CI)", y=NULL) +
  theme(axis.text.y=element_text(size=8.2),
        legend.position="right",
        legend.title=element_text(size=7.5), legend.text=element_text(size=7),
        legend.key.size=unit(9,"pt"))

# Round-3: "the equivalence forest is visually overannotated, with long embedded notes and arrows …
# move most annotation text to the legend; keep the plot as points, intervals, SESOI bands, and a
# small rare-binary marker." The two four-line callouts and their curved arrows are removed: both
# cases are ALREADY encoded by the `status` colour scale and its key ("rare binary (underpowered)",
# "inconclusive"), and the Additional file 10 legend carries the compression argument, the raw
# contrast (−0.18 SD, 95% CI −0.37 to +0.01) and τ's inconclusive verdict in full. What remains is a
# single compact tag on the rare-binary row, matching Figure 3a's treatment.
yid <- setNames(seq_len(yN), levels(t$flabel))
yc  <- yid[["genetic causal status"]]
p <- p +
  geom_text(data=data.frame(x=-0.0036, y=yc), aes(x=x, y=y), inherit.aes=FALSE,
            label="rare binary, n = 9", hjust=-0.28, vjust=-1.1, size=2.6, color="#5A2A6E")

save_figure_ckm(p, "Figure_S8", width_mm=183, height_mm=105, output_dir=FIGDIR)
cat("done S8\n")

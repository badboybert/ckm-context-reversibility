# render_S9.R — Supplementary Fig. S9: within- vs cross-tissue pairwise reversal agreement (descriptive).
# Demoted from main Fig 2a on 2026-07-21: the permutation-based variance partition is now main Fig 2a,b;
# this descriptive Welch-t view of the 153 pairwise correlations is retained here as a raw-data companion.
suppressMessages({library(ggplot2); library(scales)})
source("theme_publication_ckm.R")  # bundled alongside this script
set_publication_defaults_ckm()
SD <- "source_data"; rd <- function(p) read.csv(file.path(SD, p), stringsAsFactors = FALSE)
TIS <- c(liver="#8C6D31", adipose="#E1A100", muscle="#B2182B", blood="#2166AC", cross="grey70")
NEU <- "#333333"

d <- rd("S9_pair_correlations.csv")
d$kind <- factor(d$kind, levels=c("within","cross"), labels=c("within\ntissue","cross\ntissue"))
d$tissue <- factor(d$tissue, levels=names(TIS))
mu <- aggregate(rho~kind, d, mean); set.seed(7)
p <- ggplot(d, aes(kind, rho)) +
  geom_hline(yintercept=0, linewidth=0.4, color="grey55") +
  geom_jitter(aes(fill=tissue), width=0.16, height=0, shape=21, color="white", stroke=0.25, size=2.0, alpha=0.9) +
  geom_crossbar(data=mu, aes(y=rho, ymin=rho, ymax=rho), width=0.45, linewidth=0.5, color=NEU) +
  scale_fill_manual(values=TIS, breaks=c("liver","adipose","muscle","blood"), name=NULL) +
  scale_y_continuous(limits=c(-0.6,0.9), breaks=c(-0.5,0,0.5)) +
  annotate("text", x=1, y=0.76, label="mean +0.162", size=3.0) +
  annotate("text", x=2, y=0.58, label="mean +0.020", size=3.0) +
  annotate("text", x=1.5, y=0.88, label="Welch t = 4.4 (descriptive)", size=2.8, fontface="italic") +
  labs(x=NULL, y="Spearman ρ (panel pair)") +
  theme(legend.position="right", legend.key.size=unit(9,"pt"), legend.text=element_text(size=8))
save_figure_ckm(p, "Figure_S9", width_mm=120, height_mm=95, output_dir=".")
cat("done S9\n")

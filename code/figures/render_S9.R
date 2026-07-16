# render_S9.R — Supplementary Fig. S9 (provisional): VARIANCE PARTITION (workstream B4; review §3.5b/§4i).
# (a) unique R2 per factor with jackknife range, Type-II % and panel-permutation P — platform-null made obvious;
# (b) 2x2 tissue x intervention-family cell means (super-additive "both"). Source: results/variance_partition.tsv.
suppressMessages({library(ggplot2); library(patchwork); library(scales)})
source("C:/Users/Bert/Downloads/CKM papers/paper 3.eas/analysis/tools/theme_publication_ckm.R")
set_publication_defaults_ckm()
FIGDIR <- "C:/Users/Bert/Downloads/CKM papers/weight-loss.omics/signature-pivot/manuscript/figures"
RES    <- "C:/Users/Bert/Downloads/CKM papers/weight-loss.omics/signature-pivot/results"
SD     <- file.path(FIGDIR, "source_data")

v <- read.delim(file.path(RES, "variance_partition.tsv"), stringsAsFactors = FALSE)
f <- v[v$term %in% c("same_tissue","same_intervention_family","same_platform"),]
disp <- c(same_tissue="shared tissue", same_intervention_family="shared intervention family",
          same_platform="shared assay platform")
f$disp <- disp[f$term]
f$grp  <- ifelse(f$term=="same_platform","assay platform (control)","biological context")
# order bottom->top by unique R2: platform, tissue, intervention
f$disp <- factor(f$disp, levels=c("shared assay platform","shared tissue","shared intervention family"))
pP <- c(same_tissue="1.5×10⁻⁴", same_intervention_family="5×10⁻⁵", same_platform="0.27 (NS)")
# NOTE: typeII_pct_totalvar is a FRACTION in the TSV (0.134 = 13.4%) -> ×100 for display.
f$lab <- sprintf("R² = %.3f · Type-II %.1f%% · P = %s", f$unique_R2, 100*f$typeII_pct_totalvar, pP[f$term])

# ---- source CSV ----
write.csv(f[,c("term","disp","grp","unique_R2","typeII_pct_totalvar","perm_p_unique",
               "unique_R2_jack_min","unique_R2_jack_max")],
          file.path(SD,"S9a_variance_partition.csv"), row.names=FALSE)

grpcol <- c("biological context"="#1B7837","assay platform (control)"="grey72")
pa <- ggplot(f, aes(unique_R2, disp, fill=grp)) +
  geom_col(width=0.58, color="black", linewidth=0.3) +
  geom_errorbar(aes(xmin=unique_R2_jack_min, xmax=unique_R2_jack_max), orientation="y",
                width=0.20, linewidth=0.4, color="grey25") +
  geom_text(aes(label=lab), x=0.004, hjust=0, position=position_nudge(y=0.42),
            size=2.05, color="grey15") +
  annotate("text", x=0.035, y=1, label="platform ≈ 0  →  not a batch/assay artifact", hjust=0,
           size=2.05, fontface="italic", color="#7A1E1E") +
  scale_fill_manual(values=grpcol, guide="none") +
  scale_x_continuous(limits=c(0,0.44), breaks=c(0,0.1,0.2,0.3), expand=expansion(mult=c(0,0.02))) +
  labs(x="unique variance explained  (unique R², jackknife range)", y=NULL) +
  theme(axis.text.y=element_text(size=8.4))

# ---- Panel b: 2x2 cell means (stored as term-rows; value in beta_ols, count in n_same) ----
cval <- function(term) v$beta_ols[v$term==term]
cn   <- function(term) v$n_same[v$term==term]
cm <- data.frame(
  tissue = factor(c("different","different","same","same"), levels=c("different","same")),
  family = factor(c("different","same","different","same"), levels=c("different","same")),
  rho = c(cval("cell_diffT_diffF_mean"), cval("cell_diffT_sameF_mean"),
          cval("cell_sameT_diffF_mean"), cval("cell_sameT_sameF_mean")),
  n   = c(cn("cell_diffT_diffF_mean"), cn("cell_diffT_sameF_mean"),
          cn("cell_sameT_diffF_mean"), cn("cell_sameT_sameF_mean")))
write.csv(cm, file.path(SD,"S9b_cell_means.csv"), row.names=FALSE)
pb <- ggplot(cm, aes(family, tissue, fill=rho)) +
  geom_tile(color="white", linewidth=1.4) +
  geom_text(aes(label=sprintf("ρ = %+.3f\nn = %d", rho, n),
                color=rho>0.20), size=2.5, lineheight=0.95, show.legend=FALSE) +
  scale_color_manual(values=c(`TRUE`="white",`FALSE`="grey10")) +
  scale_fill_gradient(low="white", high="#1B7837", limits=c(0,0.36),
                      breaks=c(0,0.1,0.2,0.3), name="mean ρ") +
  labs(x="intervention family", y="tissue") +
  theme(legend.position="right", legend.key.width=unit(6,"pt"), legend.key.height=unit(11,"pt"),
        legend.title=element_text(size=7.5), legend.text=element_text(size=7))

fig <- (pa | pb) + plot_layout(widths=c(1.55,1)) +
  plot_annotation(tag_levels="a") &
  theme(plot.tag=element_text(face="bold", size=11, hjust=0, vjust=1), plot.tag.position=c(0.005,1.03),
        plot.margin=margin(12,6,7,6,"pt"))
save_figure_ckm(fig, "Figure_S9", width_mm=183, height_mm=92, output_dir=FIGDIR)
cat("done S9\n")

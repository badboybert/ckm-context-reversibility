# render_S7.R — Supplementary Fig. S7 (provisional; final number reconciled at repackage)
# DRUG-CLASS RECONCILIATION (review §3.6/§125; workstream B1).
# Two complementary metrics per powered intervention pair, side by side:
#   left  = effect-size Spearman rho over jointly-measured proteins (ref 0)
#   right = power-normalised significant-reverser overlap fold (ref 1.0)
# Colour = metric-agreement flag; faded = within-metric n.s.; orange ring = the
# empagliflozin-specific decoupling cases. Source: results/drug_class_fig5_panel.tsv.
suppressMessages({library(ggplot2); library(patchwork); library(scales)})
source("C:/Users/Bert/Downloads/CKM papers/paper 3.eas/analysis/tools/theme_publication_ckm.R")
set_publication_defaults_ckm()
FIGDIR <- "C:/Users/Bert/Downloads/CKM papers/weight-loss.omics/signature-pivot/manuscript/figures"
RES    <- "C:/Users/Bert/Downloads/CKM papers/weight-loss.omics/signature-pivot/results"
SD     <- file.path(FIGDIR, "source_data")

d <- read.delim(file.path(RES, "drug_class_fig5_panel.tsv"), stringsAsFactors = FALSE)
d$sig  <- d$significant %in% c("True","TRUE",TRUE)
d$pair <- gsub("MS_bariatric","MS", d$pair_label); d$pair <- gsub(" x "," × ", d$pair)

# class ordering (facet rows, top -> bottom) and within-class order by descending rho
clev  <- c("DIET x SURGERY","within-SURGERY","DRUG x SURGERY","DIET x DRUG","within-DRUG")
cdisp <- c("DIET x SURGERY"="diet ×\nsurgery","within-SURGERY"="within\nsurgery",
           "DRUG x SURGERY"="drug ×\nsurgery","DIET x DRUG"="diet ×\ndrug","within-DRUG"="within\ndrug")
d$classf <- factor(cdisp[d$class_pair], levels = cdisp[clev])
rho_by_pair <- setNames(d$value[d$metric=="effect_size_rho"], d$pair[d$metric=="effect_size_rho"])
ord <- unique(d[,c("pair","class_pair")]); ord$rho <- rho_by_pair[ord$pair]
ord$crank <- as.integer(factor(ord$class_pair, levels=clev))
ord <- ord[order(ord$crank, ord$rho),]           # asc rho within class => highest at top
d$pair <- factor(d$pair, levels = ord$pair)

# empagliflozin-specific decoupling set (task): metformin x empa + empa x {sema,BBS,DiRECT}
flagset <- c("metformin × empagliflozin","empagliflozin × semaglutide",
             "empagliflozin × BBS","empagliflozin × DiRECT")
d$flag  <- as.character(d$pair) %in% flagset

d$agreef <- factor(d$agree, levels=c("concordant-positive","concordant-null","discordant"))
agcol <- c("concordant-positive"="#0072B2","concordant-null"="#8C8C8C","discordant"="#C6C6C6")
d$sigf <- factor(ifelse(d$sig,"significant","n.s."), levels=c("significant","n.s."))

A  <- subset(d, metric=="effect_size_rho")
Bd <- subset(d, metric=="overlap_fold")

# ---- write source CSV ----
src <- d[order(d$classf, d$pair, d$metric),
         c("pair","class_pair","metric","value","p","sig","samedir_frac","agree","flag")]
names(src)[names(src)=="sig"] <- "within_metric_significant"
write.csv(src, file.path(SD,"S7_drug_class_reconciliation.csv"), row.names=FALSE)

base_pt <- function(p) p +
  geom_point(aes(fill=agreef, alpha=sigf), shape=21, size=2.7, color="black", stroke=0.5) +
  geom_point(data=function(x) x[x$flag,], shape=21, size=4.4, color="#D55E00", fill=NA, stroke=0.85) +
  scale_fill_manual(values=agcol, name="metric agreement", drop=FALSE) +
  scale_alpha_manual(values=c(significant=1, `n.s.`=0.38), name="within-metric test", drop=FALSE) +
  guides(fill=guide_legend(order=1, override.aes=list(alpha=1)),
         alpha=guide_legend(order=2, override.aes=list(fill="grey30")))

# ---- Panel A: effect-size rho ----
pa <- ggplot(A, aes(value, pair)) +
  geom_vline(xintercept=0, linewidth=0.4, color="grey50", linetype="22") +
  geom_segment(aes(x=0, xend=value, yend=pair), linewidth=0.45, color="grey78")
pa <- base_pt(pa) +
  facet_grid(classf ~ ., scales="free_y", space="free_y") +
  scale_x_continuous(limits=c(-0.10,0.40), breaks=c(0,0.1,0.2,0.3)) +
  geom_label(data=data.frame(classf=factor(cdisp["within-DRUG"],levels=cdisp[clev]),
                             value=0.385, pair=factor("empagliflozin × semaglutide",levels=ord$pair)),
             aes(x=value,y=pair), inherit.aes=FALSE, hjust=1, vjust=0.5, size=2.05,
             label="empagliflozin ≈ 0 vs\nweight-loss panels", fill="white", linewidth=0,
             lineheight=0.92, color="grey20") +
  labs(x="effect-size Spearman ρ", y=NULL) +
  theme(strip.text.y=element_blank(), strip.background=element_blank(),
        axis.text.y=element_text(size=8), legend.position="none")

# ---- Panel B: significant-reverser overlap fold ----
pb <- ggplot(Bd, aes(value, pair)) +
  geom_vline(xintercept=1, linewidth=0.4, color="grey50", linetype="22") +
  geom_segment(aes(x=1, xend=value, yend=pair), linewidth=0.45, color="grey78")
Bd$dircol <- ifelse(Bd$flag, "#B2182B", "grey40")   # red = empagliflozin decoupling (chance-level same-dir)
pb <- base_pt(pb) +
  facet_grid(classf ~ ., scales="free_y", space="free_y") +
  scale_x_continuous(limits=c(0.62,2.02), breaks=c(0.8,1.0,1.2,1.4,1.6,1.8)) +
  geom_text(aes(x=2.0, label=sprintf("%.0f%%",100*samedir_frac)), hjust=1, size=1.9, color=Bd$dircol) +
  labs(x="sig-reverser overlap (fold)", y=NULL) +
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(),
        strip.text.y=element_text(size=7.4, angle=0, lineheight=0.9),
        legend.position="right", legend.box="vertical",
        legend.title=element_text(size=7.5), legend.text=element_text(size=7))

fig <- (pa | pb) + plot_layout(widths=c(1,1.06)) +
  plot_annotation(
    # title removed (duplicates the Additional file 9 legend headline). Subtitle KEPT: it is the
    # panel's only colour key (orange ring / blue) and the legend does not carry that encoding.
    subtitle="15 adequately powered pairs (underpowered excluded by design). Orange ring marks the empagliflozin-specific decoupling; blue = concordant on both metrics.",
    tag_levels="a") &
  theme(plot.tag=element_text(face="bold", size=11, hjust=0, vjust=1), plot.tag.position=c(0.005,1.045),
        plot.title=element_text(size=8.8, hjust=0, margin=margin(t=2,b=2)),
        plot.subtitle=element_text(size=7, color="grey30", margin=margin(b=4)),
        plot.margin=margin(16,6,7,6,"pt"))
save_figure_ckm(fig, "Figure_S7", width_mm=183, height_mm=150, output_dir=FIGDIR)
cat("done S7\n")

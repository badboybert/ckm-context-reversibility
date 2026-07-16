# render_F1.R — Figure 1 (SPINE): the question, the data landscape, a reproducible score.
# R/ggplot, publication-figures-tables skill + locked theme. 4 panels: a schematic, b data landscape,
# c not-a-power-proxy QC, d reproducibility (split-half + independent-cohort replication).
suppressMessages({library(ggplot2); library(patchwork); library(scales)})
source("C:/Users/Bert/Downloads/CKM papers/paper 3.eas/analysis/tools/theme_publication_ckm.R")
set_publication_defaults_ckm()
RES <- "C:/Users/Bert/Downloads/CKM papers/weight-loss.omics/signature-pivot/results"
SD <- "source_data"
PRO<-"#0072B2"; TXN<-"#009E73"; MET<-"#CC79A7"; MEY<-"#D55E00"; NEU<-"#333333"
LAYERCOL <- c(proteome=PRO, transcriptome=TXN, metabolome=MET, methylation=MEY)
arr <- arrow(length=unit(4.5,"pt"), type="closed")

# a) study-design schematic (programmatic)
pa <- ggplot() + xlim(0,10) + ylim(0,10) +
  annotate("rect", xmin=1.0,xmax=9.0, ymin=8.55,ymax=9.85, fill="#E7F1ED", color=TXN, linewidth=0.4) +
  annotate("text", x=5,y=9.2, size=1.4, lineheight=1.0,
           label="40 human intervention cohorts\n4 layers × 4 tissues × 6 interventions") +
  annotate("segment", x=5,xend=5, y=8.55,yend=7.45, arrow=arr, linewidth=0.45) +
  annotate("rect", xmin=2.0,xmax=8.0, ymin=6.35,ymax=7.35, fill="#E7F1ED", color=TXN, linewidth=0.4) +
  annotate("text", x=5,y=6.85, size=2.2, label="per-mark reversibility score") +
  annotate("segment", x=5,xend=5, y=6.35,yend=5.25, arrow=arr, linewidth=0.45) +
  annotate("rect", xmin=0.6,xmax=9.4, ymin=2.15,ymax=5.15, fill="#FBEEE4", color=MEY, linewidth=0.4) +
  annotate("text", x=5,y=4.60, size=2.2, fontface="bold", lineheight=1.0,
           label="Does any mark-intrinsic\nproperty predict reversal?") +
  annotate("text", x=5,y=3.18, size=1.5, lineheight=0.92,
           label="constraint · druggability\ndisease assoc. · cis-architecture\ncausal status") +
  annotate("segment", x=5,xend=5, y=2.15,yend=1.4, arrow=arr, linewidth=0.45) +
  annotate("text", x=5,y=0.9, size=2.2, fontface="italic", color=TXN, label="→ context, not the molecule") +
  theme_void(base_family=BASE_FONT) +
  theme(plot.title=element_text(size=8.5, hjust=0, margin=margin(b=4)))

# b) data landscape (layer x tissue, bubble = paired n)
m <- read.delim(file.path(RES,"panel_manifest_full.tsv"), stringsAsFactors=FALSE)
m <- m[!grepl("EXCLUDED", toupper(m$status)),]
tmap <- c(plasma="plasma/serum", serum="plasma/serum", whole_blood="blood", PBMC="blood",
          SAT_adipose="adipose", skeletal_muscle="muscle", liver="liver", blood="blood")
m$tis <- tmap[m$tissue]; m$np <- pmax(suppressWarnings(as.numeric(m$n_pairs)),0,na.rm=FALSE)
ag <- aggregate(np~layer+tis, m, function(x) sum(x,na.rm=TRUE))
ag$layer <- factor(ag$layer, levels=c("metabolome","methylation","proteome","transcriptome"))
ag$tis   <- factor(ag$tis,   levels=c("plasma/serum","blood","adipose","liver","muscle"))
pb <- ggplot(ag, aes(tis, layer, size=np, color=layer)) +
  geom_point(alpha=0.85) +
  scale_size_area(max_size=9, breaks=c(100,1000,3000), name="paired n") +
  scale_color_manual(values=LAYERCOL, guide="none") +
  scale_x_discrete(labels=c("plasma/serum"="plasma/serum","blood"="blood","adipose"="adipose","liver"="liver","muscle"="muscle")) +
  labs(x=NULL, y=NULL) +
  theme(axis.text.x=element_text(size=8, angle=30, hjust=1),
        legend.position="right", legend.key.size=unit(8,"pt"), legend.title=element_text(size=7), legend.text=element_text(size=6.5))

# c) score is not a sample-size proxy
qc <- data.frame(layer=factor(c("proteome","transcriptome","methylation","metabolome"),
                              levels=c("proteome","transcriptome","methylation","metabolome")),
                 r=c(-0.16,-0.28,-0.24,-0.02))
pc <- ggplot(qc, aes(layer, r, fill=layer)) +
  geom_col(width=0.66, color="black", linewidth=0.2) +
  geom_hline(yintercept=0, linewidth=0.4) +
  geom_hline(yintercept=-0.3, linetype="dashed", color="grey45") +
  annotate("text", x=2.5, y=-0.315, label="power-artifact threshold (|r| = 0.3)", size=2.0, color="grey35") +
  scale_fill_manual(values=LAYERCOL, guide="none") +
  coord_cartesian(ylim=c(-0.36,0.04)) +
  labs(x=NULL, y="corr(score, log N)") +
  theme(axis.text.x=element_text(angle=35, hjust=1))

# d) replication across INDEPENDENT cohorts (CENTRAL x DIRECT-PLUS methylome)
rp <- data.frame(test=factor(c("genome-wide","q < 0.10","sign\nconcord."), levels=c("genome-wide","q < 0.10","sign\nconcord.")),
                 val=c(0.245,0.88,0.944))
pd <- ggplot(rp, aes(test, val)) +
  geom_col(width=0.62, fill=TXN, color="black", linewidth=0.2) +
  geom_text(aes(label=sprintf("%.2f", val)), vjust=-0.4, size=2.2) +
  coord_cartesian(ylim=c(0,1.04)) +
  labs(x=NULL, y="Spearman / concordance") +
  theme(axis.text.x=element_text(size=8, angle=25, hjust=1))

# e) within-cohort split-half reproducibility (internal stability; distinct from d's cross-cohort transfer)
sh <- read.csv(file.path(SD,"F1e_split_half.csv")); sh$layer <- factor(sh$layer, levels=c("transcriptome","methylome"))
pe <- ggplot(sh, aes(layer, sp_signed, fill=layer)) +
  geom_col(width=0.6, color="black", linewidth=0.2) +
  geom_errorbar(aes(ymin=sp_signed-sd, ymax=sp_signed+sd), width=0.14, linewidth=0.4) +
  geom_text(aes(label=sprintf("%.2f", sp_signed)), vjust=-1.3, size=2.4) +
  geom_text(aes(label=sprintf("rev %.1f×", rev_enr)), y=0.06, size=2.0, color="white") +
  scale_fill_manual(values=c(transcriptome=TXN, methylome=MEY), guide="none") +
  coord_cartesian(ylim=c(0,0.95)) +
  labs(x=NULL, y="signed Spearman (half vs half)") +
  theme(axis.text.x=element_text(size=8, angle=18, hjust=1))

# f) cross-ancestry generalizability + its limits (directional-only, n=12 underpowered)
an <- read.csv(file.path(SD,"F1f_ancestry.csv")); anr <- an[an$panel=="rho",]
anr$metric <- factor(anr$metric, levels=c("BLACK vs EUR","EUR vs EUR (ceiling)"))
conc <- an$value[an$panel=="concordance"][1]
pf <- ggplot(anr, aes(metric, value, fill=metric)) +
  geom_col(width=0.58, color="black", linewidth=0.2) +
  geom_hline(yintercept=0, linewidth=0.4) +
  geom_text(aes(label=sprintf("%+.2f", value)), vjust=ifelse(anr$value>=0,-0.5,1.4), size=2.4) +
  scale_fill_manual(values=c("BLACK vs EUR"="grey62","EUR vs EUR (ceiling)"=TXN), guide="none") +
  scale_x_discrete(labels=c("BLACK vs EUR"="BLACK\nvs EUR","EUR vs EUR (ceiling)"="EUR vs EUR\n(ceiling)")) +
  coord_cartesian(ylim=c(-0.16,0.42)) +
  annotate("label", x=1.5, y=0.40, size=1.72, fill="white", linewidth=0, lineheight=0.95, vjust=1,
           label=sprintf("directional concord. %.0f%% (p=0.009)\nn=12 — genome-wide underpowered", conc*100)) +
  labs(x=NULL, y="genome-wide ρ") +
  theme(axis.text.x=element_text(size=8))

fig <- wrap_plots(pa, pb, pc, pd, pe, pf, ncol=3) +
  plot_annotation(tag_levels="a") &
  theme(plot.tag=element_text(face="bold", size=11, hjust=0, vjust=1),
        plot.tag.position=c(0.01,1.02),
        plot.margin=margin(11,8,8,6,"pt"))
save_figure_ckm(fig, "Figure_1", width_mm=183, height_mm=150, output_dir=".")
cat("done\n")

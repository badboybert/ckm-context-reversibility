# render_ED6.R — Extended Data 6: intervention classes share a reversible core by signature OVERLAP,
# with NO within-class privilege (fold-enrichment, orthogonal to F5e correlation).
suppressMessages({library(ggplot2); library(patchwork); library(scales)})
source("theme_publication_ckm.R")  # bundled alongside this script
set_publication_defaults_ckm()
SD <- "source_data"; rd <- function(p) read.csv(file.path(SD, p), stringsAsFactors=FALSE)
TXN<-"#009E73"; UP<-"#B2182B"; DN<-"#2166AC"; NEU<-"#333333"

# a) class-grouped mean fold-enrichment — within-DRUG ~chance, no privilege
cm <- rd("ED6_classmeans.csv"); cm <- cm[order(cm$mean_fold),]; cm$classpair <- factor(cm$classpair, levels=cm$classpair)
cm$hl <- cm$classpair=="within-DRUG"
pa <- ggplot(cm, aes(mean_fold, classpair, fill=hl)) +
  geom_vline(xintercept=1, linetype="dashed", linewidth=0.4, color="grey45") +
  geom_col(width=0.62, color="black", linewidth=0.2) +
  geom_text(aes(label=sprintf("%.2f×", mean_fold)), hjust=-0.2, size=2.8) +
  scale_fill_manual(values=c(`TRUE`=UP, `FALSE`="grey65"), guide="none") +
  scale_x_continuous(limits=c(0,1.6), breaks=c(0,0.5,1,1.5)) +
  annotate("text", x=1, y=0.45, label="chance", size=2.8, color="grey35", fontface="italic", hjust=1, vjust=0) +
  labs(x="mean sig-set fold-enrichment", y=NULL) +
  theme(axis.text.y=element_text(size=7.5))

# b) powered pairwise folds, colored by class-pair (within-DRUG depletion highlighted)
pp <- rd("ED6_pairs.csv"); pp$lab <- paste(pp$a, pp$b, sep="×")
# Test the within-* classes FIRST: a bare grepl("SURGERY|DIET") also matches "within-SURGERY"
# and mis-coloured that within-class pair as cross-class.
pp$grp <- ifelse(pp$classpair=="within-DRUG","within-DRUG",
          ifelse(pp$classpair=="within-SURGERY","within-SURGERY",
          ifelse(grepl("SURGERY|DIET",pp$classpair),"cross-class (weight-loss)","other")))
pp <- pp[order(pp$fold),]; pp$lab <- factor(pp$lab, levels=pp$lab)
pb <- ggplot(pp, aes(fold, lab, color=grp)) +
  geom_vline(xintercept=1, linetype="dashed", linewidth=0.4, color="grey45") +
  geom_point(size=2) +
  scale_color_manual(values=c("within-DRUG"=UP, "cross-class (weight-loss)"=TXN,
                              "within-SURGERY"="grey55", "other"="grey75"), name=NULL) +
  # upper limit must cover the largest fold (empagliflozin×MS_bariatric = 1.84): the old
  # limits=c(0.7,1.6) silently dropped that pair, so only 16 of the 17 overlap-eligible pairs plotted.
  scale_x_continuous(limits=c(0.7,1.95), breaks=c(0.8,1.0,1.2,1.4,1.6,1.8)) +
  # (in-plot callout REMOVED: the legend already states "the within-DRUG pair metformin ×
  #  empagliflozin was depleted at 0.79-fold — fewer shared significant reversers than chance",
  #  the pair is named on the y-axis and coloured as within-DRUG, and the chance line is explained
  #  in the legend. Nothing here is unique to the callout. It also could not be placed cleanly:
  #  x=0.72 is the panel floor and any string long enough to carry the text reached the dashed
  #  chance line at x=1.0, which struck through it.)
  labs(x="sig-set fold-enrichment", y=NULL) +
  theme(axis.text.y=element_text(size=6.3),
        legend.position="inside", legend.position.inside=c(0.40,0.04), legend.justification=c(0,0),
        legend.key.size=unit(7,"pt"), legend.text=element_text(size=6.0),
        legend.background=element_rect(fill=alpha("white",0.7), color=NA))

fig <- wrap_plots(pa, pb, ncol=2) +
  plot_annotation(tag_levels="a",
    theme=theme(plot.background=element_rect(fill="white", color=NA))) &
  theme(plot.tag=element_text(face="bold", size=11, hjust=0, vjust=1), plot.tag.position=c(0.01,1.04),
        plot.margin=margin(13,6,8,6,"pt"))
save_figure_ckm(fig, "Figure_ED6", width_mm=178, height_mm=82, output_dir=".")
cat("done\n")

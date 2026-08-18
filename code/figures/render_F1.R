# render_F1.R — Figure 1 (SPINE): the question, the data landscape, a reproducible score.
# R/ggplot, locked publication theme. Panel a redrawn 2026-07-21 as a horizontal 3-stage
# routing schematic (atlas -> harmonized per-mark effect -> the four planned analyses); the
# declarative "context, not the molecule" endpoint is removed (it pre-announced the conclusion)
# and "40 cohorts" is corrected to "40 analysis panels". NO "prespecified" (nothing was preregistered).
suppressMessages({library(ggplot2); library(patchwork); library(scales)})
source("theme_publication_ckm.R")  # bundled alongside this script
set_publication_defaults_ckm()
RES <- "../../results"
SD <- "source_data"
PRO<-"#0072B2"; TXN<-"#009E73"; MET<-"#CC79A7"; MEY<-"#D55E00"; NEU<-"#333333"
LAYERCOL <- c(proteome=PRO, transcriptome=TXN, metabolome=MET, methylation=MEY)
arr <- arrow(length=unit(4.5,"pt"), type="closed")
BOXG<-"#F2F7F5"; BOXO<-"#FBEEE4"   # green/orange soft panel fills

# a) study-design schematic — horizontal 3-stage flow (full-width banner)
rectb <- function(xmin,xmax,ymin,ymax,fill,col,lw=0.4)
  annotate("rect", xmin=xmin,xmax=xmax, ymin=ymin,ymax=ymax, fill=fill, color=col, linewidth=lw)
txt <- function(x,y,label,size=2.9,hjust=0.5,fontface="plain",color="black",lineheight=0.95)
  annotate("text", x=x,y=y, label=label, size=size, hjust=hjust, fontface=fontface, color=color, lineheight=lineheight)

pa <- ggplot() + xlim(0,31.2) + ylim(0,10) +
  # ---- stage 1: intervention atlas ----
  rectb(0.3,9.6,0.4,9.6,BOXG,TXN) +
  txt(4.95,8.95,"INTERVENTION ATLAS",size=2.95,fontface="bold") +
  txt(1.5,7.5,"40",size=5.0,fontface="bold",color=TXN) +
  txt(3.1,7.85,"analysis",size=2.9,hjust=0) + txt(3.1,7.1,"panels",size=2.9,hjust=0) +
  txt(0.8,6.0,"within-person & randomized",size=2.7,hjust=0,fontface="italic",color="grey30") +
  txt(0.8,5.1,"Intervention families",size=2.85,hjust=0,fontface="bold") +
  txt(0.8,4.35,"diet · surgery · drugs · exercise",size=2.7,hjust=0) +
  txt(0.8,3.35,"Molecular layers",size=2.85,hjust=0,fontface="bold") +
  annotate("point",x=1.05,y=2.55,size=2.0,color=PRO)+txt(1.4,2.55,"proteome",size=2.7,hjust=0)+
  annotate("point",x=4.9,y=2.55,size=2.0,color=TXN)+txt(5.25,2.55,"transcriptome",size=2.7,hjust=0)+
  annotate("point",x=1.05,y=1.75,size=2.0,color=MEY)+txt(1.4,1.75,"methylome",size=2.7,hjust=0)+
  annotate("point",x=4.9,y=1.75,size=2.0,color=MET)+txt(5.25,1.75,"metabolome",size=2.7,hjust=0)+
  txt(0.8,1.15,"adipose·liver·muscle·blood·plasma",size=2.2,hjust=0,fontface="italic",color="grey35") +
  annotate("segment",x=9.7,xend=10.8,y=5.0,yend=5.0,arrow=arr,linewidth=0.5) +
  # ---- stage 2: harmonized per-mark effect ----
  rectb(10.9,20.1,0.4,9.6,BOXG,TXN) +
  txt(15.5,8.95,"PER-MARK EFFECT",size=2.95,fontface="bold") +
  annotate("point",x=12.6,y=7.35,size=6.2,shape=21,fill="#DCE9F5",color=PRO)+txt(12.6,7.35,"base-\nline",size=2.2)+
  annotate("point",x=15.5,y=7.35,size=6.2,shape=21,fill="#FBE9D6",color=MEY)+txt(15.5,7.35,"inter-\nvention",size=2.05)+
  annotate("point",x=18.4,y=7.35,size=6.2,shape=21,fill="#DCEEE4",color=TXN)+txt(18.4,7.35,"follow-\nup",size=2.2)+
  annotate("segment",x=13.5,xend=14.6,y=7.35,yend=7.35,arrow=arr,linewidth=0.4)+
  annotate("segment",x=16.4,xend=17.5,y=7.35,yend=7.35,arrow=arr,linewidth=0.4)+
  rectb(11.5,19.5,4.3,6.0,"white",TXN,lw=0.35)+
  txt(15.5,5.4,"signed reversal effect",size=2.85,fontface="bold",color=TXN)+
  txt(15.5,4.7,"standardized (post − pre)",size=2.7)+
  txt(15.5,3.2,"rank-standardized within panel,",size=2.6,fontface="italic",color="grey30")+
  txt(15.5,2.45,"synthesized across panels",size=2.6,fontface="italic",color="grey30")+
  txt(15.5,1.5,"→ direction, magnitude, uncertainty",size=2.55)+
  annotate("segment",x=20.2,xend=21.3,y=5.0,yend=5.0,arrow=arr,linewidth=0.5) +
  # ---- stage 3: the four planned analyses ----
  rectb(21.4,30.9,0.4,9.6,BOXO,MEY) +
  txt(26.15,8.95,"ANALYSIS FRAMEWORK",size=2.95,fontface="bold") +
  rectb(21.75,30.55,7.55,8.55,"white",MEY,lw=0.3)+txt(22.1,8.05,"1  Reproducibility",size=2.65,hjust=0)+
  rectb(21.75,30.55,4.30,7.25,"#F6E2CE",MEY,lw=0.5)+
  txt(22.1,6.75,"2  Determinants (PRIMARY)",size=2.65,hjust=0,fontface="bold")+
  txt(22.1,6.00,"constraint · druggability · GWAS burden",size=2.35,hjust=0)+
  txt(22.1,5.35,"cis-architecture · τ · causal status",size=2.35,hjust=0)+
  txt(22.1,4.70,"transcriptome-powered (k = 9)",size=2.35,hjust=0,fontface="italic",color="grey35")+
  rectb(21.75,30.55,2.55,4.05,"white",MEY,lw=0.3)+
  txt(22.1,3.55,"3  Shared context",size=2.65,hjust=0)+
  txt(22.1,2.90,"tissue × intervention family",size=2.4,hjust=0,color="grey35")+
  rectb(21.75,30.55,0.85,2.35,"white",MEY,lw=0.3)+
  txt(22.1,1.85,"4  Restoration & durability",size=2.65,hjust=0)+
  txt(22.1,1.20,"BMI-aligned lean ref · multi-year",size=2.4,hjust=0,color="grey35")+
  theme_void(base_family=BASE_FONT) +
  theme(plot.margin=margin(2,3,2,3,"pt"))

# b) data landscape (layer x tissue, bubble = panel n)
m <- read.delim(file.path(RES,"panel_manifest_full.tsv"), stringsAsFactors=FALSE)
m <- m[!grepl("EXCLUDED", toupper(m$status)),]
tmap <- c(plasma="plasma/serum", serum="plasma/serum", whole_blood="blood", PBMC="blood",
          SAT_adipose="adipose", skeletal_muscle="muscle", liver="liver", blood="blood",
          blood_neutrophil="blood")
m$tis <- tmap[m$tissue]; m$np <- pmax(suppressWarnings(as.numeric(m$n_pairs)),0,na.rm=FALSE)
# An unmapped tissue silently becomes NA and drops the panel from the bubble grid.
stopifnot(!any(is.na(m$tis)))
ag <- aggregate(np~layer+tis, m, function(x) sum(x,na.rm=TRUE))
ag$layer <- factor(ag$layer, levels=c("metabolome","methylation","proteome","transcriptome"))
ag$tis   <- factor(ag$tis,   levels=c("plasma/serum","blood","adipose","liver","muscle"))
pb <- ggplot(ag, aes(tis, layer, size=np, color=layer)) +
  geom_point(alpha=0.85) +
  scale_size_area(max_size=9, breaks=c(100,1000,3000), name="panel n") +
  scale_color_manual(values=LAYERCOL, guide="none") +
  scale_x_discrete(labels=c("plasma/serum"="plasma/serum","blood"="blood","adipose"="adipose","liver"="liver","muscle"="muscle")) +
  scale_y_discrete(labels=c("metabolome"="metabolome","methylation"="methylome","proteome"="proteome","transcriptome"="transcriptome")) +
  labs(x=NULL, y=NULL) +
  theme(axis.text.x=element_text(size=8, angle=30, hjust=1),
        legend.position="right", legend.key.size=unit(8,"pt"), legend.title=element_text(size=7.5), legend.text=element_text(size=7))

# c) score is not a sample-size proxy
qc <- data.frame(layer=factor(c("proteome","transcriptome","methylome","metabolome"),
                              levels=c("proteome","transcriptome","methylome","metabolome")),
                 r=c(-0.16,-0.28,-0.24,-0.02))
pc <- ggplot(qc, aes(layer, r, fill=layer)) +
  geom_col(width=0.66, color="black", linewidth=0.2) +
  geom_hline(yintercept=0, linewidth=0.4) +
  geom_hline(yintercept=-0.3, linetype="dashed", color="grey45") +
  annotate("text", x=2.5, y=-0.335, label="|ρ| = 0.3 bound", size=2.55, color="grey35") +
  scale_fill_manual(values=c(proteome=PRO,transcriptome=TXN,methylome=MEY,metabolome=MET), guide="none") +
  coord_cartesian(ylim=c(-0.36,0.04)) +
  labs(x=NULL, y="Spearman ρ with log N") +
  theme(axis.text.x=element_text(angle=35, hjust=1))

# d) replication across INDEPENDENT cohorts (CENTRAL x DIRECT-PLUS methylome)
rp <- data.frame(test=factor(c("genome-wide","q < 0.10","sign\nconcord."), levels=c("genome-wide","q < 0.10","sign\nconcord.")),
                 val=c(0.245,0.88,0.944))
pd <- ggplot(rp, aes(test, val)) +
  geom_col(width=0.62, fill=TXN, color="black", linewidth=0.2) +
  geom_text(aes(label=sprintf("%.2f", val)), vjust=-0.4, size=2.85) +
  coord_cartesian(ylim=c(0,1.04)) +
  labs(x=NULL, y="Spearman / concordance") +
  theme(axis.text.x=element_text(size=8, angle=25, hjust=1))

# e) within-cohort split-half reproducibility (internal stability; distinct from d's cross-cohort transfer)
sh <- read.csv(file.path(SD,"F1e_split_half.csv")); sh$layer <- factor(sh$layer, levels=c("transcriptome","methylome"))
pe <- ggplot(sh, aes(layer, sp_signed, fill=layer)) +
  geom_col(width=0.6, color="black", linewidth=0.2) +
  geom_errorbar(aes(ymin=sp_signed-sd, ymax=sp_signed+sd), width=0.14, linewidth=0.4) +
  geom_text(aes(label=sprintf("%.2f", sp_signed)), vjust=-1.3, size=2.85) +
  geom_text(aes(label=sprintf("%.1f×", rev_enr)), y=0.06, size=2.85, color="white") +
  scale_fill_manual(values=c(transcriptome=TXN, methylome=MEY), guide="none") +
  coord_cartesian(ylim=c(0,0.95)) +
  labs(x=NULL, y="split-half Spearman") +
  theme(axis.text.x=element_text(size=8, angle=18, hjust=1))

# f) cross-ancestry generalizability + its limits (directional-only, n=12 underpowered)
an <- read.csv(file.path(SD,"F1f_ancestry.csv")); anr <- an[an$panel=="rho",]
anr$metric <- factor(anr$metric, levels=c("BLACK vs EUR","EUR vs EUR (ceiling)"))
conc <- an$value[an$panel=="concordance"][1]
pf <- ggplot(anr, aes(metric, value, fill=metric)) +
  geom_col(width=0.58, color="black", linewidth=0.2) +
  geom_hline(yintercept=0, linewidth=0.4) +
  geom_text(aes(label=sprintf("%+.2f", value)), vjust=ifelse(anr$value>=0,-0.5,1.4), size=2.85) +
  scale_fill_manual(values=c("BLACK vs EUR"="grey62","EUR vs EUR (ceiling)"=TXN), guide="none") +
  scale_x_discrete(labels=c("BLACK vs EUR"="Black cohort\nvs EUR","EUR vs EUR (ceiling)"="EUR vs EUR\n(ceiling)")) +
  coord_cartesian(ylim=c(-0.16,0.46)) +
  annotate("label", x=1.5, y=0.44, size=2.5, fill="white", linewidth=0, lineheight=0.95, vjust=1,
           label=sprintf("directional %.0f%% (p=0.009)\nn=12 · underpowered", conc*100)) +
  labs(x=NULL, y="genome-wide ρ") +
  theme(axis.text.x=element_text(size=8))

layout <- "AAA\nBCD\nEFF"
fig <- wrap_plots(a=pa, b=pb, c=pc, d=pd, e=pe, f=pf, design=layout, heights=c(1.18,1,1)) +
  plot_annotation(title="Is molecular reversibility a portable property of the mark, or a property of biological context?",
                  tag_levels="a",
                  theme=theme(plot.title=element_text(size=9, face="italic", color="grey25", hjust=0, margin=margin(b=3)))) &
  theme(plot.tag=element_text(face="bold", size=11, hjust=0, vjust=1),
        plot.tag.position=c(0.01,1.02),
        plot.margin=margin(9,8,7,6,"pt"))
save_figure_ckm(fig, "Figure_1", width_mm=183, height_mm=176, output_dir=".")
cat("done\n")

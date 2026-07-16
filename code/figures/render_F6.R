# render_F6.R — Figure 6: Reversibility is a property of context, not the molecule.
# (a) 4-context rev_beta correlation matrix (tissue block); (b) same- vs cross-tissue directional concordance;
# (c) exemplar molecules x contexts — the SAME molecule reverses differently by context. LOCKED theme.
suppressMessages({library(ggplot2); library(patchwork); library(scales)})
source("C:/Users/Bert/Downloads/CKM papers/paper 3.eas/analysis/tools/theme_publication_ckm.R")
set_publication_defaults_ckm()
SD <- "source_data"; rd <- function(p) read.csv(file.path(SD, p), stringsAsFactors = FALSE, check.names = FALSE)
TXN<-"#009E73"; UP<-"#B2182B"; DN<-"#2166AC"; NEU<-"#333333"
CTXO <- c("liver","adipose-CR","adipose-diet","blood")
CTXLAB <- c("liver"="liver","adipose-CR"="adip.\nCR","adipose-diet"="adip.\ndiet","blood"="blood")
CTXLABX <- c("liver"="liver","adipose-CR"="adip. CR","adipose-diet"="adip. diet","blood"="blood")  # single-line, for angled x-axes

# ---- a) context correlation matrix (rev_beta Spearman) ----
cm <- rd("F6a_context_corr.csv")
cm$ctx_a <- factor(cm$ctx_a, levels=CTXO); cm$ctx_b <- factor(cm$ctx_b, levels=rev(CTXO))
cm$show <- ifelse(as.integer(cm$ctx_a)==as.integer(factor(cm$ctx_b,levels=CTXO)), NA, cm$rho)  # blank diagonal
pa <- ggplot(cm, aes(ctx_a, ctx_b, fill=show)) +
  geom_tile(color="white", linewidth=0.6) +
  geom_text(aes(label=ifelse(is.na(show),"",sprintf("%+.2f",show))), size=2.3,
            color=ifelse(!is.na(cm$show) & cm$show>0.35,"white","black")) +
  scale_fill_gradient2(low=DN, mid="white", high=TXN, midpoint=0, limits=c(-0.15,0.55), oob=scales::squish,
                       na.value="grey92", name="ρ", breaks=c(0,0.25,0.5)) +
  scale_x_discrete(labels=CTXLABX) + scale_y_discrete(labels=CTXLAB) +
  annotate("rect", xmin=1.5, xmax=3.5, ymin=1.5, ymax=3.5, fill=NA, color="black", linewidth=0.6) +
  labs(x=NULL, y=NULL) +
  theme(axis.text.x=element_text(size=9, angle=30, hjust=1), axis.text.y=element_text(size=9),
        legend.key.width=unit(5,"pt"), legend.key.height=unit(9,"pt"),
        legend.title=element_text(size=9), legend.text=element_text(size=8))

# ---- b) directional concordance: same- vs cross-tissue ----
cc <- rd("F6b_concordance.csv")
cc$kind <- factor(cc$kind, levels=c("same tissue","cross tissue"))
set.seed(5)
pb <- ggplot(cc, aes(kind, concordance, color=kind)) +
  geom_hline(yintercept=0.5, linetype="dashed", linewidth=0.4, color="grey45") +
  geom_jitter(width=0.12, height=0, size=2, alpha=0.95) +
  scale_color_manual(values=c("same tissue"=TXN,"cross tissue"="grey55"), guide="none") +
  scale_y_continuous(limits=c(0.3,1.04), breaks=c(0.5,0.75,1.0)) +
  scale_x_discrete(labels=c("same tissue"="same\ntissue","cross tissue"="cross\ntissue")) +
  annotate("text", x=2.42, y=0.5, label="chance", size=2.1, color="grey35", fontface="italic", hjust=1) +
  annotate("text", x=1, y=0.99, label="0.99", size=2.6, hjust=-0.35) +
  labs(x=NULL, y="directional concordance")

# ---- c) exemplar molecules x contexts (the demonstration) ----
ex <- rd("F6c_exemplars.csv")
HOMEO <- c("liver","adipose","blood")
ex$home <- factor(ex$home, levels=HOMEO)
ord <- ex[order(ex$home, -abs(apply(ex[,CTXO],1,function(z) max(abs(z),na.rm=TRUE)))),]
glev <- ord$gene
long <- do.call(rbind, lapply(CTXO, function(c) data.frame(gene=ex$gene, context=c, rev_beta=ex[[c]], home=ex$home)))
long$gene <- factor(long$gene, levels=rev(glev)); long$context <- factor(long$context, levels=CTXO)
pc <- ggplot(long, aes(context, gene, fill=rev_beta)) +
  geom_tile(color="white", linewidth=0.6) +
  geom_text(aes(label=ifelse(is.na(rev_beta),"·",sprintf("%+.1f",rev_beta))), size=1.95,
            color=ifelse(!is.na(long$rev_beta) & abs(long$rev_beta)>0.9,"white","black")) +
  scale_fill_gradient2(low=DN, mid="white", high=UP, midpoint=0, limits=c(-1.5,1.5), oob=scales::squish,
                       na.value="grey90", name="reversal\nβ", breaks=c(-1.5,0,1.5), labels=c("down","0","up")) +
  scale_x_discrete(labels=CTXLABX) +
  labs(x=NULL, y=NULL) +
  theme(axis.text.x=element_text(size=9, angle=30, hjust=1), axis.text.y=element_text(size=9, face="italic"),
        legend.key.width=unit(5,"pt"), legend.key.height=unit(9,"pt"),
        legend.title=element_text(size=9), legend.text=element_text(size=8))

# ---- d) reversible sig-set overlap: shared within tissue, not across ----
od <- rd("F6d_reversible_overlap.csv"); od$kind <- factor(od$kind, levels=c("within tissue","cross tissue"))
set.seed(4)
pd <- ggplot(od, aes(kind, fold, color=kind)) +
  geom_hline(yintercept=1, linetype="dashed", linewidth=0.4, color="grey45") +
  geom_jitter(width=0.1, height=0, size=1.9, alpha=0.95) +
  scale_color_manual(values=c("within tissue"=TXN,"cross tissue"="grey55"), guide="none") +
  scale_y_continuous(limits=c(0.9,1.7), breaks=c(1.0,1.3,1.6)) +
  annotate("text", x=1, y=1.63, label="1.6×", size=2.3, vjust=-0.4) +
  annotate("text", x=2, y=0.95, label="chance", size=2.0, color="grey35", fontface="italic") +
  scale_x_discrete(labels=c("within tissue"="within\ntissue","cross tissue"="cross\ntissue")) +
  labs(x=NULL, y="reversible overlap (fold)")

# ---- e) persistent pole: barely shares within tissue, never across (the specificity null) ----
oe <- rd("F6e_persistent_overlap.csv"); oe$kind <- factor(oe$kind, levels=c("within tissue","cross tissue"))
set.seed(4)
pe <- ggplot(oe, aes(kind, fold, color=kind)) +
  geom_hline(yintercept=1, linetype="dashed", linewidth=0.4, color="grey45") +
  geom_jitter(width=0.1, height=0, size=1.9, alpha=0.95) +
  scale_color_manual(values=c("within tissue"="#9970AB","cross tissue"="grey55"), guide="none") +
  scale_y_continuous(limits=c(0.9,1.7), breaks=c(1.0,1.3,1.6)) +
  annotate("text", x=2, y=0.95, label="chance", size=2.0, color="grey35", fontface="italic") +
  scale_x_discrete(labels=c("within tissue"="within\ntissue","cross tissue"="cross\ntissue")) +
  labs(x=NULL, y="persistent overlap (fold)")

# 3-column layout: (a matrix / b concordance) | c exemplar (prominent, tall) | (d / e overlap companions)
fig <- ((pa / pb) | pc | (pd / pe)) + plot_layout(widths=c(0.95,1.0,0.9)) +
  plot_annotation(tag_levels="a") &
  theme(plot.tag = element_text(face="bold", size=11, hjust=0, vjust=1),
        plot.tag.position = c(0.01, 1.02),
        plot.margin = margin(11, 6, 8, 6, "pt"))
save_figure_ckm(fig, "Figure_6", width_mm=183, height_mm=150, output_dir=".")
cat("done\n")

library(clonevol)
library(fishplot)
library(RColorBrewer)
library(ggplot2)
library(DESeq2)
library(dendextend)
library(lsa)
library(dplyr)
library(ggpubr)
library(multcomp)
library(rstatix)

# Figure 5a

x <- read.table("ClonEvol_DataSet.txt", header=TRUE, sep="\t")

# shorten vaf column names as they will be
vaf.col.names <- grep('.vaf', colnames(x), value=T)
var.col.names <- grep('.var', colnames(x), value=T)
ref.col.names <- grep('.ref', colnames(x), value=T)
sample.names <- gsub('.vaf', '', vaf.col.names)
x[, sample.names] <- x[, vaf.col.names] 
vaf.col.names <- sample.names
# prepare sample grouping
#sample.groups <- c("M1", "M5A", "M5B", "M5C", "M6A", "M6B", "M7", "M8A", "M8B", "M9", "M10C", "M10D", "M10E", "M11A", "M11B", "M12A", "M12B", "M15", "M16", "M17", "M18"); names(sample.groups) <- vaf.col.names
sample.groups <- c("Chest_Wall", "Chest_Wall", "M1", "Chest_Wall", "Lung_Right", "Lung_Right", "Lung_Right", "Diaphragm", "Lung_Left", "Lung_Right", "Pericardium", "Lung_Right", "Diaphragm", "Diaphragm", "Lung_Left", "Lung_Left", "Lung_Left", "Pericardium", "Lymph_Nodes", "Lymph_Nodes", "Chest_Wall"); names(sample.groups) <- vaf.col.names
# setup the order of clusters to display in various plots (later)
x <- x[order(x$cluster),]
#If you want ClonEvol to choose colors for you, simply set it to NULL, like this:
clone.colors <- c("black", "yellow3", "#56B4E9", "#009E73", "green2", "red", "gray50", "magenta" , "orange", "maroon")
#clone.colors <- NULL
#pdf('BoxPlot.pdf', width = 10, height = 30, useDingbats = FALSE, title='') 
pp <- plot.variant.clusters(x,
    cluster.col.name = 'cluster', 
    show.cluster.size = FALSE, 
    cluster.size.text.color = 'blue', 
    vaf.col.names = vaf.col.names, 
    vaf.limits = 100, 
    sample.title.size = 20,
    violin = FALSE,
    box = TRUE,
    jitter = TRUE,
    jitter.shape = 1,
    jitter.color = clone.colors, 
    jitter.size = 1,
    jitter.alpha = 1, 
    jitter.center.method = 'median', 
    jitter.center.size = 1, 
    jitter.center.color = 'darkgray', 
    jitter.center.display.value = 'none', 
    highlight = 'is.driver', 
    highlight.shape = 21,
    highlight.color = 'white', 
    highlight.fill.color = 'white', 
    highlight.note.col.name = 'gene', 
    highlight.note.size = 2, 
    order.by.total.vaf = FALSE)
#dev.off()
plot.pairwise(x, col.names = vaf.col.names, out.prefix = 'Plots/AllMets_variants.pairwise.plot',
colors = clone.colors)
pdf('FlowPlot.pdf', width=20, height=10, useDingbats=FALSE, title='') 
plot.cluster.flow(x, vaf.col.names = vaf.col.names,
                  sample.names = c("M1", "M5A", "M5B", "M5C", "M6A", "M6B", "M7", "M8A", "M8B", "M9", "M10C", "M10D", "M10E", "M11A", "M11B", "M12A", "M12B", "M15", "M16", "M17", "M18"), 
                  colors = clone.colors, y.title = "CCF (%)")
dev.off()
y = infer.clonal.models(variants = x, 
    cluster.col.name = 'cluster',
    vaf.col.names = vaf.col.names, 
    sample.groups = sample.groups,
    #sample.groups = NULL,
    cancer.initiation.model='monoclonal', 
    subclonal.test = 'bootstrap', 
    subclonal.test.model = 'non-parametric', 
    num.boots = 1000,
    founding.cluster = 1, 
    cluster.center = 'mean', 
    ignore.clusters = 2,
    clone.colors = clone.colors,
    #JBH adjusted original 0.01
    min.cluster.vaf = 0.01,
    p.value.cutoff = 0.05,
    # min probability that CCF(clone) is non-negative
    #JBH adjusted original 0.05
    sum.p = 0.05,
    # alpha level in confidence interval estimate for CCF(clone) 
    #JBH adjusted original 0.05
    alpha = 0.05,
    random.seed = 666,
    score.model.by = "probability",
    vaf.in.percent = TRUE)
y <- transfer.events.to.consensus.trees(y, 
    x[x$is.driver,],
    cluster.col.name = 'cluster', 
    event.col.name = 'gene')

y <- convert.consensus.tree.clone.to.branch(y, branch.scale = 'sqrt')

plot.clonal.models(y,
    # box plot parameters
    box.plot = FALSE,
    fancy.boxplot = TRUE,
    fancy.variant.boxplot.highlight = 'is.driver', 
    fancy.variant.boxplot.highlight.shape = 21, 
    fancy.variant.boxplot.highlight.fill.color = 'red', 
    fancy.variant.boxplot.highlight.color = 'black', 
    fancy.variant.boxplot.highlight.note.col.name = 'gene', 
    fancy.variant.boxplot.highlight.note.color = 'blue', 
    fancy.variant.boxplot.highlight.note.size = 2,
    fancy.variant.boxplot.jitter.alpha = 1, 
    fancy.variant.boxplot.jitter.center.color = 'grey50', 
    fancy.variant.boxplot.base_size = 12, 
    fancy.variant.boxplot.plot.margin = 1, 
    fancy.variant.boxplot.vaf.suffix = '.VAF',
    # bell plot parameters
    clone.shape = 'bell',
    bell.event = TRUE, 
    bell.event.label.color = 'blue', 
    bell.event.label.angle = 60, 
    clone.time.step.scale = 1, 
    bell.curve.step = 2,
    # node-based consensus tree parameters
    merged.tree.plot = TRUE, 
    tree.node.label.split.character = NULL, 
    tree.node.shape = 'circle', 
    tree.node.size = 30, 
    tree.node.text.size = 0.5, 
    merged.tree.node.size.scale = 1.25, 
    merged.tree.node.text.size.scale = 2.5, 
    merged.tree.cell.frac.ci = FALSE,
    # branch-based consensus tree parameters
    merged.tree.clone.as.branch = TRUE, 
    mtcab.event.sep.char = ',', 
    mtcab.branch.text.size = 1, 
    mtcab.branch.width = 0.75, 
    mtcab.node.size = 3, 
    mtcab.node.label.size = 1, 
    mtcab.node.text.size = 1.5,
    # cellular population parameters
    cell.plot = TRUE,
    num.cells = 100,
    cell.border.size = 0.25, 
    cell.border.color = 'black', 
    clone.grouping = 'horizontal', 
    #meta-parameters 
    scale.monoclonal.cell.frac = TRUE, 
    show.score = FALSE,
    cell.frac.ci = TRUE, 
    disable.cell.frac = FALSE, 
    # output figure parameters 
    out.dir = 'Plots', 
    out.format = 'pdf', 
    overwrite.output = TRUE,
    width = 30,
    height = 20,
    max.num.models.to.plot = 35,
    # vector of width scales for each panel from left to right 
    #panel.widths = c(3,4,2,4,2))
    panel.widths = c(3,2,2,2))

input_names = sample.names

fishplot_input <- generateFishplotInputs(y, rescale = TRUE)

vafs = clonevol:::estimate.clone.vaf(y$variants, vaf.col.names=vaf.col.names)
scales = vafs[vafs$cluster == 1, vaf.col.names]/max(vafs[vafs$cluster == 1, vaf.col.names])
for (i in 1:length(fishplot_input$cell.fractions)){
fishplot_input$cell.fractions[[i]] = (fishplot_input$cell.fractions[[i]] *
as.matrix(scales[rep(1,nrow(fishplot_input$cell.fractions[[i]])),]))
}

pdf("FishPlot.pdf", width = 8, height = 4)
for (i in 1:length(fishplot_input$cell.fractions)){
    fishplot_objects <- createFishObject(fishplot_input$cell.fractions[[i]], fishplot_input$parents[[i]], fix.missing.clones=TRUE)
    fishplot_layout = layoutClones(fishplot_objects)
    fishplot_layout = setCol(fishplot_layout, fishplot_input$clonevol.clone.colors)
    fishPlot(fishplot_layout, shape = "spline", cex.title = 0.7, 
             vlines = seq(1, length(input_names)), vlab = input_names, pad.left = 0.5)
}

dev <- dev.off()

# Fifure 5b

df <- read.csv("ELISpot_Data.txt", header = TRUE, sep = "\t")

positions <- c("MNRRPILTIIN", "NRRPILTIINT", "RRPILTIINTG", "RPILTIINTGR", "PILTIINTGRL", "ILTIINTGRLQ", "LTIINTGRLQW", "NRRPILTIIN", "RRPILTIINT", "RPILTIINTG", "PILTIINTGR", "ILTIINTGRL", "LTIINTGRLQ", "TIINTGRLQW", "RRPILTIIN", "RPILTIINT", "PILTIINTG", "ILTIINTGR", "LTIINTGRL", "TIINTGRLQ", "IINTGRLQW")

pdf('BarPlot.pdf', height = 5, width = 8, useDingbats = FALSE)

p<-ggplot(data=df, aes(x=Peptide, y=Average)) +
  geom_bar(stat="identity", fill="steelblue") + 
  theme_minimal() + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) + 
  scale_x_discrete(limits = positions) + 
  scale_y_continuous(name="IFNg spots per 5 x 104 cells", limits=c(0, 100)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  geom_hline(yintercept = 50, colour = "red", linetype = "dashed")
p

dev.off()

# Figure 5b <- Changed?

IO<-read.delim(file="IO_signatures.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE,row.names = 1)
IO_TMB<-subset(IO, Days %in% c("799","2033"))

#remove unused rows
IO_TMB<-IO_TMB[-c(1),]
IO_TMB[1,2] <-"M1"
tcga_brca<-read.table(file="brca_tcga_merged.txt", header = TRUE, sep="\t", na.strings = "NA")



#get how many tumor from each tumor type are in the pan dataset
table(tcga_brca$Type)

#Normalize mutation counts
tcga_brca$Mutation.Count.Normalized<-tcga_brca$Mutation.Count/40
means<-aggregate(tcga_brca$Mutation.Count.Normalized, by=list(tcga_brca$Type),mean, paired=FALSE)
median<-aggregate(tcga_brca$Mutation.Count.Normalized, by=list(tcga_brca$Type),median, paired=FALSE)
means
median

pdf("TMBnormalized_TCGA_Delta.pdf", useDingbats = FALSE,height = 4, width = 5)
ggplot(data = tcga_brca, mapping = aes(x = reorder(Type, Mutation.Count.Normalized,FUN = "median"), y = log10(Mutation.Count.Normalized)))+
  geom_hline(yintercept = 0.185635, colour = 'black', linetype='dotted')+
  geom_hline(yintercept = 0.4988405, colour = 'black', linetype='dotted')+
  geom_point(color = "grey")+ 
  geom_jitter(data = IO_TMB, mapping = aes(x = reorder(Tissue, TMB.norm,FUN = "median"), y = log10(TMB.norm),color = Tissue), width = 0.2, height = 0.1)+
  stat_summary(
    geom = "crossbar",
    fun = "median",
    col = "red",
    width = 0.4,
    lwd = 0.3,
  )+
  scale_x_discrete(limits=c(
    "BRCA","TNBC", "M1", "chest wall", "diaphragm","lung","mediastinal_lymph_nodes", "pericardium"), labels=c(
      "TCGA BRCA \n (primary)","TCGA TNBC \n (primary)", "M1", "chest wall", "diaphragm","lung","lymph nodes", "pericardium"))+
  scale_colour_manual(values = c(
    'BRCA' = 'grey',
    'TNBC' = 'grey',
    'lung' = '#1F78B4',
    'diaphragm' = "#E31A1C",
    'mediastinal_lymph_nodes' = "#FF7F00",
    'pericardium' = "#33A02C",
    'chest wall' = "purple4",
    'M1' = "black"))+
  scale_y_continuous(name = "Mutations per Mb", labels =c("0", "0.1", "1", "10", "100"))+ labs(x = "")+
  theme_classic()+
  theme(legend.position="none")+rotate_x_text(angle = 45)
dev.off()

## Fig 5c

# Dendograms

# Molecular Clock -> Code in SupplFig5 for MutationTimeR

# Pie charts 

pdf('Pie_Charts.pdf', width = 5, height = 5)

slices <- c(92, 14, 85, 48)
lbls <- c("clonal [early]", "clonal [NA]", "clonal [late]", "subclonal")
pie(slices, labels = lbls, col=c("gray10", "gray30", "gray50", "yellow"), main="M5A CLS")

slices <- c(128, 9, 41, 43)
lbls <- c("clonal [early]", "clonal [NA]", "clonal [late]", "subclonal")
pie(slices, labels = lbls, col=c("gray10", "gray30", "gray50", "yellow"), main="M5B CLS")

slices <- c(88, 9, 110, 93)
lbls <- c("clonal [early]", "clonal [NA]", "clonal [late]", "subclonal")
pie(slices, labels = lbls, col=c("gray10", "gray30", "gray50", "yellow"), main="M5C CLS")

slices <- c(89, 19, 109, 3)
lbls <- c("clonal [early]", "clonal [NA]", "clonal [late]", "subclonal")
pie(slices, labels = lbls, col=c("gray10", "gray30", "gray50", "yellow"), main="M6A CLS")

slices <- c(97, 24, 140, 6)
lbls <- c("clonal [early]", "clonal [NA]", "clonal [late]", "subclonal")
pie(slices, labels = lbls, col=c("gray10", "gray30", "gray50", "yellow"), main="M6B CLS")

slices <- c(129, 2, 65, 24)
lbls <- c("clonal [early]", "clonal [NA]", "clonal [late]", "subclonal")
pie(slices, labels = lbls, col=c("gray10", "gray30", "gray50", "yellow"), main="M7 CLS")

slices <- c(124, 5, 24, 21)
lbls <- c("clonal [early]", "clonal [NA]", "clonal [late]", "subclonal")
pie(slices, labels = lbls, col=c("gray10", "gray30", "gray50", "yellow"), main="M8A CLS")

slices <- c(91, 71, 26)
lbls <- c("clonal [early]", "clonal [late]", "subclonal")
pie(slices, labels = lbls, col=c("gray10", "gray50", "yellow"), main="M8B CLS")

slices <- c(103, 24, 79, 18)
lbls <- c("clonal [early]", "clonal [NA]", "clonal [late]", "subclonal")
pie(slices, labels = lbls, col=c("gray10", "gray30", "gray50", "yellow"), main="M9 CLS")

slices <- c(106, 16, 139, 14)
lbls <- c("clonal [early]", "clonal [NA]", "clonal [late]", "subclonal")
pie(slices, labels = lbls, col=c("gray10", "gray30", "gray50", "yellow"), main="M10C CLS")

slices <- c(105, 56, 88)
lbls <- c("clonal [early]", "clonal [late]", "subclonal")
pie(slices, labels = lbls, col=c("gray10", "gray50", "yellow"), main="M10D CLS")

slices <- c(106, 9, 52, 67)
lbls <- c("clonal [early]", "clonal [NA]", "clonal [late]", "subclonal")
pie(slices, labels = lbls, col=c("gray10", "gray30", "gray50", "yellow"), main="M10E CLS")

slices <- c(69, 6, 9, 3)
lbls <- c("clonal [early]", "clonal [NA]", "clonal [late]", "subclonal")
pie(slices, labels = lbls, col=c("gray10", "gray30", "gray50", "yellow"), main="M11A CLS")

slices <- c(70, 12, 34)
lbls <- c("clonal [early]", "clonal [NA]", "clonal [late]")
pie(slices, labels = lbls, col=c("gray10", "gray30", "gray50"), main="M11B CLS")

slices <- c(108, 20, 94, 51)
lbls <- c("clonal [early]", "clonal [NA]", "clonal [late]", "subclonal")
pie(slices, labels = lbls, col=c("gray10", "gray30", "gray50", "yellow"), main="M12A CLS")

slices <- c(3, 18)
lbls <- c("clonal [early]", "subclonal")
pie(slices, labels = lbls, col=c("gray10", "yellow"), main="M12B CLS")

slices <- c(98, 2, 21)
lbls <- c("clonal [early]", "clonal [NA]", "clonal [late]")
pie(slices, labels = lbls, col=c("gray10", "gray30", "gray50"), main="M15 CLS")

slices <- c(27, 3, 18)
lbls <- c("clonal [early]", "clonal [NA]", "clonal [late]")
pie(slices, labels = lbls, col=c("gray10", "gray30", "gray50"), main="M16 CLS")

slices <- c(118, 24, 62)
lbls <- c("clonal [early]", "clonal [NA]", "clonal [late]")
pie(slices, labels = lbls, col=c("gray10", "gray30", "gray50"), main="M17 CLS")

slices <- c(98, 64, 39, 7)
lbls <- c("clonal [early]", "clonal [NA]", "clonal [late]", "subclonal")
pie(slices, labels = lbls, col=c("gray10", "gray30", "gray50", "yellow"), main="M18 CLS")

dev.off()

## Metagene

load('Perou_RNAseq_raw_counts.RData')
load('Perou_Samples_data.RData')

patients <- c('A1', 'A5', 'A7', 'A11', 'A15', 'A17', 'A20', 'A23', 'A26', 'A30')
pd_f <- subset(pd_f, pd_f$Patient %in% patients)

dds <- DESeqDataSetFromMatrix(countData=raw_counts_basal, 
                              colData=pd_f, 
                              design=~Patient+site)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

dds <- DESeq(dds)
res <- results(dds)

de_genes <- as.data.frame(res)[res$padj < 0.05, ]
de_genes <- de_genes[order(abs(de_genes$log2FoldChange), decreasing = TRUE), ]

vsdata <- vst(dds, blind=FALSE)
exhaustion <- c('PDCD1','LAG3','HAVCR2','KLRG1','TIGIT','CD244','BTLA','CTLA4','ENTPD1','CD160','ID2')
de_exhaustion <- subset(de_genes, rownames(de_genes) %in% exhaustion)

vsd_exh <- assay(vsdata)[rownames(assay(vsdata)) %in% rownames(de_exhaustion), ]

vsd_f_metagene <- vsd_exh * de_exhaustion$log2FoldChange
vsd_f_metagene <- scale(t(vsd_f_metagene))
metagene_rank <- apply(vsd_f_metagene, 1, mean)

mgene_df <- as.data.frame(metagene_rank[order(metagene_rank, decreasing = FALSE)])
colnames(mgene_df)[1] <- 'rank'
mgene_df$Patient
saveRDS(mgene_df, file = 'Metagene_exhaustion_df.rds')

## Figure 5c 

neot <- read.csv('neoantigen_time.csv')
neot$class_col <- neot$class
neot$class_col <- c('green', 'blue', 'purple2')[as.factor(neot$class)]

m <- as.matrix(neot[, c(2:4)])
rownames(m) <- neot$sample

d <- as.dist(1 - cosine(t(m)))
hc <- hclust(d, method = 'ward.D')
dend <- as.dendrogram(hc)

pal = colorRampPalette(c("navy", "blue", 'lightskyblue',"firebrick3"))
neot$T1.histBeta_cols <- pal(length(unique(neot$T1.histBeta)))[as.factor(neot$T1.histBeta)]

legend_image <- as.raster(matrix(neot$T1.histBeta_cols, ncol=1))

pdf('Figure_5c.pdf')
par(mar = c(2, 1, 1, 5))

dend %>% set("branches_lwd", 2)  %>% set("leaves_pch", 19) %>% set("leaves_cex", 2) %>% set("leaves_col", value = neot$T1.histBeta_cols[hc$order]) %>% set("labels_cex", 1.27) %>% set("labels_col", value = neot$class_col[hc$order])%>% plot(horiz=TRUE)

rasterImage(legend_image, 0.5, 24, 0.47, 18)  

text(0.49, 24.7, 'Time', cex = 1.47)
text(0.45, 24, '- 1', cex = 1.17)
text(0.433, 21.3, '- 0.77', cex = 1.17)
text(0.425, 18.3, '- 0.005', cex = 1.17)

legend('topleft', legend = c('Early', 'Inter', 'Late'), fill = c('green', 'blue', 'purple2'), border = 'white', bty = 'n', cex = 1.47, title = 'Class')
dev.off()

# Figure 5d

df <- read.csv(file = 'Neoantigens_Prop.txt', sep = "\t")
df["class"] <- as.factor(df$class)

res.aov <- aov(p.late ~ class, data = df)
summary(glht(res.aov, linfct = mcp(class = "Tukey")))

group1 <- c("Inter", "Late", "Late")
group2 <- c("Early", "Early", "Inter")
p <- c(0.8799, 0.0159, 0.0276)
y.position <- c(0.28,0.82,0.40)

stat.test <- data.frame(group1, group2, p, y.position)

plot2 <- ggboxplot(df, x = "class", y = "p.late",
color = "class", palette = c("chartreuse", "blue1", "darkviolet"),
order = c("Early", "Inter", "Late"),
ylab = "p.late", xlab = "MolClocl_Timing")
plot2 <- plot2 + stat_pvalue_manual(
  stat.test,  label = "p", tip.length = 0
  )

#Plot both in one PDF

pdf('BoxPlot.pdf', width = 8, height = 4)
grid.arrange(plot1, plot2, ncol=2)
dev.off()

## Figure 5f

pdf("Fig5f.pdf", useDingbats = FALSE, height = 8, width = 8)
col2 <- colorRampPalette(c("#053061","#2166AC","#4393C3","#92C5DE","#D1E5F0", "#FFFFFF","#FDDBC7", "#F4A582","#D6604D","#B2182B", "#67001F"))
corrplot::corrplot(HLAimbalance.IO.rcorr$r, sig.level = 0.05, insig = "blank",  p.mat = HLAimbalance.IO.rcorr$P,
         type = "upper", tl.cex = 1.0, col = col2 (50))
dev.off()

## Suppl Figure 5d 

neo_es_long <- read.csv('Figure_5d.csv')

sc1 <- ggplot(data = neo_es_long[neo_es_long$variable == 'p.early', ], aes(x=rank, y=value, color=variable)) + geom_point(size=2)  + theme_ipsum() + geom_text_repel(aes(label = SampleID), show.legend = FALSE, size=6, max.overlaps = 1000) + geom_smooth(method=lm , color="red", fill="#69b3a2", se=TRUE) + ylab('Neoantigen time early') + xlab('Exhaustion score') + theme(axis.text.y = element_text(size = 14), axis.title.y=element_text(size=14,face="bold"), axis.text.x = element_text(size = 14), axis.title.x=element_text(size=14,face="bold")) +  scale_color_manual(values = c('navy')) + theme(legend.position = "none") +  annotate("text", x=-0.3, y=0.2, label= "r= -0.60", cex=7)

sc2 <- ggplot(data = neo_es_long[neo_es_long$variable == 'p.inter', ], aes(x=rank, y=value, color=variable)) + geom_point(size=2)  + theme_ipsum() + geom_text_repel(aes(label = SampleID), show.legend = FALSE, size=6, size=6, max.overlaps = 1000) + geom_smooth(method=lm , color="red", fill="#69b3a2", se=TRUE) + ylab('Neoantigen time intermediate') + xlab('Exhaustion score') + theme(axis.text.y = element_text(size = 14), axis.title.y=element_text(size=14,face="bold"), axis.text.x = element_text(size = 14), axis.title.x=element_text(size=14,face="bold")) +  scale_color_manual(values = c('forestgreen')) + theme(legend.position = "none") + theme(legend.position = "none") +  annotate("text", x=-0.3, y=0.07, label= "r= 0.22", cex=7)

sc3 <- ggplot(data = neo_es_long[neo_es_long$variable == 'p.late', ], aes(x=rank, y=value, color=variable)) + geom_point(size=2)  + theme_ipsum() + geom_text_repel(aes(label = SampleID), show.legend = FALSE, size=6, size=6, max.overlaps = 1000) + geom_smooth(method=lm , color="red", fill="#69b3a2", se=TRUE)  + ylab('Neoantigen time late') + xlab('Exhaustion score') + theme(axis.text.y = element_text(size = 14), axis.title.y=element_text(size=14,face="bold"), axis.text.x = element_text(size = 14), axis.title.x=element_text(size=14,face="bold")) +  scale_color_manual(values = c('orange')) + theme(legend.position = "none") + theme(legend.position = "none") +  annotate("text", x=-0.3, y=0.35, label= "r= 0.67", cex=7)


ggarrange(sc1, sc2, sc3, ncol = 1, labels = c('A', 'B', 'C'))

ggsave( filename = "plots_v2.pdf", 
        device = cairo_pdf,
        plot = ggarrange(sc1, sc2, sc3, ncol = 1, labels = c('A', 'B', 'C')), 
        width = 8, height = 18)

## Suppl Figure 5f

library("MutationTimeR")

chrOffset <- MutationTimeR:::chrOffset

chrOffset <- c("1" = 0, "2" = 249250621, "3" = 492449994, "4" = 690472424, "5" = 881626700, "6" = 1062541960,
       "7" = 1233657027, "8" = 1392795690, "9" = 1539159712, "10" = 1680373143, "11" = 1815907890, "12" = 1950914406,
       "13" = 2084766301, "14" = 2199936179, "15" = 2307285719, "16" = 2409817111, "17" = 2500171864, "18" = 2581367074,
       "19" = 2659444322, "20" = 2718573305, "21" = 2781598825, "22" = 222829728720, "X" =2881033286)

.histBeta <- function(bb, time="time",n.min=5, s=seq(0.005,0.995,0.01)){
    s <- pmax(0.001,pmin(0.999, s))
    cols <- paste0(time,c("",".lo",".up"))
    w <- which(bb$n.snv_mnv > n.min & !is.na(mcols(bb)[cols[1]]) & !duplicated(bb))
    if(length(w)==0) return(rep(NA, length(s)))
    d <- apply(as.matrix(mcols(bb)[w,c(cols, "n.snv_mnv")]), 1, function(x){
                beta <- .betaFromCi(x[1:3])
                beta <- (beta * x[4] + 5*c(1,1))/(x[4]+5) # "posterior" with prior B(1,1)
                dbeta(s,beta[1],beta[2])
            })
    wd <- as.numeric(width(bb)[w])
    o <- d %*% wd
}

.betaFromCi <- function(x, weight.mode=5){
    if(any(is.na(x))) return(rep(NA,2))
    f <- function(par,x) {
        beta <- exp(par)
        sum((qbeta(c(0.025,0.975), beta[1], beta[2])-x[-1])^2)+weight.mode*((beta[1]-1)/(beta[1]+beta[2]-2)-x[1])^2
    }
    tryCatch(exp(optim(c(0.1,0.1), fn=f,x=x)$par), error=c(1,1))
}

#.plotVcf <- function(vcf, bb, col = RColorBrewer::brewer.pal(9, "Set1")[c(3,4,2,1,9)], ID = meta(header(vcf))[[1]]["ID",1], legend=TRUE, lty.grid=1, col.grid="grey", xaxt=TRUE, pch=16, pch.out=pch, cex=0.66, xlim=c(0,chrOffset["MT"])) {
.plotVcf <- function(vcf, bb, col = RColorBrewer::brewer.pal(9, "Set1")[c(3,4,2,1,9)], ID = meta(header(vcf))[[1]]["ID",1], legend=TRUE, lty.grid=1, col.grid="grey", xaxt=TRUE, pch=16, pch.out=pch, cex=2, xlim=c(0,chrOffset["MT"])) {
    cls <- factor(paste(as.character(info(vcf)$CLS)), levels = c("clonal [early]","clonal [late]","clonal [NA]","subclonal" , "NA"))
    plot(NA,NA, xlab='', ylab="VAF", ylim=c(0,1), xlim=xlim, xaxt="n", cex=cex)
    abline(v = chrOffset[1:25], lty=lty.grid, col=col.grid)
    if(xaxt) mtext(side=1, line=1, at=chrOffset[1:24] + diff(chrOffset[1:25])/2, text=names(chrOffset[1:24]))
    for(i in which(!sapply(bb$timing_param, is.null))) {
        s <- start(bb)[i]
        e <- end(bb)[i]
        x <- chrOffset[as.character(seqnames(bb)[i])]
        y <- bb$timing_param[[i]][,"f"]
        l <- bb$timing_param[[i]][,"pi.s"] * bb$timing_param[[i]][,"P.m.sX"]
        l[is.na(l)] <- 0
        if(any(is.na(c(s,e,x,y,l)))) next
        segments(s+x,y,e+x,y, lwd=l*4+.1)
        #text(x=(s+e)/2 +x, y=y, paste(signif(bb$timing_param[[i]][,"m"],2),signif(bb$timing_param[[i]][,"cfi"]/purityPloidy[meta(header(vv))["ID",1],"purity"],2), sep=":"), pos=3, cex=0.5)
    }
    points(start(vcf) + chrOffset[as.character(seqnames(vcf))], getAltCount(vcf)/getTumorDepth(vcf),col=col[cls],  pch=ifelse(info(vcf)$pMutCNTail < 0.025 | info(vcf)$pMutCNTail > 0.975, pch.out , pch),  cex=cex)                
    if(legend) legend("topleft", pch=19, col=col, legend=paste(as.numeric(table(cls)), levels(cls)), bg='white')
}

.plotBB <- function(bb, ylim=c(0,max(max(bb$total_cn, na.rm=TRUE))), col=RColorBrewer::brewer.pal(4,"Set2"), type=c("lines","bars"), legend=TRUE, lty.grid=1, col.grid="grey", xaxt=TRUE, xlim=c(min(chrOffset[as.character(seqnames(bb))]+start(bb)),max(chrOffset[as.character(seqnames(bb))]+end(bb)))){
    type <- match.arg(type)
    s <- c(1:22, "X","Y")
    l <- as.numeric(width(refLengths[as.character(seqnames(refLengths)) %in% s]))
    names(l) <- s
    plot(NA,NA, ylab="Copy number",xlab="",xlim=xlim, ylim=ylim, xaxt="n")
    c <- cumsum(l)
    axis(side=1, at=c(0,c), labels=rep('', length(l)+1))
    if(xaxt) mtext(side=1, at= cumsum(l) - l/2, text=names(l), line=1)
    #abline(v=c, lty=3)
    if(type=="lines"){
        x0 <- start(bb) + cumsum(l)[as.character(seqnames(bb))] - l[as.character(seqnames(bb))]
        x1 <- end(bb) + cumsum(l)[as.character(seqnames(bb))] - l[as.character(seqnames(bb))]
        lwd <- 5* bb$clonal_frequency / max(bb$clonal_frequency)
        segments(x0=x0, bb$major_cn ,x1, bb$major_cn, col=col[1], lwd=lwd)
        segments(x0=x0, bb$minor_cn -.125,x1, bb$minor_cn-.125, col=col[2], lwd=lwd)
        segments(x0=x0, bb$total_cn+.125,x1, bb$total_cn+.125, col=1, lwd=lwd)
#   cv <- coverage(bb)
#   cv <- cv[s[s%in%names(cv)]]
#   par(xpd=NA)
#   for(n in names(cv)){
#       cc <- cv[[n]]
#       segments(start(cc) + cumsum(l)[n] - l[n] ,-runValue(cc)/2,end(cc)+ cumsum(l)[n] - l[n], -runValue(cc)/2, col=4)
#   }
    }else{
        ub <- unique(bb)
        f <- findOverlaps(ub,bb)
        m <- t(model.matrix( ~ 0 + factor(queryHits(f))))
        ub$major_cn <- m %*% mg14::na.zero(bb$major_cn * bb$clonal_frequency) / max(bb$clonal_frequency)
        ub$minor_cn <- m %*% mg14::na.zero(bb$minor_cn * bb$clonal_frequency) / max(bb$clonal_frequency)
        ub$total_cn <- ub$major_cn + ub$minor_cn
        ub$clonal_frequency <- max(bb$clonal_frequency)
        x0 <- start(ub) + cumsum(l)[as.character(seqnames(ub))] - l[as.character(seqnames(ub))]
        x1 <- end(ub) + cumsum(l)[as.character(seqnames(ub))] - l[as.character(seqnames(ub))]
        rect(x0,0,x1, ub$major_cn, col=col[2], lwd=NA)
        rect(x0,ub$major_cn,x1, ub$total_cn, col=col[1], lwd=NA)
        abline(h = 1:floor(ylim[2]), lty=lty.grid, col=col.grid)
    }
    abline(v = chrOffset[1:25], lty=lty.grid, col=col.grid)
    if(xaxt) mtext(side=1, line=1, at=chrOffset[1:24] + diff(chrOffset[1:25])/2, text=names(chrOffset[1:24]))
    if(legend){
        if(type=="lines") legend("topleft", c("Total CN","Major CN","Minor CN"), col=c("black", col[1:2]), lty=1, lwd=2, bg='white')
        else legend("topleft", c("Major CN","Minor CN"), fill=col[1:2], bg='white')
    }
}

.plotTiming <- function(bb, time=mcols(bb), col=paste0(RColorBrewer::brewer.pal(5,"Set2")[c(3:5)],"88"), legend=TRUE, col.grid='grey', lty.grid=1, xlim=c(0,chrOffset["MT"]), plot=2){
    plot(NA,NA, xlab='', ylab="Time [mutations]", ylim=c(0,1), xlim=xlim, xaxt="n")
    if(any(!is.na(bb$time)))
        tryCatch({
                    bb <- bb[!is.na(bb$time)]
                    s <- start(bb)
                    e <- end(bb)
                    x <- chrOffset[as.character(seqnames(bb))]
                    y <- time[,"time"]
                    rect(s+x,time[,"time.lo"],e+x,time[,"time.up"], border=NA, col=col[time[,"type"]], angle = ifelse(bb$time.star=="*" | is.na(bb$time.star),45,135), density=ifelse(bb$time.star == "***", -1, 72))
                    segments(s+x,y,e+x,y)
                    
                    if("time.2nd" %in% colnames(time)){ 
                        w <- !is.na(time[,"time.2nd"])
                        if(sum(w) != 0 & plot==2){
                            s <- start(bb)[w]
                            e <- end(bb)[w]
                            x <- chrOffset[as.character(seqnames(bb))][w]
                            y <- time[w,"time.2nd"]
                            rect(s+x,time[w,"time.2nd.lo"],e+x,time[w,"time.2nd.up"], border=NA, col=sub("88$","44",col)[as.numeric(time[w,"type"])], angle = ifelse(bb$time.star[w]=="*" | is.na(bb$time.star[w]),45,135), density=ifelse(bb$time.star[w] == "***", -1, 72))
                            segments(s+x,y,e+x,y)
                        }
                    }
                }, error=function(x) warning(x))
    abline(v = chrOffset[1:25], lty=lty.grid, col=col.grid)
    s <- c(1:22, "X","Y")
    l <- as.numeric(width(refLengths[as.character(seqnames(refLengths)) %in% s]))
    names(l) <- s
    c <- cumsum(l)
    axis(side=1, at=c(0,c), labels=rep('', length(l)+1))
    mtext(side=1, line=1, at=chrOffset[1:24] + diff(chrOffset[1:25])/2, text=names(chrOffset[1:24]))
    if(legend) legend("topleft", levels(time[,"type"]), fill=col, bg="white")
}
                 
plotSampleJBH <- function (vcf, cn, sv = NULL, title = "", regions = NULL, ylim.cn = c(0, 
    5), layout.height = c(4, 1.2, 3.5), y.sv = ylim.cn[2] - 1) 
{
    if (is.null(regions)) 
        #regions <- refLengths[1:24]
        regions <- rL
    p <- par()
    layout(matrix(1:3, ncol = 1), height = layout.height)
    par(mar = c(0.5, 3, 0.5, 0.5), mgp = c(2, 0.25, 0), bty = "L", 
        las = 2, tcl = -0.25, cex = 1)
    xlim = c(min(chrOffset[as.character(seqnames(regions))] + 
        start(regions)), max(chrOffset[as.character(seqnames(regions))] + 
        end(regions)))
    bbb <- cn[cn %over% regions]
    .plotVcf(vcf[vcf %over% regions], bbb, legend = FALSE, col.grid = "white", 
        xaxt = FALSE, cex = 0.33, xlim = xlim)
    mtext(line = -1, side = 3, title, las = 1)
    .plotBB(bbb, ylim = ylim.cn, legend = FALSE, type = "bar", 
        col.grid = "white", col = c("lightgrey", "darkgrey"), 
        xaxt = FALSE, xlim = xlim)
    tryCatch({
        par(xpd = NA)
        if (!is.null(sv)) 
            .plotSv(sv, y1 = y.sv, regions = regions, add = TRUE)
        par(xpd = FALSE)
    }, error = function(x) warning(x))
    par(mar = c(3, 3, 0.5, 0.5))
    .plotTiming(bbb, xlim = xlim, legend = FALSE, col.grid = NA)
    if (length(regions) == 1) 
        axis(side = 1, at = pretty(c(start(regions), end(regions))) + 
            chrOffset[as.character(seqnames(regions))], labels = sitools::f2si(pretty(c(start(regions), 
            end(regions)))))
    if (any(!is.na(cn$time))) {
        y0 <- seq(0.005, 0.995, 0.01)
        s <- .histBeta(cn)
        g <- colorRampPalette(RColorBrewer::brewer.pal(4, "Set1")[c(3, 
            2, 4)])(100)
        segments(x0 = chrOffset["MT"], y0 = y0, x1 = chrOffset["MT"] + 
            s/max(s) * 1e+08, col = g, lend = 3)
        getMode <- function(s) {
            if (all(is.na(s))) 
                return(NA)
            w <- which.max(s)
            if (w %in% c(1, length(s))) {
                m <- which(c(0, diff(s)) > 0 & c(diff(s), 0) < 
                  0)
                if (length(m) == 0) 
                  return(w)
                m <- m[which.max(s[m])]
                return(if (s[w] > 2 * s[m]) w else m)
            }
            else return(w)
        }
        abline(h = y0[getMode(s)], lty = 5)
        if ("time.2nd" %in% colnames(mcols(cn))) 
            if (any(!is.na(cn$time.2nd))) {
                s2 <- .histBeta(cn, time = "time.2nd")
                segments(x0 = chrOffset["MT"] + s/max(s) * 1e+08, 
                  y0 = y0, x1 = chrOffset["MT"] + s/max(s) * 
                    1e+08 + s2/max(s) * 1e+08, col = paste0(g, 
                    "44"), lend = 3)
                abline(h = y0[getMode(s2)], lty = 3)
            }
    }
    par(p[setdiff(names(p), c("cin", "cra", "csi", "cxy", "din", 
        "page"))])
}
                       
chrOffset <- MutationTimeR:::chrOffset

rL <- GRanges(
    seqnames = Rle(c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y"), c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)),
    ranges = IRanges(start = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1), end = c(248956422, 242193529, 198295559, 190214555, 181538259, 170805979, 159345973, 145138636, 138394717, 133797422, 135086622, 133275309, 114364328, 107043718, 101991189, 90338345, 83257441, 80373285, 58617616, 64444167, 46709983, 50818468, 156040895, 57227415)),
    strand = Rle(strand(c("*", "*", "*", "*", "*", "*", "*", "*", "*", "*", "*", "*", "*", "*", "*", "*", "*", "*", "*", "*", "*", "*", "*", "*")), c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)))
refLengths <- rL
                       
getMode <- function(s){
    if(all(is.na(s))) return(NA)
    w <- which.max(s)
    if(w %in% c(1, length(s))){
        m <- which(c(0,diff(s))>0 & c(diff(s),0)<0)
        if(length(m)==0) return(w)
        m <- m[which.max(s[m])]
        return(if(s[w] > 2*s[m]) w else m) 
    } else return(w)
}

#Plotting

vcf <- readVcf("M1_MutationTimeR.vcf")
cn <- read.csv("M1_CN.txt", sep="\t")
bb <- makeGRangesFromDataFrame(cn,
                         keep.extra.columns=TRUE,
                         ignore.strand=FALSE,
                         seqinfo=NULL,
                         seqnames.field=c("seqnames", "seqname",
                                          "chromosome", "chrom",
                                          "chr", "chromosome_name",
                                          "seqid"),
                         start.field="start",
                         end.field=c("end", "stop"),
                         strand.field="strand",
                         starts.in.df.are.0based=FALSE)
clusters <- data.frame(cluster=1:3, proportion=c(0.7308, 0.1221, 0.0103), n_ssms=c(394, 513, 264))
mt <- mutationTime(vcf, bb, clusters=clusters, n.boot=10)
vcf <- addMutTime(vcf, mt$V)
mcols(bb) <- cbind(mcols(bb),mt$T)

pdf("M1_MutationTimeR_JBH.pdf", 10, 8)
plotSampleJBH(vcf,bb)
dev.off()

T <- table(mt$V$CLS)
s <- .histBeta(bb)
s2 <- .histBeta(bb, time = "time.2nd")
y0 <- seq(0.005, 0.995, 0.01)
h = y0[getMode(s)]
h2 = y0[getMode(s2)]

T$T1.histBeta <- h  
T$T2.histBeta <- h2  
write.table(T, file='M1_Results_MutationTimeR.tsv', quote=FALSE, sep='\t', row.names = FALSE)

vcf <- readVcf("M5A_MutationTimeR.vcf")
cn <- read.csv("M5A_CN.txt", sep="\t")
bb <- makeGRangesFromDataFrame(cn,
                         keep.extra.columns=TRUE,
                         ignore.strand=FALSE,
                         seqinfo=NULL,
                         seqnames.field=c("seqnames", "seqname",
                                          "chromosome", "chrom",
                                          "chr", "chromosome_name",
                                          "seqid"),
                         start.field="start",
                         end.field=c("end", "stop"),
                         strand.field="strand",
                         starts.in.df.are.0based=FALSE)
clusters <- data.frame(cluster=1:4, proportion=c(0.9616, 0.1665, 0.1317, 0.8731), n_ssms=c(394, 260, 459,264))
mt <- mutationTime(vcf, bb, clusters=clusters, n.boot=10)
vcf <- addMutTime(vcf, mt$V)
mcols(bb) <- cbind(mcols(bb),mt$T)

pdf("M5A_MutationTimeR_JBH.pdf", 10, 8)
plotSampleJBH(vcf,bb)
dev.off()

T <- table(mt$V$CLS)
s <- .histBeta(bb)
s2 <- .histBeta(bb, time = "time.2nd")
y0 <- seq(0.005, 0.995, 0.01)
h = y0[getMode(s)]
h2 = y0[getMode(s2)]

T$T1.histBeta <- h  
T$T2.histBeta <- h2  
write.table(T, file='M5A_Results_MutationTimeR.tsv', quote=FALSE, sep='\t', row.names = FALSE)

vcf <- readVcf("M5B_MutationTimeR.vcf")
cn <- read.csv("M5B_CN.txt", sep="\t")
bb <- makeGRangesFromDataFrame(cn,
                         keep.extra.columns=TRUE,
                         ignore.strand=FALSE,
                         seqinfo=NULL,
                         seqnames.field=c("seqnames", "seqname",
                                          "chromosome", "chrom",
                                          "chr", "chromosome_name",
                                          "seqid"),
                         start.field="start",
                         end.field=c("end", "stop"),
                         strand.field="strand",
                         starts.in.df.are.0based=FALSE)
clusters <- data.frame(cluster=1:4, proportion=c(0.9691, 0.1641, 0.1918, 0.4606), n_ssms=c(394, 260, 459, 264))
mt <- mutationTime(vcf, bb, clusters=clusters, n.boot=10)
vcf <- addMutTime(vcf, mt$V)
mcols(bb) <- cbind(mcols(bb),mt$T)

pdf("M5B_MutationTimeR_JBH.pdf", 10, 8)
plotSampleJBH(vcf,bb)
dev.off()

T <- table(mt$V$CLS)
s <- .histBeta(bb)
s2 <- .histBeta(bb, time = "time.2nd")
y0 <- seq(0.005, 0.995, 0.01)
h = y0[getMode(s)]
h2 = y0[getMode(s2)]

T$T1.histBeta <- h  
T$T2.histBeta <- h2  
write.table(T, file='M5B_Results_MutationTimeR.tsv', quote=FALSE, sep='\t', row.names = FALSE)

vcf <- readVcf("M5C_MutationTimeR.vcf")
cn <- read.csv("M5C_CN.txt", sep="\t")
bb <- makeGRangesFromDataFrame(cn,
                         keep.extra.columns=TRUE,
                         ignore.strand=FALSE,
                         seqinfo=NULL,
                         seqnames.field=c("seqnames", "seqname",
                                          "chromosome", "chrom",
                                          "chr", "chromosome_name",
                                          "seqid"),
                         start.field="start",
                         end.field=c("end", "stop"),
                         strand.field="strand",
                         starts.in.df.are.0based=FALSE)
clusters <- data.frame(cluster=1:5, proportion=c(0.9609, 0.1829, 0.1957, 0.0303, 0.889), n_ssms=c(394, 260, 171, 459, 264))
mt <- mutationTime(vcf, bb, clusters=clusters, n.boot=10)
vcf <- addMutTime(vcf, mt$V)
mcols(bb) <- cbind(mcols(bb),mt$T)

pdf("M5C_MutationTimeR_JBH.pdf", 10, 8)
plotSampleJBH(vcf,bb)
dev.off()

T <- table(mt$V$CLS)
s <- .histBeta(bb)
s2 <- .histBeta(bb, time = "time.2nd")
y0 <- seq(0.005, 0.995, 0.01)
h = y0[getMode(s)]
h2 = y0[getMode(s2)]

T$T1.histBeta <- h  
T$T2.histBeta <- h2  
write.table(T, file='M5C_Results_MutationTimeR.tsv', quote=FALSE, sep='\t', row.names = FALSE)

vcf <- readVcf("M6A_MutationTimeR.vcf")
cn <- read.csv("M6A_CN.txt", sep="\t")
bb <- makeGRangesFromDataFrame(cn,
                         keep.extra.columns=TRUE,
                         ignore.strand=FALSE,
                         seqinfo=NULL,
                         seqnames.field=c("seqnames", "seqname",
                                          "chromosome", "chrom",
                                          "chr", "chromosome_name",
                                          "seqid"),
                         start.field="start",
                         end.field=c("end", "stop"),
                         strand.field="strand",
                         starts.in.df.are.0based=FALSE)
clusters <- data.frame(cluster=1:5, proportion=c(0.9999, 0.1651, 0.0202, 0.1883, 0.8928), n_ssms=c(394, 260, 513, 245, 264))
mt <- mutationTime(vcf, bb, clusters=clusters, n.boot=10)
vcf <- addMutTime(vcf, mt$V)
mcols(bb) <- cbind(mcols(bb),mt$T)

pdf("M6A_MutationTimeR_JBH.pdf", 10, 8)
plotSampleJBH(vcf,bb)
dev.off()

T <- table(mt$V$CLS)
s <- .histBeta(bb)
s2 <- .histBeta(bb, time = "time.2nd")
y0 <- seq(0.005, 0.995, 0.01)
h = y0[getMode(s)]
h2 = y0[getMode(s2)]

T$T1.histBeta <- h  
T$T2.histBeta <- h2  
write.table(T, file='M6A_Results_MutationTimeR.tsv', quote=FALSE, sep='\t', row.names = FALSE)

vcf <- readVcf("M6B_MutationTimeR.vcf")
cn <- read.csv("M6B_CN.txt", sep="\t")
bb <- makeGRangesFromDataFrame(cn,
                         keep.extra.columns=TRUE,
                         ignore.strand=FALSE,
                         seqinfo=NULL,
                         seqnames.field=c("seqnames", "seqname",
                                          "chromosome", "chrom",
                                          "chr", "chromosome_name",
                                          "seqid"),
                         start.field="start",
                         end.field=c("end", "stop"),
                         strand.field="strand",
                         starts.in.df.are.0based=FALSE)
clusters <- data.frame(cluster=1:5, proportion=c(1, 0.2717, 0.0202, 0.3392, 0.7257), n_ssms=c(394, 260, 513, 245, 264))
mt <- mutationTime(vcf, bb, clusters=clusters, n.boot=10)
vcf <- addMutTime(vcf, mt$V)
mcols(bb) <- cbind(mcols(bb),mt$T)

pdf("M6B_MutationTimeR_JBH.pdf", 10, 8)
plotSampleJBH(vcf,bb)
dev.off()

T <- table(mt$V$CLS)
s <- .histBeta(bb)
s2 <- .histBeta(bb, time = "time.2nd")
y0 <- seq(0.005, 0.995, 0.01)
h = y0[getMode(s)]
h2 = y0[getMode(s2)]

T$T1.histBeta <- h  
T$T2.histBeta <- h2  
write.table(T, file='M6B_Results_MutationTimeR.tsv', quote=FALSE, sep='\t', row.names = FALSE)

vcf <- readVcf("M7_MutationTimeR.vcf")
cn <- read.csv("M7_CN.txt", sep="\t")
bb <- makeGRangesFromDataFrame(cn,
                         keep.extra.columns=TRUE,
                         ignore.strand=FALSE,
                         seqinfo=NULL,
                         seqnames.field=c("seqnames", "seqname",
                                          "chromosome", "chrom",
                                          "chr", "chromosome_name",
                                          "seqid"),
                         start.field="start",
                         end.field=c("end", "stop"),
                         strand.field="strand",
                         starts.in.df.are.0based=FALSE)
clusters <- data.frame(cluster=1:4, proportion=c(0.9913, 0.1323, 0.0808, 0.52), n_ssms=c(394, 260, 513, 264))
mt <- mutationTime(vcf, bb, clusters=clusters, n.boot=10)
vcf <- addMutTime(vcf, mt$V)
mcols(bb) <- cbind(mcols(bb),mt$T)

pdf("M7_MutationTimeR_JBH.pdf", 10, 8)
plotSampleJBH(vcf,bb)
dev.off()

T <- table(mt$V$CLS)
s <- .histBeta(bb)
s2 <- .histBeta(bb, time = "time.2nd")
y0 <- seq(0.005, 0.995, 0.01)
h = y0[getMode(s)]
h2 = y0[getMode(s2)]

T$T1.histBeta <- h  
T$T2.histBeta <- h2  
write.table(T, file='M7_Results_MutationTimeR.tsv', quote=FALSE, sep='\t', row.names = FALSE)

vcf <- readVcf("M8A_MutationTimeR.vcf")
cn <- read.csv("M8A_CN.txt", sep="\t")
bb <- makeGRangesFromDataFrame(cn,
                         keep.extra.columns=TRUE,
                         ignore.strand=FALSE,
                         seqinfo=NULL,
                         seqnames.field=c("seqnames", "seqname",
                                          "chromosome", "chrom",
                                          "chr", "chromosome_name",
                                          "seqid"),
                         start.field="start",
                         end.field=c("end", "stop"),
                         strand.field="strand",
                         starts.in.df.are.0based=FALSE)
clusters <- data.frame(cluster=1:5, proportion=c(0.9622, 0.1508, 0.1193, 0.0101, 0.5928), n_ssms=c(394, 260, 449, 245, 264))
mt <- mutationTime(vcf, bb, clusters=clusters, n.boot=10)
vcf <- addMutTime(vcf, mt$V)
mcols(bb) <- cbind(mcols(bb),mt$T)

pdf("M8A_MutationTimeR_JBH.pdf", 10, 8)
plotSampleJBH(vcf,bb)
dev.off()

T <- table(mt$V$CLS)
s <- .histBeta(bb)
s2 <- .histBeta(bb, time = "time.2nd")
y0 <- seq(0.005, 0.995, 0.01)
h = y0[getMode(s)]
h2 = y0[getMode(s2)]

T$T1.histBeta <- h  
T$T2.histBeta <- h2  
write.table(T, file='M8A_Results_MutationTimeR.tsv', quote=FALSE, sep='\t', row.names = FALSE)

vcf <- readVcf("M8B_MutationTimeR.vcf")
cn <- read.csv("M8B_CN.txt", sep="\t")
bb <- makeGRangesFromDataFrame(cn,
                         keep.extra.columns=TRUE,
                         ignore.strand=FALSE,
                         seqinfo=NULL,
                         seqnames.field=c("seqnames", "seqname",
                                          "chromosome", "chrom",
                                          "chr", "chromosome_name",
                                          "seqid"),
                         start.field="start",
                         end.field=c("end", "stop"),
                         strand.field="strand",
                         starts.in.df.are.0based=FALSE)
clusters <- data.frame(cluster=1:4, proportion=c(0.9761, 0.134, 0.126, 0.7951), n_ssms=c(394, 260, 449, 264))
mt <- mutationTime(vcf, bb, clusters=clusters, n.boot=10)
vcf <- addMutTime(vcf, mt$V)
mcols(bb) <- cbind(mcols(bb),mt$T)

pdf("M8B_MutationTimeR_JBH.pdf", 10, 8)
plotSampleJBH(vcf,bb)
dev.off()

T <- table(mt$V$CLS)
s <- .histBeta(bb)
s2 <- .histBeta(bb, time = "time.2nd")
y0 <- seq(0.005, 0.995, 0.01)
h = y0[getMode(s)]
h2 = y0[getMode(s2)]

T$T1.histBeta <- h  
T$T2.histBeta <- h2  
write.table(T, file='M8B_Results_MutationTimeR.tsv', quote=FALSE, sep='\t', row.names = FALSE)

vcf <- readVcf("M9_MutationTimeR.vcf")
cn <- read.csv("M9_CN.txt", sep="\t")
bb <- makeGRangesFromDataFrame(cn,
                         keep.extra.columns=TRUE,
                         ignore.strand=FALSE,
                         seqinfo=NULL,
                         seqnames.field=c("seqnames", "seqname",
                                          "chromosome", "chrom",
                                          "chr", "chromosome_name",
                                          "seqid"),
                         start.field="start",
                         end.field=c("end", "stop"),
                         strand.field="strand",
                         starts.in.df.are.0based=FALSE)
clusters <- data.frame(cluster=1:4, proportion=c(0.9682, 0.149, 0.1453, 0.8544), n_ssms=c(394, 260, 459, 264))
mt <- mutationTime(vcf, bb, clusters=clusters, n.boot=10)
vcf <- addMutTime(vcf, mt$V)
mcols(bb) <- cbind(mcols(bb),mt$T)

pdf("M9_MutationTimeR_JBH.pdf", 10, 8)
plotSampleJBH(vcf,bb)
dev.off()

T <- table(mt$V$CLS)
s <- .histBeta(bb)
s2 <- .histBeta(bb, time = "time.2nd")
y0 <- seq(0.005, 0.995, 0.01)
h = y0[getMode(s)]
h2 = y0[getMode(s2)]

T$T1.histBeta <- h  
T$T2.histBeta <- h2  
write.table(T, file='M9_Results_MutationTimeR.tsv', quote=FALSE, sep='\t', row.names = FALSE)

vcf <- readVcf("M10C_MutationTimeR.vcf")
cn <- read.csv("M10C_CN.txt", sep="\t")
bb <- makeGRangesFromDataFrame(cn,
                         keep.extra.columns=TRUE,
                         ignore.strand=FALSE,
                         seqinfo=NULL,
                         seqnames.field=c("seqnames", "seqname",
                                          "chromosome", "chrom",
                                          "chr", "chromosome_name",
                                          "seqid"),
                         start.field="start",
                         end.field=c("end", "stop"),
                         strand.field="strand",
                         starts.in.df.are.0based=FALSE)
clusters <- data.frame(cluster=1:6, proportion=c(0.9803, 0.1598, 0.3422, 0.0101, 0.0101, 0.8349), n_ssms=c(394, 260, 171, 851, 245, 264))
mt <- mutationTime(vcf, bb, clusters=clusters, n.boot=10)
vcf <- addMutTime(vcf, mt$V)
mcols(bb) <- cbind(mcols(bb),mt$T)

pdf("M10C_MutationTimeR_JBH.pdf", 10, 8)
plotSampleJBH(vcf,bb)
dev.off()

T <- table(mt$V$CLS)
s <- .histBeta(bb)
s2 <- .histBeta(bb, time = "time.2nd")
y0 <- seq(0.005, 0.995, 0.01)
h = y0[getMode(s)]
h2 = y0[getMode(s2)]

T$T1.histBeta <- h  
T$T2.histBeta <- h2  
write.table(T, file='M10C_Results_MutationTimeR.tsv', quote=FALSE, sep='\t', row.names = FALSE)

vcf <- readVcf("M10D_MutationTimeR.vcf")
cn <- read.csv("M10D_CN.txt", sep="\t")
bb <- makeGRangesFromDataFrame(cn,
                         keep.extra.columns=TRUE,
                         ignore.strand=FALSE,
                         seqinfo=NULL,
                         seqnames.field=c("seqnames", "seqname",
                                          "chromosome", "chrom",
                                          "chr", "chromosome_name",
                                          "seqid"),
                         start.field="start",
                         end.field=c("end", "stop"),
                         strand.field="strand",
                         starts.in.df.are.0based=FALSE)
clusters <- data.frame(cluster=1:4, proportion=c(0.9785, 0.1308, 0.3532, 0.7935), n_ssms=c(394, 260, 171, 264))
mt <- mutationTime(vcf, bb, clusters=clusters, n.boot=10)
vcf <- addMutTime(vcf, mt$V)
mcols(bb) <- cbind(mcols(bb),mt$T)

pdf("M10D_MutationTimeR_JBH.pdf", 10, 8)
plotSampleJBH(vcf,bb)
dev.off()

T <- table(mt$V$CLS)
s <- .histBeta(bb)
s2 <- .histBeta(bb, time = "time.2nd")
y0 <- seq(0.005, 0.995, 0.01)
h = y0[getMode(s)]
h2 = y0[getMode(s2)]

T$T1.histBeta <- h  
T$T2.histBeta <- h2  
write.table(T, file='M10D_Results_MutationTimeR.tsv', quote=FALSE, sep='\t', row.names = FALSE)

vcf <- readVcf("M10E_MutationTimeR.vcf")
cn <- read.csv("M10E_CN.txt", sep="\t")
bb <- makeGRangesFromDataFrame(cn,
                         keep.extra.columns=TRUE,
                         ignore.strand=FALSE,
                         seqinfo=NULL,
                         seqnames.field=c("seqnames", "seqname",
                                          "chromosome", "chrom",
                                          "chr", "chromosome_name",
                                          "seqid"),
                         start.field="start",
                         end.field=c("end", "stop"),
                         strand.field="strand",
                         starts.in.df.are.0based=FALSE)
clusters <- data.frame(cluster=1:4, proportion=c(0.9723, 0.1358, 0.263, 0.5898), n_ssms=c(394, 260, 171, 264))
mt <- mutationTime(vcf, bb, clusters=clusters, n.boot=10)
vcf <- addMutTime(vcf, mt$V)
mcols(bb) <- cbind(mcols(bb),mt$T)

pdf("M10E_MutationTimeR_JBH.pdf", 10, 8)
plotSampleJBH(vcf,bb)
dev.off()

T <- table(mt$V$CLS)
s <- .histBeta(bb)
s2 <- .histBeta(bb, time = "time.2nd")
y0 <- seq(0.005, 0.995, 0.01)
h = y0[getMode(s)]
h2 = y0[getMode(s2)]

T$T1.histBeta <- h  
T$T2.histBeta <- h2  
write.table(T, file='M10E_Results_MutationTimeR.tsv', quote=FALSE, sep='\t', row.names = FALSE)

vcf <- readVcf("M11A_MutationTimeR.vcf")
cn <- read.csv("M11A_CN.txt", sep="\t")
bb <- makeGRangesFromDataFrame(cn,
                         keep.extra.columns=TRUE,
                         ignore.strand=FALSE,
                         seqinfo=NULL,
                         seqnames.field=c("seqnames", "seqname",
                                          "chromosome", "chrom",
                                          "chr", "chromosome_name",
                                          "seqid"),
                         start.field="start",
                         end.field=c("end", "stop"),
                         strand.field="strand",
                         starts.in.df.are.0based=FALSE)
clusters <- data.frame(cluster=1:4, proportion=c(0.9882, 0.1253, 0.0303, 0.3517), n_ssms=c(394, 260, 851, 264))
mt <- mutationTime(vcf, bb, clusters=clusters, n.boot=10)
vcf <- addMutTime(vcf, mt$V)
mcols(bb) <- cbind(mcols(bb),mt$T)

pdf("M11A_MutationTimeR_JBH.pdf", 10, 8)
plotSampleJBH(vcf,bb)
dev.off()

T <- table(mt$V$CLS)
s <- .histBeta(bb)
s2 <- .histBeta(bb, time = "time.2nd")
y0 <- seq(0.005, 0.995, 0.01)
h = y0[getMode(s)]
h2 = y0[getMode(s2)]

T$T1.histBeta <- h  
T$T2.histBeta <- h2  
write.table(T, file='M11A_Results_MutationTimeR.tsv', quote=FALSE, sep='\t', row.names = FALSE)

vcf <- readVcf("M11B_MutationTimeR.vcf")
cn <- read.csv("M11B_CN.txt", sep="\t")
bb <- makeGRangesFromDataFrame(cn,
                         keep.extra.columns=TRUE,
                         ignore.strand=FALSE,
                         seqinfo=NULL,
                         seqnames.field=c("seqnames", "seqname",
                                          "chromosome", "chrom",
                                          "chr", "chromosome_name",
                                          "seqid"),
                         start.field="start",
                         end.field=c("end", "stop"),
                         strand.field="strand",
                         starts.in.df.are.0based=FALSE)
clusters <- data.frame(cluster=1:4, proportion=c(0.9949, 0.1316, 0.101, 0.2529), n_ssms=c(394, 260, 851, 264))
mt <- mutationTime(vcf, bb, clusters=clusters, n.boot=10)
vcf <- addMutTime(vcf, mt$V)
mcols(bb) <- cbind(mcols(bb),mt$T)

pdf("M11B_MutationTimeR_JBH.pdf", 10, 8)
plotSampleJBH(vcf,bb)
dev.off()

T <- table(mt$V$CLS)
s <- .histBeta(bb)
s2 <- .histBeta(bb, time = "time.2nd")
y0 <- seq(0.005, 0.995, 0.01)
h = y0[getMode(s)]
h2 = y0[getMode(s2)]

T$T1.histBeta <- h  
T$T2.histBeta <- h2  
write.table(T, file='M11B_Results_MutationTimeR.tsv', quote=FALSE, sep='\t', row.names = FALSE)

vcf <- readVcf("M12A_MutationTimeR.vcf")
cn <- read.csv("M12A_CN.txt", sep="\t")
bb <- makeGRangesFromDataFrame(cn,
                         keep.extra.columns=TRUE,
                         ignore.strand=FALSE,
                         seqinfo=NULL,
                         seqnames.field=c("seqnames", "seqname",
                                          "chromosome", "chrom",
                                          "chr", "chromosome_name",
                                          "seqid"),
                         start.field="start",
                         end.field=c("end", "stop"),
                         strand.field="strand",
                         starts.in.df.are.0based=FALSE)
clusters <- data.frame(cluster=1:5, proportion=c(0.9833, 0.1513, 0.0404, 0.0202, 0.6716), n_ssms=c(394, 260, 851, 245, 264))
mt <- mutationTime(vcf, bb, clusters=clusters, n.boot=10)
vcf <- addMutTime(vcf, mt$V)
mcols(bb) <- cbind(mcols(bb),mt$T)

pdf("M12A_MutationTimeR_JBH.pdf", 10, 8)
plotSampleJBH(vcf,bb)
dev.off()

T <- table(mt$V$CLS)
s <- .histBeta(bb)
s2 <- .histBeta(bb, time = "time.2nd")
y0 <- seq(0.005, 0.995, 0.01)
h = y0[getMode(s)]
h2 = y0[getMode(s2)]

T$T1.histBeta <- h  
T$T2.histBeta <- h2  
write.table(T, file='M12A_Results_MutationTimeR.tsv', quote=FALSE, sep='\t', row.names = FALSE)

vcf <- readVcf("M12B_MutationTimeR.vcf")
cn <- read.csv("M12B_CN.txt", sep="\t")
bb <- makeGRangesFromDataFrame(cn,
                         keep.extra.columns=TRUE,
                         ignore.strand=FALSE,
                         seqinfo=NULL,
                         seqnames.field=c("seqnames", "seqname",
                                          "chromosome", "chrom",
                                          "chr", "chromosome_name",
                                          "seqid"),
                         start.field="start",
                         end.field=c("end", "stop"),
                         strand.field="strand",
                         starts.in.df.are.0based=FALSE)
clusters <- data.frame(cluster=1:4, proportion=c(1, 0.1624, 0.0717, 0.527), n_ssms=c(394, 260, 449, 264))
mt <- mutationTime(vcf, bb, clusters=clusters, n.boot=10)
vcf <- addMutTime(vcf, mt$V)
mcols(bb) <- cbind(mcols(bb),mt$T)

pdf("M12B_MutationTimeR_JBH.pdf", 10, 8)
plotSampleJBH(vcf,bb)
dev.off()

T <- table(mt$V$CLS)
s <- .histBeta(bb)
s2 <- .histBeta(bb, time = "time.2nd")
y0 <- seq(0.005, 0.995, 0.01)
h = y0[getMode(s)]
h2 = y0[getMode(s2)]

T$T1.histBeta <- h  
T$T2.histBeta <- h2  
write.table(T, file='M12B_Results_MutationTimeR.tsv', quote=FALSE, sep='\t', row.names = FALSE)

vcf <- readVcf("M15_MutationTimeR.vcf")
cn <- read.csv("M15_CN.txt", sep="\t")
bb <- makeGRangesFromDataFrame(cn,
                         keep.extra.columns=TRUE,
                         ignore.strand=FALSE,
                         seqinfo=NULL,
                         seqnames.field=c("seqnames", "seqname",
                                          "chromosome", "chrom",
                                          "chr", "chromosome_name",
                                          "seqid"),
                         start.field="start",
                         end.field=c("end", "stop"),
                         strand.field="strand",
                         starts.in.df.are.0based=FALSE)
clusters <- data.frame(cluster=1:4, proportion=c(0.9816, 0.1325, 0.1111, 0.1449), n_ssms=c(394, 260, 851, 264))
mt <- mutationTime(vcf, bb, clusters=clusters, n.boot=10)
vcf <- addMutTime(vcf, mt$V)
mcols(bb) <- cbind(mcols(bb),mt$T)

pdf("M15_MutationTimeR_JBH.pdf", 10, 8)
plotSampleJBH(vcf,bb)
dev.off()

T <- table(mt$V$CLS)
s <- .histBeta(bb)
s2 <- .histBeta(bb, time = "time.2nd")
y0 <- seq(0.005, 0.995, 0.01)
h = y0[getMode(s)]
h2 = y0[getMode(s2)]

T$T1.histBeta <- h  
T$T2.histBeta <- h2  
write.table(T, file='M15_Results_MutationTimeR.tsv', quote=FALSE, sep='\t', row.names = FALSE)

vcf <- readVcf("M16_MutationTimeR.vcf")
cn <- read.csv("M16_CN.txt", sep="\t")
bb <- makeGRangesFromDataFrame(cn,
                         keep.extra.columns=TRUE,
                         ignore.strand=FALSE,
                         seqinfo=NULL,
                         seqnames.field=c("seqnames", "seqname",
                                          "chromosome", "chrom",
                                          "chr", "chromosome_name",
                                          "seqid"),
                         start.field="start",
                         end.field=c("end", "stop"),
                         strand.field="strand",
                         starts.in.df.are.0based=FALSE)
clusters <- data.frame(cluster=1:4, proportion=c(0.9721, 0.1424, 0.1956, 0.2403), n_ssms=c(394, 260, 513, 264))
mt <- mutationTime(vcf, bb, clusters=clusters, n.boot=10)
vcf <- addMutTime(vcf, mt$V)
mcols(bb) <- cbind(mcols(bb),mt$T)

pdf("M16_MutationTimeR_JBH.pdf", 10, 8)
plotSampleJBH(vcf,bb)
dev.off()

T <- table(mt$V$CLS)
s <- .histBeta(bb)
s2 <- .histBeta(bb, time = "time.2nd")
y0 <- seq(0.005, 0.995, 0.01)
h = y0[getMode(s)]
h2 = y0[getMode(s2)]

T$T1.histBeta <- h  
T$T2.histBeta <- h2  
write.table(T, file='M16_Results_MutationTimeR.tsv', quote=FALSE, sep='\t', row.names = FALSE)

vcf <- readVcf("M17_MutationTimeR.vcf")
cn <- read.csv("M17_CN.txt", sep="\t")
bb <- makeGRangesFromDataFrame(cn,
                         keep.extra.columns=TRUE,
                         ignore.strand=FALSE,
                         seqinfo=NULL,
                         seqnames.field=c("seqnames", "seqname",
                                          "chromosome", "chrom",
                                          "chr", "chromosome_name",
                                          "seqid"),
                         start.field="start",
                         end.field=c("end", "stop"),
                         strand.field="strand",
                         starts.in.df.are.0based=FALSE)
clusters <- data.frame(cluster=1:3, proportion=c(0.9847, 0.0909, 0.2626), n_ssms=c(394, 260, 825))
mt <- mutationTime(vcf, bb, clusters=clusters, n.boot=10)
vcf <- addMutTime(vcf, mt$V)
mcols(bb) <- cbind(mcols(bb),mt$T)

pdf("M17_MutationTimeR_JBH.pdf", 10, 8)
plotSampleJBH(vcf,bb)
dev.off()

T <- table(mt$V$CLS)
s <- .histBeta(bb)
s2 <- .histBeta(bb, time = "time.2nd")
y0 <- seq(0.005, 0.995, 0.01)
h = y0[getMode(s)]
h2 = y0[getMode(s2)]

T$T1.histBeta <- h  
T$T2.histBeta <- h2  
write.table(T, file='M17_Results_MutationTimeR.tsv', quote=FALSE, sep='\t', row.names = FALSE)

vcf <- readVcf("M18_MutationTimeR.vcf")
cn <- read.csv("M18_CN.txt", sep="\t")
bb <- makeGRangesFromDataFrame(cn,
                         keep.extra.columns=TRUE,
                         ignore.strand=FALSE,
                         seqinfo=NULL,
                         seqnames.field=c("seqnames", "seqname",
                                          "chromosome", "chrom",
                                          "chr", "chromosome_name",
                                          "seqid"),
                         start.field="start",
                         end.field=c("end", "stop"),
                         strand.field="strand",
                         starts.in.df.are.0based=FALSE)
clusters <- data.frame(cluster=1:3, proportion=c(0.9906, 0.1112, 0.2626), n_ssms=c(394, 260, 825))
mt <- mutationTime(vcf, bb, clusters=clusters, n.boot=10)
vcf <- addMutTime(vcf, mt$V)
mcols(bb) <- cbind(mcols(bb),mt$T)

pdf("M18_MutationTimeR_JBH.pdf", 10, 8)
plotSampleJBH(vcf,bb)
dev.off()

T <- table(mt$V$CLS)
s <- .histBeta(bb)
s2 <- .histBeta(bb, time = "time.2nd")
y0 <- seq(0.005, 0.995, 0.01)
h = y0[getMode(s)]
h2 = y0[getMode(s2)]

T$T1.histBeta <- h  
T$T2.histBeta <- h2  
write.table(T, file='M18_Results_MutationTimeR.tsv', quote=FALSE, sep='\t', row.names = FALSE)



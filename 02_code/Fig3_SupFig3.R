
rm(list=ls())

library(ggplot2)
library(pheatmap)
library(ggpubr)
library(dplyr)
library(ggthemes)
library(ggupset)
library(ggrepel)

IO<-read.delim(file="IO_signatures.txt", header = TRUE, sep = "\t", row.names = 1)
IO$Date <- as.Date(IO$Date, format = "%d-%b-%y")
IO
heatmap_sig<-IO[, c('Date','Tissue', 'Response', 'Prolif',
                    'Stroma',
                    'Lymphoid',
                    'Myeloid',
                    'Endothelial.Cells',
                    'APM',
                    'MHC2',
                    'IFN.Gamma',
                    'Cytotoxicity',
                    'Immunoproteasome',
                    'Apoptosis',
                    'Inflammatory.Chemokines',
                    'Hypoxia',
                    'MAGEs',
                    'Glycolytic.Activity',
                    'IFN.Downstream',
                    'Myeloid.Inflam',
                    'B.Cells',
                    'CD45',
                    'CD8.T.Cells',
                    'Cytotoxic.Cells',
                    'DC',
                    'Exhausted.CD8',
                    'Macrophages',
                    'Mast.Cells',
                    'Neutrophils',
                    'NK.CD56dim.Cells',
                    'NK.Cells',
                    'T.Cells',
                    'Th1.Cells',
                    'Treg',
                    'TIS',
                    'ARG1',
                    'NOS2',
                    'IDO1',
                    'PDL1',
                    'CTLA4',
                    'IL10',
                    'PDL2',
                    'B7.H3',
                    'TIGIT',
                    'TGF.Beta',
                    'PD1',
                    'MMR.loss',
                    'APM.loss',
                    'JAKSTAT.loss')]

heatmap_sig$Date <- format(as.Date(heatmap_sig$Date), "%d-%b-%y")
heatmap_sig
heatmap_sig[heatmap_sig < 0] <- 0

#remove unused rows
heatmap_sig <- heatmap_sig[-c(18),]

scores.scaled <- as.data.frame(apply(heatmap_sig[,4:49], 2, scale))
rownames(scores.scaled)<-rownames(heatmap_sig)
my_sample_col <- data.frame(Timing = rep(c("Sequential", "Parallel"), c(13,19)))
annotation_row<-data.frame(my_sample_col,heatmap_sig[,2])
colnames(annotation_row)<-c("Timing", "Parallel metastases" )
rownames(annotation_row)<-rownames(scores.scaled)
scores.scaled
annotation_row

anno_col = list(
  Timing = c(Sequential = "darkorchid4", Parallel = "darkgrey"),
  "Parallel metastases" = c(
    'lung' = '#1F78B4',
    'diaphragm' = "#E31A1C",
    'mediastinal_lymph_nodes' = "#FF7F00",
    'pericardium' = "#33A02C",
    'chest wall' = "purple4"))

pdf("heatmap_IOsignatures_unsupervised.pdf", useDingbats = FALSE, height = 8, width =11)
scores.scaled<-na.omit(scores.scaled)
pheatmap(as.matrix(scores.scaled),  row_order = FALSE, fontsize_col = 8,  fontsize_row = 8, 
         breaks = seq(min(t(scores.scaled)),abs(min(t(scores.scaled))),length.out = 101),
         legend_breaks = c(-3,-2,-1,0,1,2,3), annotation_row = annotation_row, 
         annotation_colors = anno_col, reverse_colorscale = FALSE, clustering_distance_rows = "canberra", 
         cellheight = 11, color = colorRampPalette(rev(brewer.pal(n = 12, name ="RdYlBu")))(100))
dev.off()

IO<-read.delim(file="IO_signatures.txt", header = TRUE, sep = "\t", row.names = 1)
IO$Date <- as.Date(IO$Date, format = "%d-%b-%y")


RESP= read.table(file="Response_duration.txt", header = TRUE, sep = "\t")
RESP$Start <- as.Date(RESP$Start, format = "%Y-%m-%d")
RESP$End <- as.Date(RESP$End, format = "%Y-%m-%d")

# Subset chest wall data
IO_thoracic<-subset(IO, Tissue %in% "chest wall")


# Plot longitudinals
a <-ggplot(data = IO_thoracic, mapping = aes(x = Date, y = APM))+
  scale_x_date(breaks=as.Date(c("2013-11-28","2015-01-28","2017-01-17","2017-07-04", "2018-02-02","2018-06-15")), date_labels = c("d373","d799","d1,519","d1,687","d1,900", "d2,033"))+
  geom_point(aes(y = APM, colour = "APM"))+
  stat_summary(fun = mean, geom = 'line', colour = "blue", aes(y = APM, group = 1))+
  geom_point(aes(y = IFN.Gamma, color = "IFN-gamma signaling"))+
  stat_summary(fun = mean, geom = 'line', colour = "black", aes(y = IFN.Gamma, group = 1))+
  geom_point(aes(y = Prolif, color = "Proliferation"))+
  stat_summary(fun = mean, geom = 'line', colour = "red", aes(y = Prolif, group = 1))+
  scale_linetype_manual(name ="", values = c('dashed')) +
  geom_rect(data = RESP[2:9,], aes(xmin = as.Date(Start,"%Y-%m-%d"), xmax = as.Date(End,"%Y-%m-%d"), fill = Response), 
            ymin = -Inf, ymax = Inf, alpha = 0.2,
            inherit.aes = FALSE, show.legend = FALSE)+
  scale_fill_manual(values = c("green","grey", "red", "orange", "yellow" )) +
  scale_color_manual(values = c(
    'APM' = 'blue',
    'IFN-gamma signaling' = 'black',
    'Proliferation' = 'red'),
    labels = c(
      'APM' = 'APM',
      'IFN-gamma signaling' = expression(IFN~gamma~ signaling),
      'Proliferation' = 'Tumor proliferation'
    )) +
  labs(color = 'Signatures scores', y = "Signatures scores")+
  geom_hline(aes(yintercept = 6.5, linetype= "TIS, primary TNBC (TCGA)"),colour = 'grey47')+
  theme_classic()+
  theme(text = element_text(family="Helvetica"), legend.position="right",axis.title.x=element_blank(),
        axis.text.x=element_blank(),legend.text.align = 0, plot.margin=margin(t = 1, b = 1, unit = "pt"))

b<-ggplot(data = IO_thoracic, mapping = aes(x = Date, y = CD8.T.Cells)) +
  scale_x_date(breaks=as.Date(c("2013-11-28","2015-01-28","2017-01-17","2017-07-04", "2018-02-02","2018-06-15")), date_labels = c("d373","d799","d1,519","d1,687","d1,900", "d2,033"))+
  geom_point(aes(y = CD8.T.Cells, color = "CD8 T Cells"))+
  geom_point(aes(y = Th1.Cells, color = "Th1 Cells"))+
  geom_point(aes(y = NK.Cells, color = "NK Cells"))+
  stat_summary(fun = mean, geom = 'line', colour = "darkred", aes(y = CD8.T.Cells, group = 1))+
  stat_summary(fun = mean, geom = 'line', colour = "darkorange", aes(y = Th1.Cells, group = 1))+
  stat_summary(fun = mean, geom = 'line', colour = "gray47", aes(y = NK.Cells, group = 1))+
  geom_rect(data = RESP[2:9,], aes(xmin = as.Date(Start,"%Y-%m-%d"), xmax = as.Date(End,"%Y-%m-%d"), fill = Response), 
            ymin = -Inf, ymax = Inf, alpha = 0.2,
            inherit.aes = FALSE, show.legend = FALSE)+
  scale_fill_manual(values = c("green","grey", "red", "orange", "yellow" ))+
  scale_color_manual(values = c(
    'CD8 T Cells' = 'darkred',
    'Th1 Cells' = 'darkorange',
    'NK Cells' = 'gray47'))+
  labs(color = 'CTL abundance', y = "CTL abundance scores")+
  theme_classic()+
  theme(text = element_text(family="Helvetica"),axis.title.x=element_blank(),legend.text.align = 0, plot.margin=margin(t = 1, b = 1, unit = "pt"))

# Facet longitudinal plots
figure3b <- ggarrange(a, b,
                    ncol = 1, nrow = 2, align = "hv") 
pdf("longitudinal3b.pdf", useDingbats = FALSE, height = 5, width = 8)
figure3b
dev.off()

IO<-read.delim(file="IO_signatures.txt", header = TRUE, sep = "\t", row.names = 1)

fit<-lm(IFN.Gamma~APM, data=IO)

summary(fit)

pdf("Immunophenotype.pdf", useDingbats = FALSE, height = 4, width = 7)
ggplot(IO, aes(x = IFN.Gamma , y = APM))+
  geom_point(aes(color = Response)) +
  geom_smooth(method = "lm")+
  geom_hline(yintercept = 4.25, colour = 'black', linetype='dotted')+
  geom_vline(xintercept = 3.70, colour = 'black', linetype='dotted')+
  scale_color_manual(values = c("black","green", "red", "orange", "yellow")) +
  theme_bw(base_size = 12)+
  geom_text_repel(aes(label = rownames(IO)), size = 3)+
  labs(x=expression(IFN~gamma~ signaling),y="APM", color = "Best response")+
  theme(text = element_text(family="Helvetica", size=10), axis.title = element_text(size = 10),
        axis.text.x = element_text(size=9), axis.text.y = element_text(size=9), legend.position = "top", 
        panel.grid = element_blank(),plot.margin=margin(t = 1, b = 1, unit = "pt"))+
  annotate("text", x = 6.2, y = 3.6, label = "paste(italic(p) < 2.2e-16)", colour="black", size = 4, parse=TRUE)
dev.off()

## Suppl Figure 3b

ihc_long <- read.csv('IHC_data.csv')

bp2 <- ggplot( aes(x=Status, y=value, fill=variable), data = ihc_long[ihc_long$variable == 'CD8_count' | ihc_long$variable == 'CD4_count', ]) +
  geom_boxplot() + theme_ipsum() +
  scale_fill_viridis(discrete = TRUE, alpha=0.6, labels=c('CD8', 'CD4'), option = 'magma') +
  geom_point(aes(fill = variable), position = position_jitterdodge()) + theme(
    axis.title.x =  element_text(size=16),
    axis.text.x  =  element_text(size=16),
    axis.title.y =  element_text(size=16),
    axis.text.y  =  element_text(size=16),
    legend.text = element_text(size = 16),
    legend.title = element_text(size=16)
  ) + xlab('Status') + ylab('Counts') + labs(fill = "Cell")
bp2

## Suppl Figure 3e

prot_f <- read.csv('Protein.csv')

ggplot(data = prot_f, aes(x=INF_gamma, y=MIC_AB, color=Response)) + geom_point(size=2)  + theme_ipsum() + geom_text_repel(aes(label = sample_id), show.legend = FALSE, size=10) + ylab('MIC-A/B') + xlab('INF-gamma') + scale_color_manual(values = col_vector[c(16,11:15)], name = "Response") + theme(axis.text.y = element_text(size = 18), axis.title.y=element_text(size=18,face="bold"), axis.text.x = element_text(size = 18), axis.title.x=element_text(size=18,face="bold"), legend.title = element_text(size=14), legend.text = element_text(size=14))

## Suppl Figure 3f

prot_final <- read.csv('Protein_final.csv')

ggplot( aes(y=INF_gamma, x=Status, fill=Status), data = prot_final) +
  geom_boxplot() + scale_fill_viridis(discrete = TRUE, alpha=0.6, option = 'magma') + theme_ipsum() +
  geom_jitter(color="black", size=1, alpha=0.9) +
  theme(
    legend.position="none",
    axis.title.x =  element_text(size=16),
    axis.text.x  =  element_text(size=16),
    axis.title.y =  element_text(size=16),
    axis.text.y  =  element_text(size=16),
    legend.text = element_text(size = 16),
    legend.title = element_text(size=16)
  ) + xlab('Status') + ylab('INF gamma')


```

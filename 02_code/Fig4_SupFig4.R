library(RColorBrewer)
library(ggplot2)
library(DESeq2)
library(dendextend)
library(lsa)

## Fig 4a

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

pd <- read.csv('ex_score_tnbc.csv')

pdf('Figure_4a.pdf', width = 16, height = 8, fillOddEven = TRUE)
barplot(-1*pd$rank[order(-1*pd$rank, decreasing = T)], col = pd$site_cols[order(-1*pd$rank, decreasing = T)], ylim = c(-1.5,2.5), space=0, border=NA, names=pd$SampleID[order(-1*pd$rank, decreasing = T)], las=2, cex.axis = 1.5, cex.names = 0.7, ylab = 'Exhaustion score')
legend('topright', legend = c('Primary', 'Metastases Siegal', '302-M1', 'Metastases 302'), col = c('purple', 'gray', 'firebrick3', 'gray37'), pch=19, bty='n')
dev.off()

## Fig 4b

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

## SupFig 4a

df_merged <- read.csv('cyt_score_rank.csv')

ggplot(data = df_merged, aes(x=rank*-1, y=cyt_score, size=tCD8*100)) + geom_point(aes(color=site_cols)) + ylab('Cytolitic score') + xlab('Metagene rank') + labs(size='CD8+ cells (%)', color='Samples') + scale_color_manual(values=c("purple","gray", "firebrick3", "gray37"))  + theme_bw(base_size = 16) + geom_hline(yintercept=mean(df_merged$cyt_score),linetype=2, color='red') + geom_vline(xintercept=0,linetype=2, color='red') + annotate('text', x=c(-1, 0.5, -1, 0.5), y=c(6.5, 6.5, 9.5, 9.5), label=c('Desert', 'Exhausted', 'Inflamed', 'Inflamed/pre-exhausted'), size=5, hjust = 0, color=c('blue','blue','red','red')) + guides(color = guide_legend(override.aes = list(size = 4)))

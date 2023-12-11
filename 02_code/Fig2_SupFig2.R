library(Seurat) 
library(ggplot2) 
library(dplyr)
library(data.table)
library(grid)
library(RColorBrewer)

#### 1. Load the samples. ####

data <- Read10X(data.dir = "~/filtered_feature_bc_matrix/")
Seurat.obj <- CreateSeuratObject(counts = data, project = "")
Seurat.obj
#Add metadata
mito.genes <- grep(pattern = "^MT-", x = rownames(x = Seurat.obj), value = TRUE)
percent.counts <- Matrix::colSums(Seurat.obj@assays$RNA@data[mito.genes,])
percent.mito <- Matrix::colSums(Seurat.obj@assays$RNA@data[mito.genes,])/Matrix::colSums(Seurat.obj@assays$RNA@data)
# ribo <- fread("ribosomal_genes.tab")
ribo <- grep(pattern = "^RP[SL][[:digit:]]|^RP[[:digit:]]|^RPSA",
             rownames(Seurat.obj),
             value=TRUE, invert=F) #get list of non-ribosomal genes 
ribo.genes <- rownames(Seurat.obj)[which(rownames(Seurat.obj) %in% ribo)]
percent.ribo <- Matrix::colSums(Seurat.obj@assays$RNA@data[ribo.genes,])/Matrix::colSums(Seurat.obj@assays$RNA@data)
Seurat.obj <- AddMetaData(object = Seurat.obj, metadata = percent.ribo, col.name = "percent.ribo")
Seurat.obj <- AddMetaData(object = Seurat.obj, metadata = percent.mito, col.name = "percent.mito")
Seurat.obj <- AddMetaData(object = Seurat.obj, metadata = percent.counts, col.name = "percent.counts")
Seurat.obj$log10_nCount_RNA <- log10(Seurat.obj$nCount_RNA)
Seurat.obj$Sample_id <- "Seurat.obj"

#### -> Supplementary Figure 2a <- ####

VlnPlot(object = Seurat.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mito", "percent.ribo"), ncol = 4)
upper_bound <- mean(Seurat.obj$log10_nCount_RNA) + 2*sd(Seurat.obj$log10_nCount_RNA)
lower_bound <- mean(Seurat.obj$log10_nCount_RNA) - 2*sd(Seurat.obj$log10_nCount_RNA)
p1 <- FeatureScatter(object = Seurat.obj, feature1 = "log10_nCount_RNA", feature2 = "percent.mito")
p2 <- FeatureScatter(object = Seurat.obj, feature1 = "log10_nCount_RNA", feature2 = "nFeature_RNA")
p1 + geom_vline(xintercept=lower_bound) + geom_vline(xintercept=3) + geom_vline(xintercept=4.6) + geom_hline(yintercept=0.2)  + geom_hline(yintercept=0.1)
p2 + geom_vline(xintercept=lower_bound) + geom_hline(yintercept=200) + geom_hline(yintercept=5000)

# We filter out cells that have unique gene counts less than 200 and cells too expressed (could be doublets)
Seurat.obj <- subset(x = Seurat.obj, subset = nCount_RNA > 1000 & nCount_RNA < 40000 & nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mito < 0.1)
Seurat.obj

Seurat.obj$nGenes <- Seurat.obj$nFeature_RNA
Seurat.obj$nUMI <- Seurat.obj$nCount_RNA

#### 2. Normalize the data ####

Seurat.obj <- NormalizeData(object = Seurat.obj, normalization.method = "LogNormalize", scale.factor = 1e4)
Seurat.obj <- FindVariableFeatures(Seurat.obj, selection.method = "vst", nfeatures = 2000)
Seurat.obj <- ScaleData(Seurat.obj,vars.to.regress = c("nCount_RNA","percent.mito"), features = VariableFeatures(Seurat.obj))

#### 3. PCA and markers ####

Seurat.obj <- RunPCA(Seurat.obj, features = VariableFeatures(Seurat.obj), npcs = 40)
ElbowPlot(Seurat.obj,ndims = 40)

#Use first 9 PC
n_PC = 9
Seurat.obj <- FindNeighbors(Seurat.obj, dims = 1:n_PC)

Seurat.obj <- FindClusters(Seurat.obj, resolution = 0.5)
head(Idents(Seurat.obj), 5)

#PCA
DimPlot(Seurat.obj, reduction = "pca")

#### -> Figure 2c <- ###

#TSNE
n_cells = ncol(Seurat.obj)
grob <- grobTree(textGrob(paste0("n = ",n_cells), x=0.05,  y=0.95, hjust=0,gp=gpar(fontsize=22)))
getPalette = colorRampPalette(brewer.pal(11, "Set1"))
colors = getPalette(11)

#UMAP
Seurat.obj <- RunUMAP(Seurat.obj, dims = 1:n_PC)
# Plot
DimPlot(Seurat.obj, reduction = "umap",  pt.size = 1.2) 

# Get the markers from the unsup clusters
markers_all_cells <- FindAllMarkers(Seurat.obj, only.pos = TRUE, logfc.threshold = 0.25)

#### -> Figure 2e <- ###

#### 4. Incoroprate the TCR info ####

# Load the data generated of the TCR
TCR_dict <- Dict$new()
Clonotype_size_dict <- Dict$new()
Clonotype_id_dict <- Dict$new()
TCR_UMIs_dict <- Dict$new()

# Control
Control_TCR <- read.csv(file =  "~/filtered_contig_annotations.csv")
#Take valid cells with productive contigs
Control_TCR_valid <- Control_TCR[which(Control_TCR$full_length=="True" & Control_TCR$productive=="True"),]

#Nº cells with TRA or TRB
Control_TCR_valid_n <- length(unique(Control_TCR_valid$barcode))
Control_TCR_clonotype <- read.csv(file="~/clonotypes.csv")
Control_TCR_valid_clonotype <- merge(Control_TCR_valid,Control_TCR_clonotype,by.x="raw_clonotype_id",by.y="clonotype_id")
Control_TCR_valid_unique <- as.data.frame(table(Control_TCR_valid$barcode,Control_TCR_valid$chain))
Control_TCR_valid_unique2 <- Control_TCR_valid_unique[which(Control_TCR_valid_unique$Freq!=0),]

full_cont = 0
i <- 1
for(i in 1:nrow(Control_TCR_valid_unique2)){
  barcode_id <- paste0(substr(as.character(Control_TCR_valid_unique2[i,1]),1,nchar(as.character(Control_TCR_valid_unique2[i,1]))-2),"-1")
  tcr_info <- as.character(Control_TCR_valid_unique2[i,2])
  clone_size_info <- unique(Control_TCR_valid_clonotype$frequency[which(as.character(Control_TCR_valid_clonotype$barcode)%in%as.character(Control_TCR_valid_unique2[i,1]))])
  clonotype_id <- unique(Control_TCR_valid_clonotype$cdr3s_aa[which(as.character(Control_TCR_valid_clonotype$barcode)%in%as.character(Control_TCR_valid_unique2[i,1]))])
  UMIs <- unique(Control_TCR_valid_clonotype$umis[which(as.character(Control_TCR_valid_clonotype$barcode)%in%as.character(Control_TCR_valid_unique2[i,1]))])
  #Update the info related with this barcode
  #Check if the barcode is already in the dict
  if(TCR_dict$has(barcode_id)){
    TCR_dict$set(barcode_id,paste0(TCR_dict$peek(barcode_id),"/",tcr_info))
    full_cont = full_cont + 1
  } else{
    TCR_dict$add(barcode_id,tcr_info)
  }
  if(Clonotype_size_dict$has(barcode_id)){
    Clonotype_size_dict$set(barcode_id,clone_size_info)
  } else{
    Clonotype_size_dict$add(barcode_id,clone_size_info)
  }
  if(Clonotype_id_dict$has(barcode_id)){
    Clonotype_id_dict$set(barcode_id,clonotype_id)
  } else{
    Clonotype_id_dict$add(barcode_id,clonotype_id)
  }
  if(TCR_UMIs_dict$has(barcode_id)){
    TCR_UMIs_dict$set(barcode_id,UMIs)
  } else{
    TCR_UMIs_dict$add(barcode_id,UMIs)
  }
}
paste0("Number of cells with TCR (alpha or beta) reconstruction: ",i)
paste0("Number of cells with full TCR (alpha and beta) reconstruction: ",full_cont)

#Add the info to the Seurat object
n_cells <- length(colnames(Seurat.obj))

#For each barcode, check its status in TCR_dict
TCR_list <- NULL
Clono_size_list <- NULL
Clono_id_list <- NULL
TCR_UMIs_list <- NULL
i <- 1
cont <- 0
for(i in 1:ncol(Seurat.obj)){
  id <- colnames(Seurat.obj)[i]
  if(is.null(TCR_dict$peek(id))){
    TCR_list <- c(TCR_list,"No TCR") 
    Clono_size_list <- c(Clono_size_list,0)  
    Clono_id_list <- c(Clono_id_list,"No clonotype")  
    TCR_UMIs_list <- c(TCR_UMIs_list,0) 
  } else{
    # print(i)
    Clono_size <- Clonotype_size_dict$get(id)
    Clono_size_list <- c(Clono_size_list,Clono_size)
    Clono_id <- Clonotype_id_dict$get(id)
    Clono_id_list <- c(Clono_id_list,as.character(Clono_id))
    #Update TCR status according to Clono_id
    TCR_status <- ifelse(grepl("TRA",Clono_id),
                         ifelse(grepl("TRB",Clono_id),
                                "TRA/TRB",
                                "TRA"),
                         ifelse(grepl("TRB",Clono_id),
                                "TRB",
                                "No TCR"))
    # TCR_status <- TCR_dict$get(id)
    TCR_list <- c(TCR_list,TCR_status)
    TCR_UMIs_list <- c(TCR_UMIs_list,round(mean(TCR_UMIs_dict$get(id)),0)) 
    cont = cont + 1
  }
}

#Add this info to the seurat object
names(TCR_list) <- names(Seurat.obj$orig.ident)
Seurat.obj <- AddMetaData(object = Seurat.obj, metadata = TCR_list, col.name = "TCR_info")
names(Clono_size_list) <- names(Seurat.obj$orig.ident)
Seurat.obj <- AddMetaData(object = Seurat.obj, metadata = Clono_size_list, col.name = "Clono_size")
names(Clono_id_list) <- names(Seurat.obj$orig.ident)
Seurat.obj <- AddMetaData(object = Seurat.obj, metadata = Clono_id_list, col.name = "Clono_id")
names(TCR_UMIs_list) <- names(Seurat.obj$orig.ident)
Seurat.obj <- AddMetaData(object = Seurat.obj, metadata = TCR_UMIs_list, col.name = "TCR_UMIs")

actual_sizes <- table(Clono_id_list)
#Create another object with the actual sizes of the clonotypes present
Clono_size_list_updated <- Clono_size_list
# i <- 1
for(i in 1:length(Clono_size_list)){
  # print(i)
  id_cell <- names(Clono_size_list_updated[i])
  #Take the clonotype_id 
  clonotype_id <- Clono_id_list[which(names(Clono_id_list)%in%id_cell)]
  #Take the real size of each cell
  if(clonotype_id=="No clonotype"){
    new_size <- 0
  } else{
    new_size <- actual_sizes[which(names(actual_sizes)%in%clonotype_id)]
  }
  #Update the size
  Clono_size_list_updated[i] <- new_size
}

names(Clono_size_list_updated) <- names(Seurat.obj$orig.ident)
Seurat.obj <- AddMetaData(object = Seurat.obj, metadata = Clono_size_list_updated, col.name = "Clono_size_updated")

#### -> Figure 2f <- ###

# How many expanded T cells are from each sample?
n_cells = ncol(Seurat.obj)
Expanded <- rep("Non-TCR",n_cells)
Expanded[which(Seurat.obj$Clono_size>1)] <- "Expanded"
Expanded[which(Seurat.obj$Clono_size==1)] <- "Singleton"
table(Expanded)
names(Expanded) <- colnames(Seurat.obj)
Seurat.obj <- AddMetaData(object = Seurat.obj, metadata = Expanded, col.name = "Expanded")
as.matrix(table(Seurat.obj$Expanded,Idents(Seurat.obj)))
aux_table <- as.matrix(table(Seurat.obj$Expanded,Seurat.obj$Sample_id))

grob <- grobTree(textGrob(paste0("n = ",n_cells), x=0.05,  y=0.95, hjust=0,gp=gpar(fontsize=22)))
number_clusters <- length(unique(Idents(Seurat.obj)))
getPalette = colorRampPalette(brewer.pal(number_clusters, "Set1"))
colors = getPalette(number_clusters)
p1 <- DimPlot(Seurat.obj, reduction = "umap",  pt.size = 1.2) + annotation_custom(grob) + scale_color_manual(values = colors) + ggtitle("Clusters")
colors2 <- c("#b8bab9", '#0020f0', '#b152ff', '#ff3d17')
pt_sizes <- rep(0.6,n_cells)
pt_sizes[which(Seurat.obj$TCR_info!="No TCR")] <- 2
n_cells2 = length(which(Seurat.obj$TCR_info!="No TCR"))
grob2 <- grobTree(textGrob(paste0("n = ",n_cells2), x=0.05,  y=0.95, hjust=0,gp=gpar(fontsize=22)))
p2 <- DimPlot(Seurat.obj, reduction = "umap",  pt.size = pt_sizes, group.by = "TCR_info") + annotation_custom(grob2) + scale_color_manual(values = colors2)
# p3 <- DimPlot(Seurat.obj, reduction = "umap",  pt.size = pt_sizes, group.by = "Sample_id") + annotation_custom(grob)
getPalette = colorRampPalette(brewer.pal(32, "Set1"))
colors3 = getPalette(32)
colors3[1] <- "#e0e0e0"
pt_sizes <- rep(0.4,n_cells)
pt_sizes[which(Seurat.obj$Clono_size==1)] <- 0.6
pt_sizes[which(Seurat.obj$Clono_size>=2)] <- 0.8
pt_sizes[which(Seurat.obj$Clono_size>=10)] <- 1
pt_sizes[which(Seurat.obj$Clono_size>=50)] <- 1.2
pt_sizes[which(Seurat.obj$Clono_size>=100)] <- 1.4
pt_sizes[which(Seurat.obj$Clono_size>=200)] <- 1.6
p4 <- DimPlot(Seurat.obj, reduction = "umap",  group.by = "Clono_size",pt.size = pt_sizes) + annotation_custom(grob) + scale_color_manual(values = colors3) + ggtitle("Clone size")

#Create another plot in which its colored by clonotype (only with the top expanded (clonal size >= 10))
#From all the clono ids, only assign colour to those ones which clono szie >= 10
#Associate a color to each clono id
getPalette = colorRampPalette(brewer.pal(19, "Paired"))
myColors <- getPalette(31)
Clono_id_aux <- Seurat.obj$Clono_id
Clono_id_aux[which(Seurat.obj$Clono_size<10)] <- "No clonotype"
Seurat.obj <- AddMetaData(object = Seurat.obj, metadata = Clono_id_aux, col.name = "Clono_id_aux")
pt_sizes <- rep(0.4,n_cells)
pt_sizes[which(Seurat.obj$Clono_size>=10)] <- 2
p5 <- DimPlot(Seurat.obj, reduction = "umap",  group.by = "Clono_id_aux",  pt.size = pt_sizes) + 
  # annotation_custom(grob) + 
  scale_color_manual(values = c("#c9c9c9",myColors)) +
  # theme(legend.position = "none") +
  ggtitle("Expanded Clones (>10)")

p1; p2; p4; p5;

#### -> Figure 2c <- ###

#### 5. Rename the clusters ####
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
colors = getPalette(12)

new.cluster.ids <- c("1 CD8 Tem", "2 CD8 gd", "3 CD4 Naive", "4 CD8 Naive", "5 CD8 Temra", "6", "7 CD4 CTL", "8 CD8 Interm c1-2", "9 CD8 Degraded")
names(new.cluster.ids) <- levels(Seurat.obj)
Seurat.obj <- RenameIdents(Seurat.obj, new.cluster.ids)
DimPlot(Seurat.obj, reduction = "umap",  pt.size = 1.2) + scale_color_manual(values = colors[c(1:6,10:12)])

#### Get the markers for all the new clusters
markers <- FindAllMarkers(Seurat.obj, only.pos = TRUE, logfc.threshold = 0.25)
table(markers$cluster)
top8 <- markers %>% group_by(cluster) %>% top_n(n = 8, wt = avg_log2FC)

#### -> Figure 2d <- ####

DoHeatmap(Seurat.obj, features = top8$gene, group.colors = colors[c(1:6,10:12)]) #+ NoLegend()

### -> Supplementary Figure 2c <- ###

DotPlot(Seurat.obj, features = unique(top8$gene)) + 
  RotatedAxis() +
  scale_colour_gradient2(low = "blue", mid = "#e3e3e3", high = "red") 

#### Take cluster 6 and recluster
Seurat.obj.6 <- subset(Seurat.obj, idents="6")
Seurat.obj.6 <- NormalizeData(object = Seurat.obj.6, normalization.method = "LogNormalize", scale.factor = 1e4)
Seurat.obj.6 <- FindVariableFeatures(Seurat.obj.6, selection.method = "vst", nfeatures = 2000)
Seurat.obj.6 <- ScaleData(Seurat.obj.6,vars.to.regress = c("nCount_RNA","percent.mito"))
Seurat.obj.6 <- RunPCA(Seurat.obj.6, features = VariableFeatures(Seurat.obj.6), npcs = 40)
ElbowPlot(Seurat.obj.6,ndims = 40)

#Use first 9 PC
n_PC = 9
Seurat.obj.6 <- FindNeighbors(Seurat.obj.6, dims = 1:n_PC)
Seurat.obj.6 <- FindClusters(Seurat.obj.6, resolution = 0.5)

#PCA
DimPlot(Seurat.obj.6, reduction = "pca")

n_cells = ncol(Seurat.obj.6)
grob <- grobTree(textGrob(paste0("n = ",n_cells), x=0.05,  y=0.95, hjust=0,gp=gpar(fontsize=22)))

#UMAP
Seurat.obj.6 <- RunUMAP(Seurat.obj.6, dims = 1:n_PC)
# Plot
DimPlot(Seurat.obj.6, reduction = "umap",  pt.size = 1.2) + annotation_custom(grob) + scale_color_manual(values = colors[c(6:9)])
DimPlot(object = Seurat.obj.6, reduction = 'umap', pt.size = 1.2, group.by = "Phase")

markers_6 <- FindAllMarkers(Seurat.obj.6, only.pos = TRUE, logfc.threshold = 0.25)

new.cluster.ids <- c("6.0 CD4 Interm N-Eff", "6.1 Treg", "6.2 CD8 early Tem", "6.3 CD8 Prolif")
names(new.cluster.ids) <- levels(Seurat.obj.6)
Seurat.obj.6 <- RenameIdents(Seurat.obj.6, new.cluster.ids)

#### -> Figure 2d <- ####

DoHeatmap(Seurat.obj.6, features = top10$gene, group.colors = colors[c(6:9)]) #+ NoLegend()

DotPlot(Seurat.obj.6, features = unique(top10$gene)) + 
  RotatedAxis() +
  scale_colour_gradient2(low = "blue", mid = "#e3e3e3", high = "red") 

#### Give this clustering to the original dataset
old_clustering <- as.character(Idents(Seurat.obj))
old_clustering[which(names(Idents(Seurat.obj))%in%names(Idents(Seurat.obj.6)))] <- as.character(Idents(Seurat.obj.6))
names(old_clustering) <- names(Idents(Seurat.obj))
Seurat.obj <- AddMetaData(object = Seurat.obj, metadata = old_clustering, col.name = "new_clustering")

DimPlot(Seurat.obj, reduction = "umap",  pt.size = 1.2, group.by = "new_clustering", label = T) + scale_color_manual(values = colors) + NoAxes()

#### 6. Cross the scrna-seq info with the immunoSEQ one ####

rearrangementDetails <- read.csv(file="./immunoSEQ/RearrangementDetails_10-27-2021_8-41-20_AM.tsv",sep = "\t")
combinedRearrangements <- read.csv(file="./immunoSEQ/CombinedRearrangements_11-05-2021_10-00-55_AM.tsv",sep = "\t")

#Format the name ids
rearrangementDetails$sample_name2 <- rearrangementDetails$sample_name
rearrangementDetails$sample_name2[rearrangementDetails$sample_name =="TNBC302_M0_TCRB"]<-"01_M0"
rearrangementDetails$sample_name2[rearrangementDetails$sample_name =="TNBC302_M1_TCRB"]<-"02_M1"
rearrangementDetails$sample_name2[rearrangementDetails$sample_name =="TNBC302_M2_TCRB"]<-"03_M2"
rearrangementDetails$sample_name2[rearrangementDetails$sample_name =="TNBC302_M3_TCRB"]<-"04_M3"
rearrangementDetails$sample_name2[rearrangementDetails$sample_name =="TNBC302_M5A_TCRB"]<-"06_M5A"
rearrangementDetails$sample_name2[rearrangementDetails$sample_name =="TNBC302_M10E_TCRB"]<-"07_M10E"
rearrangementDetails$sample_name2[rearrangementDetails$sample_name =="TNBC302_M12A_TCRB"]<-"08_M12A"
rearrangementDetails$sample_name2[rearrangementDetails$sample_name =="TNBC302_M15_TCRB"]<-"09_M15"
rearrangementDetails$sample_name2[rearrangementDetails$sample_name =="TNBC302_M18_TCRB"]<-"10_M18"
colnames(combinedRearrangements)[4:12] <- c("01_M0","02_M1","03_M2","04_M3","06_M5A","07_M10E","08_M12A","09_M15","10_M18")

#Pile all the samples in common for each sample
aux2 <- rearrangementDetails %>% 
  group_by(amino_acid) %>% 
  mutate(Sample = paste0(unique(sample_name2), collapse = ","))
clones_samples_association <- unique(data.frame(cloneAA=aux2$amino_acid,Samples_associated=aux2$Sample))
#Merge the info back
combinedRearrangements2 <- merge(combinedRearrangements,clones_samples_association,by.x="Amino.Acid",by.y="cloneAA")

Idents(Seurat.obj) <- Seurat.obj$new_clustering

# Load the clonotypes info from the sc 
Control_TCR <- read.csv(file =  "~/filtered_contig_annotations.csv")

#Take valid cells with productive contigs
Control_TCR_valid <- Control_TCR[which(Control_TCR$full_length=="True" & Control_TCR$productive=="True"),]
#Separate between alpha and beta regions
Alpha_TCR_valid <- Control_TCR_valid[which(Control_TCR_valid$chain=="TRA"),]
Beta_TCR_valid <- Control_TCR_valid[which(Control_TCR_valid$chain=="TRB"),]

#Beta
beta_merge <- merge(Beta_TCR_valid,combinedRearrangements2,by.x="cdr3",by.y="Amino.Acid")
#How many clones?
cat(nrow(beta_merge)," total cells sharing beta region between sc and bulk")
cat(length(unique(beta_merge$cdr3))," unique clones")
#UMI distribution
ggplot(beta_merge,aes(umis)) +
  geom_bar(fill="darkgrey", color="black") +
  geom_vline(xintercept = 5, linetype='dashed')

#There are some cells with two productive TRB that mapps to two different TCR. Collapse this info
length(unique(beta_merge$barcode))
table_aux <- table(beta_merge$barcode)
# table_aux[which(table_aux>1)]

aux3 <- beta_merge %>% 
  group_by(barcode) %>% 
  mutate(new_associated_samples = paste0(unique(Samples_associated), collapse = ","))
aux4 <- unique(data.frame(barcode=aux3$barcode,Samples_associated=aux3$new_associated_samples))
#Remove repeated elements on Samples_associated
aux4$Samples_associated2 <- unlist(lapply(aux4$Samples_associated,function(x)paste0(sort(unique(strsplit(x,",")[[1]])), collapse=",")))

#How many cells in our Seurat.obj has beta region coincidence in the immunoSEQ data?
match2 <- rep("No match",length(colnames(Seurat.obj)))
names(match2) <- colnames(Seurat.obj)
match2 <- as.data.frame(match2)
merge_f <- merge(match2,aux4,by.x="row.names",by.y="barcode",all.x=TRUE)
merge_f$Samples_associated[which(is.na(merge_f$Samples_associated))] <- "No match"
merge_f$Samples_associated2[which(is.na(merge_f$Samples_associated2))] <- "No match"
match3 <- merge_f$Samples_associated2
names(match3) <- merge_f$Row.names
Seurat.obj <- AddMetaData(object = Seurat.obj, metadata = match3, col.name = "match_immunoSEQ")
merge_f$Samples_associated2[which(grepl("\\,",merge_f$Samples_associated2))] <- "Multisample"
match4 <- merge_f$Samples_associated2
names(match4) <- merge_f$Row.names
Seurat.obj <- AddMetaData(object = Seurat.obj, metadata = match4, col.name = "match_immunoSEQ2")

#### -> Figure 2e <- ####

#Plot where these dots are in the UMAP
getPalette = colorRampPalette(brewer.pal(11, "Set1"))
n_cells = ncol(Seurat.obj)
pt_sizes <- rep(0.6,n_cells)
pt_sizes[which(Seurat.obj$match_immunoSEQ2!="No match")] <- 2
n_cells2 = length(which(Seurat.obj$match_immunoSEQ2!="No match"))
colors3 = getPalette(length(unique(Seurat.obj$match_immunoSEQ2)))
DimPlot(Seurat.obj, reduction = "umap",  pt.size = pt_sizes, group.by = "match_immunoSEQ2") + scale_color_manual(values = colors3) 
colors4 = getPalette(length(unique(Seurat.obj$match_immunoSEQ)))
DimPlot(Seurat.obj, reduction = "umap",  pt.size = pt_sizes, group.by = "match_immunoSEQ") + scale_color_manual(values = colors4) + NoLegend() 

#### -> Figure 2f <- ####

#Label them according to if they are longitudinal or autopsy
pos_match <- Seurat.obj$match_immunoSEQ
pos_long <- which(grepl("01_M0|02_M1|03_M2|04_M3",Seurat.obj$match_immunoSEQ))
pos_autopsy <- which(grepl("06_M5A|07_M10E|08_M12A|09_M15|10_M18",Seurat.obj$match_immunoSEQ))
pos_both <- intersect(pos_long,pos_autopsy)
pos_match[pos_long] <- "Longitudinal"
pos_match[pos_autopsy] <- "Autopsy"
pos_match[pos_both] <- "Both"
table(pos_match)
names(pos_match) <- colnames(Seurat.obj)
Seurat.obj <- AddMetaData(object = Seurat.obj, metadata = pos_match, col.name = "match_immunoSEQ3")
getPalette = colorRampPalette(brewer.pal(10, "Set1"))
colors = getPalette(4)
DimPlot(Seurat.obj, reduction = "umap",  pt.size = 1, group.by = "match_immunoSEQ3") + scale_color_manual(values = colors)

# Create a barplot from this information
getPalette = colorRampPalette(brewer.pal(10, "Set1"))
colors = getPalette(4)
meta.data <- Seurat.obj@meta.data
ggplot(meta.data, aes(x=new_clustering,fill=match_immunoSEQ3)) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  theme_classic() + 
  scale_fill_manual(values = colors) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

#### 7. Generate plots comparing the immunoSEQ data from pre and posmortem metastasis ####

#Format the name ids
rearrangementDetails$sample_name2 <- rearrangementDetails$sample_name
rearrangementDetails$sample_name2[rearrangementDetails$sample_name =="TNBC302_M0_TCRB"]<-"01_M0"
rearrangementDetails$sample_name2[rearrangementDetails$sample_name =="TNBC302_M1_TCRB"]<-"02_M1"
rearrangementDetails$sample_name2[rearrangementDetails$sample_name =="TNBC302_M2_TCRB"]<-"03_M2"
rearrangementDetails$sample_name2[rearrangementDetails$sample_name =="TNBC302_M3_TCRB"]<-"04_M3"
rearrangementDetails$sample_name2[rearrangementDetails$sample_name =="TNBC302_M5A_TCRB"]<-"06_M5A"
rearrangementDetails$sample_name2[rearrangementDetails$sample_name =="TNBC302_M10E_TCRB"]<-"07_M10E"
rearrangementDetails$sample_name2[rearrangementDetails$sample_name =="TNBC302_M12A_TCRB"]<-"08_M12A"
rearrangementDetails$sample_name2[rearrangementDetails$sample_name =="TNBC302_M15_TCRB"]<-"09_M15"
rearrangementDetails$sample_name2[rearrangementDetails$sample_name =="TNBC302_M18_TCRB"]<-"10_M18"
colnames(combinedRearrangements)[4:12] <- c("01_M0","02_M1","03_M2","04_M3","06_M5A","07_M10E","08_M12A","09_M15","10_M18")

#PIle all the samples in common for each sample
aux2 <- rearrangementDetails %>% 
  group_by(amino_acid) %>% 
  mutate(Sample = paste0(unique(sample_name2), collapse = ","))
clones_samples_association <- unique(data.frame(cloneAA=aux2$amino_acid,Samples_associated=aux2$Sample))
#Merge the info back
combinedRearrangements2 <- merge(combinedRearrangements,clones_samples_association,by.x="Amino.Acid",by.y="cloneAA")

# 3. Get back to the immunoSEQ info. Add the sc info
#Control_TCR <- read.csv(file =  "filtered_contig_annotations.csv")
#Take valid cells with productive contigs
Control_TCR_valid <- Control_TCR[which(Control_TCR$full_length=="True" & Control_TCR$productive=="True"),]
#Separate the beta regions
Beta_TCR_valid <- Control_TCR_valid[which(Control_TCR_valid$chain=="TRB"),]
#Get the size of each clone
cdr3_table <- as.data.frame(table(Beta_TCR_valid$cdr3))
colnames(cdr3_table) <- c("cdr3","clonal_size")
Beta_TCR_valid_merge <- merge(Beta_TCR_valid,cdr3_table,by="cdr3")
aux <- Beta_TCR_valid_merge %>%
  group_by(cdr3) %>%
  summarise(n = n()) %>%
  mutate(Proportion = n / sum(n))
Beta_TCR_valid_merge <- merge(Beta_TCR_valid_merge,aux,by="cdr3")

#Bud the table in the expected format for immunarch
scmetadata <- unique(data.frame(Clones=Beta_TCR_valid_merge$clonal_size,Proportion=Beta_TCR_valid_merge$Proportion,CDR3.nt=Beta_TCR_valid_merge$cdr3_nt,CDR3.aa=Beta_TCR_valid_merge$cdr3,V.name=Beta_TCR_valid_merge$v_gene,D.name=Beta_TCR_valid_merge$d_gene,J.name=Beta_TCR_valid_merge$j_gene,V.end=NA,D.start=NA,D.end=NA,J.start=NA,VJ.ins=NA,VD.ins=NA,DJ.ins=NA,Sequence=Beta_TCR_valid_merge$cdr3))
scmetadata <- scmetadata[order(scmetadata$Clones,decreasing = TRUE),]
rownames(scmetadata) <- NULL

#Load the immunoSEQ data
myData <- repLoad("~/immunoSEQ/samples/")
myData$data$sc <- scmetadata
myData$meta <- rbind(myData$meta,c("sc","sc","blood","2018","2018_sc"))

#Downsample the number of clonotypes to the smallest one 
myData_downsampled <- list()
myData_downsampled$data <- repSample(myData$data,  "sample")
# myData_downsampled$data <- repSample(myData$data,  "sample", .n = 1000)
myData_downsampled$meta <- myData$meta

#Nº of clones
exp_vol <- repExplore(myData_downsampled$data, .method = "volume")
vis(exp_vol, .by = c("Sample"), .meta = myData_downsampled$meta, .test = F, .points = FALSE)
# vis(exp_vol, .by = c("Status_Year"), .meta = myData_downsampled$meta, .test = F)

# Clonality
imm_top <- repClonality(myData_downsampled$data, .method = "top", .head = c(10, 100, 1000, 3000, 10000))
vis(imm_top, .meta = myData_downsampled$meta) 

# Overlap
imm_ov2 <- repOverlap(myData_downsampled$data, .method = "morisita", .verbose = F, .downsample=TRUE)
p2 <- vis(imm_ov2, .text.size = 4, .signif.digits = 1)
p2

#### -> Supplementary Figure 2f <- ####

#Sum all Mets into single one
Mets <- rbind(myData_downsampled$data$"06_M5A",myData_downsampled$data$"07_M10E",myData_downsampled$data$"08_M12A",myData_downsampled$data$"09_M15",myData_downsampled$data$"10_M18")
myData_downsampled2 <- myData_downsampled
myData_downsampled2$data$"Mets" <- Mets
myData_downsampled2$data[[5]] <- NULL
myData_downsampled2$data[[5]] <- NULL
myData_downsampled2$data[[5]] <- NULL
myData_downsampled2$data[[5]] <- NULL
myData_downsampled2$data[[5]] <- NULL
imm_ov3 <- repOverlap(myData_downsampled2$data, .method = "morisita", .verbose = F, .downsample=TRUE)
p3 <- vis(imm_ov3, .text.size = 4, .signif.digits = 1)
p3

vis(imm_ov2, "heatmap2") 

ov <- repOverlap(myData_downsampled$data)
vis(ov, .plot = "circos")

ov2 <- repOverlap(myData_downsampled2$data)
vis(ov2, .plot = "circos")

# Gene usage analysis
imm_gu <- geneUsage(myData_downsampled$data, "hs.trbv", .norm = T)
imm_gu <- imm_gu[!is.na(imm_gu$Names),]

imm_gu_js <- geneUsageAnalysis(imm_gu, .method = "js", .verbose = F)
imm_gu_cor <- geneUsageAnalysis(imm_gu, .method = "cor", .verbose = F)

p1 <- vis(imm_gu_js, .title = "Gene usage JS-divergence", .leg.title = "JS", .text.size = 1.5)
p2 <- vis(imm_gu_cor, .title = "Gene usage correlation", .leg.title = "Cor", .text.size = 1.5)
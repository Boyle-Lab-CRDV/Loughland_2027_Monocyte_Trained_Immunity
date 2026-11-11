## 13th of Nov 2023

## IBSM40 monocytes analysis
## Author: Zuly P.


#.libPaths("/Users/zuleip/miniconda3/envs/scrnaseq_2/lib/R/library")


library(glmGamPoi)
library(tidyverse)
library(Seurat) 
library(patchwork)
library(sctransform)
library(limma)
library(data.table)
library(ape)
library(SC3)
library(clustree)
library(MAST)
library(UpSetR)

##Integrations using SCTtransform normalised data
##workflow https://hbctraining.github.io/scRNA-seq/lessons/06_SC_SCT_and_integration.html

#### Reading the rds file produce in the data_exploration_rds0.R script
ibsm.mono.F <- readRDS(file = "data/IBSM_mono_filtered_210823.rds")

#### NORMALISATION #########

#adjust the limit for allowable object sizes within R 
options(future.globals.maxSize = 4000 * 1024^2)

# splitting the RNA measurements into layers for each group. Day and donor samples
split_seurat <- SplitObject(ibsm.mono.F, split.by = "day_donor")

split_seurat <- split_seurat[c("day0_donor0", "day0_donor1", "day0_donor2", "day0_donor3","day8_donor0", "day8_donor1","day8_donor2", "day8_donor3","day36_donor0", "day36_donor1", "day36_donor2", "day36_donor3")]

for (i in 1:length(split_seurat)) {
   split_seurat[[i]] <- SCTransform(split_seurat[[i]], vars.to.regress = c("percent.mt"))
}

# Check which assays are stored in objects
split_seurat$day0_donor0@assays

# Select the most variable features to use for integration
integ_features <- SelectIntegrationFeatures(object.list = split_seurat, 
                                            nfeatures = 3000) 

# Prepare the SCT list object for integration
split_seurat <- PrepSCTIntegration(object.list = split_seurat, 
                                   anchor.features = integ_features)

# Find best buddies - can take a while to run
integ_anchors <- FindIntegrationAnchors(object.list = split_seurat, 
                                        normalization.method = "SCT", 
                                        anchor.features = integ_features)
####INTEGRATION###

# Integrate across conditions
immune.combined <- IntegrateData(anchorset = integ_anchors, 
                                   normalization.method = "SCT")

# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(immune.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)
ElbowPlot(immune.combined, ndims = 30)## 17 is good but we'll go with 20.
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:20)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:20)

##Fixing levels of some variables
day_levels <- c("day0", "day8", "day36")
immune.combined@meta.data$day <- factor(x = immune.combined@meta.data$day, levels = day_levels)
day_donor_levels <- c("day0_donor0", "day0_donor1", "day0_donor2", "day0_donor3","day8_donor0", "day8_donor1","day8_donor2", "day8_donor3","day36_donor0", "day36_donor1", "day36_donor2", "day36_donor3", "NA")
immune.combined@meta.data$day_donor <- factor(x = immune.combined@meta.data$day_donor, levels = day_donor_levels)
donor_levels <- c("donor0", "donor1","donor2","donor3")
immune.combined@meta.data$donor_id <- factor(x = immune.combined@meta.data$donor_id, levels = donor_levels)


####Elbow plot in numbers.####
# Determine percent of variation associated with each PC
pct <- immune.combined[["pca"]]@stdev / sum(immune.combined[["pca"]]@stdev) * 100
# Calculate cumulative percents for each PC
cumu <- cumsum(pct)
# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]
co1
# Determine the difference between variation of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
# last point where change of % of variation is more than 0.1%.
co2
# Minimum of the two calculation
pcs <- min(co1, co2)
pcs
# Create a dataframe with values
plot_df <- data.frame(pct = pct, 
                      cumu = cumu, 
                      rank = 1:length(pct))
# Elbow plot to visualize 
a <- ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) + 
  geom_text() + 
  geom_vline(xintercept = 90, color = "grey") + 
  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  theme_bw()
a

##Number of PC=17 going with 20 it's ok.
# Visualization
p1 <- DimPlot(immune.combined, reduction = "umap", group.by = "day_donor")
p2 <- DimPlot(immune.combined, reduction = "umap", label = TRUE, repel = TRUE)
p1 
p2

##Clusters
##Graph-based clustering approach using the distance metric previously identified (i.e. Number of PCs using Elbowplot)
##Algorithm for modularity optimization (1 = original Louvain algorithm; 2 = Louvain algorithm with multilevel refinement; 3 = SLM ##algorithm; 4 = Leiden algorithm). Leiden requires the leidenalg python.

#We inspected resolutions -0.1-0.9 finding 0.8 as optimal based on cluster tree.

#unsupervised clustering
immune.combined <- FindClusters(immune.combined, algorithm = 3, resolution = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1), verbose = TRUE)

# DimPlot(immune.combined, label = TRUE, pt.size = 2, group.by = "integrated_snn_res.0.1") + NoLegend() 
# DimPlot(immune.combined, label = TRUE, pt.size = 2, group.by = "integrated_snn_res.0.2") + NoLegend() 
# DimPlot(immune.combined, label = TRUE, pt.size = 2, group.by = "integrated_snn_res.0.3") + NoLegend() 
# DimPlot(immune.combined, label = TRUE, pt.size = 2, group.by = "integrated_snn_res.0.4") + NoLegend() 
# DimPlot(immune.combined, label = TRUE, pt.size = 2, group.by = "integrated_snn_res.0.5") + NoLegend() 
# DimPlot(immune.combined, label = TRUE, pt.size = 2, group.by = "integrated_snn_res.0.6") + NoLegend() 
# DimPlot(immune.combined, label = TRUE, pt.size = 2, group.by = "integrated_snn_res.0.7") + NoLegend() 
# DimPlot(immune.combined, label = TRUE, pt.size = 2, group.by = "integrated_snn_res.0.8") + NoLegend() 
# DimPlot(immune.combined, label = TRUE, pt.size = 2, group.by = "integrated_snn_res.0.9") + NoLegend() 
# DimPlot(immune.combined, label = TRUE, pt.size = 2, group.by = "integrated_snn_res.1") + NoLegend() 
# DimPlot(immune.combined, label = TRUE, pt.size = 2, group.by = "integrated_snn_res.1.9") + NoLegend() 

# ClusTree of reference editors. Each level of the tree (from top to bottom) corresponds to the k used, and each node on that level corresponds to a cluster of that k. The edges (arrows) of the tree represent the editors that "move" from one cluster in level k to another cluster in level k+1. The legend of the graph displays 4 scales, from top to bottom: (1) the transparency level of arrows (in_prop) shows the proportion editors from one group that end up in another group. (2) The arrow colour (count) shows the number of editors that "move" from one cluster to another, (3) the node size is proportional to the number
pdf("graphs/clustersSCTintdaydonor/clustree_SCTintDayDonor_Nov23.pdf", h=20, w=10)
clustree(immune.combined, prefix = "integrated_snn_res.", layout = "sugiyama")
dev.off()

pdf("graphs/clustersSCTintdaydonor/clustreeStability_SCTintDayDonor_Nov23.pdf", h=20, w=10)
clustree(immune.combined, prefix = "integrated_snn_res.", node_colour = "sc3_stability")
dev.off()
##Based on the tree we will explore resolution 0.8 (17 clusters) 


####STOP here and assess clustering.
# Select clustering resolution
Idents(immune.combined) <- "integrated_snn_res.0.8"

###QC clusters
pdf("graphs/clustersSCTintdaydonor/VlnPlot_clustersQC_Res0.8_Nov23.pdf", h=15, w=10)
VlnPlot(immune.combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), group.by = "integrated_snn_res.0.8", pt.size = 0, ncol = 4)
dev.off()
##How clusters look
pdf("graphs/clustersSCTintdaydonor/Dimplot_clusters_Res0.8_Nov23.pdf", h=15, w=10)
DimPlot(immune.combined, label = TRUE, pt.size = 2, group.by = "integrated_snn_res.0.8") + NoLegend() 
dev.off()
##How clusters look per day
pdf("graphs/clustersSCTintdaydonor/Dimplot_clustersbyDay_Res0.8_Nov23.pdf", h=15, w=10)
DimPlot(immune.combined, label = TRUE, pt.size = 2, group.by = "integrated_snn_res.0.8", split.by = "day") + NoLegend() 
dev.off()
##How clusters look per donor
pdf("graphs/clustersSCTintdaydonor/Dimplot_clustersbyDonor_Res0.8_Nov23.pdf", h=15, w=10)
DimPlot(immune.combined, label = TRUE, pt.size = 2, group.by = "integrated_snn_res.0.8", split.by = "donor_id") + NoLegend() 
dev.off()

meta <- immune.combined@meta.data
p <-  ggplot(data = meta, aes(x = integrated_snn_res.0.8, fill = day)) +
  geom_bar(position = "fill") +
  ylab("Proportion_day") + xlab("Cluster")
p + theme_light()
pdf("graphs/clustersSCTintdaydonor/Barplot_ClustersperDay_Res0.8_Nov23.pdf", h=15, w=10)
print(p)
dev.off()

d <-  ggplot(data = meta, aes(x = integrated_snn_res.0.8, fill = donor_id)) +
  geom_bar(position = "fill") +
  ylab("Proportion_donor") + xlab("Cluster")
d + theme_light()
pdf("graphs/clustersSCTintdaydonor/Barplot_ClustersperDonor_Res0.8_Nov23.pdf", h=15, w=10)
print(d)
dev.off()

##Data looks better for days and donor. Clusters are all present across days and 
#donors

##Cluster ID
# Run the standard workflow for visualization and clustering
#Following this workflow https://satijalab.org/seurat/articles/sctransform_v2_vignette
# We set up the default assay as "SCT" and run standard workflow
DefaultAssay(immune.combined) <- "SCT"

immune.combined <- SCTransform(immune.combined, vars.to.regress = "percent.mt", verbose = FALSE)
immune.combined <- RunPCA(immune.combined, features = VariableFeatures(object = immune.combined), npcs = 30, verbose = FALSE)
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:20)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:20)
ElbowPlot(immune.combined, ndims = 30)## we'll go with 30.

pdf("graphs/clustersSCTintdaydonor/Dotplot_monomarkers_Res0.8_Nov23.pdf", h=15, w=10)
DotPlot(immune.combined, assay = "SCT", features=c("CD14", "FCGR3A", "CD3E","CD4", "CCL3", "RETN", "CD8A", "CCR7", "TCF7", "CD27", "IL7R", "GZMA", "GZMB","CCL5", 
                                                    "FOXP3", "CTLA4", "HLA-DRA", "PPBP", "CLEC4C", "CLIC3","IL1B", "CD1C"),
        cols=c("white", "red3"))+theme_bw(base_size=14)+theme(axis.text.x = element_text(angle = 45,  vjust = 1, hjust=1))+theme (panel.grid.major = element_blank (), panel.grid.minor = element_blank ())
dev.off()

pdf("graphs/clustersSCTintdaydonor/Dotplot_monomarkers2_Res0.8_Nov23.pdf", h=15, w=10)
  DotPlot(immune.combined, assay = "SCT", features=c("CD14", "FCGR3A", "CX3CR1", "CDKN1C", "CD1C", "FCER1A", "CLEC4C", "NRP1","XCR1", "CLEC9A","CADM1","IL3RA"),
          cols=c("white", "red3"))+theme_bw(base_size=14)+theme(axis.text.x = element_text(angle = 45,  vjust = 1, hjust=1))+theme (panel.grid.major = element_blank (), panel.grid.minor = element_blank ())
dev.off()

##Main clusters IDs

##Clusters 0, 1, 2, 6, 7, 9, 10, 11, 12,13 and 14 classical monocytes
##Cluster 3 intermediate monocytes (CD14+, FCGR3A+, IL3RA)
##Cluster 4 and 8 non-classical monocytes (FCGR3A, CX3CR1, CDKN1C)
##Cluster 5 DCs (Dendritic cells 2) (CD1C+,FCER1A+)
##Cluster 15 DCs1 (Dendritic cells 1) (XCR1+,CLEC9A+,CADM1+)
##Cluster 16 pDCs (Plasmacytoids cells) (CLEC4C+, NRP1+,IL3RA+)


#This step takes a long time. To make it run faster just read the file
#cluster_markers <- FindAllMarkers(immune.combined)
cluster_markers <- read.csv(file = "graphs/clustersSCTintdaydonor/clustersmarkers.csv")
tapply(cluster_markers$gene, cluster_markers$cluster, head)
#write.csv(cluster_markers, file = "graphs/clustersSCTintdaydonor/clustersmarkers.csv")


##Classical monocytes 

#Cluster 0 - Classical monocyte - IL1B
#Cluster 1 - Classical monocyte - COQ7(Modulation of the amount of β2-integrins on the surface of blood monocytes)/HNRNPU (with longncRNA, enhance stability of several inflammatory mRNAs)
#Cluster 2 - Classical monocyte - CDKN1A
#Cluster 6 - Classical monocyte - DUSP1 (involved in the regulation of the cell cycle and apoptosis.)
#Cluster 7 - Classical monocyte - MALAT1 (regulates antigen uptake, processing, and phagocytosis in macrophages)
#Cluster 9 - Classical monocyte - S100A8
#Cluster 10 - Classical monocyte - ISG15 (controls cellular metabolism and mitochondrial activity.)
#Cluster 11 - Classical monocyte - HLA-DPB1/RPS18 Ribosomal proteins associated with synthesis pro-inflammatory
#Cluster 12 - Classical monocyte - DCs2-RPS18
#Cluster 13 - Classical monocyte - longnoncoding genes and HLa gens
#Cluster 14 - Classical monocyte - IL1RN

##Plotting similar genes between clusters using upset plot
upset0 <-  cluster_markers
upset0$clusterid <- factor(upset0$cluster,levels=c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12","13", "14", "15", "16"),
                           labels=c("cMonocyte-IL1B_0","cMonocyte-COQ7_1",
                                    "cMonocyte-CDKN1A_2","intMonocyte_3",
                                    "ncMonocyte-CD16_4", "DCs2_5", 
                                    "cMonocyte-DUSP1_6", "cMonocyte-MALAT_7",
                                    "ncMonocyte-CD16_8", "cMonocyte-S100A8_9",
                                    "cMonocyte-ISG15_10","DCs2-RPS18_11",
                                    "cMonocyte-mixed_12","cMonocyte-HLA_13",
                                    "cMonocyte-IL1RN_14", "DCs1_15", "pDCs_16"))

colnames(upset0)
upset1 <- upset0
##Creating new grouping variables
upset1$dir <- NA
upset1$dir[upset1$avg_log2FC>0] <- "up"
upset1$dir[upset1$avg_log2FC<0] <- "dw"
tapply(upset1$avg_log2FC, upset1$dir, summary)

##selecting only up genes
upset1 <- upset1[(upset1$dir[]=="up"),]
tapply(upset1$avg_log2FC, upset1$dir, summary)

upset1$cluster_dir <- NA
upset1$cluster_dir <- paste0(upset1$clusterid,"_",upset1$dir)

head(upset1)
colnames(upset1)
##creating lists of the sets
upset2 <- upset1[,c(8,11)]
head(upset2)

upset2Split = split(upset2[,1], f = upset2$cluster_dir)
head(upset2Split)

##Upset plot ordering by set frequency
pdf("graphs/clustersSCTintdaydonor/UpsetR_allclusters0.8_UPorderbyFreq.pdf", h=8, w =20)
upset(fromList(upset2Split),
      nsets = 17,
      nintersects = 100,
      order.by = "freq",
      decreasing = T,
      mb.ratio = c(0.5, 0.5),
      number.angles = 0,
      text.scale = 1.1,
      point.size = 2.8,
      line.size = 1)
dev.off()
##Note: there are more variations of the upset plot here "notused/cluster_id_postsctIntegration.R

##Labelling clusters
immune.combined@meta.data$subtype <- NA
immune.combined@meta.data$subtype[which(str_detect(immune.combined@meta.data$integrated_snn_res.0.8, "0"))] <- "cMonocyte-IL1B_0"
immune.combined@meta.data$subtype[which(str_detect(immune.combined@meta.data$integrated_snn_res.0.8, "1"))] <- "cMonocyte-COQ7_1"
immune.combined@meta.data$subtype[which(str_detect(immune.combined@meta.data$integrated_snn_res.0.8, "2"))] <- "cMonocyte-CDKN1A_2"
immune.combined@meta.data$subtype[which(str_detect(immune.combined@meta.data$integrated_snn_res.0.8, "3"))] <- "intMonocyte_3"
immune.combined@meta.data$subtype[which(str_detect(immune.combined@meta.data$integrated_snn_res.0.8, "4"))] <- "ncMonocyte-CD16_4"
immune.combined@meta.data$subtype[which(str_detect(immune.combined@meta.data$integrated_snn_res.0.8, "5"))] <- "DCs2_5"
immune.combined@meta.data$subtype[which(str_detect(immune.combined@meta.data$integrated_snn_res.0.8, "6"))] <- "cMonocyte-DUSP1_6"
immune.combined@meta.data$subtype[which(str_detect(immune.combined@meta.data$integrated_snn_res.0.8, "7"))] <- "cMonocyte-MALAT_7"
immune.combined@meta.data$subtype[which(str_detect(immune.combined@meta.data$integrated_snn_res.0.8, "8"))] <- "ncMonocyte-CD16_8"
immune.combined@meta.data$subtype[which(str_detect(immune.combined@meta.data$integrated_snn_res.0.8, "9"))] <- "cMonocyte-S100A8_9"
immune.combined@meta.data$subtype[which(str_detect(immune.combined@meta.data$integrated_snn_res.0.8, "10"))] <- "cMonocyte-ISG15_10"
immune.combined@meta.data$subtype[which(str_detect(immune.combined@meta.data$integrated_snn_res.0.8, "11"))] <- "DCs2-RPS18_11"
immune.combined@meta.data$subtype[which(str_detect(immune.combined@meta.data$integrated_snn_res.0.8, "12"))] <- "cMonocyte-_12"
immune.combined@meta.data$subtype[which(str_detect(immune.combined@meta.data$integrated_snn_res.0.8, "13"))] <- "cMonocyte-HLA_13"
immune.combined@meta.data$subtype[which(str_detect(immune.combined@meta.data$integrated_snn_res.0.8, "14"))] <- "cMonocyte-IL1RN_14"
immune.combined@meta.data$subtype[which(str_detect(immune.combined@meta.data$integrated_snn_res.0.8, "15"))] <- "DCs1_15"
immune.combined@meta.data$subtype[which(str_detect(immune.combined@meta.data$integrated_snn_res.0.8, "16"))] <- "pDCs_16"

pdf("graphs/clustersSCTintdaydonor/Dimplot_subtype_Res0.8_Nov23.pdf", h=15, w=10)
DimPlot(immune.combined, label = TRUE, pt.size = 2, group.by = "subtype") + NoLegend() 
dev.off()

##Regrouping cells based on the genes and directions that they shared (upsetplot) and the origin of the cluster (clustree)
immune.combined@meta.data$celltype <- NA
immune.combined@meta.data$celltype[which(str_detect(immune.combined@meta.data$integrated_snn_res.0.8, "0"))] <- "cMonocyte-IL_02614"
immune.combined@meta.data$celltype[which(str_detect(immune.combined@meta.data$integrated_snn_res.0.8, "1"))] <- "cMonocyte-COQ7_11213"
immune.combined@meta.data$celltype[which(str_detect(immune.combined@meta.data$integrated_snn_res.0.8, "2"))] <- "cMonocyte-IL_02614"
immune.combined@meta.data$celltype[which(str_detect(immune.combined@meta.data$integrated_snn_res.0.8, "3"))] <- "intMonocyte_3"
immune.combined@meta.data$celltype[which(str_detect(immune.combined@meta.data$integrated_snn_res.0.8, "4"))] <- "ncMonocyte-CD16_48"
immune.combined@meta.data$celltype[which(str_detect(immune.combined@meta.data$integrated_snn_res.0.8, "5"))] <- "DCs2_511"
immune.combined@meta.data$celltype[which(str_detect(immune.combined@meta.data$integrated_snn_res.0.8, "6"))] <- "cMonocyte-IL_02614"
immune.combined@meta.data$celltype[which(str_detect(immune.combined@meta.data$integrated_snn_res.0.8, "7"))] <- "cMonocyte-S100A8_79"
immune.combined@meta.data$celltype[which(str_detect(immune.combined@meta.data$integrated_snn_res.0.8, "8"))] <- "ncMonocyte-CD16_48"
immune.combined@meta.data$celltype[which(str_detect(immune.combined@meta.data$integrated_snn_res.0.8, "9"))] <- "cMonocyte-S100A8_79"
immune.combined@meta.data$celltype[which(str_detect(immune.combined@meta.data$integrated_snn_res.0.8, "10"))] <- "cMonocyte-ISG15_10"
immune.combined@meta.data$celltype[which(str_detect(immune.combined@meta.data$integrated_snn_res.0.8, "11"))] <- "DCs2_511"
immune.combined@meta.data$celltype[which(str_detect(immune.combined@meta.data$integrated_snn_res.0.8, "12"))] <- "cMonocyte-COQ7_11213"
immune.combined@meta.data$celltype[which(str_detect(immune.combined@meta.data$integrated_snn_res.0.8, "13"))] <- "cMonocyte-COQ7_11213"
immune.combined@meta.data$celltype[which(str_detect(immune.combined@meta.data$integrated_snn_res.0.8, "14"))] <- "cMonocyte-IL_02614"
immune.combined@meta.data$celltype[which(str_detect(immune.combined@meta.data$integrated_snn_res.0.8, "15"))] <- "DCs1_15"
immune.combined@meta.data$celltype[which(str_detect(immune.combined@meta.data$integrated_snn_res.0.8, "16"))] <- "pDCs_16"

# now save 
saveRDS(immune.combined, file = "data/IBSM_monoSCT_integbyDayDonor_220823.rds")

#Checking new clustering
DefaultAssay(immune.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
#immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)
ElbowPlot(immune.combined, ndims = 30)## 17 is good but we'll go with 20.
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:30)

Idents(immune.combined) <- "celltype"

pdf("graphs/clustersSCTintdaydonor/Dimplot_celltype_Res0.8_Nov23.pdf", h=15, w=10)
DimPlot(immune.combined, label = TRUE, pt.size = 2, group.by = "celltype") + NoLegend() 
dev.off()





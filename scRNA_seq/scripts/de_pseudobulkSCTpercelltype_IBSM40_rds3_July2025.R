#Author: Zuleima P
#Pseudobulk DE analysis using edgeR
#https://github.com/hbc/knowledgebase/blob/master/scrnaseq/pseudobulkDE_edgeR.md
#This scripts includes differential analysis using the pseudobluk DE analysis from edger

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("edgeR")

#.libPaths("/Users/zuleip/miniconda3/envs/scrnaseq_2/lib/R/library")

#library(scater)
library(Seurat)
library(tidyverse)
library(cowplot)
#library(Matrix.utils)#not installed
library(edgeR)
library(Matrix)
library(reshape2)
library(S4Vectors)
library(SingleCellExperiment)
library(pheatmap)
library(png)
library(RColorBrewer)
library(limma)
library(magrittr)
library(gridExtra)
library(knitr)
library(ggrepel)
library(data.table)
library(grr)

##Converting Seurat object to edgeR input file. 
##Note; We're going to perform two analysis the current one using major celltypes.
##The second one using subtypes. script de_pseudobulkSCT_IBSM40_rds2.R.


#### Day 0 vs Day8 #####
seurat <- readRDS( "/working_groups/boylelab/shared/ZulyPava/Training_paper_scRNAseq_R4.4/input/IBSM_monoSCT_integbyDayDonor_220823.rds")

#Creating a variable for de
seurat@meta.data$cluster_day <- paste(seurat@meta.data$day, seurat@meta.data$celltype, sep = "_")

#Set up identity to filter two days comparison
Idents(object = seurat) <- "day"

seurat_s <-subset(x = seurat, idents = c("day0", "day8"))

# Define an order of cluster identities
day_levels <- c("day0", "day8")
day_donor_levels <- c("day0_donor0", "day0_donor1", "day0_donor2", "day0_donor3","day8_donor0", "day8_donor1","day8_donor2", "day8_donor3")
# Re-level object@ident
seurat_s@meta.data$day <- factor(x = seurat_s@meta.data$day, levels = day_levels)
seurat_s@meta.data$day_donor <- factor(x = seurat_s@meta.data$day_donor, levels = day_donor_levels)

########Extract Counts and Metadata###############
# Extract raw counts and metadata to create SingleCellExperiment object
counts <- GetAssayData(object = seurat_s, slot = "counts", assay="SCT")
metadata <- seurat_s@meta.data

#Set up metadata as desired for aggregation and DE analysis
Idents(object = seurat_s) <- "celltype"
metadata$cluster_id <- factor(seurat_s@active.ident)

# Create single cell experiment object
sce <- SingleCellExperiment(assays = list(counts = counts),
                            colData = metadata)

## Explore the cellular metadata for the dataset
### (number of cells) * (number of meta columns)
dim(colData(sce)) 
#13672    33

#####Explore the counts data##########

#Explore the raw counts for the dataset
## Check the assays present (only counts)
assays(sce)

#Explore the raw counts for the dataset
### (genes) * (cells)
dim(counts(sce))  
# 19088 13672
counts(sce)[1:6, 1:6]

#####Additional QC filtering##########
## Remove lowly expressed genes which have less than 10 cells with any counts
#sce <- sce[rowSums(counts(sce) > 1) >= 10, ] A cleanup step was already done.

# (genes) x (cells)
dim(sce)

#####Preparing data for count aggregation##########
# Named vector of cluster names =
## 2
kids <- purrr::set_names(levels(sce$cluster_id))

# Total number of clusters
## 7
nk <- length(kids)
nk
# Named vector of sample names
sids <- purrr::set_names(levels(as.factor(sce$day_donor)))
sids
# Total number of samples = donor_id 
## 16
ns <- length(sids)
ns
##8

# Generate sample level metadata

## Determine the number of cells per sample
table(sce$day_donor)

#day0_donor0 day0_donor1 day0_donor2 day0_donor3 day8_donor0 day8_donor1 day8_donor2 day8_donor3 
#2425        1630        1122        1463        2794        2072        1160        1006 

## Turn class "table" into a named vector of cells per sample
n_cells <- table(sce$day_donor) %>%  as.vector()
names(n_cells) <- names(table(sce$day_donor))

## Match the named vector with metadata to combine it
m <- match(names(n_cells), sce$day_donor)

## Create the sample level metadata by selecting specific columns
ei <- data.frame(colData(sce)[m, ], 
                 n_cells, row.names = NULL) %>% 
  dplyr::select("day_donor", "day", "donor_id", "n_cells")
kable(ei)

#######Count aggregation ##########

# Aggregate the counts per sample_id and cluster_id

# Subset metadata to only include the cluster and sample IDs to aggregate across
groups <- colData(sce)[, c("cluster_id", "day_donor")]
groups$day_donor <- factor(groups$day_donor)

# Aggregate across cluster-sample groups
# Each row corresponds to aggregate counts for a cluster-sample combo
#Run the script agregatematrixfunction before continue with this part.
source("scripts/agregatematrixfunction.R")

pb <- aggregate.Matrix(t(counts(sce)), 
                       groupings = groups, fun = "sum") 

# class(pb)
# dim(pb)
pb[1:8, 1:8]

####Split/subsetting###########
##Here, we split the aggregated matrix into a matrix for each cluster. 
##In the table below, we report the total number of cells in each sample correposnding to each cluster

# create a vector that represents how to split samples
splitf <- sapply(stringr::str_split(rownames(pb), 
                                    pattern = "_day",n = 2), `[`, 1)

# Split data and turn into a list
# Each component corresponds to a cluster; storing associated expression matrix (counts)
# Transform data i.e, so rows are genes and columns are samples 
pb <- split.data.frame(pb,factor(splitf)) %>%
  lapply(function(u) 
    set_colnames(t(u), gsub(".*_da", "", rownames(u))))

# Explore the different components of list
class(pb)
str(pb)

###Print cluster-sample table

options(width = 100)
kable(table(sce$cluster_id, sce$day_donor))

# |                     | day0_donor0| day0_donor1| day0_donor2| day0_donor3| day8_donor0| day8_donor1| day8_donor2| day8_donor3|
#   |:--------------------|-----------:|-----------:|-----------:|-----------:|-----------:|-----------:|-----------:|-----------:|
#   |ncMonocyte-CD16_48   |         278|         198|         159|         180|         412|         223|         144|         165|
#   |DCs2_511             |         163|         181|         113|         178|         237|         207|          66|         151|
#   |cMonocyte-IL_02614   |        1171|         679|         467|         476|         935|         770|         331|         204|
#   |intMonocyte_3        |         234|          92|          76|         108|         310|         149|         166|          90|
#   |cMonocyte-ISG15_10   |          86|          41|          35|          30|          90|          55|          41|          20|
#   |cMonocyte-S100A8_79  |         148|         121|          90|         151|         222|         201|         137|         112|
#   |cMonocyte-COQ7_11213 |         326|         294|         154|         313|         557|         433|         267|         251|
#   |DCs1_15              |          13|          18|          22|          21|          27|          26|           5|           8|
#   |pDCs_16              |           6|           6|           6|           6|           4|           8|           3|           5|


######Starting EdgeR Analysis##########

######Subset Clusters ###############

### construct SCE of pseudo-bulk counts for only select clusters
# If you are interested in all clusters AND you have the same samples represented in each cluster you can just use pb

## Select the cluster you're interested?? no option to use pb directly without subsetting
# Check that every cluster is present in each of the samples
# > dim(pb$`cMonocyte-COQ7_11213`)
# [1] 19088     8
# > dim(pb$`cMonocyte-IL_02614`)
# [1] 19088     8
# > dim(pb$`cMonocyte-ISG15_10`)
# [1] 19088     8
# > dim(pb$`cMonocyte-S100A8_79`)
# [1] 19088     8
# > dim(pb$DCs1_15)
# [1] 19088     8
# > dim(pb$DCs2_511)
# dim(pb$pDCs_16)
# [1] 19088     8
# [1] 19088     8
# > dim(pb$intMonocyte_3)
# [1] 19088     8
# > dim(pb$`ncMonocyte-CD16_48`)
# [1] 19088     8

# If one of the cluster is absent from one of the samples then:
# Create a character vector of the clusters to use for DE
keepClusters <- as.character(c("cMonocyte-COQ7_11213", "cMonocyte-IL_02614", "cMonocyte-ISG15_10", "cMonocyte-S100A8_79", "pDCs_16", "DCs1_15", "DCs2_511", "intMonocyte_3", "ncMonocyte-CD16_48"))
# # Subset the sce object
pbc <- SingleCellExperiment(assays = pb[keepClusters])

######## MDS plots #######
# Similar to a PCA plot, the MDS allows us to visulize the distances between samples 
# for each cluster. Ideally we'd like "het" samples to segregate together and "ko" together.

# compute MDS coordinates
mds <-  lapply(as.list(assays(pbc)), function(a){
  DGEList(a, remove.zeros = TRUE) %>% 
    calcNormFactors %>% 
    plotMDS.DGEList(plot = FALSE)
})  
mds

# Add cluster names
cnames <- paste("Cluster", keepClusters)
for (m in 1:length(mds)){
  mds[[m]]$cluster <- cnames[m]
}

# prep. data.frame for plotting
plots <- lapply(mds, function(m){
  gg_df <<- data.frame(m[c("x", "y")],
                      sample_id = sids,
                      group_id = ei$day,
                      cluster_id = rep(m$cluster, length(m$x)))})

# Create a plotting function
plotFunc <- function(x) {
  ggplot(x, aes(x, y, col = group_id)) + 
    geom_point(size = 3, alpha = 0.8) +
    labs(x = "MDS dim. 1", y = "MDS dim. 2") + 
    ggtitle(unique(x$cluster_id)) + 
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(hjust = 0.5)) +
    coord_fixed() 
}

# Plot all MDS plots
mds_day_donor <- do.call(grid.arrange,c(lapply(plots, plotFunc)))
mds_day<- do.call(grid.arrange,c(lapply(plots, plotFunc)))

ggsave(plot = mds_day, filename = "output_for_paper/edger_graphs/MDSplots_Day8perday_07032024.pdf",
       width = 420, height= 297, units = "mm", device = "pdf" )

ggsave(plot = mds_day_donor, filename = "output_for_paper/edger_graphs/MDSplots_Day8_perday_donor_07032024.pdf",
       width = 420, height= 297, units = "mm", device = "pdf" )

##### Run differential expression analysis ########
#First set up the experiment design and a contrast matrix.
set.seed(12345)

# Build fixed effect design matrix (only for 'day')
design <- model.matrix(~ ei$day)
colnames(design) <- make.names(colnames(design))  # ensure valid names

# Define contrast (day8 vs day0)
contrast <- makeContrasts(ei.dayday8, levels = design)

# Run voom-limma workflow with donor as random effect
res <- lapply(keepClusters, function(k) {
  y <- assays(pbc)[[k]]
  y <- DGEList(y, remove.zeros = TRUE)
  y <- calcNormFactors(y)
  
  # Apply voom transformation
  v <- voom(y, design, plot = FALSE)
  
  # Estimate correlation due to repeated measures per donor
  dupcor <- duplicateCorrelation(v, design, block = ei$donor_id)
  
  # Fit linear model with correlation
  fit <- lmFit(v, design, block = ei$donor_id, correlation = dupcor$consensus.correlation)
  
  # Apply contrast
  fit2 <- contrasts.fit(fit, contrast)
  fit2 <- eBayes(fit2)
  
  # Extract results
  topTable(fit2, number = Inf, sort.by = "none") %>%
    dplyr::mutate(gene = rownames(.), cluster_id = k) %>%
    dplyr::rename(p_val = P.Value, p_adj = adj.P.Val)
})
# Results filtering & overview

# filter FDR < 0.05, |logFC| > 1 & sort by FDR
res_fil <- lapply(res, 
                  function(u)  u %>% 
                    dplyr::filter(p_adj < 0.05, abs(logFC) > 1) %>% 
                    dplyr::arrange(p_adj))

#######Significant genes ###########
# We filter the results with the criteria FDR < 0.05, and |logFC| > 1. 
# For each cluster we report the number of differentially expressed genes and what percentage that represents (of total genes).

## Count the number of differential findings by cluster.
# nb. & % of DE genes per cluster
n_de <- vapply(res_fil, nrow, numeric(1))
cbind(cluster=keepClusters, numDE_genes=n_de, 
      percentage = round(n_de / nrow(pbc) * 100, digits =2)) %>%  kable()

# |cluster              |numDE_genes |percentage |
#   |:--------------------|:-----------|:----------|
#   |cMonocyte-COQ7_11213 |133         |0.7        |
#   |cMonocyte-IL_02614   |154         |0.81       |
#   |cMonocyte-ISG15_10   |20          |0.1        |
#   |cMonocyte-S100A8_79  |49          |0.26       |
#   |pDCs_16              |9           |0.05       |
#   |DCs1_15              |4           |0.02       |
#   |DCs2_511             |41          |0.21       |
#   |intMonocyte_3        |64          |0.34       |
#   |ncMonocyte-CD16_48   |48          |0.25       |


##Write results to file

for(cluster in 1:length(keepClusters)){
  # Full results
  filePath <- paste0("output_for_paper/edger/", keepClusters[cluster])
  out <- res[[cluster]][,c("gene", "logFC", "AveExpr", "B", "p_val", "p_adj")]
  write.csv(out, file = paste0(filePath, "_", "D0D36_edgeR_pseudobulk_Allresults.csv"), quote=F, row.names = F)
  
  # Sig genes
  filePath <- paste0("output_for_paper/edger/", keepClusters[cluster])
  out <- res_fil[[cluster]][,c("gene", "logFC", "AveExpr", "B", "p_val", "p_adj")]
  write.csv(out, file = paste0(filePath, "_", "D0D36_edgeR_pseudobulk_sigGenes.csv"), quote=F, row.names = F)
  
}

rm(list = ls())
####Day 0 vs Day36 analysis ####

seurat <- readRDS( "/working_groups/boylelab/shared/ZulyPava/Training_paper_scRNAseq_R4.4/input/IBSM_monoSCT_integbyDayDonor_220823.rds")

#Creating a variable for de
seurat@meta.data$cluster_day <- paste(seurat@meta.data$day, seurat@meta.data$celltype, sep = "_")

#Set up identity to filter two days comparison
Idents(object = seurat) <- "day"

seurat_s <-subset(x = seurat, idents = c("day0", "day36"))

# Define an order of cluster identities
day_levels <- c("day0", "day36")
day_donor_levels <- c("day0_donor0", "day0_donor1", "day0_donor2", "day0_donor3","day36_donor0", "day36_donor1","day36_donor2", "day36_donor3")
# Re-level object@ident
seurat_s@meta.data$day <- factor(x = seurat_s@meta.data$day, levels = day_levels)
seurat_s@meta.data$day_donor <- factor(x = seurat_s@meta.data$day_donor, levels = day_donor_levels)

########Extract Counts and Metadata###############
# Extract raw counts and metadata to create SingleCellExperiment object
counts <- GetAssayData(object = seurat_s, slot = "counts", assay="SCT")
metadata <- seurat_s@meta.data

#Set up metadata as desired for aggregation and DE analysis
Idents(object = seurat_s) <- "celltype"
metadata$cluster_id <- factor(seurat_s@active.ident)

# Create single cell experiment object
sce <- SingleCellExperiment(assays = list(counts = counts),
                            colData = metadata)

## Explore the cellular metadata for the dataset
### (number of cells) * (number of meta columns)
dim(colData(sce)) 

##8603 33

#####Explore the counts data##########

#Explore the raw counts for the dataset
## Check the assays present (only counts)
assays(sce)

#Explore the raw counts for the dataset
### (genes) * (cells)
dim(counts(sce))  
#19088 8603

counts(sce)[1:6, 1:6]

#####Additional QC filtering##########
## Remove lowly expressed genes which have less than 10 cells with any counts
#sce <- sce[rowSums(counts(sce) > 1) >= 10, ] A cleanup step was already done.

# (genes) x (cells)
dim(sce)
#19088 8603

#####Preparing data for count aggregation##########
# Named vector of cluster names =
kids <- purrr::set_names(levels(sce$cluster_id))

# Total number of clusters
nk <- length(kids)
#9
# Named vector of sample names
sids <- purrr::set_names(levels(as.factor(sce$day_donor)))

# Total number of samples = donor_id 
ns <- length(sids)
#8

# Generate sample level metadata

## Determine the number of cells per sample
table(sce$day_donor)

# day0_donor0  day0_donor1  day0_donor2  day0_donor3 day36_donor0 day36_donor1 day36_donor2 
# 2425         1630         1122         1463          625          477          381 
# day36_donor3 
# 480 

## Turn class "table" into a named vector of cells per sample
n_cells <- table(sce$day_donor) %>%  as.vector()
names(n_cells) <- names(table(sce$day_donor))

## Match the named vector with metadata to combine it
m <- match(names(n_cells), sce$day_donor)

## Create the sample level metadata by selecting specific columns
ei <- data.frame(colData(sce)[m, ], 
                 n_cells, row.names = NULL) %>% 
  dplyr::select("day_donor", "day","donor_id", "n_cells")
kable(ei)

#######Count aggregation ##########

# Aggregate the counts per sample_id and cluster_id

# Subset metadata to only include the cluster and sample IDs to aggregate across
groups <- colData(sce)[, c("cluster_id", "day_donor")]
groups$day_donor <- factor(groups$day_donor)

# Aggregate across cluster-sample groups
# Each row corresponds to aggregate counts for a cluster-sample combo
#Run the script agregatematrixfunction before continue with this part.
source("scripts/agregatematrixfunction.R")

pb <- aggregate.Matrix(t(counts(sce)), 
                       groupings = groups, fun = "sum") 

# class(pb)
# dim(pb)
pb[1:8, 1:8]

####Split/subsetting###########
##Here, we split the aggregated matrix into a matrix for each cluster. 
##In the table below, we report the total number of cells in each sample correposnding to each cluster

# create a vector that represents how to split samples
splitf <- sapply(stringr::str_split(rownames(pb), 
                                    pattern = "_day",n = 2), `[`, 1)

# Split data and turn into a list
# Each component corresponds to a cluster; storing associated expression matrix (counts)
# Transform data i.e, so rows are genes and columns are samples 
pb <- split.data.frame(pb,factor(splitf)) %>%
  lapply(function(u) 
    set_colnames(t(u), gsub(".*_da", "", rownames(u))))

# Explore the different components of list
class(pb)
str(pb)

###Print cluster-sample table

options(width = 100)
kable(table(sce$cluster_id, sce$day_donor))


######Starting EdgeR Analysis##########

######Subset Clusters ###############

### construct SCE of pseudo-bulk counts for only select clusters
# If you are interested in all clusters AND you have the same samples represented in each cluster you can just use pb

#Note: We want all clusters and already remove cluster 8 at the beginning.

## Select the cluster you're interested?? no option to use pb directly without subsetting
# > dim(pb$`cMonocyte-COQ7_11213`)
# [1] 19088     8
# > dim(pb$`cMonocyte-IL_02614`)
# [1] 19088     8
# > dim(pb$`cMonocyte-ISG15_10`)
# [1] 19088     8
# > dim(pb$`cMonocyte-S100A8_79`)
# [1] 19088     8
# > dim(pb$DCs1_15)
# [1] 19088     8
# > dim(pb$DCs2_511)
# [1] 19088     8
# > dim(pb$pDCs_16)
# [1] 19088     8
# > dim(pb$intMonocyte_3)
# [1] 19088    8
# > dim(pb$`ncMonocyte-CD16_48`)
# [1] 19088     8

# Create a character vector of the clusters to use for DE
keepClusters <- as.character(c("cMonocyte-COQ7_11213", "cMonocyte-IL_02614", "cMonocyte-ISG15_10", "cMonocyte-S100A8_79", "pDCs_16", "DCs1_15", "DCs2_511", "intMonocyte_3", "ncMonocyte-CD16_48"))
# Subset the sce object
pbc <- SingleCellExperiment(assays = pb[keepClusters])

######## MDS plots #######
# Similar to a PCA plot, the MDS allows us to visulize the distances between samples 
# for each cluster. Ideally we'd like "het" samples to segregate together and "ko" together.

# compute MDS coordinates
mds <-  lapply(as.list(assays(pbc)), function(a){
  DGEList(a, remove.zeros = TRUE) %>% 
    calcNormFactors %>% 
    plotMDS.DGEList(plot = FALSE)
})  

# Add cluster names
cnames <- paste("Cluster", keepClusters)
for (m in 1:length(mds)){
  mds[[m]]$cluster <- cnames[m]
}

# prep. data.frame for plotting
plots <- lapply(mds, function(m){
  gg_df <<- data.frame(m[c("x", "y")],
                       sample_id = sids,
                       group_id = ei$day,
                       cluster_id = rep(m$cluster, length(m$x)))})

# Create a plotting function
plotFunc <- function(x) {
  ggplot(x, aes(x, y, col = group_id)) + 
    geom_point(size = 3, alpha = 0.8) +
    labs(x = "MDS dim. 1", y = "MDS dim. 2") + 
    ggtitle(unique(x$cluster_id)) + 
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(hjust = 0.5)) +
    coord_fixed() 
}

# Plot all MDS plots
mds_d36 <- do.call(grid.arrange,c(lapply(plots, plotFunc)))

ggsave(plot = mds_d36, filename = "output_for_paper/edger_graphs/MDSplots_d36_day_donor_070324.pdf",
       width = 420, height= 297, units = "mm", device = "pdf" )
ggsave(plot = mds_d36, filename = "output_for_paper/edger_graphs/MDSplots_d36_day_070324.pdf",
       width = 420, height= 297, units = "mm", device = "pdf" )

##### Run differential expression analysis ########
#First set up the experiment design and a contrast matrix.
set.seed(12345)

# Build fixed effect design matrix (only for 'day')
design <- model.matrix(~ ei$day)
colnames(design) <- make.names(colnames(design))  # ensure valid names

# Define contrast (day36 vs day0)
contrast <- makeContrasts(ei.dayday36, levels = design)

# Run voom-limma workflow with donor as random effect
res <- lapply(keepClusters, function(k) {
  y <- assays(pbc)[[k]]
  y <- DGEList(y, remove.zeros = TRUE)
  y <- calcNormFactors(y)
  
  # Apply voom transformation
  v <- voom(y, design, plot = FALSE)
  
  # Estimate correlation due to repeated measures per donor
  dupcor <- duplicateCorrelation(v, design, block = ei$donor_id)
  
  # Fit linear model with correlation
  fit <- lmFit(v, design, block = ei$donor_id, correlation = dupcor$consensus.correlation)
  
  # Apply contrast
  fit2 <- contrasts.fit(fit, contrast)
  fit2 <- eBayes(fit2)
  
  # Extract results
  topTable(fit2, number = Inf, sort.by = "none") %>%
    dplyr::mutate(gene = rownames(.), cluster_id = k) %>%
    dplyr::rename(p_val = P.Value, p_adj = adj.P.Val)
})
# Results filtering & overview

# filter FDR < 0.05, |logFC| > 1 & sort by FDR
res_fil <- lapply(res, 
                  function(u)  u %>% 
                    dplyr::filter(p_adj < 0.05, abs(logFC) > 1) %>% 
                    dplyr::arrange(p_adj))

#######Significant genes ###########
# We filter the results with the criteria FDR < 0.05, and |logFC| > 1. 
# For each cluster we report the number of differentially expressed genes and what percentage that represents (of total genes).

## Count the number of differential findings by cluster.
# nb. & % of DE genes per cluster
n_de <- vapply(res_fil, nrow, numeric(1))
cbind(cluster=keepClusters, numDE_genes=n_de, 
      percentage = round(n_de / nrow(pbc) * 100, digits =2)) %>%  kable()

#   |cluster              |numDE_genes |percentage |
#   |:--------------------|:-----------|:----------|
#   |cMonocyte-COQ7_11213 |266         |2.53       |
#   |cMonocyte-IL_02614   |381         |3.62       |
#   |cMonocyte-ISG15_10   |148         |1.41       |
#   |cMonocyte-S100A8_79  |165         |1.57       |
#   |pDCs_16              |52          |0.49       |
#   |DCs1_15              |36          |0.34       |
#   |DCs2_511             |97          |0.92       |
#   |intMonocyte_3        |143         |1.36       |
#   |ncMonocyte-CD16_48   |34          |0.32       |

#with random effects
# |cluster              |numDE_genes |percentage |
#   |:--------------------|:-----------|:----------|
#   |cMonocyte-COQ7_11213 |711         |3.72       |
#   |cMonocyte-IL_02614   |1560        |8.17       |
#   |cMonocyte-ISG15_10   |181         |0.95       |
#   |cMonocyte-S100A8_79  |234         |1.23       |
#   |pDCs_16              |57          |0.3        |
#   |DCs1_15              |333         |1.74       |
#   |DCs2_511             |752         |3.94       |
#   |intMonocyte_3        |332         |1.74       |
#   |ncMonocyte-CD16_48   |426         |2.23       |

##Write results to file

for(cluster in 1:length(keepClusters)){
  # Full results
  filePath <- paste0("output_for_paper/edger/", keepClusters[cluster])
  out <- res[[cluster]][,c("gene", "logFC", "logCPM", "F", "p_val", "p_adj")]
  write.csv(out, file = paste0(filePath, "_", "D0D36_edgeR_pseudobulk_Allresults.csv"), quote=F, row.names = F)
  
  # Sig genes
  filePath <- paste0("output_for_paper/edger/", keepClusters[cluster])
  out <- res_fil[[cluster]][,c("gene", "logFC", "logCPM", "F", "p_val", "p_adj")]
  write.csv(out, file = paste0(filePath, "_", "D0D36_edgeR_pseudobulk_sigGenes.csv"), quote=F, row.names = F)
  
}

##Write results to file

for(cluster in 1:length(keepClusters)){
  # Full results
  filePath <- paste0("output_for_paper/edger/", keepClusters[cluster])
  out <- res[[cluster]][,c("gene", "logFC", "AveExpr", "B", "p_val", "p_adj")]
  write.csv(out, file = paste0(filePath, "_", "D0D36_edgeR_pseudobulk_Allresults.csv"), quote=F, row.names = F)
  
  # Sig genes
  filePath <- paste0("output_for_paper/edger/", keepClusters[cluster])
  out <- res_fil[[cluster]][,c("gene", "logFC", "AveExpr", "B", "p_val", "p_adj")]
  write.csv(out, file = paste0(filePath, "_", "D0D36_edgeR_pseudobulk_sigGenes.csv"), quote=F, row.names = F)
  
}
rm(list = ls())
#### Day8vsDay36 analysis ####

seurat <- readRDS( "/working_groups/boylelab/shared/ZulyPava/Training_paper_scRNAseq_R4.4/input/IBSM_monoSCT_integbyDayDonor_220823.rds")

#Creating a variable for de
seurat@meta.data$cluster_day <- paste(seurat@meta.data$day, seurat@meta.data$celltype, sep = "_")

#Set up identity to filter two days comparison
Idents(object = seurat) <- "day"

seurat_s <-subset(x = seurat, idents = c("day8", "day36"))

# Define an order of cluster identities
day_levels <- c("day8", "day36")
day_donor_levels <- c("day8_donor0", "day8_donor1", "day8_donor2", "day8_donor3","day36_donor0", "day36_donor1","day36_donor2", "day36_donor3")
# Re-level object@ident
seurat_s@meta.data$day <- factor(x = seurat_s@meta.data$day, levels = day_levels)
seurat_s@meta.data$day_donor <- factor(x = seurat_s@meta.data$day_donor, levels = day_donor_levels)

########Extract Counts and Metadata###############
# Extract raw counts and metadata to create SingleCellExperiment object
counts <- GetAssayData(object = seurat_s, slot = "counts", assay="SCT")
metadata <- seurat_s@meta.data

#Set up metadata as desired for aggregation and DE analysis
Idents(object = seurat_s) <- "celltype"
metadata$cluster_id <- factor(seurat_s@active.ident)

# Create single cell experiment object
sce <- SingleCellExperiment(assays = list(counts = counts),
                            colData = metadata)

## Explore the cellular metadata for the dataset
### (number of cells) * (number of meta columns)
dim(colData(sce)) 

##8995 33

#####Explore the counts data##########

#Explore the raw counts for the dataset
## Check the assays present (only counts)
assays(sce)

#Explore the raw counts for the dataset
### (genes) * (cells)
dim(counts(sce))  
#19088 8995

counts(sce)[1:6, 1:6]

#####Additional QC filtering##########
## Remove lowly expressed genes which have less than 10 cells with any counts
#sce <- sce[rowSums(counts(sce) > 1) >= 10, ]

# (genes) x (cells)
dim(sce)
#19088  8995

#####Preparing data for count aggregation##########
# Named vector of cluster names =
## 10
kids <- purrr::set_names(levels(sce$cluster_id))

# Total number of clusters
## 10
nk <- length(kids)

# Named vector of sample names
sids <- purrr::set_names(levels(as.factor(sce$day_donor)))

# Total number of samples = donor_id 
## 12
ns <- length(sids)

# Generate sample level metadata

## Determine the number of cells per sample
table(sce$day_donor)

## Turn class "table" into a named vector of cells per sample
n_cells <- table(sce$day_donor) %>%  as.vector()
names(n_cells) <- names(table(sce$day_donor))

## Match the named vector with metadata to combine it
m <- match(names(n_cells), sce$day_donor)

## Create the sample level metadata by selecting specific columns
ei <- data.frame(colData(sce)[m, ], 
                 n_cells, row.names = NULL) %>% 
  dplyr::select("day_donor", "day", "donor_id", "n_cells")
kable(ei)

#######Count aggregation ##########

# Aggregate the counts per sample_id and cluster_id

# Subset metadata to only include the cluster and sample IDs to aggregate across
groups <- colData(sce)[, c("cluster_id", "day_donor")]
groups$day_donor <- factor(groups$day_donor)

# Aggregate across cluster-sample groups
# Each row corresponds to aggregate counts for a cluster-sample combo
#Run the script agregatematrixfunction before continue with this part.
source("scripts/agregatematrixfunction.R")

pb <- aggregate.Matrix(t(counts(sce)), 
                       groupings = groups, fun = "sum") 

# class(pb)
# dim(pb)
pb[1:8, 1:8]

####Split/subsetting###########
##Here, we split the aggregated matrix into a matrix for each cluster. 
##In the table below, we report the total number of cells in each sample correposnding to each cluster

# create a vector that represents how to split samples
splitf <- sapply(stringr::str_split(rownames(pb), 
                                    pattern = "_day",n = 2), `[`, 1)

# Split data and turn into a list
# Each component corresponds to a cluster; storing associated expression matrix (counts)
# Transform data i.e, so rows are genes and columns are samples 
pb <- split.data.frame(pb,factor(splitf)) %>%
  lapply(function(u) 
    set_colnames(t(u), gsub(".*_da", "", rownames(u))))

# Explore the different components of list
class(pb)
str(pb)

###Print cluster-sample table

options(width = 100)
kable(table(sce$cluster_id, sce$day_donor))


######Starting EdgeR Analysis##########

######Subset Clusters ###############

### construct SCE of pseudo-bulk counts for only select clusters
# If you are interested in all clusters AND you have the same samples represented in each cluster you can just use pb

#Note: We want all clusters and already remove cluster 8 at the beginning.

## Select the cluster you're interested?? no option to use pb directly without subsetting
# > dim(pb$`cMonocyte-COQ7_11213`)
# [1] 10165     8
# > dim(pb$`cMonocyte-IL_02614`)
# [1] 19088     8
# > dim(pb$`cMonocyte-ISG15_10`)
# [1] 19088     8
# > dim(pb$`cMonocyte-S100A8_79`)
# [1] 19088     8
# > dim(pb$DCs1_15)
# [1] 19088     8
# > dim(pb$DCs2_511)
# [1] 19088     8
# > dim(pb$pDCs_16)
# [1] 19088     8
# > dim(pb$intMonocyte_3)
# [1] 19088    8
# > dim(pb$`ncMonocyte-CD16_48`)
# [1] 19088     8

# Create a character vector of the clusters to use for DE
keepClusters <- as.character(c("cMonocyte-COQ7_11213", "cMonocyte-IL_02614", "cMonocyte-ISG15_10", "cMonocyte-S100A8_79", "pDCs_16", "DCs1_15", "DCs2_511", "intMonocyte_3", "ncMonocyte-CD16_48"))
# Subset the sce object
pbc <- SingleCellExperiment(assays = pb[keepClusters])

######## MDS plots #######
# Similar to a PCA plot, the MDS allows us to visulize the distances between samples 
# for each cluster. Ideally we'd like "het" samples to segregate together and "ko" together.

# compute MDS coordinates
mds <-  lapply(as.list(assays(pbc)), function(a){
  DGEList(a, remove.zeros = TRUE) %>% 
    calcNormFactors %>% 
    plotMDS.DGEList(plot = FALSE)
})  

# Add cluster names
cnames <- paste("Cluster", keepClusters)
for (m in 1:length(mds)){
  mds[[m]]$cluster <- cnames[m]
}

# prep. data.frame for plotting
plots <- lapply(mds, function(m){
  gg_df <<- data.frame(m[c("x", "y")],
                       sample_id = sids,
                       group_id = ei$day_donor,
                       cluster_id = rep(m$cluster, length(m$x)))})

# Create a plotting function
plotFunc <- function(x) {
  ggplot(x, aes(x, y, col = group_id)) + 
    geom_point(size = 3, alpha = 0.8) +
    labs(x = "MDS dim. 1", y = "MDS dim. 2") + 
    ggtitle(unique(x$cluster_id)) + 
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(hjust = 0.5)) +
    coord_fixed() 
}

# Plot all MDS plots
mds_d36 <- do.call(grid.arrange,c(lapply(plots, plotFunc)))

ggsave(plot = mds_d36, filename = "/working_groups/boylelab/shared/ZulyPava/Training_paper_scRNAseq_R4.4/output_for_paper/edger_graphs/MDSplots_d8d36_day_donor_070324.pdf",
       width = 420, height= 297, units = "mm", device = "pdf" )
ggsave(plot = mds_d36, filename = "/working_groups/boylelab/shared/ZulyPava/Training_paper_scRNAseq_R4.4/output_for_paper/edger_graphs/MDSplots_d8d36_day_070324.pdf",
       width = 420, height= 297, units = "mm", device = "pdf" )

##### Run differential expression analysis ########
#First set up the experiment design and a contrast matrix.
set.seed(12345)

# Ensure ei$day is a factor and set day36 as the reference
ei$day <- factor(ei$day, levels = c("day36", "day8"))  # include day36 first

# Build design matrix with day36 as reference
design <- model.matrix(~ei$day)
colnames(design) <- make.names(colnames(design))  # make valid column names

# To compare day8 vs day36 (i.e., day8 - day36)
contrast <- makeContrasts(ei.dayday8, levels = design)

# Run voom-limma workflow with donor as random effect
res <- lapply(keepClusters, function(k) {
  y <- assays(pbc)[[k]]
  y <- DGEList(y, remove.zeros = TRUE)
  y <- calcNormFactors(y)
  
  # Apply voom transformation
  v <- voom(y, design, plot = FALSE)
  
  # Estimate correlation due to repeated measures per donor
  dupcor <- duplicateCorrelation(v, design, block = ei$donor_id)
  
  # Fit linear model with correlation
  fit <- lmFit(v, design, block = ei$donor_id, correlation = dupcor$consensus.correlation)
  
  # Apply contrast
  fit2 <- contrasts.fit(fit, contrast)
  fit2 <- eBayes(fit2)
  
  # Extract results
  topTable(fit2, number = Inf, sort.by = "none") %>%
    dplyr::mutate(gene = rownames(.), cluster_id = k) %>%
    dplyr::rename(p_val = P.Value, p_adj = adj.P.Val)
})

# Results filtering & overview

# filter FDR < 0.05, |logFC| > 1 & sort by FDR
res_fil <- lapply(res, 
                  function(u)  u %>% 
                    dplyr::filter(p_adj < 0.05, abs(logFC) > 1) %>% 
                    dplyr::arrange(p_adj))

#######Significant genes ###########
# We filter the results with the criteria FDR < 0.05, and |logFC| > 1. 
# For each cluster we report the number of differentially expressed genes and what percentage that represents (of total genes).

## Count the number of differential findings by cluster.
# nb. & % of DE genes per cluster
n_de <- vapply(res_fil, nrow, numeric(1))
cbind(cluster=keepClusters, numDE_genes=n_de, 
      percentage = round(n_de / nrow(pbc) * 100, digits =2)) %>%  kable()

#   |cluster              |numDE_genes |percentage |
#   |:--------------------|:-----------|:----------|
#   |cMonocyte-COQ7_11213 |266         |2.53       |
#   |cMonocyte-IL_02614   |381         |3.62       |
#   |cMonocyte-ISG15_10   |148         |1.41       |
#   |cMonocyte-S100A8_79  |165         |1.57       |
#   |pDCs_16              |52          |0.49       |
#   |DCs1_15              |36          |0.34       |
#   |DCs2_511             |97          |0.92       |
#   |intMonocyte_3        |143         |1.36       |
#   |ncMonocyte-CD16_48   |34          |0.32       |

#with random effects
# |cluster              |numDE_genes |percentage |
#   |:--------------------|:-----------|:----------|
#   |cMonocyte-COQ7_11213 |621         |3.25       |
#   |cMonocyte-IL_02614   |343         |1.8        |
#   |cMonocyte-ISG15_10   |54          |0.28       |
#   |cMonocyte-S100A8_79  |144         |0.75       |
#   |pDCs_16              |61          |0.32       |
#   |DCs1_15              |8           |0.04       |
#   |DCs2_511             |81          |0.42       |
#   |intMonocyte_3        |423         |2.22       |
#   |ncMonocyte-CD16_48   |419         |2.2        |

##Write results to file

for(cluster in 1:length(keepClusters)){
  # Full results
  filePath <- paste0("output_for_paper/edger/", keepClusters[cluster])
  out <- res[[cluster]][,c("gene", "logFC", "AveExpr", "B", "p_val", "p_adj")]
  write.csv(out, file = paste0(filePath, "_", "D36D8_edgeR_pseudobulk_Allresults.csv"), quote=F, row.names = F)
  
  # Sig genes
  filePath <- paste0("output_for_paper/edger/", keepClusters[cluster])
  out <- res_fil[[cluster]][,c("gene", "logFC", "AveExpr", "B", "p_val", "p_adj")]
  write.csv(out, file = paste0(filePath, "_", "D36D8_edgeR_pseudobulk_sigGenes.csv"), quote=F, row.names = F)
  
}

##Plotting results
##Combining all csv for significant genes into one
#dir.create("/Seurat_analysis_test/IBSM_40/output_for_paper/edger/")
setwd("/working_groups/boylelab/shared/ZulyPava/Training_paper_scRNAseq_R4.4/output_for_paper/edger/")
files <- list.files(pattern = "_sigGenes")
pbulkDe <- rbindlist(sapply(files, fread, simplify = FALSE), idcol = 'filename')

##Making new variables
separate(pbulkDe, filename, sep="_", into = c("cluster_name", "cluster_no", "day", "V1", "V2", "V3"), remove = FALSE) -> pbulkDe1

##Day variable
pbulkDe1$day <- ifelse (pbulkDe1$day == c("D0D8"), c("day8"), pbulkDe1$day)
pbulkDe1$day <- ifelse (pbulkDe1$day == c("D0D36"), c("day36"), pbulkDe1$day)
pbulkDe1$day <- ifelse (pbulkDe1$day == c("D36D8"), c("day368"), pbulkDe1$day)
#removing extra column
pbulkDe2 <- pbulkDe1[, c(2,4, 8:13)]

##saving to csv.
write.csv (pbulkDe2, "psbulk_sigGens__all_DayClusters.csv")


##Combining all csv for all genes into one
#dir.create("D:/'OneDrive - Burnet Institute'/Seurat_analysis_test/IBSM_40/data/edger/sigGenes/")
files <- list.files(pattern = "_Allresults")
#setwd("D:/'OneDrive - Burnet Institute'/Seurat_analysis_test/IBSM_40/data/edger/")
pbulkDeall <- rbindlist(sapply(files, fread, simplify = FALSE), idcol = 'filename')

##Making new variables
separate(pbulkDeall, filename, sep="_", into = c("cluster_name", "cluster_no", "day", "V1", "V2", "V3"), remove = FALSE) -> pbulkDeAll1
#983 7
##Day variable
pbulkDeAll1$day <- ifelse (pbulkDeAll1$day == c("D0D8"), c("day8"), pbulkDeAll1$day)
pbulkDeAll1$day <- ifelse (pbulkDeAll1$day == c("D0D36"), c("day36"), pbulkDeAll1$day)
pbulkDeAll1$day <- ifelse (pbulkDeAll1$day == c("D36D8"), c("day368"), pbulkDeAll1$day)

#removing extra column
pbulkDeAll2 <- pbulkDeAll1[, c(2,4, 8:13)]

##saving to csv.
#setwd("~/Seurat_analysis_test/IBSM_40/")
write.csv (pbulkDeAll2, "psbulk_allgenes_DayClusters.csv")









#Author: Zuleima P
#Pseudobulk DE analysis using edgeR
#https://github.com/hbc/knowledgebase/blob/master/scrnaseq/pseudobulkDE_edgeR.md
#This scripts includes differential analysis using the pseudobluk DE analysis from edger

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("edgeR")


#library(scater)#not installed
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
seurat <- readRDS( "data/IBSM_monoSCT_integbyDayDonor_220823.rds")

#Creating a variable for de
seurat@meta.data$cluster_day <- paste(seurat@meta.data$day, seurat@meta.data$subtype, sep = "_")

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
Idents(object = seurat_s) <- "subtype"
metadata$cluster_id <- factor(seurat_s@active.ident)

# Create single cell experiment object
sce <- SingleCellExperiment(assays = list(counts = counts),
                            colData = metadata)

## Explore the cellular metadata for the dataset
### (number of cells) * (number of meta columns)
dim(colData(sce)) 
##13672    33

#####Explore the counts data##########

#Explore the raw counts for the dataset
## Check the assays present (only counts)
assays(sce)

#Explore the raw counts for the dataset
### (genes) * (cells)
dim(counts(sce))  
## 19088 13672

counts(sce)[1:6, 1:6]

#####Additional QC filtering##########
## Remove lowly expressed genes which have less than 10 cells with any counts
sce <- sce[rowSums(counts(sce) > 1) >= 10, ] 

# (genes) x (cells)
dim(sce)

#10749 13672

#####Preparing data for count aggregation##########
# Named vector of cluster names =
## 2
kids <- purrr::set_names(levels(sce$cluster_id))
kids
# Total number of clusters
## 17
nk <- length(kids)
nk
# Named vector of sample names
sids <- purrr::set_names(levels(as.factor(sce$day_donor)))
sids
# Total number of samples = donor_id 
## 16
ns <- length(sids)

# Generate sample level metadata

## Determine the number of cells per sample
table(sce$day_donor)

# day0_donor0 day0_donor1 day0_donor2 day0_donor3 day8_donor0 day8_donor1 day8_donor2 day8_donor3 
# 2425        1630        1122        1463        2794        2072        1160        1006 


## Turn class "table" into a named vector of cells per sample
n_cells <- table(sce$day_donor) %>%  as.vector()
names(n_cells) <- names(table(sce$day_donor))

## Match the named vector with metadata to combine it
m <- match(names(n_cells), sce$day_donor)

## Create the sample level metadata by selecting specific columns
ei <- data.frame(colData(sce)[m, ], 
                 n_cells, row.names = NULL) %>% 
  dplyr::select("day_donor", "day", "n_cells")
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

class(pb)
dim(pb)
# 136 10749
pb[1:6, 1:6]

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
# Check that every cluster is present in each of the samples
dim(pb$`cMonocyte-IL1B_0`)
#10749     8
dim(pb$`cMonocyte-COQ7_1`)
#10749     8
dim(pb$`cMonocyte-CDKN1A_2`)
#10749     8
dim(pb$`intMonocyte_3`)
#10749     8
dim(pb$`ncMonocyte-CD16_4`)
#10749     8
dim(pb$`DCs2_5`)
#10749     8
dim(pb$`cMonocyte-DUSP1_6`)
#10749     8
dim(pb$`cMonocyte-MALAT_7`)
#10749     8
dim(pb$`ncMonocyte-CD16_8`)
#10749     8
dim(pb$`cMonocyte-S100A8_9`)
#10749     8
dim(pb$`cMonocyte-ISG15_10`)
#10749     8
dim(pb$`DCs2-RPS18_11`)
#10749     8
dim(pb$`cMonocyte-_12`)
#10749     8
dim(pb$`cMonocyte-HLA_13`)
#10749     8
dim(pb$`cMonocyte-IL1RN_14`)
#10749     8
dim(pb$`DCs1_15`)
#10749     8
dim(pb$`pDCs_16`)
#10749     8

##All clusters are present accross all the samples. 
##but we only select the interesting ones
# If one of the cluster is absent from one of the samples then:
# Create a character vector of the clusters to use for DE
keepClusters <- as.character(c("cMonocyte-IL1B_0","cMonocyte-COQ7_1", 
                                "cMonocyte-CDKN1A_2","intMonocyte_3", 
                                "ncMonocyte-CD16_4", "DCs2_5", 
                                "cMonocyte-DUSP1_6", "cMonocyte-MALAT_7", 
                                "ncMonocyte-CD16_8", "cMonocyte-S100A8_9",
                                "cMonocyte-ISG15_10","DCs2-RPS18_11",
                                "cMonocyte-_12","cMonocyte-HLA_13",
                                "cMonocyte-IL1RN_14", "DCs1_15", "pDCs_16"))
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
                      group_id = ei$day_donor,#change to day to generate the other plot
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

ggsave(plot = mds_day, filename = "D:/OneDrive - Burnet Institute/Seurat_analysis_test/IBSM_40/data/edger/graphs/MDSplots_perday_Nov23.pdf",
       width = 620, height= 497, units = "mm", device = "pdf" )

ggsave(plot = mds_day_donor, filename = "D:/OneDrive - Burnet Institute/Seurat_analysis_test/IBSM_40/data/edger/graphs/MDSplots_perday_donor_Nov23.pdf",
       width = 620, height= 497, units = "mm", device = "pdf" )



##### Run differential expression analysis ########
#First set up the experiment design and a contrast matrix.
set.seed(1234)

# construct design & contrast matrix
(design <- model.matrix(~ 0 + ei$day) %>% 
   set_rownames(ei$day_donor) %>% 
   set_colnames(levels(factor(ei$day))))

# A positive FC is increased expression in the ko compared to het
(contrast <- makeContrasts("day8-day0", levels = design))


# for ea. cluster, run edgeR w/ default parameters
res <- lapply(keepClusters, function(k) {
  y <- assays(pbc)[[k]]
  y <- DGEList(y, remove.zeros = TRUE)
  y <- calcNormFactors(y)
  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design)
  fit <- glmQLFTest(fit, contrast = contrast)
  topTags(fit, n = Inf, sort.by = "none")$table %>% 
    dplyr::mutate(gene = rownames(.), cluster_id = k) %>% 
    dplyr::rename(p_val = PValue, p_adj = FDR)
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

##Write results to file

for(cluster in 1:length(keepClusters)){
  # Full results
  filePath <- paste0("D:/OneDrive - Burnet Institute/Seurat_analysis_test/IBSM_40/data/edger/", keepClusters[cluster])
  out <- res[[cluster]][,c("gene", "logFC", "logCPM", "F", "p_val", "p_adj")]
  write.csv(out, file = paste0(filePath, "D0D8_edgeR_pseudobulk_Allresults.csv"), quote=F, row.names = F)
  
  # Sig genes
  filePath <- paste0("D:/OneDrive - Burnet Institute/Seurat_analysis_test/IBSM_40/data/edger/", keepClusters[cluster])
  out <- res_fil[[cluster]][,c("gene", "logFC", "logCPM", "F", "p_val", "p_adj")]
  write.csv(out, file = paste0(filePath, "_", "D0D8_edgeR_pseudobulk_sigGenes.csv"), quote=F, row.names = F)
  
}

rm(list = ls())
####Day 0 vs Day 36 analysis ####

seurat <- readRDS( "data/IBSM_monoSCT_integbyDayDonor_220823.rds")

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
Idents(object = seurat_s) <- "subtype"
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
sce <- sce[rowSums(counts(sce) > 1) >= 10, ]

# (genes) x (cells)
dim(sce)
#10513 8603

#####Preparing data for count aggregation##########
# Named vector of cluster names = 
kids <- purrr::set_names(levels(sce$cluster_id))

# Total number of clusters
## 17
nk <- length(kids)

# Named vector of sample names
sids <- purrr::set_names(levels(as.factor(sce$day_donor)))
sids
# Total number of samples = day_donor_id 
## 8
ns <- length(sids)

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
  dplyr::select("day_donor", "day", "n_cells")
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
dim(pb$`cMonocyte-IL1B_0`)
dim(pb$`cMonocyte-COQ7_1`)
dim(pb$`cMonocyte-CDKN1A_2`)
dim(pb$`intMonocyte_3`)
dim(pb$`ncMonocyte-CD16_4`)
dim(pb$`DCs2_5`)
dim(pb$`cMonocyte-DUSP1_6`)
dim(pb$`cMonocyte-MALAT_7`)
dim(pb$`ncMonocyte-CD16_8`)
dim(pb$`cMonocyte-S100A8_9`)
dim(pb$`cMonocyte-ISG15_10`)
dim(pb$`DCs2-RPS18_11`)
dim(pb$`cMonocyte-_12`)
dim(pb$`cMonocyte-HLA_13`)
dim(pb$`cMonocyte-IL1RN_14`)
dim(pb$`DCs1_15`)
dim(pb$`pDCs_16`)

##Clusters "CD14-FOLR3" was excluded because there are not present for all donors 
# Create a character vector of the clusters to use for DE
keepClusters <- as.character(c("cMonocyte-IL1B_0","cMonocyte-COQ7_1", 
                               "cMonocyte-CDKN1A_2","intMonocyte_3", 
                               "ncMonocyte-CD16_4", "DCs2_5", 
                               "cMonocyte-DUSP1_6", "cMonocyte-MALAT_7", 
                               "ncMonocyte-CD16_8", "cMonocyte-S100A8_9",
                               "cMonocyte-ISG15_10","DCs2-RPS18_11",
                               "cMonocyte-_12","cMonocyte-HLA_13",
                               "cMonocyte-IL1RN_14", "DCs1_15", "pDCs_16"))
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
                       group_id = ei$day, ##change to day_donor to generate the second plot
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

ggsave(plot = mds_d36, filename = "D:/OneDrive - Burnet Institute/Seurat_analysis_test/IBSM_40/data/edger/graphs/MDSplots_d36_perday_250923.pdf",
       width = 620, height= 497, units = "mm", device = "pdf" )


##### Run differential expression analysis ########
#First set up the experiment design and a contrast matrix.
set.seed(12345)
# construct design & contrast matrix

(design <- model.matrix(~ 0 + ei$day) %>% 
   set_rownames(ei$day_donor) %>% 
   set_colnames(levels(factor(ei$day))))

# A positive FC is increased expression in the ko compared to het
(contrast <- makeContrasts("day36-day0", levels = design))


# for ea. cluster, run edgeR w/ default parameters
res <- lapply(keepClusters, function(k) {
  y <- assays(pbc)[[k]]
  y <- DGEList(y, remove.zeros = TRUE)
  y <- calcNormFactors(y)
  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design)
  fit <- glmQLFTest(fit, contrast = contrast)
  topTags(fit, n = Inf, sort.by = "none")$table %>% 
    dplyr::mutate(gene = rownames(.), cluster_id = k) %>% 
    dplyr::rename(p_val = PValue, p_adj = FDR)
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

##Write results to file

for(cluster in 1:length(keepClusters)){
  # Full results
  filePath <- paste0("D:/OneDrive - Burnet Institute/Seurat_analysis_test/IBSM_40/data/edger/", keepClusters[cluster])
  out <- res[[cluster]][,c("gene", "logFC", "logCPM", "F", "p_val", "p_adj")]
  write.csv(out, file = paste0(filePath, "D0D36_edgeR_pseudobulk_Allresults.csv"), quote=F, row.names = F)
  
  # Sig genes
  filePath <- paste0("D:/OneDrive - Burnet Institute/Seurat_analysis_test/IBSM_40/data/edger/", keepClusters[cluster])
  out <- res_fil[[cluster]][,c("gene", "logFC", "logCPM", "F", "p_val", "p_adj")]
  write.csv(out, file = paste0(filePath, "_", "D0D36_edgeR_pseudobulk_sigGenes.csv"), quote=F, row.names = F)
  
}

##Plotting results

##Combining all csv for significant genes into one
dir.create("D:/OneDrive - Burnet Institute/Seurat_analysis_test/IBSM_40/data/edger/sigGenes/")
files <- list.files(path = "D:/OneDrive - Burnet Institute/Seurat_analysis_test/IBSM_40/data/edger/", pattern = "_sigGenes")
setwd("D:/OneDrive - Burnet Institute/Seurat_analysis_test/IBSM_40/data/edger/")
pbulkDe <- rbindlist(sapply(files, fread, simplify = FALSE), idcol = 'filename')
head(pbulkDe)
##Making new variables
separate(pbulkDe, filename, sep="_", into = c("cluster_name", "cluster_no", "day", "V1", "V2", "V3"), remove = FALSE) -> pbulkDe1

head(pbulkDe1)
#983 7
##Day variable
pbulkDe1$day <- ifelse (pbulkDe1$day == c("D0D8"), c("day8"), pbulkDe1$day)
pbulkDe1$day <- ifelse (pbulkDe1$day == c("D0D36"), c("day36"), pbulkDe1$day)

#removing extra columns
colnames(pbulkDe1)
pbulkDe2 <- pbulkDe1[, c(2:4, 8:13)]

##saving to csv.
write.csv (pbulkDe2, "D:/OneDrive - Burnet Institute/Seurat_analysis_test/IBSM_40/data/edger/psbulk_sigGens_all_DayClusters.csv")


#######Significant genes ###########
# We filter the results with the criteria FDR < 0.05, and |logFC| > 1. 
# For each cluster we report the number of differentially expressed genes and what percentage that represents (of total genes).

## Count the number of differential findings by cluster.
# nb. & % of DE genes per cluster
n_de <- vapply(res_fil, nrow, numeric(1))
cbind(cluster=keepClusters, numDE_genes=n_de, 
      percentage = round(n_de / nrow(pbc) * 100, digits =2)) %>%  kable()

##Write results to file

for(cluster in 1:length(keepClusters)){
  # Full results
  filePath <- paste0("output_for_paper/edger/", keepClusters[cluster])
  out <- res[[cluster]][,c("gene", "logFC", "logCPM", "F", "p_val", "p_adj")]
  write.csv(out, file = paste0(filePath, "_", "D8D36_edgeR_pseudobulk_Allresults.csv"), quote=F, row.names = F)
  
  # Sig genes
  filePath <- paste0("output_for_paper/edger/", keepClusters[cluster])
  out <- res_fil[[cluster]][,c("gene", "logFC", "logCPM", "F", "p_val", "p_adj")]
  write.csv(out, file = paste0(filePath, "_", "D8D36_edgeR_pseudobulk_sigGenes.csv"), quote=F, row.names = F)
  
}


##Plotting results
##Combining all csv for significant genes into one
#dir.create("/Seurat_analysis_test/IBSM_40/output_for_paper/edger/")
setwd("output_for_paper/edger/")
files <- list.files(pattern = "_sigGenes")
pbulkDe <- rbindlist(sapply(files, fread, simplify = FALSE), idcol = 'filename')

##Making new variables
separate(pbulkDe, filename, sep="_", into = c("cluster_name", "cluster_no", "day", "V1", "V2", "V3"), remove = FALSE) -> pbulkDe1

##Day variable
pbulkDe1$day <- ifelse (pbulkDe1$day == c("D0D8"), c("day8"), pbulkDe1$day)
pbulkDe1$day <- ifelse (pbulkDe1$day == c("D0D36"), c("day36"), pbulkDe1$day)
pbulkDe1$day <- ifelse (pbulkDe1$day == c("D8D36"), c("day836"), pbulkDe1$day)
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
pbulkDeAll1$day <- ifelse (pbulkDeAll1$day == c("D8D36"), c("day836"), pbulkDeAll1$day)

#removing extra column
pbulkDeAll2 <- pbulkDeAll1[, c(2,4, 8:13)]

##saving to csv.
#setwd("~/Seurat_analysis_test/IBSM_40/")
write.csv (pbulkDeAll2, "psbulk_allgenes_DayClusters.csv")










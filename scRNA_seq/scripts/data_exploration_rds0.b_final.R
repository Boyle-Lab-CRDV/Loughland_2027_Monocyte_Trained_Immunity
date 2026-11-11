###Revised:30 July 2024
### 5 July 2022
#Author: Zuleima Pava
#Notes: This scripts integrates the donor ID information obtained from vireo.
#It follows the SCTransform approach to normalize data. Following script 5B 
#Tests normalization via NormalizeData function use in the classical seurat pipeline.
#Integration analysis via NormalizeData is described in the script Mono_5C_IntegratedNormAnalysis


#####IMPORTANT####
#######  cellranger aggr  has already normalised read depth of individual samples 
# https://github.com/satijalab/seurat/issues/672
# https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger
# use the individual counts generated just using cellranger count without cellranger normalisation and just run the normal pipeline! #####


## This data is using indiviudal counts  
### Use this vignette: https://satijalab.org/seurat/articles/merge_vignette.html

library(here)
library(tidyverse)
library(Seurat)  ##version 4.4.0
library(patchwork)
library(sctransform)#
library(BiocManager)
library(limma)

####COMBINING DATA & CREATING A FEW VARIABLES####
day0mono.data <- Read10X(data.dir = here("input/rawdata/Monocytes_Day0/filtered_feature_bc_matrix"))
day0mono <- CreateSeuratObject(counts = day0mono.data, project = "IBSM.mono", min.cells = 0, min.features = 200)

day8mono.data <- Read10X(data.dir = here("input/rawdata/Monocytes_Day8/filtered_feature_bc_matrix"))
day8mono <- CreateSeuratObject(counts = day8mono.data, project = "IBSM.mono", min.cells = 0, min.features = 200)

day36mono.data <- Read10X(data.dir = here("input/rawdata/Monocytes_Day36/filtered_feature_bc_matrix"))
day36mono <- CreateSeuratObject(counts = day36mono.data, project = "IBSM.mono", min.cells = 0, min.features = 200)

#Due to poor quality day 16 was remove from the analysis.
#For more information check Monocytes_with_D16data/

day0mono@meta.data$day <- c("day0")
View(day0mono@meta.data)

day8mono@meta.data$day <- c("day8")
View(day8mono@meta.data)

day36mono@meta.data$day <- c("day36")
View(day36mono@meta.data)

# ##Adding cell names from row names.
# day0mono@meta.data$cellname <- rownames(day0mono@meta.data)
# day8mono@meta.data$cellname <- rownames(day8mono@meta.data)
# day36mono@meta.data$cellname <- rownames(day36mono@meta.data)
# 
# ### change suffix so that they are the unique per cell type -> 5-D0, 6-D8, 8-D36.Note this step was done late, after receiving donor deconvolution
# day0mono@meta.data$cellname<- gsub("*-1", "-5", day0mono@meta.data$cellname)
# day8mono@meta.data$cellname<- gsub("*-1", "-6", day8mono@meta.data$cellname)
# day36mono@meta.data$cellname<- gsub("*-1", "-8", day36mono@meta.data$cellname)
# 
# ## changing the cell names in the rows
# rownames(day0mono@meta.data) <- day0mono@meta.data$cellname
# rownames(day8mono@meta.data) <- day8mono@meta.data$cellname
# rownames(day36mono@meta.data) <- day36mono@meta.data$cellname
# 
# ### change suffix so that they are the unique per cell type -> 5-D0, 6-D8, 8-D36.Note this step was done late, after receiving donor deconvolution
# 
# day0mono@meta.data$colname <- colnames(day0mono)
# day8mono@meta.data$colname <- colnames(day8mono)
# day36mono@meta.data$colname <- colnames(day36mono)
# 
# day0mono@meta.data$colname<- gsub("*-1", "-5", day0mono@meta.data$colname)
# day8mono@meta.data$colname<- gsub("*-1", "-6", day8mono@meta.data$colname)
# day36mono@meta.data$colname<- gsub("*-1", "-8", day36mono@meta.data$colname)
# 
# vectorname0 <- day0mono@meta.data$colname
# vectorname8 <- day8mono@meta.data$colname
# vectorname36 <- day36mono@meta.data$colname
# 
# names(day0mono@active.ident) <- day0mono@meta.data$colname
# names(day8mono@active.ident) <- day8mono@meta.data$colname
# names(day36mono@active.ident) <- day36mono@meta.data$colname
# 
# ##This part doesn't work in Seurat version
# day0mono <- RenameCells(day0mono, new.names = vectorname0)
# day8mono <- RenameCells(day8mono, new.names = vectorname8)
# day36mono <- RenameCells(day36mono, new.names = vectorname36)

ibsm.mono <- merge(day0mono, y = c(day8mono, day36mono), add.cell.ids = c("day0", "day8", "day36"), project = "IBSM.mono")

##Adding suffix 
#getting cellnames
ibsm.mono@meta.data$cellname <- rownames(ibsm.mono@meta.data)
##changing suffix
ibsm.mono <- ibsm.mono %>%
  AddMetaData(
    metadata = ibsm.mono@meta.data %>%
                mutate(
    cellname_s = case_when(
      day =="day0" ~ gsub("-1$", "-5", cellname),
      day =="day8" ~ gsub("-1$", "-6", cellname),
      day =="day36" ~ gsub("-1$", "-8", cellname),
      TRUE ~ cellname
    )
  ))

## changing the cell names in the rows
rownames(ibsm.mono@meta.data) <- ibsm.mono@meta.data$cellname_s
## changing the cell names in the columns
colnames(ibsm.mono) <- ibsm.mono@meta.data$cellname_s

##merging the layers
ibsm.mono[["RNA"]] <- JoinLayers(ibsm.mono)

##Check if suffix is not an issue
table(check <- day0mono@meta.data$day %in% day8mono@meta.data$day)

saveRDS(ibsm.mono, file = "/Users/zuleip/Documents/Seurat_analysis_test/IBSM_40/data/IBSM_mono_COUNTS_190522.rds")

##adding celltype for future merging
ibsm.mono <- readRDS(file = here("input/IBSM_mono_COUNTS_190522.rds"))

ibsm.mono@meta.data$celltype <- c("Mono")

##check that all the cell names are unique and match the original input
head(colnames(ibsm.mono))
tail(colnames(ibsm.mono))

##how many cells do you have?
table(ibsm.mono@meta.data$day)

#day0 day36 day8 
#7521 2181  8758

##Adding donor information

##import the demultiplexing
donor.ids <- read.csv("input/donor_ids.csv")

##Subsetting only cellnames for monocytes only
donor_ids1 <- donor.ids %>% filter(grepl(pattern = ".*-[568]$", x = cell))

##Changing the colname in donor_ids
donor_ids <- donor.ids %>%
  plyr::rename(c("cell" = "cellname"))

##now join the donor_id information into the Meta.data.
##when you do this it overwrites the rownames in the meta.data.  These are needed to map!!
#Preparing a new seurat object.
ibsm.mono.d <- ibsm.mono

ibsm.mono.d@meta.data <- plyr::join(ibsm.mono.d@meta.data, donor_ids, by = "cellname")
head(ibsm.mono.d@meta.data)

##add the rownames back in from colnames
rownames(ibsm.mono.d@meta.data) <- colnames(ibsm.mono.d)

head(ibsm.mono.d@meta.data)
##rownames are back

#merging
ibsm.mono.D <- merge(ibsm.mono, y = ibsm.mono.d, add.cells.ids = NULL, project = "IBSM40")
# Warning message:
#   In CheckDuplicateCellNames(object.list = objects) :
#   Some cell names are duplicated across objects provided. Renaming to enforce unique cell names.

##removing duplicates and douplets
ibsm.mono.D <- subset(ibsm.mono.D, subset = donor_id != "unassigned")
ibsm.mono.D <- subset(ibsm.mono.D, subset = donor_id != "doublet")
table(ibsm.mono.D@meta.data$donor_id)

##Creating more variables with factors

##make a slot that has both the day and the donor id (from the demultiplexing)
ibsm.mono.D@meta.data$day_donor <- paste(ibsm.mono.D@meta.data$day, ibsm.mono.D@meta.data$donor_id, sep = "_")
donor_levels <- c("day0_donor0", "day0_donor1", "day0_donor2", "day0_donor3","day8_donor0", "day8_donor1","day8_donor2", "day8_donor3","day36_donor0", "day36_donor1", "day36_donor2", "day36_donor3", "NA")
ibsm.mono.D@meta.data$day_donor <- factor(x = ibsm.mono.D@meta.data$day_donor, levels = donor_levels)

# Define an order of identities
day_levels <- c("day0", "day8", "day36")
ibsm.mono.D@meta.data$day <- factor(x = ibsm.mono.D@meta.data$day, levels = day_levels)

# add other details 
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
ibsm.mono.D[["percent.mt"]] <- PercentageFeatureSet(ibsm.mono.D, pattern = "^MT-")

# Add number of genes per UMI for each cell to metadata
ibsm.mono.D$log10GenesPerUMI <- log10(ibsm.mono.D$nFeature_RNA) / log10(ibsm.mono.D$nCount_RNA)

# now save 
saveRDS(ibsm.mono.D, file = "/Users/zuleip/Documents/Seurat_analysis_test/IBSM_40/data/IBSM_mono_prefiltered_210823.rds")

#######################               START QUALITY CHECKS                     #######################

ibsm.mono.D <- readRDS(file = "/Users/zuleip/Documents/Seurat_analysis_test/IBSM_40/data/IBSM_mono_prefiltered_210823.rds")
#filtering # select the parameters especially percent.mt off what looks reasonable from the QC
##Day36 have more genes than day0 or day8. Thus threshold kept at 7000.

#visualise the number of cells per day and relative donor composition
metadata <- ibsm.mono.D@meta.data
metadata %>%
  ggplot(aes(x=day, fill=donor_id)) +
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells per donor_day")
ggsave("graphs/NCells per donor_day.pdf")

##cells look evenly distributed among donors. There are fewer cells on day 36 compared to D0 or D8.

# Visualize the number UMIs/transcripts per cell
#Day8_donor3 has a bimodal distribution of umis vs expression
metadata %>%
  ggplot(aes(color=day_donor[]=="day8_donor3", x=nCount_RNA, fill= day_donor)) +
  geom_density(alpha = 0.2) +
  scale_x_log10() +
  theme_classic() +
  ylab("log10 Cell density") +
  geom_vline(xintercept = 500) +
  ggtitle("UMI counts(transcripts per cell)")

ggsave("graphs/raw_number of UMIs or transcripts per day.pdf")

#Scatterplot number UMIs vs transcripts per cell
FeatureScatter(ibsm.mono.D, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
FeatureScatter(ibsm.mono.D, feature1 = "nCount_RNA", feature2 = "percent.mt")

# Visualize the number mito ratio/transcripts per cell
metadata %>%
  ggplot(aes(color=day, x=percent.mt, fill= day)) +
  geom_density(alpha = 0.2) +
  scale_x_log10() +
  theme_classic() +
  ylab("log10 Cell density") +
  geom_vline(xintercept = 500) +
  ggtitle("percent.mt(transcripts per cell)")
ggsave("graphs/raw_number of percent.mt per day.pdf")

metadata %>%
  ggplot(aes(x=log10GenesPerUMI, color = day, fill=day)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.8)+
  ggtitle("Complexity (genes detected per UMI)")
ggsave("graphs/raw_complexity, genes detected per UMI.pdf")


##Remove low quality cells based on mitochondrial content and gene counts
pdf("graphs/IBSM_mono_VlnPlot_counts_mt_prefilter_210823.pdf", h=4, w=6)
VlnPlot(ibsm.mono.D, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), group.by = "day", pt.size = 0, ncol = 4)
dev.off()

table(ibsm.mono.D@meta.data$day)
# day0 day36  day8 
# 7169  2122  8044

#### LOW QUALITY CELLS FILTERING ####
ibsm.mono.f <- subset(ibsm.mono.D, subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & percent.mt < 20 & log10GenesPerUMI > 0.8)  
table(ibsm.mono.f@meta.data$day)

# day0  day8 day36 
# 6640  7032  1963 

pdf("graphs/IBSM_mono_VlnPlot_counts_mt_postfilter_210823.pdf", h=4, w=6)
VlnPlot(ibsm.mono.f, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), group.by = "day", pt.size = 0, ncol = 4)
dev.off()

# Output a logical vector for every gene on whether the more than zero counts per cell
# > dim(ibsm.mono.f@assays$RNA)
# [1] 36601 15635

# Extract counts
counts0 <- GetAssayData(object = ibsm.mono.f, slot = "counts", assay = "RNA")

# Output a logical vector for every gene on whether the more than zero counts per cell
nonzero <- counts0 > 0

# Sums all TRUE values and returns TRUE if more than 10 TRUE values per gene
keep_genes <- Matrix::rowSums(nonzero) >= 10

# Only keeping those genes expressed in more than 10 cells
counts1 <- counts0[keep_genes, ]

###Add an additional step to remove residual RBC 
#Checking that we don't have hemoglobin genes
#gene names from matrix rows
rcounts1 <- row.names(counts1)
#Check if the gene names are present
hb_present <- c("HBB","HBA1","HBA2","HBG2","HBG1","HBD") %in% rcounts1
hb_present
#They are not present

# Reassign to filtered Seurat object
ibsm.mono.F <- CreateSeuratObject(counts1, meta.data = ibsm.mono.f@meta.data)

options(future.globals.maxSize = 4000 * 1024^4)

##Chceking data distribution after clean up.
#Scatterplot number UMIs vs transcripts per cell
Idents(ibsm.mono.F) <- ibsm.mono.F@meta.data$day_donor

FeatureScatter(ibsm.mono.F, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
FeatureScatter(ibsm.mono.F, feature1 = "nCount_RNA", feature2 = "percent.mt")


#cell cycle sorting
seurat_phase <- NormalizeData(ibsm.mono.F)

# Load cell cycle markers
load("/Users/zuleip/Documents/Seurat_analysis_test/IBSM_40/data/cycle.rda")

# Score cells for cell cycle
seurat_phase <- CellCycleScoring(seurat_phase,
                                 g2m.features = g2m_genes,
                                 s.features = s_genes, set.ident = TRUE)
seurat_phase <- FindVariableFeatures(seurat_phase,
                                        selection.method = "vst",
                                        nfeatures = 2000,
                                        verbose = FALSE)
# Scale the counts
seurat_phase <- ScaleData(seurat_phase)
seurat_phase <- RunPCA(seurat_phase, features = c(s_genes, g2m_genes), npcs = 24)

DimPlot(seurat_phase,
        reduction = "pca",
        group.by= "Phase")
seurat_phase@meta.data$Phase_day <- paste(seurat_phase@meta.data$Phase, seurat_phase@meta.data$day, sep = "_")
DimPlot(seurat_phase, reduction = "pca", group.by = "Phase_day")

##Cell cycle do not seem to be affecting the clustering.
# now save 
saveRDS(ibsm.mono.F, file = "data/IBSM_mono_filtered_210823.rds")


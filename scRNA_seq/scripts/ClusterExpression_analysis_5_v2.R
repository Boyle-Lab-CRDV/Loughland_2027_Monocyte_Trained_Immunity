library(tidyverse)
library(Seurat)
library(here)
library(data.table)
library(UpSetR)

##Read in data
immune.combined <- readRDS( "data/IBSM_monoSCT_integbyDayDonor_220823.rds")

###FIGURE 1a.
##Clustering

#Checking new clustering

##Changing clusters names
immune.combined@meta.data$celltype2 <- NA
immune.combined@meta.data$celltype2[which(str_detect(immune.combined@meta.data$celltype, "cMonocyte-COQ7_11213"))] <- "cMonocyte_COQ7"
immune.combined@meta.data$celltype2[which(str_detect(immune.combined@meta.data$celltype, "cMonocyte-IL_02614"))] <- "cMonocyte_IL1B"
immune.combined@meta.data$celltype2[which(str_detect(immune.combined@meta.data$celltype, "cMonocyte-ISG15_10"))] <- "cMonocyte_ISG15"
immune.combined@meta.data$celltype2[which(str_detect(immune.combined@meta.data$celltype, "cMonocyte-S100A8_79"))] <- "cMonocyte_S100A8"
immune.combined@meta.data$celltype2[which(str_detect(immune.combined@meta.data$celltype, "DCs1_15"))] <- "DCs1"
immune.combined@meta.data$celltype2[which(str_detect(immune.combined@meta.data$celltype, "DCs2_511"))] <- "DCs2"
immune.combined@meta.data$celltype2[which(str_detect(immune.combined@meta.data$celltype, "intMonocyte_3"))] <- "int_Monocyte"
immune.combined@meta.data$celltype2[which(str_detect(immune.combined@meta.data$celltype, "ncMonocyte-CD16_48"))] <- "nc_Monocyte"
immune.combined@meta.data$celltype2[which(str_detect(immune.combined@meta.data$celltype, "pDCs_16"))] <- "pDCs"

##Fixing levels of some variables
celltype_levels <- c("cMonocyte_IL1B","cMonocyte_S100A8",
                     "cMonocyte_COQ7","cMonocyte_ISG15",
                     "int_Monocyte", "nc_Monocyte", 
                     "DCs1", "DCs2", "pDCs")

immune.combined@meta.data$celltype2 <- factor(x = immune.combined@meta.data$celltype2, levels = celltype_levels)

Idents(immune.combined) <- "celltype2"
Idents(immune.combined) <- "integrated_snn_res.0.3"

Tfh_colors_list <- c("#9e0142", "#f4a582", "#5e4fa2", "#66c2a5","#74add1",
                     "#8c510a", "#bf812d", "#fee090", "#01665e", "#b2abd2", "#f8ac0e", "#007b40")



DefaultAssay(immune.combined) <- "SCT"

immune.combined <- SCTransform(immune.combined, vars.to.regress = "percent.mt", verbose = FALSE)
immune.combined <- RunPCA(immune.combined, features = VariableFeatures(object = immune.combined), npcs = 30, verbose = FALSE)
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:20)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:20)
ElbowPlot(immune.combined, ndims = 30)## we'll go with 30.


##How classical markers are driven by donors RETN is one of them!
Idents(immune.combined) <- "day_donor"
pdf("output_for_paper/clusters_expression/Dotplot_monomarkers_per_daydonor_Res0.3_Feb25.pdf", h=15, w=10)
DotPlot(immune.combined, assay = "SCT", 
        features=c("CD14", "FCGR3A", "CD3E","CD4", "CCL3", "RETN", "CD8A", 
                   "CCR7", "TCF7", "CD27", "IL7R", "GZMA", "GZMB","CCL5", 
                   "FOXP3", "CTLA4", "HLA-DRA", "PPBP", "CLEC4C", "CLIC3","IL1B", "CD1C"),
        cols=c("white", "#9e0142")) +
        theme_bw(base_size=14) + 
        theme(axis.text.x = element_text(angle = 45,  vjust = 1, hjust=1)) + 
  theme (panel.grid.major = element_blank (), panel.grid.minor = element_blank ())
dev.off()


##Creating cluster day variable
immune.combined@meta.data$cluster_day <- paste0(immune.combined@meta.data$celltype2, "_", immune.combined@meta.data$day)

Idents(object = immune.combined) <- immune.combined@meta.data$'cluster_day'
pdf("output_for_paper/clusters_expression/Dotplot_monomarkers_per_cluster_day_Feb25.pdf", h=15, w=10)
DotPlot(immune.combined, assay = "SCT", 
        features=c("CD14", "FCGR3A", "CD3E","CD4", "CCL3", "RETN", "CD8A", 
                   "CCR7", "TCF7", "CD27", "IL7R", "GZMA", "GZMB","CCL5", 
                   "FOXP3", "CTLA4", "HLA-DRA", "PPBP", "CLEC4C", "CLIC3","IL1B", "CD1C"),
        cols=c("#74add1","#9e0142",  "#b2abd2"), split.by = "day") +
  theme_bw(base_size=14) + 
  theme(axis.text.x = element_text(angle = 45,  vjust = 1, hjust=1)) + 
  theme (panel.grid.major = element_blank (), panel.grid.minor = element_blank ())+ coord_flip()
dev.off()

Idents(object = immune.combined) <- immune.combined@meta.data$'celltype2'
pdf("output_for_paper/clusters_expression/Dotplot_monomarkers_per_cluster_day_seq_Feb25.pdf", h=15, w=10)
DotPlot(immune.combined, assay = "SCT", 
        features=c("CD14","IL1B","LRRK2","COQ7", "ISG15", "FCGR3A", "CDKN1C", "CX3CR1",
                   "CLEC9A","CADM1", "XCR1",
                   "CD1C","FCER1A","CLEC10A",
                   "CLIC3","NRP1","GZMB"),
        cols=c("#74add1","#9e0142",  "#b2abd2"), split.by = "day") +
  theme_bw(base_size=14) + 
  theme(axis.text.x = element_text(angle = 45,  vjust = 1, hjust=1)) + 
  theme (panel.grid.major = element_blank (), panel.grid.minor = element_blank ())+ coord_flip()
dev.off()

###Differential expression analysis Clusters were identified via cannonical expression of the previous markers
##Diff expression among classical monocytes
Idents(object = immune.combined) <- immune.combined@meta.data$'celltype2'
cMonocyte_IL1B <-  FindMarkers(immune.combined, ident.1 = "cMonocyte_IL1B", ident.2 = c("cMonocyte_S100A8", "cMonocyte_COQ7", "cMonocyte_ISG15", "int_Monocyte", "nc_Monocyte"), test.use= "MAST", random.seed = 1234, min.pct = 0.25, grouping.var = "day")
cMonocyte_IL1B$cluster <- "cMonocyte_IL1B"
cMonocyte_IL1B$Gene <- row.names(cMonocyte_IL1B)
head(cMonocyte_IL1B)

cMonocyte_S100A8 <-  FindMarkers(immune.combined, ident.1 = "cMonocyte_S100A8", ident.2 = c("cMonocyte_IL1B", "cMonocyte_COQ7", "cMonocyte_ISG15", "int_Monocyte", "nc_Monocyte"), test.use= "MAST", random.seed = 1234, min.pct = 0.25, grouping.var = "day")
cMonocyte_S100A8$cluster <- "cMonocyte_S100A8"
cMonocyte_S100A8$Gene <- row.names(cMonocyte_S100A8)
head(cMonocyte_S100A8)

cMonocyte_COQ7 <-  FindMarkers(immune.combined, ident.1 = "cMonocyte_COQ7", ident.2 = c("cMonocyte_IL1B", "cMonocyte_S100A8", "cMonocyte_ISG15", "int_Monocyte", "nc_Monocyte"), test.use= "MAST", random.seed = 1234, min.pct = 0.25, grouping.var = "day")
cMonocyte_COQ7$cluster <- "cMonocyte_COQ7"
cMonocyte_COQ7$Gene <- row.names(cMonocyte_COQ7)
head(cMonocyte_COQ7)

cMonocyte_ISG15 <-  FindMarkers(immune.combined, ident.1 = "cMonocyte_ISG15", ident.2 = c("cMonocyte_IL1B", "cMonocyte_S100A8", "cMonocyte_COQ7", "int_Monocyte", "nc_Monocyte"), test.use= "MAST", random.seed = 1234, min.pct = 0.25, grouping.var = "day")
cMonocyte_ISG15$cluster <- "cMonocyte_ISG15"
cMonocyte_ISG15$Gene <- row.names(cMonocyte_ISG15)
head(cMonocyte_ISG15)

int_Monocyte <-  FindMarkers(immune.combined, ident.1 = "int_Monocyte", ident.2 = c("cMonocyte_IL1B", "cMonocyte_S100A8", "cMonocyte_COQ7", "cMonocyte_ISG15", "nc_Monocyte"), test.use= "MAST", random.seed = 1234, min.pct = 0.25, grouping.var = "day")
int_Monocyte$cluster <- "int_Monocyte"
int_Monocyte$Gene <- row.names(int_Monocyte)
head(int_Monocyte)

nc_Monocyte <-  FindMarkers(immune.combined, ident.1 = "nc_Monocyte", ident.2 = c("cMonocyte_IL1B", "cMonocyte_S100A8", "cMonocyte_COQ7", "cMonocyte_ISG15", "int_Monocyte"), test.use= "MAST", random.seed = 1234, min.pct = 0.25, grouping.var = "day")
nc_Monocyte$cluster <- "nc_Monocyte"
nc_Monocyte$Gene <- row.names(nc_Monocyte)
head(nc_Monocyte)

#combined
cMono_FC_FDR <- rbind(cMonocyte_IL1B,cMonocyte_S100A8,cMonocyte_COQ7,cMonocyte_ISG15,int_Monocyte,nc_Monocyte)
rownames(cMono_FC_FDR)
head(cMono_FC_FDR)

write.csv(cMono_FC_FDR, here("output_for_paper/clusters_expression/FCfdr_cMonosDE_int_0.3_190225.csv"))
cMono_FC_FDR_f <- cMono_FC_FDR %>% 
  filter(p_val_adj<=0.05)
write.csv(cMono_FC_FDR_f, here("output_for_paper/clusters_expression/FCfdr_cMonosDE_Filtered_int_0.3_190225.csv"))


##There is a test for discriminatory markers AUC, which quantifies how accurate a marker is.
cluster0.markers <- FindMarkers(immune.combined, ident.1 = 0, logfc.threshold = 0.25, min.pct = 0.25, test.use = "roc", only.pos = TRUE)
head(cluster0.markers, n = 5)
###Differential expression analysis Clusters were identified via cannonical expression of the previous markers
##Diff expression among classical monocytes
Idents(object = immune.combined) <- immune.combined@meta.data$'celltype2'
cMonocyte_IL1B_AUC <-  FindMarkers(immune.combined, ident.1 = "cMonocyte_IL1B", ident.2 = c("cMonocyte_S100A8", "cMonocyte_COQ7", "cMonocyte_ISG15", "int_Monocyte", "nc_Monocyte"), random.seed=1234, logfc.threshold = 0.25, min.pct = 0.25, test.use = "roc", only.pos = TRUE)
cMonocyte_IL1B_AUC$cluster <- "cMonocyte_IL1B"
cMonocyte_IL1B_AUC$Gene <- row.names(cMonocyte_IL1B_AUC)
head(cMonocyte_IL1B_AUC)

cMonocyte_S100A8_AUC <-  FindMarkers(immune.combined, ident.1 = "cMonocyte_S100A8", ident.2 = c("cMonocyte_IL1B", "cMonocyte_COQ7", "cMonocyte_ISG15", "int_Monocyte", "nc_Monocyte"), random.seed=1234, logfc.threshold = 0.25, min.pct = 0.25, test.use = "roc", only.pos = TRUE)
cMonocyte_S100A8_AUC$cluster <- "cMonocyte_S100A8"
cMonocyte_S100A8_AUC$Gene <- row.names(cMonocyte_S100A8_AUC)
head(cMonocyte_S100A8_AUC)

cMonocyte_COQ7_AUC <-  FindMarkers(immune.combined, ident.1 = "cMonocyte_COQ7", ident.2 = c("cMonocyte_IL1B", "cMonocyte_S100A8", "cMonocyte_ISG15", "int_Monocyte", "nc_Monocyte"), random.seed=1234, logfc.threshold = 0.25, min.pct = 0.25, test.use = "roc", only.pos = TRUE)
cMonocyte_COQ7_AUC$cluster <- "cMonocyte_COQ7"
cMonocyte_COQ7_AUC$Gene <- row.names(cMonocyte_COQ7_AUC)
head(cMonocyte_COQ7_AUC)

cMonocyte_ISG15_AUC <-  FindMarkers(immune.combined, ident.1 = "cMonocyte_ISG15", ident.2 = c("cMonocyte_IL1B", "cMonocyte_S100A8", "cMonocyte_COQ7", "int_Monocyte", "nc_Monocyte"), random.seed=1234, logfc.threshold = 0.25, min.pct = 0.25, test.use = "roc", only.pos = TRUE)
cMonocyte_ISG15_AUC$cluster <- "cMonocyte_ISG15"
cMonocyte_ISG15_AUC$Gene <- row.names(cMonocyte_ISG15_AUC)
head(cMonocyte_ISG15_AUC)

int_Monocyte_AUC <-  FindMarkers(immune.combined, ident.1 = "int_Monocyte", ident.2 = c("cMonocyte_IL1B", "cMonocyte_S100A8", "cMonocyte_COQ7", "cMonocyte_ISG15", "nc_Monocyte"), random.seed=1234, logfc.threshold = 0.25, min.pct = 0.25, test.use = "roc", only.pos = TRUE)
int_Monocyte_AUC$cluster <- "int_Monocyte"
int_Monocyte_AUC$Gene <- row.names(int_Monocyte_AUC)
head(int_Monocyte_AUC)

nc_Monocyte_AUC <-  FindMarkers(immune.combined, ident.1 = "nc_Monocyte", ident.2 = c("cMonocyte_IL1B", "cMonocyte_S100A8", "cMonocyte_COQ7", "cMonocyte_ISG15", "int_Monocyte"), random.seed=1234, logfc.threshold = 0.25, min.pct = 0.25, test.use = "roc", only.pos = TRUE)
nc_Monocyte_AUC$cluster <- "nc_Monocyte"
nc_Monocyte_AUC$Gene <- row.names(nc_Monocyte_AUC)
head(nc_Monocyte_AUC)

#combined
cMono_AUC <- rbind(cMonocyte_IL1B_AUC,cMonocyte_S100A8_AUC,cMonocyte_COQ7_AUC,cMonocyte_ISG15_AUC,int_Monocyte_AUC,nc_Monocyte_AUC)
rownames(cMono_AUC)
head(cMono_AUC)

write.csv(cMono_AUC, here("output_for_paper/clusters_expression/AUC_cMonosDE_int_0.3_190225.csv"))
cMono_AUC_f <- cMono_AUC %>% 
  filter(myAUC>=0.6)
write.csv(cMono_AUC_f, here("output_for_paper/clusters_expression/AUC_cMonosDE_filtered_int_0.3_190225.csv"))

####Dotplot FC fdr clusters markers ####
##Let's do plots of the 10 more significant genes based on fdr

clust_fdr <- cMono_FC_FDR %>% 
  group_by(cluster) %>%
  arrange(p_val_adj) %>%
  slice_head(n=10)
head(clust_fdr)
### markers fdr and fold change
##Creating celltype factors
clust_fact<- c("cMonocyte_IL1B", "cMonocyte_S100A8", "cMonocyte_COQ7", "cMonocyte_ISG15", "int_Monocyte", "nc_Monocyte")
##Creating day factors
clust_fdr$cluster <- factor(clust_fdr$cluster, levels = clust_fact)

##Nicks plot
clust_fdr$log_p <- -log10(clust_fdr$p_val_adj)
clust_fdr$log_p_cat <- NA
clust_fdr$log_p_cat[clust_fdr$log_p <200] <- "<200"
clust_fdr$log_p_cat[clust_fdr$log_p >200 & clust_fdr$log_p<=300] <- "200-300"
clust_fdr$log_p_cat[clust_fdr$log_p >300 ] <- ">300"
clust_fdr$log_p_cat <- factor(clust_fdr$log_p_cat, levels = c("<200","200-300", ">300"))

Tfh_colors_list <- c("#9e0142", "#f4a582", "#5e4fa2", "#66c2a5","#74add1",
                     "#8c510a", "#bf812d", "#fee090", "#01665e", "#b2abd2", "#f8ac0e", "#007b40")

clust_fdr <- clust_fdr %>% 
  filter(avg_log2FC >0) %>%
  group_by(cluster) %>%
  arrange(cluster, Gene)  # Changed to arrange by cluster first, then Gene

# Create dotplot with modified reordering
clust_fdr_dotplot <- ggplot(clust_fdr, aes(reorder(Gene, desc(cluster)), cluster)) + 
  geom_point(aes(colour=avg_log2FC, size=log_p_cat))+
  geom_point(aes(size=log_p_cat), shape=1, colour="black")+
  scale_colour_gradient2(low ="#5e4fa2", mid = "white", high = "#9e0142", midpoint = 0)+
  theme_bw()+
  theme(axis.text.x = element_text(angle=90, vjust = 0.10, hjust=1)) + coord_flip()
clust_fdr_dotplot
ggsave(clust_fdr_dotplot, filename = here("output_for_paper/clusters_expression/FC_fdrmarkers_btw_cMonos_190225.pdf"), h = 10,w =10, dpi = 400,limitsize = FALSE)

####Dotplot AUC clusters markers ####
##Let's do plots of the 10 more significant genes based on fdr

clust_AUC <- cMono_AUC %>%
  filter(myAUC>=0.6)%>%
  group_by(cluster) %>%
  arrange(desc(myAUC)) %>%
  slice_head(n=10)
head(clust_AUC)

### markers fdr and fold change
##Creating celltype factors
clust_fact<- c("cMonocyte_IL1B", "cMonocyte_S100A8", "cMonocyte_COQ7", "cMonocyte_ISG15", "int_Monocyte", "nc_Monocyte")
##Creating day factors
clust_AUC$cluster <- factor(clust_AUC$cluster, levels = clust_fact)

##Nicks plot
clust_AUC$AUC_cat  <- NA
clust_AUC$AUC_cat[clust_AUC$myAUC  <0.7] <- "<0.7"
clust_AUC$AUC_cat[clust_AUC$myAUC >0.7 & clust_AUC$myAUC<=0.8] <- "0.7-0.8"
clust_AUC$AUC_cat[clust_AUC$myAUC >0.8 & clust_AUC$myAUC<=0.9] <- "0.8-0.9"
clust_AUC$AUC_cat[clust_AUC$myAUC >0.9 ] <- ">0.9"
clust_AUC$AUC_cat <- factor(clust_AUC$AUC_cat, levels = c("<0.7","0.7-0.8", "0.8-0.9", ">0.9"))

Tfh_colors_list <- c("#9e0142", "#f4a582", "#5e4fa2", "#66c2a5","#74add1",
                     "#8c510a", "#bf812d", "#fee090", "#01665e", "#b2abd2", "#f8ac0e", "#007b40")

clust_AUC <- clust_AUC %>% 
  group_by(cluster) %>%
  arrange(cluster, Gene)  # Changed to arrange by cluster first, then Gene

# Create dotplot with modified reordering
clust_AUC_dotplot <- ggplot(clust_AUC, aes(reorder(Gene, desc(cluster)), cluster)) + 
  geom_point(aes(colour=avg_log2FC, size=AUC_cat))+
  geom_point(aes(size=AUC_cat), shape=1, colour="black")+
  scale_colour_gradient2(low ="#5e4fa2", mid = "white", high = "#9e0142", midpoint = 0)+
  theme_bw()+
  theme(axis.text.x = element_text(angle=90, vjust = 0.10, hjust=1)) + coord_flip()
clust_AUC_dotplot

ggsave(clust_AUC_dotplot, filename = here("output_for_paper/clusters_expression/AUC_markers_btw_cMonos_190225.pdf"), h = 12,w =6, dpi = 400,limitsize = FALSE)


####Plots for AUC genes filtered to equal or higher than 0.6 analysis using Pantherdb####
##Getting results from https://geneontology.org/ 

path <- (here("output_for_paper/clusters_expression/pantherdb/AUC/filtered_morethan06/"))

file.names <- list.files(pattern = 'AUC.csv', recursive = TRUE)
#Adding names
names <- c("cMonosIL1B", "cMonosISG15", "cMonosS100A8","int_Monocyte", "nc_Monocyte")

file_seq <- 1:length(file.names)


# Read CSV files and immediately name the list
if(any(grepl("*.csv", file.names))) {
  panther_AUC <- setNames(
    lapply(file.names, function(i) {
      fread(input = i)
    }),
    names
  )
}
head(panther_AUC)
#setNames(panther_AUC, names)
# Add a column with the list name to each nested data frame
##normalise colnames
panther_AUC_named <- lapply(names(panther_AUC), function(name) {
  df <- panther_AUC[[name]]
  df$cluster <- name
  colnames(df)[1] <- "GO_biological_proces"
  colnames(df)[2] <- "Homo_sapiens_REFLIST"
  colnames(df)[3] <- "input_genes"
  colnames(df)[4] <- "expected_genes"
  colnames(df)[5] <- "over_under"
  colnames(df)[6] <- "fold_enrichment"
  colnames(df)[7] <- "raw_pval"
  colnames(df)[8] <- "FDR"
  return(df)
})

##Upset plot
##Input  
## Long format: gene, pval, logFC, cluster, day
## From pseudobulk analysis subset of genes with p_adj <0.005 and
## LogFC>abs(1)

panther_AUD_db = as.data.frame(do.call(rbind, panther_AUC_named))
head(panther_AUD_db)

##creating lists of the sets
panther_AUD_db2 <- panther_AUD_db[,c(1,9)]
upset2Split = split(panther_AUD_db2[,1], f = panther_AUD_db2$cluster)
head(upset2Split)

pdf("output_for_paper/clusters_expression/UpsetR_sharedPantherAUC.pdf", h=8, w =14)
upset(fromList(upset2Split),
      nintersects = 15,
      nsets = 15,
      order.by = "freq",
      decreasing = T,
      mb.ratio = c(0.5, 0.5),
      number.angles = 0,
      text.scale = 1.1,
      point.size = 2.8,
      line.size = 1)
dev.off()

##Making log fold change continuous
panther_AUD_db1 = panther_AUD_db%>%
  mutate(
    fold_enrichment_n = 
      case_when(fold_enrichment == "> 100" ~ 100,
                over_under == "-" ~ -1,
                TRUE ~ as.numeric(fold_enrichment)
      ))

##Selecting the top 10 pathways
clust_AUCpath <- panther_AUD_db1 %>%
  group_by(cluster) %>%
  arrange(desc(fold_enrichment_n),FDR) %>%
  slice_head(n=10)

##Changing FDR values for plotting
clust_AUCpath$log_p <- -log10(clust_AUCpath$FDR)
hist(clust_AUCpath$log_p)
clust_AUCpath$log_p_cat <- NA
clust_AUCpath$log_p_cat[clust_AUCpath$log_p <=10] <- "<=10"
clust_AUCpath$log_p_cat[clust_AUCpath$log_p >10 & clust_AUCpath$log_p<=20] <- ">10-20"
clust_AUCpath$log_p_cat[clust_AUCpath$log_p >20 & clust_AUCpath$log_p<=30] <- ">20-30"
clust_AUCpath$log_p_cat[clust_AUCpath$log_p >30 ] <- ">30"
clust_AUCpath$log_p_cat <- factor(clust_AUCpath$log_p_cat, levels = c("<=10",">10-20",">20-30",">30"))

head(clust_AUCpath)
# Create dotplot with modified reordering
clust_AUCpath_dotplot <- ggplot(clust_AUCpath, aes(reorder(GO_biological_proces, desc(cluster)), cluster)) + 
  geom_point(aes(colour=fold_enrichment_n, size=log_p_cat))+
  geom_point(aes(size=log_p_cat), shape=1, colour="black")+
  scale_colour_gradient2(low ="#5e4fa2", mid = "white", high = "#9e0142", midpoint = 0)+
  theme_bw()+
  theme(axis.text.x = element_text(angle=90, vjust = 0.10, hjust=1)) + coord_flip()
clust_AUCpath_dotplot

ggsave(clust_AUCpath_dotplot, filename = here("output_for_paper/clusters_expression/Dotplot_AUCpanthers_GO_btw_cMonos_200225.pdf"), h = 14,w =14, dpi = 400,limitsize = FALSE)

####Plots for AUC all genes analysis using Pantherdb####
##Getting results from https://geneontology.org/ 

path <- (here("output_for_paper/clusters_expression/pantherdb/AUC/unfiltered/"))

file.names <- list.files(pattern = 'AUC_unfilt.csv', recursive = TRUE)
#Adding names
names <- c("cMonosCOQ7", "cMonosIL1B", "cMonosISG15", "cMonosS100A8","int_Monocyte", "nc_Monocyte")

file_seq <- 1:length(file.names)

# if(any(grepl("*.csv", file.names))==TRUE){
#   panther_AUC <- file.names %>% 
#     lapply(., function(i){
#       fread(input = i) 
#     }) 
# }

# Read CSV files and immediately name the list
if(any(grepl("*.csv", file.names))) {
  panther_AUC <- setNames(
    lapply(file.names, function(i) {
      fread(input = i)
    }),
    names
  )
}
head(panther_AUC)
#setNames(panther_AUC, names)
# Add a column with the list name to each nested data frame
##normalise colnames
panther_AUC_named <- lapply(names(panther_AUC), function(name) {
  df <- panther_AUC[[name]]
  df$cluster <- name
  colnames(df)[1] <- "GO_biological_proces"
  colnames(df)[2] <- "Homo_sapiens_REFLIST"
  colnames(df)[3] <- "input_genes"
  colnames(df)[4] <- "expected_genes"
  colnames(df)[5] <- "over_under"
  colnames(df)[6] <- "fold_enrichment"
  colnames(df)[7] <- "raw_pval"
  colnames(df)[8] <- "FDR"
  return(df)
})

##Upset plot
##Input  
## Long format: gene, pval, logFC, cluster, day
## From pseudobulk analysis subset of genes with p_adj <0.005 and
## LogFC>abs(1)

panther_AUD_db = as.data.frame(do.call(rbind, panther_AUC_named))
head(panther_AUD_db)

##creating lists of the sets
panther_AUD_db2 <- panther_AUD_db[,c(1,9)]
upset2Split = split(panther_AUD_db2[,1], f = panther_AUD_db2$cluster)
head(upset2Split)

pdf("output_for_paper/clusters_expression/UpsetR_sharedPantherAUC_unfilt.pdf", h=8, w =14)
upset(fromList(upset2Split),
      nintersects = 15,
      nsets = 15,
      order.by = "freq",
      decreasing = T,
      mb.ratio = c(0.5, 0.5),
      number.angles = 0,
      text.scale = 1.1,
      point.size = 2.8,
      line.size = 1)
dev.off()

##Selecting the top 10 pathways
panther_AUD_db1 = panther_AUD_db%>%
  mutate(
    fold_enrichment_n = 
      case_when(fold_enrichment == "> 100" ~ 100,
                over_under == "-" ~ -1,
                TRUE ~ as.numeric(fold_enrichment)
      ))



clust_AUCpath <- panther_AUD_db1 %>%
  group_by(cluster) %>%
  arrange(desc(fold_enrichment_n), FDR) %>%
  slice_head(n=10)

##Changing FDR values for plotting
clust_AUCpath$log_p <- -log10(clust_AUCpath$FDR)
clust_AUCpath$log_p_cat <- NA
clust_AUCpath$log_p_cat[clust_AUCpath$log_p <=2] <- "<=2"
clust_AUCpath$log_p_cat[clust_AUCpath$log_p >2 & clust_AUCpath$log_p<=6] <- ">2-6"
clust_AUCpath$log_p_cat[clust_AUCpath$log_p >6 & clust_AUCpath$log_p<=10] <- ">6-10"
clust_AUCpath$log_p_cat[clust_AUCpath$log_p >10 ] <- ">10"
clust_AUCpath$log_p_cat <- factor(clust_AUCpath$log_p_cat, levels = c("<=2",">2-6",">6-10",">10"))
#clust_AUCpath$fold_enrichment_n <- NA
# clust_AUCpath = clust_AUCpath%>%
#   mutate(
#     fold_enrichment_n = 
#       case_when(fold_enrichment == "> 100" ~ 100,
#                 TRUE ~ as.numeric(fold_enrichment)
#       ))

head(clust_AUCpath)
# Create dotplot with modified reordering
clust_AUCpath_dotplot <- ggplot(clust_AUCpath, aes(reorder(GO_biological_proces, desc(cluster)), cluster)) + 
  geom_point(aes(colour=fold_enrichment_n, size=log_p_cat))+
  geom_point(aes(size=log_p_cat), shape=1, colour="black")+
  scale_colour_gradient2(low ="#5e4fa2", mid = "white", high = "#9e0142", midpoint = 0)+
  theme_bw()+
  theme(axis.text.x = element_text(angle=90, vjust = 0.10, hjust=1)) + coord_flip()
clust_AUCpath_dotplot

ggsave(clust_AUCpath_dotplot, filename = here("output_for_paper/clusters_expression/Dotplot_AUC_unfilt_panthedb_GO_btw_cMonos_200225.pdf"), h = 14,w =14, dpi = 400,limitsize = FALSE)


####Plots for FDR analysis using Pantherdb####
##Getting results from https://geneontology.org/ 

fpath <- (here("output_for_paper/clusters_expression/pantherdb/FDR/"))

ffile.names <- list.files(pattern = 'FDR.csv', recursive = TRUE)
#Adding names
clust_name <- c("cMonosCOQ7","cMonosIL1B", "cMonosISG15", "cMonosS100A8","int_Monocyte", "nc_Monocyte")


file_seq <- 1:length(ffile.names)

# Read CSV files and immediately name the list
if(any(grepl("*.csv", ffile.names))) {
  panther_FDR <- setNames(
    lapply(ffile.names, function(i) {
      fread(input = i)
    }),
    clust_name
  )
}
head(panther_FDR)
# Add a column with the list name to each nested data frame
##normalise colnames
panther_FDR_named <- lapply(names(panther_FDR), function(name) {
  df <- panther_FDR[[name]]
  df$cluster <- name
  colnames(df)[1] <- "GO_biological_proces"
  colnames(df)[2] <- "Homo_sapiens_REFLIST"
  colnames(df)[3] <- "input_genes"
  colnames(df)[4] <- "expected_genes"
  colnames(df)[5] <- "over_under"
  colnames(df)[6] <- "fold_enrichment"
  colnames(df)[7] <- "raw_pval"
  colnames(df)[8] <- "FDR"
  return(df)
})

##Upset plot
##Input  
## Long format: gene, pval, logFC, cluster, day
## From pseudobulk analysis subset of genes with p_adj <0.005 and
## LogFC>abs(1)

panther_FDR_db = as.data.frame(do.call(rbind, panther_FDR_named))
head(panther_FDR_db)

##creating lists of the sets
panther_FDR_db2 <- panther_FDR_db[,c(1,9)]
upset2Split = split(panther_FDR_db2[,1], f = panther_FDR_db2$cluster)
head(upset2Split)

pdf("output_for_paper/clusters_expression/UpsetR_sharedPanther_FDR.pdf", h=8, w =14)
upset(fromList(upset2Split),
      nintersects = 20,
      nsets = 20,
      order.by = "freq",
      decreasing = T,
      mb.ratio = c(0.5, 0.5),
      number.angles = 0,
      text.scale = 1.1,
      point.size = 2.8,
      line.size = 1)
dev.off()

#Converting fold_enrichment to numeric
panther_FDR_db2 = panther_FDR_db%>%
  mutate(
    fold_enrichment_n = 
      case_when(fold_enrichment == "> 100" ~ 100,
                over_under == "-" ~ -1,
                TRUE ~ as.numeric(fold_enrichment)
      ))

##Selecting the top 10 pathways based on FDR order
# clust_FDRpath <- panther_FDR_db2 %>%
#   group_by(cluster) %>%
#   arrange(FDR) %>%
#   slice_head(n=10)

##Changing FDR values for plotting based on FDR
# clust_FDRpath$log_p <- -log10(clust_FDRpath$FDR)
# hist(clust_FDRpath$log_p)
# clust_FDRpath$log_p_cat <- NA
# clust_FDRpath$log_p_cat[clust_FDRpath$log_p <=40] <- "<=40"
# clust_FDRpath$log_p_cat[clust_FDRpath$log_p >40 & clust_FDRpath$log_p<=80] <- ">40-80"
# clust_FDRpath$log_p_cat[clust_FDRpath$log_p >80 & clust_FDRpath$log_p<=120] <- ">80-120"
# clust_FDRpath$log_p_cat[clust_FDRpath$log_p >120 ] <- ">120"
# clust_FDRpath$log_p_cat <- factor(clust_FDRpath$log_p_cat, levels = c("<=40",">40-80",">80-120",">120"))


##Selecting the top 10 pathways based on Fold enrichment
clust_FDRpath <- panther_FDR_db2 %>%
  group_by(cluster) %>%
  arrange(desc(fold_enrichment_n), FDR) %>%
  slice_head(n=10)

##Changing FDR values for plotting based on Fold enrichment
clust_FDRpath$log_p <- -log10(clust_FDRpath$FDR)
hist(clust_FDRpath$log_p)
clust_FDRpath$log_p_cat <- NA
clust_FDRpath$log_p_cat[clust_FDRpath$log_p <=2] <- "<=2"
clust_FDRpath$log_p_cat[clust_FDRpath$log_p >2 & clust_FDRpath$log_p<=4] <- ">2-4"
clust_FDRpath$log_p_cat[clust_FDRpath$log_p >4 & clust_FDRpath$log_p<=6] <- ">4-6"
clust_FDRpath$log_p_cat[clust_FDRpath$log_p >6 ] <- ">6"
clust_FDRpath$log_p_cat <- factor(clust_FDRpath$log_p_cat, levels = c("<=2",">2-4",">4-6",">6"))


head(clust_FDRpath)
# Create dotplot with modified reordering
clust_FDRpath_dotplot <- ggplot(clust_FDRpath, aes(reorder(GO_biological_proces, dplyr::desc(cluster)), cluster)) + 
  geom_point(aes(colour=fold_enrichment_n, size=log_p_cat))+
  geom_point(aes(size=log_p_cat), shape=1, colour="black")+
  scale_colour_gradient2(low ="#5e4fa2", mid = "white", high = "#9e0142", midpoint = 0)+
  theme_bw()+
  theme(axis.text.x = element_text(angle=90, vjust = 0.10, hjust=1)) + coord_flip()
clust_FDRpath_dotplot

ggsave(clust_FDRpath_dotplot, filename = here("output_for_paper/clusters_expression/Dotplot_FDRpanthers_GO_btw_cMonos_200225.pdf"), h = 14,w =14, dpi = 400,limitsize = FALSE)


####Plots for FDR analysis using reactome db####
##Getting results from https://geneontology.org/ 

rpath <- (here("output_for_paper/clusters_expression/reactome/all_sig_N_upreg_sorted/"))
##only upregulated genes and FDR was sorted in the page.
rfile.names <- list.files(pattern = 'FDRsorted_UP.txt', recursive = TRUE)
#Adding names
clust_name <- c("cMonosCOQ7","cMonosIL1B", "cMonosISG15", "cMonosS100A8","int_Monocyte", "nc_Monocyte")


file_seq <- 1:length(rfile.names)

# Read CSV files and immediately name the list
if(any(grepl("*.txt", rfile.names))) {
  reactome_FDR <- setNames(
    lapply(ffile.names, function(i) {
      fread(input = i)
    }),
    clust_name
  )
}
head(reactome_FDR)
# Add a column with the list name to each nested data frame
##normalise colnames
reactome_FDR_named <- lapply(names(reactome_FDR), function(name) {
  df <- reactome_FDR[[name]]
  df$cluster <- name
  colnames(df)[1] <- "GO_biological_proces"
  colnames(df)[2] <- "Homo_sapiens_REFLIST"
  colnames(df)[3] <- "input_genes"
  colnames(df)[4] <- "expected_genes"
  colnames(df)[5] <- "over_under"
  colnames(df)[6] <- "fold_enrichment"
  colnames(df)[7] <- "raw_pval"
  colnames(df)[8] <- "FDR"
  return(df)
})

##Upset plot
##Input  
## Long format: gene, pval, logFC, cluster, day
## From pseudobulk analysis subset of genes with p_adj <0.005 and
## LogFC>abs(1)

reactome_FDR_db = as.data.frame(do.call(rbind, reactome_FDR_named))
head(reactome_FDR_db)

##creating lists of the sets
reactome_FDR_db2 <- reactome_FDR_db[,c(1,9)]
upset2Split = split(reactome_FDR_db2[,1], f = reactome_FDR_db2$cluster)
head(upset2Split)

pdf("output_for_paper/clusters_expression/UpsetR_sharedreactome_FDR.pdf", h=8, w =14)
upset(fromList(upset2Split),
      nintersects = 20,
      nsets = 20,
      order.by = "freq",
      decreasing = T,
      mb.ratio = c(0.5, 0.5),
      number.angles = 0,
      text.scale = 1.1,
      point.size = 2.8,
      line.size = 1)
dev.off()

##Selecting the top 10 pathways
clust_reactome_FDRpath <- reactome_FDR_db %>%
  group_by(cluster) %>%
  arrange(FDR) %>%
  slice_head(n=10)

##Changing FDR values for plotting
clust_reactome_FDRpath$log_p <- -log10(clust_FDRpath$FDR)
clust_reactome_FDRpath$log_p_cat <- NA
clust_reactome_FDRpath$log_p_cat[clust_reactome_FDRpath$log_p <=40] <- "<=40"
clust_reactome_FDRpath$log_p_cat[clust_reactome_FDRpath$log_p >40 & clust_reactome_FDRpath$log_p<=80] <- ">40-80"
clust_reactome_FDRpath$log_p_cat[clust_reactome_FDRpath$log_p >80 & clust_reactome_FDRpath$log_p<=120] <- ">80-120"
clust_reactome_FDRpath$log_p_cat[clust_reactome_FDRpath$log_p >120 ] <- ">120"
clust_reactome_FDRpath$log_p_cat <- factor(clust_reactome_FDRpath$log_p_cat, levels = c("<=40",">40-80",">80-120",">120"))
clust_reactome_FDRpath$fold_enrichment_n <- NA
clust_reactome_FDRpath = clust_reactome_FDRpath%>%
  mutate(
    fold_enrichment_n = 
      case_when(fold_enrichment == "> 100" ~ 100,
                TRUE ~ as.numeric(fold_enrichment)
      ))

head(clust_reactome_FDRpath)
# Create dotplot with modified reordering
clust_reactome_FDRpath_dotplot <- ggplot(clust_reactome_FDRpath, aes(reorder(GO_biological_proces, desc(cluster)), cluster)) + 
  geom_point(aes(colour=fold_enrichment_n, size=log_p_cat))+
  geom_point(aes(size=log_p_cat), shape=1, colour="black")+
  scale_colour_gradient2(low ="#5e4fa2", mid = "white", high = "#9e0142", midpoint = 0)+
  theme_bw()+
  theme(axis.text.x = element_text(angle=90, vjust = 0.10, hjust=1)) + coord_flip()
clust_reactome_FDRpath_dotplot

ggsave(clust_reactome_FDRpath_dotplot, filename = here("output_for_paper/clusters_expression/Dotplot_reactome_FDRUPsorted_GO_btw_cMonos_200225.pdf"), h = 14,w =14, dpi = 400,limitsize = FALSE)


####Without sorting ####
rpath <- (here("output_for_paper/clusters_expression/reactome/all_sig_N_upreg/"))
##only upregulated genes and FDR was sorted in the page.
rfile.names <- list.files(pattern = 'FDRsorted_UP.txt', recursive = TRUE)
#Adding names
clust_name <- c("cMonosCOQ7","cMonosIL1B", "cMonosISG15", "cMonosS100A8","int_Monocyte", "nc_Monocyte")


file_seq <- 1:length(rfile.names)

# Read CSV files and immediately name the list
if(any(grepl("*.txt", rfile.names))) {
  reactome_FDR <- setNames(
    lapply(ffile.names, function(i) {
      fread(input = i)
    }),
    clust_name
  )
}
head(reactome_FDR)
# Add a column with the list name to each nested data frame
##normalise colnames
reactome_FDR_named <- lapply(names(reactome_FDR), function(name) {
  df <- reactome_FDR[[name]]
  df$cluster <- name
  colnames(df)[1] <- "GO_biological_proces"
  colnames(df)[2] <- "Homo_sapiens_REFLIST"
  colnames(df)[3] <- "input_genes"
  colnames(df)[4] <- "expected_genes"
  colnames(df)[5] <- "over_under"
  colnames(df)[6] <- "fold_enrichment"
  colnames(df)[7] <- "raw_pval"
  colnames(df)[8] <- "FDR"
  return(df)
})

##Upset plot
##Input  
## Long format: gene, pval, logFC, cluster, day
## From pseudobulk analysis subset of genes with p_adj <0.005 and
## LogFC>abs(1)

reactome_FDR_db = as.data.frame(do.call(rbind, reactome_FDR_named))
head(reactome_FDR_db)

##creating lists of the sets
reactome_FDR_db2 <- reactome_FDR_db[,c(1,9)]
upset2Split = split(reactome_FDR_db2[,1], f = reactome_FDR_db2$cluster)
head(upset2Split)

pdf("output_for_paper/clusters_expression/UpsetR_sharedreactome_FDRnotsorted.pdf", h=8, w =14)
upset(fromList(upset2Split),
      nintersects = 20,
      nsets = 20,
      order.by = "freq",
      decreasing = T,
      mb.ratio = c(0.5, 0.5),
      number.angles = 0,
      text.scale = 1.1,
      point.size = 2.8,
      line.size = 1)
dev.off()

##Selecting the top 10 pathways
clust_reactome_FDRpath <- reactome_FDR_db %>%
  group_by(cluster) %>%
  arrange(FDR) %>%
  slice_head(n=10)

##Changing FDR values for plotting
clust_reactome_FDRpath$log_p <- -log10(clust_FDRpath$FDR)
clust_reactome_FDRpath$log_p_cat <- NA
clust_reactome_FDRpath$log_p_cat[clust_reactome_FDRpath$log_p <=40] <- "<=40"
clust_reactome_FDRpath$log_p_cat[clust_reactome_FDRpath$log_p >40 & clust_reactome_FDRpath$log_p<=80] <- ">40-80"
clust_reactome_FDRpath$log_p_cat[clust_reactome_FDRpath$log_p >80 & clust_reactome_FDRpath$log_p<=120] <- ">80-120"
clust_reactome_FDRpath$log_p_cat[clust_reactome_FDRpath$log_p >120 ] <- ">120"
clust_reactome_FDRpath$log_p_cat <- factor(clust_reactome_FDRpath$log_p_cat, levels = c("<=40",">40-80",">80-120",">120"))
clust_reactome_FDRpath$fold_enrichment_n <- NA
clust_reactome_FDRpath = clust_reactome_FDRpath%>%
  mutate(
    fold_enrichment_n = 
      case_when(fold_enrichment == "> 100" ~ 100,
                TRUE ~ as.numeric(fold_enrichment)
      ))

head(clust_FDRpath)
# Create dotplot with modified reordering
clust_reactome_FDRpath_dotplot <- ggplot(clust_reactome_FDRpath, aes(reorder(GO_biological_proces, desc(cluster)), cluster)) + 
  geom_point(aes(colour=fold_enrichment_n, size=log_p_cat))+
  geom_point(aes(size=log_p_cat), shape=1, colour="black")+
  scale_colour_gradient2(low ="#5e4fa2", mid = "white", high = "#9e0142", midpoint = 0)+
  theme_bw()+
  theme(axis.text.x = element_text(angle=90, vjust = 0.10, hjust=1)) + coord_flip()
clust_reactome_FDRpath_dotplot

ggsave(clust_reactome_FDRpath_dotplot, filename = here("output_for_paper/clusters_expression/Dotplot_reactome_FDRUPNOTsorted_GO_btw_cMonos_200225.pdf"), h = 14,w =14, dpi = 400,limitsize = FALSE)


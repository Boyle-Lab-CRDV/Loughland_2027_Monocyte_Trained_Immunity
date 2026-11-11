# SCRIPT TO RUN ENRICHMENT ANALYSIS
# Author:         Damian Oyong
# Date created:   7 November 2024
# Date modified:  12 March 2025
# Requirements:   - Set your annotation list as a list object containing a vector list of genes/proteins

# Descriptions: There are two methods to choose here, do_gsea and do_enrich.
# do_gsea uses a ranked based enrichment analysis, whereby you need to input entire ranked list of genes
# do_enrich uses over-representation analysis, where you just need to input a set of genes. Normally these genes have been pre-filtered by p-values
# Original workflow: https://biostatsquid.com/fgsea-tutorial-gsea/ 

# Libraries -------------------------------------------------------------------------------------------------------

# if(!require(BloodGen3Module)){
#   BiocManager::install("BloodGen3Module")
#   library(BloodGen3Module)
# }
# if(!require(biomaRt)){
#   BiocManager::install("biomaRt")
#   library(biomaRt)
# }
# if(!require(fgsea)){
#   BiocManager::install("fgsea")
#   library(fgsea)
# }
# if(!require(dplyr)){
#   package.install("dplyr")
#   library(dplyr)
# }
# if(!require(tidyr)){
#   package.install("tidyr")
#   library(tidyr)
# }
# if(!require(msigdbr)){
#   package.install("msigdbr")
#   library(msigdbr)
# }
# if(!require(data.table)){
#   package.install("data.table")
#   library(data.table)
# }
# if(!require(enrichR)){
#   install.packages("enrichR")
#   library(enrichR)
# }

library(ggplot2)
library(purrr)
library(BloodGen3Module)
library (biomaRt)
library(fgsea)
library(dplyr)
library(tidyr)
library(msigdbr)
library(data.table)
library(enrichR)
#library(ggrepel)
library(ComplexHeatmap)
library(circlize)
library(grid)
#Set seed
# Set seed
set.seed(12345)

# Gene sets -------------------------------------------------------------------------------------------------------
## 1. Preparing gene sets
# from msigdb https://www.gsea-msigdb.org/gsea/msigdb/human/collections.jsp
# Hallmark collection (MSigDB category "H"): ~50 curated, non-redundant gene sets
GO_HL <- msigdbr(species = "Homo sapiens", category = "H")
GO_HL <- split(GO_HL$gene_symbol, GO_HL$gs_name)  # list: gs_name -> vector of gene symbols

# Gene Ontology: Biological Process (GO:BP) (MSigDB category C5, subcategory GO:BP)
GO_BP <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP")
GO_BP <- split(GO_BP$gene_symbol, GO_BP$gs_name)

# Gene Ontology: Cellular Component (GO:CC) (MSigDB category C5, subcategory GO:CC)
GO_CC <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:CC")
GO_CC <- split(GO_CC$gene_symbol, GO_CC$gs_name)

# Gene Ontology: Molecular Function (GO:MF) (MSigDB category C5, subcategory GO:MF)
GO_MF <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:MF")
GO_MF <- split(GO_MF$gene_symbol, GO_MF$gs_name)

# Human Phenotype Ontology (HPO) gene sets (MSigDB category C5, subcategory HPO)
GO_HPO <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "HPO")
GO_HPO <- split(GO_HPO$gene_symbol, GO_HPO$gs_name)

# Canonical Pathways: Reactome (MSigDB category C2, subcategory CP:REACTOME)
CP_RC <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME")
CP_RC <- split(CP_RC$gene_symbol, CP_RC$gs_name)

# Immunologic Signatures (ImmuneSigDB) (MSigDB category C7)
immdb <- msigdbr(species = "Homo sapiens", category = "C7")
immdb <- split(immdb$gene_symbol, immdb$gs_name)

# MSigDB category C8: Cell type signatures (often cell-state/cell-type marker sets)
C8 <- msigdbr(species = "Homo sapiens", category = "C8")
C8 <- split(C8$gene_symbol, C8$gs_name)

##Selecting hallmark database
GO_HL = msigdbr(species = "Homo sapiens", category = "H")
#Formatting for fsea pathway/genesymbol
GO_HL = split(GO_HL$gene_symbol, GO_HL$gs_name)
#Selecting Gene Ontology Biological Process (GO:BP)
GO_BP = msigdbr(species = "Homo sapiens", subcategory = "GO:BP")
GO_BP = split(GO_BP$gene_symbol, GO_BP$gs_name)

#BTM3 https://academic.oup.com/bioinformatics/article/37/16/2382/6149124?login=true
btm3 <- BloodGen3Module:::Module_listGen3 %>% 
  mutate(module_names = paste0(Module, "_", Function))  # a new naming column
# Create the list by splitting the 'Gene' column by the new 'Module_Function' column
btm3 <- split(btm3$Gene, btm3$module_names)
# Let's delete all "TBD" modules (To Be Determined)
btm3 <- btm3[!grepl("TBD", names(btm3))]

# Functions -------------------------------------------------------------------------------------------------------

do_gsea <- function(gsets, ranked_df){
  res <- apply(ranked_df, 2, function(ranks) {
    # Remove NA values
    ranks <- ranks[!is.na(ranks)]
    ranks <- sort(ranks, decreasing = TRUE)
    # Run fgsea on the ranked_df
    fgsea(pathways = gsets, stats = ranks)
  })
  return(res)
}

do_enrich <- function(df, dbs){
  dbs = dbs
  for(i in 2:ncol(df)){
    col_name <- colnames(df)[i]
    print(paste0("Analysing ",col_name))
    
    df_pos <- df %>% 
      filter(., !is.na(.data[[col_name]])) %>% 
      filter(., .data[[col_name]] > 0)
    
    if (nrow(df_pos) == 0) {
      print("No genes after filtering, proceeding with next lists")
      enrich_pos <- NULL  # Set as NULL to indicate no data
    } else {
      enrich_pos <- enrichr(df_pos[[1]], databases = dbs)
      # Concatenating the dataframes with an additional column of dataframe names
      enrich_pos <- bind_rows(
        map(enrich_pos, ~ mutate(.x, database = cur_data_all())), 
        .id = "database") %>% 
        mutate(direction = "Upregulated")
    }
    
    
    df_neg <- df %>% 
      filter(., !is.na(.data[[col_name]])) %>% 
      filter(., .data[[col_name]] < 0)
    
    if (nrow(df_neg) == 0) {
      print("No genes after filtering, proceeding with next lists")
      enrich_neg <- NULL  # Set as NULL to indicate no data
    } else {
      enrich_neg <- enrichr(df_neg[[1]], databases = dbs)
      # Concatenating the dataframes with an additional column of dataframe names
      enrich_neg <- bind_rows(
        map(enrich_neg, ~ mutate(.x, database = cur_data_all())), 
        .id = "database") %>% 
        mutate(direction = "Downregulated")
    }
    
    enrich_comb <- bind_rows(enrich_pos, enrich_neg) %>% 
      mutate(Variable = col_name)
    
    if(i == 2){
      res_comb = enrich_comb
    }else{
      res_comb = rbind(res_comb,enrich_comb)
    }
    
  }
  return(res_comb)
}

# Main Script -----------------------------------------------------------------------------------------------------
## 2. Reading input and formatting
db = read.csv("input/psbulk_allgenes_DayClusters.csv", row.names = 1)
db
##Input require
db.f = db%>%
  filter(db$day %in% c("day8", "day36"))%>%
  mutate(
  cluster_name = factor(cluster_name,
                                  labels = c("cMonocyte-IL"="CM1",
                                             "cMonocyte-S100A8"="CM2",
                                             "cMonocyte-COQ7"="CM3",
                                             "cMonocyte-ISG15"="CM4",
                                             "DCs1" ="DCs1",
                                             "DCs2" = "DCs2",
                                             "pDCs" = "pDCs",
                                             "intMonocyte" = "intMonocyte",
                                             "ncMonocyte-CD16" ="ncMonocyte")),
  cluster_day = paste0(cluster_name,"_",day),
  rankings = round(sign(logFC)*(-log10(p_adj)),4))# we will use the signed p values from spatial DGE as ranking


# Some genes have such low p values that the signed pval is +- inf,
# we need to change it to the maximum * constant to avoid problems with fgsea
# max_ranking <- max(rankings[is.finite(db.f$rankings)])
# min_ranking <- min(rankings[is.finite(db.f$rankings)])
# rankings <- replace(rankings, rankings > max_ranking, max_ranking * 10)
# rankings <- replace(rankings, rankings < min_ranking, min_ranking * 10)



## 3. Prepare ranked list of genes -----------------------------------------------
## FGSEA needs a named vector with the rankings and the genes as names
## We will use  signed p value as ranking: sign(df$log2fc)*(-log10(df$pval)))
# Note that you can also just use the log2FC only - but it will not take into account 
# genes with a large FC but non-significant

# names(rankings) <- db.f$gene # genes as names#
# rankings <- sort(rankings, decreasing = TRUE) # sort genes by ranking
# plot(rankings)

db.f.0 = db.f%>%
  select(cluster_day, gene, rankings)%>%
  pivot_wider(., names_from = cluster_day, values_from = rankings, values_fill= 0)%>% 
  tibble::column_to_rownames(., var = "gene") 

gsets = c(GO_BP)
test <- do_gsea(gsets = gsets, ranked_df = db.f.0)
#how many sig pathways we got

lapply(test, function(x) sum(x$padj < 0.01, na.rm = TRUE))
# $CM1_day36
# [1] 233
# 
# $CM1_day8
# [1] 211
# 
# $CM2_day36
# [1] 407
# 
# $CM2_day8
# [1] 162
# 
# $CM3_day36
# [1] 169
# 
# $CM3_day8
# [1] 6
# 
# $CM4_day36
# [1] 174
# 
# $CM4_day8
# [1] 96
# 
# $DCs1_day36
# [1] 0
# 
# $DCs1_day8
# [1] 11
# 
# $DCs2_day36
# [1] 8
# 
# $DCs2_day8
# [1] 132
# 
# $pDCs_day36
# [1] 197
# 
# $pDCs_day8
# [1] 88
# 
# $intMonocyte_day36
# [1] 84
# 
# $intMonocyte_day8
# [1] 137
# 
# $ncMonocyte_day36
# [1] 51
# 
# $ncMonocyte_day8
# [1] 3


test.db <- map2(test, names(test), ~ mutate(.x, cluster = .y)) %>% 
  bind_rows() %>%
  filter(padj < 0.01) %>%
  separate(.,cluster, into=c("cluster_name", "day"), remove = FALSE)

##Converting the leadingEdge variable from list to character
# Convert list elements to a single comma-separated string
test.db$leadingEdge <- vapply(test.db$leadingEdge, paste, collapse = ", ", character(1L))

write.csv(test.db, "output/data/gsea_ibsm40_monos_170326.csv")
### Plotting heatmap ####

#Id the top ten pathways per cluster
CM1_d36_p <- test.db %>%
  filter(cluster %in% "CM1_day36")%>%
  arrange(padj) %>%
  slice_head(.,n=15) %>%
  pull(pathway)
CM2_d36_p <- test.db %>%
  filter(cluster %in% "CM2_day36")%>%
  arrange(padj) %>%
  slice_head(.,n=15) %>%
  pull(pathway)
CM3_d36_p <- test.db %>%
  filter(cluster %in% "CM3_day36")%>%
  arrange(padj) %>%
  slice_head(.,n=15) %>%
  pull(pathway)
CM4_d36_p <- test.db %>%
  filter(cluster %in% "CM4_day36")%>%
  arrange(padj) %>%
  slice_head(.,n=15) %>%
  pull(pathway)

CM_d36_p <- Reduce(union, list(CM1_d36_p, CM2_d36_p, CM3_d36_p, CM4_d36_p))

CM1_d8_p <- test.db %>%
  filter(cluster %in% "CM1_day8")%>%
  arrange(padj) %>%
  slice_head(.,n=15) %>%
  pull(pathway)
CM2_d8_p <- test.db %>%
  filter(cluster %in% "CM2_day8")%>%
  arrange(padj) %>%
  slice_head(.,n=15) %>%
  pull(pathway)
CM3_d8_p <- test.db %>%
  filter(cluster %in% "CM3_day8")%>%
  arrange(padj) %>%
  slice_head(.,n=15) %>%
  pull(pathway)
CM4_d8_p <- test.db %>%
  filter(cluster %in% "CM4_day8")%>%
  arrange(padj) %>%
  slice_head(.,n=15) %>%
  pull(pathway)
  
CM_d8_p <- Reduce(union, list(CM1_d8_p, CM2_d8_p, CM3_d8_p, CM4_d8_p))

# The captioned code produces heatmaps ordered alphabetically.                     
# test.plot.d8 <- test.db[test.db$pathway %in% CM_d8_p,] %>%
#     select(pathway, cluster, NES) %>%
#     filter(grepl("CM", cluster))%>%
#     filter(grepl("day8", cluster), !is.na(NES)) %>%
#     pivot_wider(names_from = pathway, values_from = NES) %>%
#     mutate(cluster = factor(cluster, levels = c("CM1_day8", "CM2_day8", "CM3_day8", "CM4_day8"),
#                             labels=c("CM1_day8"="CM1", "CM2_day8"="CM2", "CM3_day8"="CM3", "CM4_day8"="CM4"))) %>%
#     arrange(cluster) %>%
#     tibble::column_to_rownames(., var = "cluster") 
# 
# test.plot.d36 <- test.db[test.db$pathway %in% CM_d36_p,]%>%
#   select(pathway, cluster, NES) %>%
#   filter(grepl("CM", cluster))%>%
#   filter(grepl("day36", cluster), !is.na(NES)) %>%
#   pivot_wider(names_from = pathway, values_from = NES) %>%
#   mutate(cluster = factor(cluster, levels = c("CM1_day36", "CM2_day36", "CM3_day36", "CM4_day36"),
#                           labels=c("CM1_day36"="CM1", "CM2_day36"="CM2", "CM3_day36"="CM3", "CM4_day36"="CM4"))) %>%
#   arrange(cluster) %>%
#   tibble::column_to_rownames(., var = "cluster") 
# 
# 
# #mat_d8 <- t(test.plot.d8[1:25])
# mat_d8 <- t(test.plot.d8)
# col_fun <- colorRamp2(
#   c(-2, 0, 2),
#   c("#2166AC", "#FFFFFF", "#B2182B")   # blue - white - red (swap colors if you want)
# )
# 
# pdf("output/heatmap_nex_ibsm40_monos_clasicMd8_top25_padj.pdf", height = 15, width = 15)
# d8 <- ComplexHeatmap::Heatmap(mat_d8,
#                                na_col = "grey",
#                                cluster_rows = F,
#                                cluster_columns = F,
#                                column_names_rot = 0,
#                                column_names_centered = TRUE,
#                                column_names_side = "top",
#                                column_names_gp = gpar(fontsize = 8, face = 'bold'),
#                                row_names_gp = gpar(fontsize = 8, face = 'bold'),
#                                heatmap_legend_param = list(
#                                  title = "NES",
#                                  legend_direction = "horizontal",
#                                  legend_width = unit(6, "cm")),
#                                col = col_fun,
#                                # Control heatmap size
#                                width = unit(ncol(mat_d8) * 10, "mm"), 
#                                height = unit(nrow(mat_d8) * 10, "mm"))
# draw(d8, heatmap_legend_side = "bottom")
# 
# dev.off()
# 
# #mat_d36 <- t(test.plot.d36[1:25])
# mat_d36 <- t(test.plot.d36)
# ##Adjusting colos
# col_fun <- colorRamp2(
#   c(-2, 0, 2),
#   c("#2166AC", "#FFFFFF", "#B2182B")   # blue - white - red (swap colors if you want)
# )
# 
# pdf("output/heatmap_nex_ibsm40_monos_clasicMd36_top25_padj.pdf", height = 15, width = 15)
# d36 <- ComplexHeatmap::Heatmap(mat_d36,
#                         na_col = "grey",
#                         cluster_rows = F,
#                         cluster_columns = F,
#                         column_names_rot = 0,
#                         column_names_centered = TRUE,
#                         column_names_side = "top",
#                         column_names_gp = gpar(fontsize = 8, face = 'bold'),
#                         row_names_gp = gpar(fontsize = 8, face = 'bold'),
#                         heatmap_legend_param = list(
#                           title = "NES",
#                           legend_direction = "horizontal",
#                           legend_width = unit(6, "cm")),
#                         col = col_fun,
#                         # Control heatmap size
#                         width = unit(ncol(mat_d36) * 10, "mm"), 
#                         height = unit(nrow(mat_d36) * 10, "mm"))
# draw(d36, heatmap_legend_side = "bottom")
# 
# dev.off()


#### Sig pathways accross cluster really ordered by padj value####

###Day 8####
# 3. Determine pathway order
# Smallest padj observed across CM1-CM4 comes first

pathway_order_d8 <- test.db %>%
  filter(
    pathway %in% CM_d8_p,
    cluster %in% paste0("CM", 1:4, "_day8")
  ) %>%
  group_by(pathway) %>%
  summarise(
    min_padj = min(padj, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(min_padj) %>%
  pull(pathway)


# 4. Create pathway x cell-type table containing NES

test.plot.d8 <- test.db %>%
  filter(
    pathway %in% CM_d8_p,
    cluster %in% paste0("CM", 1:4, "_day8"),
    !is.na(NES)
  ) %>%
  select(pathway, cluster, NES) %>%
  pivot_wider(
    names_from = cluster,
    values_from = NES
  )


# 5. Convert to matrix
# Pathways = rows
# CM subsets = columns

mat_d8 <- test.plot.d8 %>%
  tibble::column_to_rownames("pathway") %>%
  as.matrix()


# 6. Explicitly order columns CM1 -> CM4

mat_d8 <- mat_d8[
  ,
  paste0("CM", 1:4, "_day8"),
  drop = FALSE
]

# Simplify column labels
colnames(mat_d8) <- paste0("CM", 1:4)


# 7. Explicitly order pathways by minimum padj

mat_d8 <- mat_d8[
  match(pathway_order_d8, rownames(mat_d8)),
  ,
  drop = FALSE
]

####Heatmaps####

col_fun <- colorRamp2(
  c(-2, 2),
  c("#FFFFFF", "#B2182B")   # blue - white - red (swap colors if you want)
)

pdf("output/heatmap_nex_ibsm40_monos_clasicMd8_top25_padj_fix2.pdf", height = 15, width = 15)
d8 <- ComplexHeatmap::Heatmap(mat_d8,
                              na_col = "grey",
                              cluster_rows = F,
                              cluster_columns = F,
                              column_names_rot = 0,
                              column_names_centered = TRUE,
                              column_names_side = "top",
                              column_names_gp = gpar(fontsize = 8, face = 'bold'),
                              row_names_gp = gpar(fontsize = 8, face = 'bold'),
                              heatmap_legend_param = list(
                                title = "NES",
                                legend_direction = "horizontal",
                                legend_width = unit(6, "cm")),
                              col = col_fun,
                              # Control heatmap size
                              width = unit(ncol(mat_d8) * 10, "mm"), 
                              height = unit(nrow(mat_d8) * 10, "mm"))
draw(d8, heatmap_legend_side = "bottom")

dev.off()

####Day36####
# 3. Determine pathway order
# For each selected pathway, find its smallest padj across CM1-CM4.
# Pathways with the smallest padj will appear at the TOP of the heatmap.

pathway_order_d36 <- test.db %>%
  filter(
    pathway %in% CM_d36_p,
    cluster %in% paste0("CM", 1:4, "_day36")
  ) %>%
  group_by(pathway) %>%
  summarise(
    min_padj = min(padj, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(min_padj) %>%
  pull(pathway)


# 4. Create pathway x cell-type table containing NES

test.plot.d36 <- test.db %>%
  filter(
    pathway %in% CM_d36_p,
    cluster %in% paste0("CM", 1:4, "_day36"),
    !is.na(NES)
  ) %>%
  select(pathway, cluster, NES) %>%
  pivot_wider(
    names_from = cluster,
    values_from = NES
  )


# 5. Convert to matrix
# Pathways = rows
# CM subsets = columns

mat_d36 <- test.plot.d36 %>%
  tibble::column_to_rownames("pathway") %>%
  as.matrix()


# 6. Explicitly order columns CM1 -> CM4

mat_d36 <- mat_d36[
  ,
  paste0("CM", 1:4, "_day36"),
  drop = FALSE
]

# Simplify column labels
colnames(mat_d36) <- paste0("CM", 1:4)


# 7. Explicitly order pathways by minimum padj
# match() ensures the matrix follows pathway_order_d36,
# regardless of the order produced by pivot_wider()

mat_d36 <- mat_d36[
  match(pathway_order_d36, rownames(mat_d36)),
  ,
  drop = FALSE
]

# 
# #mat_d36 <- t(test.plot.d36[1:25])
# mat_d36 <- t(test.plot.d36)
##Adjusting colos
col_fun <- colorRamp2(
  c(-2,2),
  c("#FFFFFF", "#B2182B")   # blue - white - red (swap colors if you want)
)

pdf("output/heatmap_nex_ibsm40_monos_clasicMd36_top25_padj_fix2.pdf", height = 15, width = 15)
d36 <- ComplexHeatmap::Heatmap(mat_d36,
                               na_col = "grey",
                               cluster_rows = F,
                               cluster_columns = F,
                               column_names_rot = 0,
                               column_names_centered = TRUE,
                               column_names_side = "top",
                               column_names_gp = gpar(fontsize = 8, face = 'bold'),
                               row_names_gp = gpar(fontsize = 8, face = 'bold'),
                               heatmap_legend_param = list(
                                 title = "NES",
                                 legend_direction = "horizontal",
                                 legend_width = unit(6, "cm")),
                               col = col_fun,
                               # Control heatmap size
                               width = unit(ncol(mat_d36) * 10, "mm"), 
                               height = unit(nrow(mat_d36) * 10, "mm"))
draw(d36, heatmap_legend_side = "bottom")

dev.off()

###Dotplot day 36

test.dotplot.d36 <- test.db %>%
  filter(padj < 0.01) %>%
  arrange(padj) %>%
  slice_head(n = 25) %>%
  mutate(enrichment_score = ES)%>%
dplyr::select(pathway, NES,enrichment_score, padj, size
  )

ggplot(test.dotplot.d36,
  aes(x = enrichment_score, y = pathway)) +
  geom_point(
    aes(size = size, fill = NES),
    shape = 21,
    colour = "grey30",
    stroke = 0.2
  ) +
  scale_fill_gradient(
    low = "white",
    high = "#B2182B",
    name = "NES"
  ) +
  scale_size_continuous(
    range = c(2, 10),
    name = "Gene count"
  ) +
  labs(
    x = "Fold enrichment",
    y = NULL,
    title = "Pathway enrichment analysis"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.y = element_text(size = 9),
    axis.text.x = element_text(size = 10),
    legend.position = "right"
  )

pdf("output/dotplot_nex_ibsm40_monos_clasicMd36_top25_padj.pdf", height = 15, width = 15)
dotplot_d36
dev.off()

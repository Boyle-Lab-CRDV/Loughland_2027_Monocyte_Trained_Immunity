library(readr)
library(dplyr)
library(tidyr)
library(tibble)
library(ComplexHeatmap)
library(GenomicRanges)
library(ggplot2)
library(stringr)

##We did the pathwayenrinchement analysis in the web interphase, to run the rscript version
##follow this workflow
#https://jokergoo.github.io/rGREAT/articles/local-GREAT.html#get-enrichment-table 

##prepare list of DAA
db_s1 = read.csv("output/sabah_1_m/Dream_fitmm3_group_batch_1sv_sample_results200456.csv")%>%
  filter(Coefficient %in%"groupd0"&adj.P.Val <= 0.05)%>%
  mutate(
    Comparison = gsub("group", "", Coefficient),
    studyID = paste0("s1","_",Comparison),
    studygeneID = paste0("s1","_",Comparison,"_",GeneID),
    rankings = round(sign(logFC)*(-log10(adj.P.Val)),4)# we will use the signed p values from spatial DGE as ranking
  ) %>% dplyr::select(GeneID, studygeneID, genes, logFC, adj.P.Val, studyID, rankings)


db_s2 = read.csv("output/sabah_2_m/Dream_fitmm3_group_2sv__ci95_160426.csv")%>%
  filter(Coefficient %in%"grouppostinfx"&adj.P.Val <= 0.05)%>%
  mutate(
    Comparison = case_when(Coefficient =="grouppostinfx" ~ "d28"),
    studyID = paste0("s2","_",Comparison),
    studygeneID = paste0("s2","_",Comparison,"_",GeneID),
    rankings = round(sign(logFC)*(-log10(adj.P.Val)),4)# we will use the signed p values from spatial DGE as ranking
  )%>% dplyr::select(GeneID, studygeneID, genes, logFC, adj.P.Val, studyID, rankings)

db_s3 = read.csv("output/sabah_3_m/Dream_fitmm3_group_2sv__batchsample_results150426.csv")%>%
  filter(Coefficient %in% c("groupd28","groupd0")&adj.P.Val <= 0.05)%>%
  mutate(
    Comparison = gsub("group", "", Coefficient),
    studyID = paste0("s3","_",Comparison),
    studygeneID = paste0("s3","_",Comparison,"_",GeneID),
    rankings = round(sign(logFC)*(-log10(adj.P.Val)),4)# we will use the signed p values from spatial DGE as ranking
  )%>% dplyr::select(GeneID, studygeneID, genes, logFC, adj.P.Val, studyID, rankings)

db_s4 = read.csv("output/sabah_4_m/Dream_fitmm3_group_3sv__batch_ci95_200426.csv")%>%
  filter(Coefficient %in% c("groupqem_d0", "groupqem_d28", "groupqec_neg")&adj.P.Val <= 0.05)%>%
  mutate(
    Comparison = gsub("group", "", Coefficient),
    studyID = paste0("s4","_",Comparison),
    studygeneID = paste0("s4","_",Comparison,"_",GeneID),
    rankings = round(sign(logFC)*(-log10(adj.P.Val)),4)# we will use the signed p values from spatial DGE as ranking
  )%>% dplyr::select(GeneID, studygeneID, genes, logFC, adj.P.Val, studyID, rankings)

#chmi_1_m doesn't have sig genes

db_s5 = read.csv("output/chmi_2_m/Dream_fitmm3_group_3sv__batch_ci95_150426.csv")%>%
  filter(Coefficient %in% c("groupendemic")&adj.P.Val <= 0.05)%>%
  mutate(
    Comparison = gsub("group", "", Coefficient),
    studyID = paste0("c2","_",Comparison),
    studygeneID = paste0("c2","_",Comparison,"_",GeneID),
    rankings = round(sign(logFC)*(-log10(adj.P.Val)),4)# we will use the signed p values from spatial DGE as ranking
  )%>% dplyr::select(GeneID, studygeneID, genes, logFC, adj.P.Val, studyID, rankings)



##Input require
#db.f = rbind(db_s1,db_s2,db_s3)%>% filter(!is.na(genes))
##Anno peaks
ann_peaks_s1 = read.csv("output/sabah_1_m/annotated_filtered_peaks.csv")
ann_peaks_s2 = read.csv("output/sabah_2_m/annotated_filtered_peaks.csv")
ann_peaks_s3 = read.csv("output/sabah_3_m/annotated_filtered_peaks.csv")
ann_peaks_s4 = read.csv("output/sabah_4_m/annotated_filtered_peaks.csv")
ann_peaks_s5 = read.csv("output/chmi_2_m/annotated_filtered_peaks.csv")

##Tables for chip-enrich analysis
table_s1_0 = left_join(db_s1, ann_peaks_s1, by = "GeneID")
table_s2_0 = left_join(db_s2, ann_peaks_s2, by = "GeneID")
table_s3_0 = left_join(db_s3, ann_peaks_s3, by = "GeneID")
table_s4_0 = left_join(db_s4, ann_peaks_s4, by = "GeneID")
table_s5_0 = left_join(db_s5, ann_peaks_s5, by = "GeneID")

####Preparing input file for GREAT standford enrichment tool####
###http://great.stanford.edu/public/html/

#selecting relevant variables
sabah_1_m_enrichgU = table_s1_0 %>% filter(logFC>0) %>% dplyr::select(seqnames,  start, end, distanceToTSS)%>%dplyr::rename(chr = seqnames)
write.table(
  sabah_1_m_enrichgU[, c("chr", "start", "end", "distanceToTSS")],
  file = "/working_groups/boylelab/shared/ZulyPava/formating4hpc/ATACseq_analysis/pathway_analysis/sabah_1_m_input_greatenrixh_sigDAA_annotation_220426_Up.bed",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

sabah_1_m_enrichgD = table_s1_0 %>% filter(logFC<0)%>% dplyr::select(seqnames,  start, end, distanceToTSS)%>%dplyr::rename(chr = seqnames)
write.table(
  sabah_1_m_enrichgD[, c("chr", "start", "end","distanceToTSS")],
  file = "/working_groups/boylelab/shared/ZulyPava/formating4hpc/ATACseq_analysis/pathway_analysis/sabah_1_m_input_greatenrixh_sigDAA_annotation_220426Down.bed",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)


sabah_2_m_enrichgU = table_s2_0 %>% filter(logFC>0) %>% dplyr::select(seqnames,  start, end, distanceToTSS)%>%dplyr::rename(chr = seqnames)
write.table(
  sabah_2_m_enrichgU[, c("chr", "start", "end", "distanceToTSS")],
  file = "/working_groups/boylelab/shared/ZulyPava/formating4hpc/ATACseq_analysis/pathway_analysis/sabah_2_m_input_greatenrixh_sigDAA_annotation_220426_Up.bed",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

sabah_2_m_enrichgD = table_s2_0 %>% filter(logFC<0)%>% dplyr::select(seqnames,  start, end, distanceToTSS)%>%dplyr::rename(chr = seqnames)

write.table(
  sabah_2_m_enrichgD[, c("chr", "start", "end", "distanceToTSS")],
  file = "/working_groups/boylelab/shared/ZulyPava/formating4hpc/ATACseq_analysis/pathway_analysis/sabah_2_m_input_greatenrixh_sigDAA_annotation_220426_Down.bed",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

sabah_3_m_enrichgU = table_s3_0 %>% filter(logFC>0) %>% dplyr::select(seqnames,  start, end, distanceToTSS)%>%dplyr::rename(chr = seqnames)
write.table(
  sabah_3_m_enrichgU[, c("chr", "start", "end", "distanceToTSS")],
  file = "/working_groups/boylelab/shared/ZulyPava/formating4hpc/ATACseq_analysis/pathway_analysis/sabah_3_m_input_greatenrixh_sigDAA_annotation_220426_Up.bed",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

sabah_3_m_enrichgD = table_s3_0 %>% filter(logFC<0)%>% dplyr::select(seqnames,  start, end, distanceToTSS)%>%dplyr::rename(chr = seqnames)

write.table(
  sabah_3_m_enrichgD[, c("chr", "start", "end", "distanceToTSS")],
  file = "/working_groups/boylelab/shared/ZulyPava/formating4hpc/ATACseq_analysis/pathway_analysis/sabah_3_m_input_greatenrixh_sigDAA_annotation_220426_Down.bed",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

sabah_4_m_enrichgU = table_s4_0 %>% filter(logFC>0) %>% dplyr::select(seqnames,  start, end, distanceToTSS)%>%dplyr::rename(chr = seqnames)
write.table(
  sabah_4_m_enrichgU[, c("chr", "start", "end", "distanceToTSS")],
  file = "/working_groups/boylelab/shared/ZulyPava/formating4hpc/ATACseq_analysis/pathway_analysis/sabah_4_m_input_greatenrixh_sigDAA_annotation_220426_Up.bed",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

sabah_4_m_enrichgD = table_s4_0 %>% filter(logFC<0)%>% dplyr::select(seqnames,  start, end, distanceToTSS)%>%dplyr::rename(chr = seqnames)

write.table(
  sabah_4_m_enrichgD[, c("chr", "start", "end", "distanceToTSS")],
  file = "/working_groups/boylelab/shared/ZulyPava/formating4hpc/ATACseq_analysis/pathway_analysis/sabah_4_m_input_greatenrixh_sigDAA_annotation_220426_Down.bed",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

chmi_2_m_enrichgU = table_s5_0 %>% filter(logFC>0) %>% dplyr::select(seqnames,  start, end, distanceToTSS)%>%dplyr::rename(chr = seqnames)
write.table(
  chmi_2_m_enrichgU[, c("chr", "start", "end", "distanceToTSS")],
  file = "/working_groups/boylelab/shared/ZulyPava/formating4hpc/ATACseq_analysis/pathway_analysis/chmi_2_m_input_greatenrixh_sigDAA_annotation_220426_Up.bed",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

chmi_2_m_enrichgD = table_s5_0 %>% filter(logFC<0)%>% dplyr::select(seqnames,  start, end, distanceToTSS)%>%dplyr::rename(chr = seqnames)

write.table(
  chmi_2_m_enrichgD[, c("chr", "start", "end", "distanceToTSS")],
  file = "/working_groups/boylelab/shared/ZulyPava/formating4hpc/ATACseq_analysis/pathway_analysis/chmi_2_m_input_greatenrixh_sigDAA_annotation_220426_Down.bed",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)


####Reading and formating input from great enrichment website

# #sabah_1_m_up
# 
# lines <- readLines("/working_groups/boylelab/shared/ZulyPava/formating4hpc/ATACseq_analysis/pathway_analysis/sabah_1_m_up_shown-GOBiologicalProcess.tsv")
# 
# # find the real header line
# header_idx <- grep("^# Term Name", lines)
# 
# # keep only real table
# lines_clean <- lines[header_idx:length(lines)]
# 
# sabah_1_paths <- read_delim(
#   I(lines_clean),
#   delim = "\t",
#   quote = "",
#   col_names = TRUE,
#   trim_ws = TRUE
# )
# 
# sabah_1_paths <- sabah_1_paths %>%
#   separate(
#     col = ncol(.),
#     into = c(
#       "Binom_Raw_P",
#       "Binom_FDR",
#       "Binom_Fold_Enrichment",
#       "Observed_Genes",
#       "Expected_Genes"
#     ),
#     sep = "\t"
#   )
# 
# write.csv(sabah_1_paths,"/working_groups/boylelab/shared/ZulyPava/formating4hpc/ATACseq_analysis/pathway_analysis/sabah_1_m_up_shown-GOBiologicalProcess.csv")
# 
# #sabah_1_m_down
# 
# lines <- readLines("/working_groups/boylelab/shared/ZulyPava/formating4hpc/ATACseq_analysis/pathway_analysis/sabah_1_m_down_shown-GOBiologicalProcess.tsv")
# 
# # find the real header line
# header_idx <- grep("^# Term Name", lines)
# 
# # keep only real table
# lines_clean <- lines[header_idx:length(lines)]
# 
# sabah_1_paths_d <- read_delim(
#   I(lines_clean),
#   delim = "\t",
#   quote = "",
#   col_names = TRUE,
#   trim_ws = TRUE
# )
# 
# sabah_1_paths_d <- sabah_1_paths_d %>%
#   separate(
#     col = ncol(.),
#     into = c(
#       "Binom_Raw_P",
#       "Binom_FDR",
#       "Binom_Fold_Enrichment",
#       "Observed_Genes",
#       "Expected_Genes"
#     ),
#     sep = "\t"
#   )
# 
# write.csv(sabah_1_paths_d,"/working_groups/boylelab/shared/ZulyPava/formating4hpc/ATACseq_analysis/pathway_analysis/sabah_1_m_down_shown-GOBiologicalProcess.csv")
# 
# 
# #sabah_2_m_up
# 
# lines <- readLines("/working_groups/boylelab/shared/ZulyPava/formating4hpc/ATACseq_analysis/pathway_analysis/sabah_2_m_up_shown-GOBiologicalProcess.tsv")
# 
# # find the real header line
# header_idx <- grep("^# Term Name", lines)
# 
# # keep only real table
# lines_clean <- lines[header_idx:length(lines)]
# 
# sabah_2_paths <- read_delim(
#   I(lines_clean),
#   delim = "\t",
#   quote = "",
#   col_names = TRUE,
#   trim_ws = TRUE
# )
# 
# sabah_2_paths <- sabah_2_paths %>%
#   separate(
#     col = ncol(.),
#     into = c(
#       "Binom_Raw_P",
#       "Binom_FDR",
#       "Binom_Fold_Enrichment",
#       "Observed_Genes",
#       "Expected_Genes"
#     ),
#     sep = "\t"
#   )
# 
# write.csv(sabah_2_paths,"/working_groups/boylelab/shared/ZulyPava/formating4hpc/ATACseq_analysis/pathway_analysis/sabah_2_m_up_shown-GOBiologicalProcess.csv")
# 
# head(sabah_2_paths)
# 
# #sabah_2_m_down
# 
# lines <- readLines("/working_groups/boylelab/shared/ZulyPava/formating4hpc/ATACseq_analysis/pathway_analysis/sabah_2_m_down_shown-GOBiologicalProcess.tsv")
# 
# # find the real header line
# header_idx <- grep("^# Term Name", lines)
# 
# # keep only real table
# lines_clean <- lines[header_idx:length(lines)]
# 
# sabah_2_pathsD <- read_delim(
#   I(lines_clean),
#   delim = "\t",
#   quote = "",
#   col_names = TRUE,
#   trim_ws = TRUE
# )
# 
# sabah_2_pathsD <- sabah_2_pathsD %>%
#   separate(
#     col = ncol(.),
#     into = c(
#       "Binom_Raw_P",
#       "Binom_FDR",
#       "Binom_Fold_Enrichment",
#       "Observed_Genes",
#       "Expected_Genes"
#     ),
#     sep = "\t"
#   )
# 
# write.csv(sabah_2_pathsD,"/working_groups/boylelab/shared/ZulyPava/formating4hpc/ATACseq_analysis/pathway_analysis/sabah_2_m_down_shown-GOBiologicalProcess.csv")
# 
# 
# #sabah_3_m_up
# 
# lines <- readLines("/working_groups/boylelab/shared/ZulyPava/formating4hpc/ATACseq_analysis/pathway_analysis/sabah_3_m_up_shown-GOBiologicalProcess.tsv")
# 
# # find the real header line
# header_idx <- grep("^# Term Name", lines)
# 
# # keep only real table
# lines_clean <- lines[header_idx:length(lines)]
# 
# sabah_3_paths <- read_delim(
#   I(lines_clean),
#   delim = "\t",
#   quote = "",
#   col_names = TRUE,
#   trim_ws = TRUE
# )
# 
# sabah_3_paths <- sabah_3_paths %>%
#   separate(
#     col = ncol(.),
#     into = c(
#       "Binom_Raw_P",
#       "Binom_FDR",
#       "Binom_Fold_Enrichment",
#       "Observed_Genes",
#       "Expected_Genes"
#     ),
#     sep = "\t"
#   )
# 
# write.csv(sabah_3_paths,"/working_groups/boylelab/shared/ZulyPava/formating4hpc/ATACseq_analysis/pathway_analysis/sabah_3_m_up_shown-GOBiologicalProcess.csv")
# 
# head(sabah_3_paths)
# 
# #sabah_3_m_down
# 
# lines <- readLines("/working_groups/boylelab/shared/ZulyPava/formating4hpc/ATACseq_analysis/pathway_analysis/sabah_3_m_down_shown-GOBiologicalProcess.tsv")
# 
# # find the real header line
# header_idx <- grep("^# Term Name", lines)
# 
# # keep only real table
# lines_clean <- lines[header_idx:length(lines)]
# 
# sabah_3_pathsD <- read_delim(
#   I(lines_clean),
#   delim = "\t",
#   quote = "",
#   col_names = TRUE,
#   trim_ws = TRUE
# )
# 
# sabah_3_pathsD <- sabah_3_paths %>%
#   separate(
#     col = ncol(.),
#     into = c(
#       "Binom_Raw_P",
#       "Binom_FDR",
#       "Binom_Fold_Enrichment",
#       "Observed_Genes",
#       "Expected_Genes"
#     ),
#     sep = "\t"
#   )
# 
# write.csv(sabah_3_pathsD,"/working_groups/boylelab/shared/ZulyPava/formating4hpc/ATACseq_analysis/pathway_analysis/sabah_3_m_up_shown-GOBiologicalProcess.csv")
# 
# head(sabah_3_pathsD)

# ####Plotting####
# sabah_1_paths_plot = sabah_1_paths%>%
#   rename(`# Term Name` = "Term_name", "Binom Fold Enrichment" = "Fold_Enrichment")
# 
# pdf("/working_groups/boylelab/shared/ZulyPava/formating4hpc/ATACseq_analysis/pathway_analysis/up/sabah_1_m_up_heatmap.pdf", width = 8, height = 10)
# ggplot(sabah_1_paths_plot, aes(y = Term_name, x = 1, fill = Fold_Enrichment)) +
#   geom_tile(width = 0.5) +
#   scale_fill_gradient(low = "white", high = "#B2182B", limits = c(0, NA))+
#   theme_bw()
# dev.off()
# 
# sabah_1_paths_plotd = sabah_1_paths_d%>%
#   rename(`# Term Name` = "Term_name", "Binom Fold Enrichment" = "Fold_Enrichment")
# 
# pdf("/working_groups/boylelab/shared/ZulyPava/formating4hpc/ATACseq_analysis/pathway_analysis/down/sabah_1_m_down_heatmap.pdf", width = 8, height = 10)
# ggplot(sabah_1_paths_plotd, aes(y = Term_name, x = 1, fill = Fold_Enrichment)) +
#   geom_tile(width = 0.5) +
#   scale_fill_gradient(low = "white", high = "dodgerblue4", limits = c(0, NA))+
#   theme_bw()
# dev.off()
# 
# 
# ##Sabah2
# 
# sabah_2_paths_plot = sabah_2_paths%>%
#   rename(`# Term Name` = "Term_name", "Binom Fold Enrichment" = "Fold_Enrichment")
# 
# pdf("/working_groups/boylelab/shared/ZulyPava/formating4hpc/ATACseq_analysis/pathway_analysis/up/sabah_2_m_up_heatmap.pdf", width = 8, height = 8)
# ggplot(sabah_2_paths_plot, aes(y = Term_name, x = 1, fill = Fold_Enrichment)) +
#   geom_tile() +
#   scale_fill_gradient(low = "white", high = "#B2182B", limits = c(0, NA))+
#   theme_bw()
# dev.off()
# 
# sabah_2_paths_plotd = sabah_2_pathsD%>%
#   rename(`# Term Name` = "Term_name", "Binom Fold Enrichment" = "Fold_Enrichment")
# 
# pdf("/working_groups/boylelab/shared/ZulyPava/formating4hpc/ATACseq_analysis/pathway_analysis/down/sabah_2_m_down_heatmap.pdf", width =10, height = 8)
# ggplot(sabah_2_paths_plotd, aes(y = Term_name, x = 1, fill = Fold_Enrichment)) +
#   geom_tile() +
#   scale_fill_gradient(low = "white", high = "dodgerblue4", limits = c(0, NA))+
#   theme_bw()
# dev.off()

# #UPregulated
# files_up <- list.files(path = "/working_groups/boylelab/shared/ZulyPava/formating4hpc/ATACseq_analysis/pathway_analysis/up", pattern = "*.csv", full.names = TRUE)
# 
# for (f in files_up) {
#   tmp = read.csv(f)
#   fname = basename(f)
#   outname = gsub(pattern = "csv", replacement = "pdf", x = fname)
#   d = tmp %>%
#     select(`Term Name` = `X..Term.Name`, `Fold Enrichment` = `Binom.Fold.Enrichment`) %>%
#     mutate(`Fold Enrichment` = as.integer(`Fold Enrichment`))
#   print(typeof(d$`Fold Enrichment`))
#   plot = ggplot(d, aes(y = `Term Name`, x = 1, fill = `Fold Enrichment`)) + geom_tile() +
#     scale_fill_gradient(low = "white", high = "#B2182B", limits = c(0, NA))
#   pdf(
#     plot = plot,
#     filename = file.path(
#       "/working_groups/boylelab/shared/ZulyPava/formating4hpc/ATACseq_analysis/pathway_analysis/up",
#       outname
#     )
#   )
# }
# 
# ##Down regulated
# 
# #UPregulated
# 
# files_down <- list.files(path = "/working_groups/boylelab/shared/ZulyPava/formating4hpc/ATACseq_analysis/pathway_analysis/down", pattern = "*.csv", full.names = TRUE)
# 
# for (f in files_down) {
#   tmp = read.csv(f)
#   fname = basename(f)
#   outname = gsub(pattern = "csv", replacement = "pdf", x = fname)
#   d = tmp %>%
#     select(`Term Name` = `X..Term.Name`, `Fold Enrichment` = `Binom.Fold.Enrichment`) %>%
#     mutate(`Fold Enrichment` = as.integer(`Fold Enrichment`))
#   print(typeof(d$`Fold Enrichment`))
#   plot = ggplot(d, aes(y = `Term Name`, x = 1, fill = `Fold Enrichment`)) + geom_tile() +
#     scale_fill_gradient(low = "white", high = "dodgerblue4", limits = c(0, NA))
#   ggsave(
#     plot = plot,
#     filename = file.path(
#       "/working_groups/boylelab/shared/ZulyPava/formating4hpc/ATACseq_analysis/pathway_analysis/down",
#       outname
#     )
#   )
# }


# #All included
# table_s2 <- table_s1_0 %>%
#   filter(Coefficient %in% c("groupd8", "groupd15", "groupd27"))%>%
#   mutate(Annotation =
#            case_when(grepl("Intron", annotation) ~ "Intron",
#                      grepl("Exon", annotation) ~ "Exon",
#                      TRUE ~ annotation),
#          Comparison = 
#            case_when(Coefficient == "groupd8" ~ "Day0vsDay8",
#                            Coefficient == "groupd15" ~ "Day0vsDay15",
#                            Coefficient == "groupd27" ~ "Day0vsDay27"),
#          Peakstart = start,
#          Peakend = end,
#          Peakwidth = width) %>%
#   select(SYMBOL, GENENAME, logFC, CI.L, CI.R, AveExpr,t, P.Value, adj.P.Val, B, GeneID, Comparison,
#          seqnames,  Peakstart, Peakend, Peakwidth, strand, Annotation, 
#          geneChr,	geneStart,	geneEnd,	geneLength,	geneStrand,	geneId,	transcriptId,	distanceToTSS,	ENSEMBL)
# write.csv(table_s2, "output/chmi_1_t/plots_for_paper/table_allDAA_annotation_120325.csv")
# 
# ##
# df_summary <- table_s1_sig %>%
#   mutate(annotation_clean =
#            case_when(grepl("Intron", annotation) ~ "Intron",
#                      grepl("Exon", annotation) ~ "Exon",
#                      TRUE ~ annotation)) %>%
#   group_by(Coefficient, annotation_clean) %>%
#   summarise(n = n(), .groups = "drop") %>%
#   filter(n > 0)
# 
# df_summary$Coefficient <- factor(df_summary$Coefficient, levels = c("groupd8","groupd27"))
# 
# pdf("output/chmi_1_t/plots_for_paper/Barplot_daa_anno_030326")
# ggplot(df_summary, aes(x = annotation_clean, y = n, fill = Coefficient)) +
#   geom_bar(stat = "identity", position = "dodge") +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
# dev.off()
# 
# table_s1_sig_w <- table_s1_sig %>%
#   mutate(
#     annotation_clean =
#            case_when(grepl("Intron", annotation) ~ "Intron",
#                      grepl("Exon", annotation) ~ "Exon",
#                      TRUE ~ annotation),
#     SYMBOL2=
#            case_when(duplicated(SYMBOL)~paste0(SYMBOL,".1"),
#                            TRUE~SYMBOL),
#     value=1) %>%
#   select(SYMBOL2, annotation_clean,value)%>%
#   pivot_wider(., names_from = annotation_clean,
#               id_cols = SYMBOL2,
#               values_from = value,
#               values_fn = ~1,
#               values_fill = 0)
# 
# 
# table_s1_sig_w_x <- apply(table_s1_sig_w[2:7],2, FUN = as.numeric)%>%cbind(table_s1_sig_w[,1])
# rownames(table_s1_sig_w_x) <- table_s1_sig_w_x$SYMBOL2
# table_s1_sig_w_x <- table_s1_sig_w_x%>%select(-SYMBOL2)
# table_s1_sig_w_x <- as.matrix(table_s1_sig_w_x)
# head(table_s1_sig_w_x)
# 
# heatmap(table_s1_sig_w_x)
# 
# ####Drawing a heatmap with all the information. ####
# 
# sigDGEs <- read.csv("/working_groups/boylelab/shared/ZulyPava/ATACseq_R4.4/ATAC_seq/output/chmi_1_t/plots_for_paper/table_sigDAA_annotation_030325.csv", row.names=1)
# 
# #Selecting the variables I want to plot
# sigDGEs_sub <- sigDGEs%>%
#   select(genes, logFC, annotation, distanceToTSS, Coefficient)%>%
#   mutate(
#     annotation_clean =
#       case_when(grepl("Intron", annotation) ~ "Intron",
#                 grepl("Exon", annotation) ~ "Exon",
#                 TRUE ~ annotation),
#     SYMBOL2=
#       case_when(duplicated(SYMBOL)~paste0(SYMBOL,".1"),
#                 TRUE~SYMBOL),
#     value=1) %>%
#   select(SYMBOL2, annotation_clean,value)
# 
# mat_ <- sigDGEs_sub %>%
#   select(genes, annotation_clean, logFC)%>%
#   pivot_wider(., names_from = annotation_clean,
#               id_cols = SYMBOL2,
#               values_from = value,
#               values_fn = ~1,
#               values_fill = 0)
# 
# 
# sigDGEs<- sigDGEs[sigDGEs$day %in% c("day36","day8"),]
# 
# # # Extract gene columns and merge
# # merged_genes2 <- unique(c(sort_genes_d8$gene, sort_genes_d36$gene))
# 
# # Ensure only the top 80 unique genes
# # merged_genes2 <- unique(merged_genes2[1:80])
# 
# merged_genes1 <- sigDGEs %>%
#   mutate(clust_day = paste0(sigDGEs$cluster_name, "_", sigDGEs$day))%>%
#   filter(set %in% c("d8_only", "d36_only")) %>%
#   group_by(gene)%>%
#   mutate(n_clupgene = length(cluster_name)) %>%
#   ungroup() %>%
#   filter(n_clupgene >2) %>%
#   arrange(set, dplyr::desc(logFC)) 
# 
# merged_genes1$gene = gsub("-", ".", merged_genes1$gene)
# merged_genes2 <- unique(merged_genes1$gene)
# 
# #### subsetting our selected genes
# dfsub_d816<- as.matrix(Average.sub[, colnames(Average.sub) %in% merged_genes2])
# 
# dfsub_d816 <- as.matrix(dfsub_d816)
# dfsub_d816 = scale(dfsub_d816)
# 
# col = list(Cluster = c("CM1" = "#9e0142",
#                        "CM2" = "#f4a582",
#                        "CM3" = "#5e4fa2",
#                        "CM4" = "#66c2a5"),
#            Donor = c("donor0" = "#cb7a88",
#                      "donor1" = "#ccab7b",
#                      "donor2" = "#a3cb7a",
#                      "donor3" = "#7cbecc"),
#            Day = c("day0" = "#bababa",
#                    "day8" = "#fdae61",
#                    "day36" = "blue"))
# 
# # ##We got high medium and low expression, we want 3 colors
# # rng <- max(abs(dfsub_d816), na.rm = TRUE)  # or use abs(expr_z) if you built it
# # rng <- if (is.finite(rng) && rng > 0) rng else 1  # fallback
# 
# # 3-color diverging palette: negative / zero / positive
# col_fun <- colorRamp2(
#   c(-2, 0, 4),
#   c("#2166AC", "#FFFFFF", "#B2182B")   # blue - white - red (swap colors if you want)
# )
# 
# ha <- HeatmapAnnotation(Cluster=Average.sub$cluster,
#                         Donor=Average.sub$donor,
#                         Day=Average.sub$day,
#                         col = col)
# 
# col_split <- factor(Average.sub$cluster, levels = c("CM1","CM2","CM3","CM4"))
# 
# row_split0 <- merged_genes1 %>% 
#   filter(gene %in% merged_genes2) %>%           # keep only genes in merged_genes2
#   select(gene, set) %>% 
#   distinct(gene, .keep_all = TRUE) %>%          # one row per gene
#   mutate(set = factor(set, levels = c("d8_only", "d36_only")))
# 
# # Align 'set' to the heatmap row order (genes)
# gene_order <- rownames(t(dfsub_d816))                       # rows of the heatmap
# set_map    <- setNames(as.character(row_split0$set), row_split0$gene)
# row_set    <- factor(set_map[gene_order], levels = c("d8_only", "d36_only"))
# 
# # Row annotation for 'set'
# row_ann <- rowAnnotation(
#   Set = row_set,
#   col = list(Set = c("d8_only" = "#fdae61", "d36_only" = "#0072B2")),
#   annotation_name_side = "top"
# )
# 
# pdf(here("output/heatmap_uniqueD8orD36genes_avexpdaydonorcluster_nrev1.pdf"), width = 8, height =10)
# Heatmap(
#   t(dfsub_d816),
#   name = "AvExpr", col = col_fun,
#   top_annotation = ha,
#   left_annotation = row_ann,                 # <<< add row annotation
#   cluster_rows = TRUE,                       # keep clustering order
#   show_row_dend = FALSE,                     # <<< hide row dendrogram graphic
#   cluster_columns = FALSE,
#   show_column_dend = FALSE,                  # (explicit; no column dendrogram)
#   column_split = col_split,
#   column_gap = unit(4, "mm"),
#   row_split = row_set,                       # <<< split rows by set (d8_only / d36_only)
#   row_gap = unit(1, "mm"),
#   column_title = "Genes",
#   row_title = "Day-unique genes",
#   row_names_gp = gpar(fontsize = 6),
#   column_names_gp = gpar(fontsize = 4)
# )
# 
# dev.off()


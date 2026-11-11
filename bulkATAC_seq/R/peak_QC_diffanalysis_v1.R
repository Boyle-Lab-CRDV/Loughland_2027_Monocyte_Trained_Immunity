
library(ATACseqQC)
library(DiffBind)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(GenomicAlignments)
library(Rsamtools)

#Creating my index files (Mark didn't transfer)
# For multiple BAMs
bam_files <- list.files("/working_groups/boylelab/shared/ZulyPava/formating4hpc/ATACseq_analysis/chmi_1_m/", pattern = "\\.bam$", full.names = TRUE)
lapply(bam_files, indexBam) #takes a while

##Calculating the TSS enrichment scores
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
#Turning this one into GenomicAlignments object
tx_gr <- transcripts(txdb)

# Output: named vector of TSS enrichment scores
bamfiles <- list.files("/working_groups/boylelab/shared/ZulyPava/formating4hpc/ATACseq_analysis/chmi_1_m", pattern = "\\.bam$", full.names = TRUE)

# Compute TSS enrichment scores using readBamFile first
tss_scores <- sapply(bamfiles, function(bam) {
  message("Processing: ", bam)
  
  # Read in BAM using readBamFile (for large BAMs)
  bam_obj <- readBamFile(bam, asMates = FALSE, bigFile = TRUE)#it's suppose to tell single from paired end but if TRUE in turns output into GAlignmentsList
  
  
  # Now compute TSS enrichment
  score <- TSSEscore(bam_obj, tx_gr)
  score$TSSEscore
})

# Review results
tss_scores <- setNames(tss_scores, basename(bamfiles))

write.csv(tss_scores,"/working_groups/boylelab/shared/ZulyPava/ATACseq_R4.4/ATAC_seq/input/chmi_1_m/tss_scores.csv")

##Getting the information ready to create the sample sheet
rawnames <- read.csv("input/chmi_1_m/all_fastqnames_sampleid.csv")
all_meta <- readxl::read_excel("input/chmi_1_m/Alignment and QC_indexsequence.xlsx", sheet = "Sheet1" )
file_manifest <- read.csv("input/chmi_1_m/file_manifest.csv")

##subsetting samples for this project
sample_sheet0 = file_manifest%>%filter(destination == "chmi_1_m")
head(sample_sheet0)

sample_sheet1 = sample_sheet0 %>%
  mutate(
    file_type = case_when(
      str_detect(filename, "\\.bam$") ~ "bam_file",
      str_detect(filename, "\\.gz$") ~ "narrow_peak",
      TRUE ~ "other"
    )
  ) %>% 
  tidyr::separate(., col= filename, into= c("rawname", "V1"), sep = ".t", extra = "merge", remove = FALSE) %>%
  dplyr::select(-V1)


head(sample_sheet1)

sample_sheet2 = pivot_wider(sample_sheet1, id_cols = rawname, names_from = file_type, values_from = filename )
head(sample_sheet2)

##getting the samples name
colnames(rawnames)[4] <- "fastq_file"
rawnames1 = rawnames %>% select(sample_name, fastq_file) %>%
  separate(col=fastq_file, into=c("rawname"), sep=".fastq.gz", extra="drop")

head(rawnames1)
sample_sheet3 = left_join(x=sample_sheet2, y=rawnames1, by= "rawname")
dim(sample_sheet3)
head(sample_sheet3)

##Getting the rest of the metadata
head(all_meta)
all_meta_db <- as.data.frame(all_meta)
head(all_meta_db)
colnames(all_meta_db)[1] = "sample_name"
sample_sheet4 = left_join(x=sample_sheet3, y=all_meta_db%>%select(sample_name, Batch, Study, "Part ID", Timepoint), by= "sample_name")
head(sample_sheet4)

##creating the key 
sample_sheet4$colcountnames = paste0("X.working_groups.boylelab.shared.ZulyPava.formating4hpc.ATACseq_analysis.chmi_1_m.", sample_sheet4$bam_file)
sample_sheet = sample_sheet4%>%
  mutate(old_id = gsub("-", ".", colcountnames),
         Replicate = gsub("R", "", `Part ID`),
         day = 
           case_when(Timepoint=="0" ~ "d0",
                     Timepoint=="8" ~ "d8",
                     Timepoint=="9" ~ "d8",
                     Timepoint=="15" ~ "d15",
                     Timepoint=="16" ~ "d15",
                     Timepoint=="27" ~ "d27",
                     Timepoint=="29" ~ "d27",
                     Timepoint=="35" ~ "d27"),
         New_id = paste0(`Part ID`, "_", day)
         )

head(sample_sheet)
key = sample_sheet%>%select(New_id,old_id)

write.csv(key, "/working_groups/boylelab/shared/ZulyPava/ATACseq_R4.4/ATAC_seq/output/chmi_1_m/sample_fastqname_key.csv")

##Adding tss score 
tss_scoresdb <- read.csv("/working_groups/boylelab/shared/ZulyPava/ATACseq_R4.4/ATAC_seq/input/chmi_1_m/tss_scores.csv", header = TRUE)
colnames(tss_scoresdb)[1] <- "bam_names"
colnames(tss_scoresdb)[2] <- "TSS_scores"
key$bam_names = key$old_id
key$bam_names = gsub(pattern = "X.working_groups.boylelab.shared.ZulyPava.formating4hpc.ATACseq_analysis.chmi_1_m.", replacement = "",key$bam_names )
key$bam_names = sub("^([^.]+)\\.([^.]+)\\.([^.]+)", "\\1-\\2-\\3", key$bam_names)


key0 = left_join(key, tss_scoresdb, by = "bam_names")
head(key0)
##Adding metadata
sample_sheet0 = left_join(sample_sheet, key0%>%select(New_id, TSS_scores), by = "New_id")
colnames(sample_sheet0)

#sample_sheet= sample_sheet0 %>% select(uid, lib_id, condition,sample_id, batch, trial, old_id, bed_files, TSS_scores)
#write.csv(sample_sheet,"output/chmi_1_m/metadata_wotss.csv",row.names = F)
#write.csv(sample_sheet,"output/chmi_1_m/metadata_wtss.csv",row.names = F)

sample_sheet = sample_sheet0 %>% select(New_id, sample_name, Part.ID, day, Batch, Study, bam_file, narrow_peak, Replicate,TSS_scores)
#uid -> SampleID
colnames(sample_sheet)[1]="SampleID"
colnames(sample_sheet)[2]="lib_id"
colnames(sample_sheet)[3]="sample"
colnames(sample_sheet)[4]="Factor"
colnames(sample_sheet)[5]="Condition"
colnames(sample_sheet)[6]="trial"
colnames(sample_sheet)[7]="bamReads"
colnames(sample_sheet)[8]="Peaks"
colnames(sample_sheet)[9]="Replicate"

##adding the path
sample_sheet5 = sample_sheet%>%
  mutate(bamReads=paste0("/working_groups/boylelab/shared/ZulyPava/formating4hpc/ATACseq_analysis/chmi_1_m/", bamReads),
         Peaks=paste0("/working_groups/boylelab/shared/ZulyPava/formating4hpc/ATACseq_analysis/chmi_1_m/", Peaks))
head(sample_sheet5)
write.csv(sample_sheet5,"output/chmi_1_m/sample_sheet_DiffBind.csv",row.names = FALSE)


###Continue 12/08/25

##Creating dba object | plain
# dba_0 = dba(sampleSheet = "output/sample_sheet_DiffBind.csv",
#             peakCaller = "narrow", config = "output/config_file_DiffBind.csv")

dba_0 = dba(sampleSheet = "output/chmi_1_m/sample_sheet_DiffBind.csv",
            peakCaller = "narrow")

#Making sure comparisons are againts day 0
dba_0$samples$Factor <- relevel(as.factor(dba_0$samples$Factor), ref = "d0")

##Save files
dba.save(DBA=dba_0, file="DBA_prefilt_precounts_monos_120825", dir="/working_groups/boylelab/shared/ZulyPava/ATACseq_R4.4/ATAC_seq/output/chmi_1_m/", pre="dba_", ext="RData", 
         bRemoveAnalysis=FALSE, bRemoveBackground=FALSE,
         bCompress=FALSE)


##We're going to calculate the counts only taking peaks present in at least 4 samples.


# # === 1. Annotate peaks ===
peak_file="/working_groups/boylelab/shared/ZulyPava/formating4hpc/ATACseq_analysis/featurecounts/chmi_1_m_featurecounts/all_peaks.bed"
# ## We can use peaks already filtered to be present in at least 4 samples. 
# 
# peak_file="/working_groups/boylelab/shared/ZulyPava/formating4hpc/ATACseq_analysis/featurecounts/chmi_1_t_featurecounts/filtered_peaks.bed"
peaks <- readPeakFile(peak_file)

####Assessing the reproducibility of peaks across replicates####
## Counting using the merged peak file
dbObj_counts <- dba.count(
  DBA=dba_0, 
  peaks=peaks, 
  minOverlap=FALSE, #Already done bedtools.Only include peaks found in at least 4 samples. Prevents single-sample artifacts.
  score=DBA_SCORE_FDR,#Use normalized read counts (e.g., RPKM-like values) or DBA_SCORE_READS raw counts.Will do TMM norm inside edger
  summits = FALSE, #Resizes the peak to 100 bp width (±100 bp from the summit).#based on plots
  filter=0.4, #Remove low-signal peaks
  bRemoveDuplicates=FALSE,#already done in encode
  bScaleControl=FALSE,#only for ChIPseq
  bSubControl = FALSE, #tells DiffBind to subtract greylist regions from control samples when counting reads.
  mapQCth=FALSE, minCount=0,
  bLog=FALSE,
  bUseSummarizeOverlaps=TRUE, # tells DiffBind to use paired-end logic
  readFormat=DBA_READS_DEFAULT,
  bParallel=FALSE)

#reading to calculate frip 170826
dbObj_counts <- dba.load(file='DBA_counts_monos_120825', dir='/working_groups/boylelab/shared/ZulyPava/ATACseq_R4.4/ATAC_seq/output/chmi_1_m/', ext='RData')

FRiP <- dba.show(dbObj_counts)

write.csv(FRiP,"/working_groups/boylelab/shared/ZulyPava/ATACseq_R4.4/ATAC_seq/output/chmi_1_m/chmi_1_m_FRiP_scores_filtpeaks.csv") 

dba.save(DBA=dbObj_counts, file="DBA_counts_monos_120825", dir="/working_groups/boylelab/shared/ZulyPava/ATACseq_R4.4/ATAC_seq/output/chmi_1_m/", pre="dba_", ext="RData", 
         bRemoveAnalysis=FALSE, bRemoveBackground=FALSE,
         bCompress=FALSE)



# 2. Set up the contrast by Factor, block by Condition
dbObj_DE <- dba.contrast(
  dbObj_counts,
  design = FALSE,
  categories = DBA_FACTOR,
  block = DBA_CONDITION  # block by Condition (e.g., batch)
)

# Optional: check contrast and design
# dba.show(dbObj_DE, bDesign = TRUE)
# dba.show(dbObj_DE, bContrast = TRUE)

# 3. Analyze using the model
dbObj_analyzed <- dba.analyze(
  dbObj_DE,
  design = "~Factor",
  method = DBA_ALL_METHODS,
  bParallel=FALSE
    
)

dba.save(DBA=dbObj_analyzed, file="DBA_DAAres_withreps_120825", dir="/working_groups/boylelab/shared/ZulyPava/ATACseq_R4.4/ATAC_seq/output/chmi_1_m/", pre="dba_", ext="RData", 
         bRemoveAnalysis=FALSE, bRemoveBackground=FALSE,
         bCompress=FALSE)

# 4. Extract results
dba.show(dbObj_analyzed, bContrast = TRUE)

# Group Samples Group2 Samples2 Block1 Blk1Samps Block2 Blk2Samps DB.edgeR DB.edgeR.block DB.DESeq2
# 1    d0       6     d8        6      2         6      3         6        0              0         0
# 2    d0       6    d15        6      2         6      3         6        0              0         0
# 3    d0       6    d27        5      2         5      3         6        0              0         0
# 4    d8       6    d15        6      2         6      3         6        0              1         0
# 5    d8       6    d27        5      2         5      3         6        0              0         0
# 6   d15       6    d27        5      2         5      3         6        0              0         0
# DB.DESeq2.block
# 1               0
# 2               0
# 3               0
# 4               0
# 5               0
# 6               0

##Doing it without a block!
# 2. Set up the contrast by Factor, block by sample
dbObj_DEsa <- dba.contrast(
  dbObj_counts,
  design = FALSE,
  categories = DBA_FACTOR,
  block = DBA_REPLICATE  # block by Condition (e.g., batch)
)

# Optional: check contrast and design
# dba.show(dbObj_DE, bDesign = TRUE)
# dba.show(dbObj_DE, bContrast = TRUE)

# 3. Analyze using the model
dbObj_analyzednb <- dba.analyze(
  dbObj_DEsa,
  design = "~Factor",
  method = DBA_ALL_METHODS,
  bParallel=FALSE
  
)

dba.save(DBA=dbObj_analyzednb, file="DBA_DAAres_bloSamp_120825", dir="/working_groups/boylelab/shared/ZulyPava/ATACseq_R4.4/ATAC_seq/output/chmi_1_m/", pre="dba_", ext="RData", 
         bRemoveAnalysis=FALSE, bRemoveBackground=FALSE,
         bCompress=FALSE)


dba.show(dbObj_analyzednb, bContrast = TRUE)

dbObj_DEnoblock <- dba.contrast(
  dbObj_counts,
  design = FALSE,
  categories = DBA_FACTOR
)

dbObj_analyzed_noblock <- dba.analyze(
  dbObj_DEnoblock,
  design = "~Factor",
  method = DBA_ALL_METHODS,
  bParallel=FALSE
  
)

dba.show(dbObj_analyzed_noblock, bContrast = TRUE)

# Group Samples Group2 Samples2 DB.edgeR DB.DESeq2
# 1    d0       6     d8        6        0         0
# 2    d0       6    d15        6        0         0
# 3    d0       6    d27        5        0         0
# 4    d8       6    d15        6        0         0
# 5    d8       6    d27        5        0         0
# 6   d15       6    d27        5        0         0

all8  <- dba.report(dbObj_analyzed_noblock, contrast=1, th=1) # all stats day8 vs day0
all15 <- dba.report(dbObj_analyzed_noblock, contrast=2, th=1)
all27 <- dba.report(dbObj_analyzed_noblock, contrast=3, th=1)

sapply(list(all8,all15,all27), length)    # how many peaks were *tested*?
# 1305570 1291304 1296399

dba.save(DBA=dbObj_analyzed_noblock, file="DBA_DAAres_noblock_120825", dir="/working_groups/boylelab/shared/ZulyPava/ATACseq_R4.4/ATAC_seq/output/chmi_1_m/", pre="dba_", ext="RData", 
         bRemoveAnalysis=FALSE, bRemoveBackground=FALSE,
         bCompress=FALSE)
###Volcano plot
##day 0 vs day 15
pdf("output/chmi_1_m/plots/DiffBind_volcano_day0vsday8_130825_deseq.pdf", width = 5, height = 5)
dba.plotVolcano(dbObj_counts_conf_an2, contrast=1, method=DBA_DESEQ2,
                th=0.9, bUsePval=FALSE,
                fold=0, factor="", bFlip=FALSE,
                bLabels=TRUE, maxLabels=50, dotSize=1)
dev.off()

pdf("output/chmi_1_m/plots/DiffBind_volcano_day0vsday15_130825_deseq.pdf", width = 5, height = 5)
dba.plotVolcano(dbObj_counts_conf_an2, contrast=2, method=DBA_DESEQ2,
                th=0.05, bUsePval=FALSE,
                fold=0, factor="", bFlip=FALSE,
                bLabels=TRUE, maxLabels=50, dotSize=1)
dev.off()

pdf("output/chmi_1_m/plots/DiffBind_volcano_day0vsday27_130825_deseq.pdf", width = 5, height = 5)
dba.plotVolcano(dbObj_counts_conf_an2, contrast=3, method=DBA_DESEQ2,
                th=0.05, bUsePval=TRUE,
                fold=0, factor="", bFlip=FALSE,
                bLabels=TRUE, maxLabels=50, dotSize=1)
dev.off()
##Extracting the results

# Set FDR threshold
fdr_cutoff <- 0.05

# Initialize list to store results
all_results <- list()

# Get total number of contrasts
n_contrasts <- length(dbObj_counts_conf_an2$contrasts)

# Loop through all contrasts
for (i in 1:n_contrasts) {
  cat("\n--- Processing contrast", i, "---\n")
  
  # Initialize sublist for both methods
  contrast_result <- list()
  
  # Get DESeq2 results
  res_deseq2 <- tryCatch({
    dba.report(dbObj_counts_conf_an2, contrast=i, method=DBA_DESEQ2, th=fdr_cutoff)
  }, error = function(e) NULL)
  
  # Get edgeR results
  res_edger <- tryCatch({
    dba.report(dbObj_counts_conf_an2, contrast=i, method=DBA_EDGER, th=fdr_cutoff)
  }, error = function(e) NULL)
  
  # Save results if any were found
  if (!is.null(res_deseq2)) {
    cat("DESeq2 peaks:", length(res_deseq2), "\n")
    contrast_result$DESeq2 <- res_deseq2
  } else {
    cat("No DESeq2 peaks found.\n")
  }
  
  if (!is.null(res_edger)) {
    cat("edgeR peaks:", length(res_edger), "\n")
    contrast_result$edgeR <- res_edger
  } else {
    cat("No edgeR peaks found.\n")
  }
  
  # Store this contrast's results
  all_results[[paste0("contrast_", i)]] <- contrast_result
}


##Annotate the peaks using Chipseeker.
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# Annotate peaks
#contrast_result_anno <- annotatePeak(contrast_result, TxDb=txdb, annoDb="org.Hs.eg.db")

# List to hold annotated results
all_results_anno <- list()

# Loop through all results
for (contrast_name in names(all_results)) {
  contrast_data <- all_results[[contrast_name]]
  contrast_anno <- list()  # Sublist for annotated peaks
  
  cat("\nAnnotating", contrast_name, "\n")
  
  # Annotate DESeq2 result if available
  if (!is.null(contrast_data$DESeq2)) {
    cat(" - DESeq2 peaks: ", length(contrast_data$DESeq2), "\n")
    contrast_anno$DESeq2 <- annotatePeak(contrast_data$DESeq2,
                                         TxDb = txdb,
                                         annoDb = "org.Hs.eg.db")
  }
  
  # Annotate edgeR result if available
  if (!is.null(contrast_data$edgeR)) {
    cat(" - edgeR peaks: ", length(contrast_data$edgeR), "\n")
    contrast_anno$edgeR <- annotatePeak(contrast_data$edgeR,
                                        TxDb = txdb,
                                        annoDb = "org.Hs.eg.db")
  }
  
  # Store annotated results
  all_results_anno[[contrast_name]] <- contrast_anno
}

# Convert to data.frame
df <- as.data.frame(all_results_anno$contrast_3$DESeq2)
head(df[, c("seqnames", "start", "end", "annotation", "SYMBOL")])



#dba_0 <- dba.load(file='prefilt_precounts_300625', dir='/working_groups/boylelab/shared/ZulyPava/ATACseq_R4.4/ATAC_seq/output', ext='RData')
#dbObj_counts<- dba.load(file='filt_counts_300625', dir='/working_groups/boylelab/shared/ZulyPava/ATACseq_R4.4/ATAC_seq/output', ext='RData')

#### Plotting ####
#We can verify that there is no obvious batch effect with a PCA plot:
pdf("output/chmi_1_m/plots/QC/DiffBind_PCA_batch_120825.pdf", width = 5, height = 5)
dba.plotPCA(DBA=dbObj_counts, attributes = DBA_CONDITION)
dev.off()
pdf("output/chmi_1_m/plots/QC/DiffBind_PCA_day_120825.pdf", width = 5, height = 5)
dba.plotPCA(DBA=dbObj_counts, attributes = DBA_FACTOR)
dev.off()

pdf("output/chmi_1_m/plots/QC/DiffBind_PCA_sample_120825.pdf", width = 5, height = 5)
dba.plotPCA(DBA=dbObj_counts, attributes = DBA_REPLICATE)
dev.off()

colors_df = c(R11003_d0 = "#990F0F", R11003_d15 = "#CC5151", R11003_d29= "#FFB2B2",
              R102_d0 = "#99540F", R102_d8= "#B26F2C", R102_d15="#E5B17E",
              R108_d0= "#6B990F",  R108_d15= "#85B22C", R108_d27="#E5FFB2",
              R11001_d0="#0F6B99", R11001_d8= "#2C85B2",  R11001_d15= "#7EC3E5" ,R11001_d29= "#B2E5FF",
              R31002_d0= "#260F99",R31002_d9="#6551CC", R31002_d16="#8F7EE5",  R31002_d35="#BFB2FF",
              R105_d0="#F100F1",  R105_d8="#CC5800", R105_d15="#FFAD66", R105_d27= "#FFAD66")

pdf("output/chmi_1_m/plots/QC/DiffBind_PCA_day_sample_120825.pdf", width = 5, height = 5)
dba.plotPCA(DBA=dbObj_counts, attributes = DBA_ID, vColors = colors_df)
dev.off()
##organising samples first
dbObj_counts$samples[, c("SampleID", "DBA_ID")]

dba.plotHeatmap(
  dbObj_counts,
  attributes = c(DBA_ID, DBA_FACTOR),
  bUsePval = FALSE,
  score = DBA_SCORE_NORMALIZED)
  
pdf("output/chmi_1_m/plots/QC/DiffBind_corr_heatmap_allsamples_normalised_120825.pdf", width = 10, height = 20)
dba.plotHeatmap(dbObj_counts, attributes = DBA_ID, bUsePval = FALSE, score = DBA_SCORE_NORMALIZED)
dev.off()


# #DBA_SCORE_RPKM
# #DBA_SCORE_NORMALIZED
# #A correlation heatmap (in which correlations values between each pair of 
# #columns in the binding matrix is calculated)
# dbObj_counts$samples$Factor <- relevel(as.factor(dbObj_counts$samples$Factor), ref = "d0")
# 
# pdf("output/plots/DiffBind_heatmap_rpkm_corr_300625.pdf", width = 5, height = 5)
# dba.plotHeatmap(dbObj_counts, ColAttributes=DBA_FACTOR, correlations=TRUE, score= DBA_SCORE_RPKM)
# dev.off()
# 
# ##Adding config file
# 
# dbObj_counts_conf <- dba(dbObj_counts,
#                   config=data.frame(AnalysisMethod=DBA_EDGER,th=0.05,
#                   DataType=DBA_DATA_GRANGES, RunParallel=FALSE, 
#                   minQCth=15, fragmentSize=125, 
#                   bCorPlot=FALSE, reportInit="DBA", 
#                   bUsePval=FALSE, design=,design="~Replicate+Factor",
#                   doBlacklist=TRUE, doGreylist=TRUE)
#                         )
# 
# ###Contrast
# dbObj_DE <- dba.contrast(dbObj_counts_conf, categories=DBA_FACTOR, minMembers = 2, )
# dba.show(dbObj_DE, bContrast=TRUE)
# 
# #actual DE
# dbObj_counts_conf_an2 <- dba.analyze(dbObj_DE, 
#                                     method=DBA_ALL_METHODS,
#                                     bBlacklist =TRUE,
#                                     bGreylist =FALSE,
#                                     bRetrieveAnalysis=FALSE,
#                                     bParallel=FALSE)
# 
# dba.show(dbObj_counts_conf_an2, bContrast=TRUE)
# 
# # Get differentially bound sites
# diff_peaks <- dba.report(dbObj_counts_conf_an2, th=0.05)  # FDR threshold
# 
# ##Annotate the peaks using Chipseeker.
# txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
# 
# # Annotate peaks
# peak_anno <- annotatePeak(diff_peaks, TxDb=txdb, annoDb="org.Hs.eg.db")
# 
# #plot genomic function
# 
# 
# ##xtracting counts for annotations
# # meta_TSS <- as.data.frame(dbObj_counts_conf_an2$samples)
# # peak_counts <- as.data.frame(dbObj_counts_conf_an2$binding)
# # write.csv(peak_counts, "output/diff_bind_peak_count_1000summit_300625.csv")
# # write.csv(meta_TSS, "output/diff_bind_meta_TSS_300625.csv")
# 
# dba.plotVenn(dbObj_counts_conf_an2,contrast=6,
#              method=DBA_EDGER,
#              th=0.1, bUsePval=FALSE,
#              bDB=TRUE,bNotDB=TRUE)
# 
# #use this to see the contrast numbers
# dba.show(dbObj_counts_conf_an2, bContrast=TRUE)
# 
# #Plot MA
# ##day 0 vs day 15
# dba.plotMA(dbObj_counts_conf_an2, contrast = 1, method = DBA_EDGER, bUsePval=TRUE, bNormalized = TRUE)
# 
# ##day 0 vs day 27
# dba.plotMA(dbObj_counts_conf_an, contrast = 2, method = DBA_EDGER, bUsePval=TRUE, bNormalized = TRUE)
# 
# ##day 0 vs day 8
# dba.plotMA(dbObj_counts_conf_an, contrast = 3, method = DBA_EDGER, bUsePval=TRUE, bNormalized = TRUE)
# 
# ##day 15 vs day 27
# dba.plotMA(dbObj_counts_conf_an, contrast = 4, method = DBA_EDGER, bUsePval=TRUE, bNormalized = TRUE)
# 
# ##day 15 vs day 8
# dba.plotMA(dbObj_counts_conf_an, contrast = 5, method = DBA_EDGER, bUsePval=TRUE, bNormalized = TRUE)
# 
# ##day 8 vs day 27
# dba.plotMA(dbObj_counts_conf_an, contrast = 6, method = DBA_EDGER, bUsePval=TRUE, bNormalized = TRUE)
# 
# ###Volcano plot
# ##day 0 vs day 15
# dba.plotVolcano(dbObj_counts_conf_an2, contrast=1, method=DBA_EDGER, 
#                 th=0.05, bUsePval=FALSE, 
#                 fold=0, factor="", bFlip=FALSE, 
#                 bLabels=FALSE, maxLabels=50, dotSize=1)
# 
# ##day 0 vs day 27
# dba.plotVolcano(dbObj_counts_conf_an, contrast=2, method=DBA_EDGER, 
#                 th=0.05, bUsePval=FALSE, 
#                 fold=0, factor="", bFlip=FALSE, 
#                 bLabels=FALSE, maxLabels=50, dotSize=1)
# ##day 0 vs day 8
# dba.plotVolcano(dbObj_counts_conf_an, contrast=3, method=DBA_EDGER, 
#                 th=0.05, bUsePval=TRUE, 
#                 fold=0, factor="", bFlip=FALSE, 
#                 bLabels=FALSE, maxLabels=50, dotSize=1)
# ##day 15 vs day 27
# dba.plotVolcano(dbObj_counts_conf_an, contrast=4, method=DBA_EDGER, 
#                 th=0.05, bUsePval=TRUE, 
#                 fold=0, factor="", bFlip=FALSE, 
#                 bLabels=FALSE, maxLabels=50, dotSize=1)
# ##day 15 vs day 8
# dba.plotVolcano(dbObj_counts_conf_an, contrast=5, method=DBA_EDGER, 
#                 th=0.05, bUsePval=TRUE, 
#                 fold=0, factor="", bFlip=FALSE, 
#                 bLabels=FALSE, maxLabels=50, dotSize=1)
# ##day 8 vs day 27
# dba.plotVolcano(dbObj_counts_conf_an, contrast=6, method=DBA_EDGER, 
#                 th=0.05, bUsePval=TRUE, 
#                 fold=0, factor="", bFlip=FALSE, 
#                 bLabels=FALSE, maxLabels=50, dotSize=1)
# 
# ##Box plots to se how the reads are distributed amongst the different classes
# ## of differentially bound sites and sample groups, we can use a boxplot:
# 
# 
# ##day 0 vs day 15
# pdf("output/plots/boxplots_counts_bindSites_d0_vs_d15_pVal0.05.pdf")
# pvals_d0V15 <- dba.plotBox(dbObj_counts_conf_an, contrast = 1, method = DBA_EDGER, th=0.05, bUsePval = TRUE)
# dev.off
# 
# ##day 0 vs day 27
# pvals_d0V27 <- dba.plotBox(dbObj_counts_conf_an, contrast = 2, method = DBA_EDGER, th=0.05, bUsePval = TRUE)
# 
# ##day 0 vs day 8
# pvals_d0V8 <- dba.plotBox(dbObj_counts_conf_an, contrast = 3, method = DBA_EDGER, th=0.05, bUsePval = TRUE)
# 
# ##day 15 vs day 27
# pvals_d15V27 <- dba.plotBox(dbObj_counts_conf_an, contrast = 4, method = DBA_EDGER, th=0.05, bUsePval = TRUE)
# 
# ##day 15 vs day 8
# pvals_d15V27 <- dba.plotBox(dbObj_counts_conf_an, contrast = 5, method = DBA_EDGER, th=0.05, bUsePval = TRUE)
# 
# ##day 8 vs day 27
# pvals_d15V27 <- dba.plotBox(dbObj_counts_conf_an, contrast = 6, method = DBA_EDGER, th=0.05, bUsePval = TRUE)
# 
# 
# 
# res_deseq <- dba.report(dbObj, method=DBA_DESEQ2, contrast = 1, th=1)
# 
# db_report <- dba.report(dbObj_counts_conf_an, contrast = c("1","2","3","4","5","6"))
# 
# 
# #### Plot heatmap of correlations
# ##day 0 vs day 15
# pdf("output/plots/heatmap_correlation_d0_vs_d15_pVal0.05.pdf")
# plot(dbObj_counts_conf_an2, contrast=1, method = DBA_EDGER, th=0.05, bUsePval = TRUE)
# dev.off()
# ##day 0 vs day 27
# pdf("output/plots/heatmap_correlation_d0_vs_d27_pVal0.05.pdf")
# plot(dbObj_counts_conf_an2, contrast=2, method = DBA_EDGER, th=0.05, bUsePval = TRUE)
# dev.off()
# ##day 0 vs day 8
# pdf("output/plots/heatmap_correlation_d0_vs_d8_pVal0.05.pdf")
# plot(dbObj_counts_conf_an2, contrast=3, method = DBA_EDGER, th=0.05, bUsePval = TRUE)
# dev.off()
# ##day 15 vs day 27
# pdf("output/plots/heatmap_correlation_d15_vs_d27_pVal0.05.pdf")
# plot(dbObj_counts_conf_an2, contrast=4, method = DBA_EDGER, th=0.05, bUsePval = TRUE)
# dev.off()
# ##day 15 vs day 8
# pdf("output/plots/heatmap_correlation_d15_vs_d8_pVal0.05.pdf")
# plot(dbObj_counts_conf_an2, contrast=5, method = DBA_EDGER, th=0.05, bUsePval = TRUE)
# dev.off()
# ##day 8 vs day 27
# pdf("output/plots/heatmap_correlation_d8_vs_d27_pVal0.05.pdf")
# plot(dbObj_counts_conf_an2, contrast=6, method = DBA_EDGER, th=0.05, bUsePval = TRUE)
# dev.off()
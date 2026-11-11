# ---- Auto-save workspace if an error occurs ----
options(error = function() {
  cat("ERROR detected — saving workspace...\n")
  save.image(file = "/working_groups/boylelab/shared/ZulyPava/formating4hpc/ATACseq_analysis/R_qc/workspace_on_error.RData")
  q(status = 1)
})


suppressPackageStartupMessages({
library(ATACseqQC)
library(DiffBind)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(dplyr)
library(tibble)
library(stringr)
library(GenomicAlignments)
library(Rsamtools)
library(cowplot)
})

# -----------------------------
# 1) List processed BAMs
# -----------------------------
processed_dir <- "/working_groups/boylelab/shared/ZulyPava/formating4hpc/ATACseq_analysis/chmi_1_t"
processed_bams <- list.files(processed_dir, pattern="\\.bam$", full.names=TRUE)

processed_base <- basename(processed_bams)

# Extract sample ID prefix: 25-01137-80_S80
processed_id <- sub("_L01.*", "", processed_base)

# -----------------------------
# 2) List unprocessed BAMs
# -----------------------------
unprocessed_dir <- "/working_groups/boylelab/shared/ZulyPava/formating4hpc/ATACseq_analysis/QC_160226/bam"
all_bams <- list.files(unprocessed_dir, pattern="\\.bam$", full.names=TRUE)

all_base <- basename(all_bams)

# Remove .bam extension
unprocessed_id <- sub("\\.sorted.bam$", "", all_base)

# -----------------------------
# 3) Match by ID
# -----------------------------
idx <- match(processed_id, unprocessed_id)

# Report missing matches
missing <- processed_id[is.na(idx)]
if (length(missing) > 0) {
  cat("These samples were NOT found in unprocessed folder:\n")
  cat(paste0("  - ", missing), sep="\n")
  cat("\n")
}

# Keep only matched
keep <- !is.na(idx)
bam_files <- all_bams[idx[keep]]
bamindex_files <- paste0(bam_files, ".bai")

# 7) Build BamFileList (safe)

# bamfiles <-BamFileList(path= bam_files, index= bamindex_files)
bamfilesP <-BamFileList(bam_files, bamindex_files)

#Reading the reference and making it into transcripts
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
#Turning this one into GenomicAlignments object
tx_gr <- transcripts(txdb)


####Calculating NFR score ####
# Nucleosome Free Regions (NFR) score
# NFR score is a ratio between cut signal adjacent to TSS
# and that flanking the corresponding TSS. Each TSS window of 400 bp is 
# first divided into 3 sub-regions: the most upstream 150 bp (n1), 
# the most downstream of 150 bp (n2), and the middle 100 bp (nf). 
# Then the number of fragments with 5’ ends overlapping each region 
# are calculated for each TSS. The NFR score for each TSS is calculated as NFR-score = log2(nf) - log2((n1+n2)/2). A plot can be generated with the NFR scores as Y-axis and the average signals of 400 bp window as X-axis, very like a MA plot for gene expression data.

####Preparing tags ####
## bamfile tags to be read in
possibleTag <- list("integer"=c("AM", "AS", "CM", "CP", "FI", "H0", "H1", "H2", 
                                "HI", "IH", "MQ", "NH", "NM", "OP", "PQ", "SM",
                                "TC", "UQ"), 
                    "character"=c("BC", "BQ", "BZ", "CB", "CC", "CO", "CQ", "CR",
                                  "CS", "CT", "CY", "E2", "FS", "LB", "MC", "MD",
                                  "MI", "OA", "OC", "OQ", "OX", "PG", "PT", "PU",
                                  "Q2", "QT", "QX", "R2", "RG", "RX", "SA", "TS",
                                  "U2"))
##Doing for all the bam files
tags <- lapply(seq_along(bamfilesP), function(i){
  message("Processing: ", i)
  # Getting indiv bams
  bam <- bamfilesP[[i]]
  # Read in BAM using readBamFile (for large BAMs)
  bam_obj <- scanBam(bam, yieldSize = 100, param = ScanBamParam(tag=unlist(possibleTag)))[[1]]$tag
  # Filtering reads
  names(bam_obj)[lengths(bam_obj)>0]

})


## files will be output into outPath
outPath <- "/working_groups/boylelab/shared/ZulyPava/formating4hpc/ATACseq_analysis/R_qc/output/chmi_1_t/shifted_bams"
 
# ##Doing the shifting for all 
si <- seqinfo(txdb)
seqlev <- paste0("chr", 1:22)  # consider adding X/Y/M
seqlens <- seqlengths(si)[seqlev]

which.all <- GRanges(seqnames = seqlev,
                     ranges = IRanges(start = 1, end = seqlens))

##We do the shifting for each bamfiles and store each of them in outpath
gal1 <- lapply(seq_along(bamfilesP), function(i) {
  nm <- names(bamfilesP)[i] %||% paste0("sample", i)  # if names missing
  message("Processing: ", nm)
  
  bam  <- bamfilesP[[i]]$path
  tagi <- tags[[i]]
  
  gal <- readBamFile(bam, tag = tagi, which = which.all,
                     asMates = TRUE, bigFile = TRUE)
  
  shiftedBamfile <- file.path(outPath, paste0(nm, ".shifted.bam"))
  shiftGAlignmentsList(gal, outbam = shiftedBamfile)
  
  shiftedBamfile   # return path
})
##This step took 48 hours so, I'm going to save the rds
#saveRDS(gal1, "/working_groups/boylelab/shared/ZulyPava/formating4hpc/ATACseq_analysis/R_qc/output/chmi_1_t/shifted_bams/chmi_1_t_gal1_090226.RDS")
saveRDS(gal1, "/working_groups/boylelab/shared/ZulyPava/formating4hpc/ATACseq_analysis/R_qc/output/chmi_1_t/shifted_bams/all_fastq_gal1_230226.RDS")


shifted_bams <- unlist(gal1, use.names = TRUE)

names(shifted_bams) <- names(bamfilesP)

nfr_scores2 <- lapply(names(shifted_bams), function(nm) {
  message("Processing: ", nm)
  
  # Read shifted BAM into a GAlignments
  ga <- readGAlignments(
    shifted_bams[[nm]],
    param = ScanBamParam(which = which.all)  # uses your existing regions
  )
  
  # Compute NFR score
  NFRscore(ga, tx_gr)
})

nfr_scores3 <- nfr_scores2
names(nfr_scores3) <- names(shifted_bams)

nfr_median <- names(nfr_scores3)

##Getting the median and the median coverage
nfr_median_all <- sapply(seq_along(nfr_scores3), function(x){
  message("Processing: ", x)
  nfr_median[x] <- median(nfr_scores3[[x]]$NFR_score, na.rm = TRUE)
  })

names(nfr_median_all) <- names(shifted_bams)

median_coverage <- sapply(seq_along(nfr_scores3), function(x){
  message("Processing: ", x)
  median(nfr_scores3[[x]]$log2meanCoverage, na.rm = TRUE)
})
names(median_coverage) <- names(shifted_bams)

nfr_cov <- cbind(nfr_median_all, median_coverage)

write.csv(nfr_cov,"/working_groups/boylelab/shared/ZulyPava/formating4hpc/ATACseq_analysis/R_qc/output/chmi_1_t/chmi_1_t_nfr_cov.csv")

##Checking the PCR duplication rate, and mit reads contamination
pcrdup <- lapply(seq_along(bamfilesP), function(i) {
  nm <- names(bamfilesP)[i] %||% paste0("sample", i)  # if names missing
  message("Processing: ", nm)
  
  bam  <- bamfilesP[[i]]$path
  
  bamQC(bam, outPath = NULL)
  
})

## name pcrdup robustly
sample_names <- vapply(seq_along(bamfilesP), function(i) {
  names(bamfilesP)[i] %||% paste0("sample", i)
}, FUN.VALUE = character(1))
names(pcrdup) <- sample_names

## initialize numeric vectors
nrf   <- setNames(rep(NA_real_, length(sample_names)), sample_names)
pcrC1 <- setNames(rep(NA_real_, length(sample_names)), sample_names)
pcrC2 <- setNames(rep(NA_real_, length(sample_names)), sample_names)

## extract
for (nm in sample_names) {
  message("Extracting: ", nm)

  nrf[nm]   <- pcrdup[[nm]]$nonRedundantFraction
  pcrC1[nm] <- pcrdup[[nm]]$PCRbottleneckCoefficient_1
  pcrC2[nm] <- pcrdup[[nm]]$PCRbottleneckCoefficient_2
}

## write once (recommended)
out_dir <- "/working_groups/boylelab/shared/ZulyPava/formating4hpc/ATACseq_analysis/R_qc/output/chmi_1_t"

write.csv(data.frame(sample = names(nrf), NRF = as.numeric(nrf)),
          file.path(out_dir, "chmi_1_t_NRF.csv"), row.names = FALSE)

write.csv(data.frame(sample = names(pcrC1), PCR_C1 = as.numeric(pcrC1)),
          file.path(out_dir, "chmi_1_t_PCR_C1.csv"), row.names = FALSE)

write.csv(data.frame(sample = names(pcrC2), PCR_C2 = as.numeric(pcrC2)),
          file.path(out_dir, "chmi_1_t_PCR_C2.csv"), row.names = FALSE)


#Compute TSS enrichment scores using readBamFile first !!!Takes a while only ran it once!!!

# Make sure bamfilesP has names (optional but helpful)
if (is.null(names(bamfilesP)) || all(names(bamfilesP) == "")) {
  names(bamfilesP) <- basename(path(bamfilesP))
}

tss_scores <- numeric(length(bamfilesP))
names(tss_scores) <- names(bamfilesP)

for (i in seq_along(bamfilesP)) {
  nm <- names(bamfilesP)[i]
  if (is.null(nm) || is.na(nm) || nm == "") nm <- paste0("sample", i)

  message("Processing: ", nm)

  bam <- path(bamfilesP[i])  # safe path extraction

  bam_obj <- readBamFile(bam, asMates = FALSE, bigFile = TRUE)
  score <- TSSEscore(bam_obj, tx_gr)
  tss_scores[i] <- score$TSSEscore
}

out_file <- "/working_groups/boylelab/shared/ZulyPava/formating4hpc/ATACseq_analysis/R_qc/output/chmi_1_t/chmi_1_t_tss_scores.csv"
dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)

write.csv(data.frame(sample = names(tss_scores), TSSEscore = tss_scores),
          out_file, row.names = FALSE)

cat("Script completed successfully — saving workspace...\n")

save.image(
  file = "/working_groups/boylelab/shared/ZulyPava/formating4hpc/ATACseq_analysis/R_qc/workspace_final.RData"
)

cat("Workspace saved.\n")
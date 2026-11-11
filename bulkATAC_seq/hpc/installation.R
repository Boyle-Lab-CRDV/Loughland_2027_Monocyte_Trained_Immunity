# Ensure BiocManager is installed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Install Bioconductor packages
BiocManager::install(c(
  "ChIPseeker",
  "TxDb.Hsapiens.UCSC.hg38.knownGene",
  "org.Hs.eg.db",
  "clusterProfiler"
))

# Optional: install CRAN packages used for visualization/data handling
install.packages(c("ggplot2", "dplyr", "tibble"))

BiocManager::install("DiffBind")
BiocManager::install('EnhancedVolcano')

devtools::install_github("ge11232002/TFMPvalue")##only for r versions 4.5, 4.4 wouldnt work
install.packages( "https://cran.r-project.org/src/contrib/Archive/TFMPvalue/TFMPvalue_0.0.9.tar.gz", repos = NULL, type = "source" )

BiocManager::install("ATACseqQC")
BiocManager::install("sva")
install.packages("factoextra")

install.packages("/working_groups/boylelab/shared/ZulyPava/formating4hpc/ATACseq_analysis/R_qc/package/ChIPseeker_1.46.1.tar.gz", repos = NULL, type="source")
install.packages("/working_groups/boylelab/shared/ZulyPava/formating4hpc/ATACseq_analysis/R_qc/package/enrichplot_1.30.4.tar.gz", repos = NULL, type="source")
install.packages("/working_groups/boylelab/shared/ZulyPava/formating4hpc/ATACseq_analysis/R_qc/package/treeio_1.34.0.tar.gz", repos = NULL, type="source")
install.packages("/working_groups/boylelab/shared/ZulyPava/formating4hpc/ATACseq_analysis/R_qc/package/ggtree_4.0.4.tar.gz", repos = NULL, type="source")



# install.packages("/working_groups/boylelab/shared/ZulyPava/formating4hpc/ATACseq_analysis/R_qc/package/ATACseqQC_1.34.0.tar.gz", repos = NULL, type="source")

# install.packages("/working_groups/boylelab/shared/ZulyPava/formating4hpc/ATACseq_analysis/R_qc/package/motifStack_1.54.0.tar.gz", repos = NULL, type="source")
# install.packages("/working_groups/boylelab/shared/ZulyPava/formating4hpc/ATACseq_analysis/R_qc/package/TFBSTools_1.48.0.tar.gz", repos = NULL, type="source")

#install.packages("/working_groups/boylelab/shared/ZulyPava/formating4hpc/ATACseq_analysis/R_qc/package/Seqinfo_1.0.0.tar.gz", repos = NULL, type="source")


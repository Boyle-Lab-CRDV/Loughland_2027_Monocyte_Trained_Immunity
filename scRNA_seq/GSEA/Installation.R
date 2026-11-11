# Install BiocManager if missing
# if (!requireNamespace("BiocManager", quietly = TRUE)) {
#   install.packages("BiocManager")
# }

# Install Bioconductor packages
bioc_pkgs <- c("BloodGen3Module", "biomaRt")
for (pkg in bioc_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
}

# Install CRAN packages
cran_pkgs <- c( "viridis", "dplyr", "tidyr", "msigdbr", "data.table", "enrichR", "ggrepel")
for (pkg in cran_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

install.packages("/working_groups/boylelab/shared/ZulyPava/GSEA_damian_R4.3/renv/fgsea_1.36.2.tar.gz", repos = NULL, type = "source")

# install.packages("devtools")
# 
# install_github("alserglab/fgsea")

library(BloodGen3Module)
library (biomaRt)
library(fgsea)
library(dplyr)
library(tidyr)
library(msigdbr)
library(data.table)
library(enrichR)

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
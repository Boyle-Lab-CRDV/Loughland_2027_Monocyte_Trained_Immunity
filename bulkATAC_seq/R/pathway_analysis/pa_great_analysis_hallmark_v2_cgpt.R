library(rGREAT)
library(ATACseqQC)
library(DiffBind)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(dplyr)
library(ggplot2)
library(tidyverse)

set.seed(12345)

source("script/tables_monocytes_GREATv2.R") #creates and format the dataset
#list of databases
db_ls <- list(
  sabah_1_m_enrichgU = sabah_1_m_enrichgU,
  sabah_1_m_enrichgD = sabah_1_m_enrichgD,
  sabah_2_m_enrichgU = sabah_2_m_enrichgU,
  sabah_2_m_enrichgD = sabah_2_m_enrichgD,
  sabah_3_m_enrichgU = sabah_3_m_enrichgU,
  sabah_3_m_enrichgD = sabah_3_m_enrichgD,
  sabah_4_m_enrichgU = sabah_4_m_enrichgU,
  sabah_4_m_enrichgD = sabah_4_m_enrichgD,
  chmi_2_m_enrichgU  = chmi_2_m_enrichgU,
  chmi_2_m_enrichgD  = chmi_2_m_enrichgD
)

# pre-allocate result containers, same names as db_ls
db_gr_ls     <- vector("list", length(db_ls)); names(db_gr_ls)     <- names(db_ls)
great_res_ls <- vector("list", length(db_ls)); names(great_res_ls) <- names(db_ls)
tb_ls_res    <- vector("list", length(db_ls)); names(tb_ls_res)    <- names(db_ls)

#loop in around
for (db in names(db_ls)) {
  db_gr_ls[[db]]     <- GRanges(db_ls[[db]])
  great_res_ls[[db]] <- great(db_gr_ls[[db]], "MSigDB:H", "txdb:hg38")
  tb_ls_res[[db]]    <- getEnrichmentTable(great_res_ls[[db]])
}
##saving output
output_path <- "/working_groups/boylelab/shared/ZulyPava/ATACseq_R4.4/ATAC_seq/output/pathway_ana_enrich/hallmark/"
# make sure the directory exists (safe to run even if it already does)

if (!dir.exists(output_path)) {
  dir.create(output_path, recursive = TRUE)
}

for (nm in names(tb_ls_res)) {
  file_name <- paste0("pa_", nm, "_greatHa_290726.csv")
  write.csv(tb_ls_res[[nm]], file = file.path(output_path, file_name), row.names = FALSE)
}


###Reading saved files again to re-plot
output_path <- "/working_groups/boylelab/shared/ZulyPava/ATACseq_R4.4/ATAC_seq/output/pathway_ana_enrich/hallmark/"

files <- list.files(
  output_path,
  pattern = "^pa_.*_greatHa_290726\\.csv$",
  full.names = TRUE
)

tb_ls_res <- lapply(files, function(f) {
  tibble::as_tibble(read.csv(f, check.names = FALSE))
})

names(tb_ls_res) <- sub(
  "^pa_(.*)_greatHa_290726\\.csv$",
  "\\1",
  basename(files)
)


#### Plotting dotplot loopies ####
# 1. Define common scales across all plots
# Use the same colour limits for upregulated and downregulated pathways.

all_padj <- unlist(lapply(tb_ls_res, function(df) df$p_adjust))

all_neg_log_padj <- -log10(all_padj)

# Avoid Inf when adjusted p-values are zero
all_neg_log_padj <- all_neg_log_padj[is.finite(all_neg_log_padj)]

color_limits <- c(0, ceiling(max(all_neg_log_padj, na.rm = TRUE)))

# Common gene-count scale
all_gene_counts <- unlist(
  lapply(tb_ls_res, function(df) df$observed_gene_hits)
)

size_limits <- c(0, max(all_gene_counts, na.rm = TRUE))


# 2. Plot all datasets using one loop

for (nm in names(tb_ls_res)) {
  
  df <- tb_ls_res[[nm]]
  
  direction <- if (grepl("U$", nm)) {
    "up"
  } else if (grepl("D$", nm)) {
    "down"
  } else {
    warning(paste("Skipping", nm, "- cannot determine direction"))
    next
  }
  
  base_name <- sub("(U|D)$", "", nm)
  
  # Direction-specific colours
  if (direction == "up") {
    low_col  <- "#FCE8E8"
    high_col <- "#B2182B"
  } else {
    low_col  <- "lightblue1"
    high_col <- "dodgerblue4"
  }
  
  # Select top 10 pathways
  df_top <- df %>%
    arrange(desc(fold_enrichment)) %>%
    slice_head(n = 10) %>%
    mutate(
      description = factor(id, levels = rev(id)),
      neg_log_padj = -log10(p_adjust),
      neg_log_padj = ifelse(
        is.infinite(neg_log_padj),
        color_limits[2],
        neg_log_padj
      )
    )
  
  # Output directory
  out_dir <- file.path(output_path, direction)
  
  dir.create(
    out_dir,
    showWarnings = FALSE,
    recursive = TRUE
  )
  
  pdf_file <- file.path(
    out_dir,
    paste0(base_name, "_", direction, "_dotplot.pdf")
  )
  
  # Plot
  p <- ggplot(
    df_top,
    aes(x = fold_enrichment, y = description)
  ) +
    
    geom_point(
      aes(
        size = observed_gene_hits,
        color = neg_log_padj
      )
    ) +
    
    scale_color_gradient(
      low = low_col,
      high = high_col,
      limits = color_limits,
      breaks = pretty(color_limits, n = 5),
      name = expression(-log[10]~(p.adjust)),
      oob = scales::squish
    ) +
    
    scale_size_continuous(
      name = "Gene count",
      limits = size_limits,
      range = c(2, 8)
    ) +
    
    scale_y_discrete(
      labels = function(x) str_wrap(x, width = 35),
      expand = expansion(add = 0.6)
    ) +
    
    scale_x_continuous(
      breaks = seq(1.5, 4.0, by = 0.5),
      expand = expansion(mult = 0)
    ) +
    
    coord_cartesian(xlim = c(1, 4.0)) +
    
    labs(
      x = "Fold enrichment",
      y = NULL,
      title = paste0(base_name, " (", direction, ")")
    ) +
    
    theme_bw(base_size = 12) +
    
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      
      axis.text.y = element_text(
        size = 9,
        color = "black"
      ),
      
      axis.text.x = element_text(color = "black"),
      
      plot.title = element_text(
        face = "bold",
        size = 11
      ),
      
      legend.position = "right",
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 8),
      legend.key.size = unit(0.4, "cm"),
      
      plot.margin = margin(
        t = 5, r = 5, b = 5, l = 5
      )
    )
  
  # Dynamic plot height
  n_rows <- nrow(df_top)
  
  plot_height <- max(
    3,
    0.45 * n_rows + 1.2
  )
  
  # Save PDF
  pdf(
    pdf_file,
    width = 7,
    height = plot_height
  )
  
  print(p)
  dev.off()
}

#### Plotting Heatmap loopies####
#base_out <- "/working_groups/boylelab/shared/ZulyPava/formating4hpc/ATACseq_analysis/pathway_analysis"
for (nm in names(tb_ls_res)) {
  df <- tb_ls_res[[nm]]
  
  # direction from the U/D suffix in the name
  direction <- if (grepl("U$", nm)) "up" else if (grepl("D$", nm)) "down" else {
    warning(paste("Skipping", nm, "- can't tell if Up or Down")); next
  }
  
  fill_color <- if (direction == "up") "#B2182B" else "dodgerblue4"
  base_name  <- sub("(U|D)$", "", nm)  # strip trailing U/D, e.g. "sabah_1_m_enrichg"
  
  # sort by p_adjust ascending, keep top 10, and lock in the plotting order
  df_top <- df %>%
    arrange(p_adjust) %>%
    slice_head(n = 10) %>%
    mutate(description = factor(id, levels = rev(id)))
  
  out_dir  <- file.path(output_path, direction)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  pdf_file <- file.path(out_dir, paste0(base_name, "_", direction, "_heatmapHa.pdf"))
  
  pdf(pdf_file, width = 8, height = 10)
  print(
    ggplot(df_top, aes(x = 1, y = description, fill = fold_enrichment)) +
      geom_tile(width = 1, height = 1) +
      scale_fill_gradient(low = "white", high = fill_color, limits = c(0, NA)) +
      scale_y_discrete(labels = function(x) str_wrap(x, width = 30)) +
      scale_x_continuous(limits = c(0.5, 1.5), breaks = NULL) +
      coord_fixed(ratio = 1) +
      theme_bw() +
      theme(
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 9)
      )
  )
  dev.off()
}


####Plotting individuals####
sabah_1_paths_plot = sabah_1_paths%>%
  rename(`# Term Name` = "Term_name", "Binom Fold Enrichment" = "Fold_Enrichment")

pdf("/working_groups/boylelab/shared/ZulyPava/formating4hpc/ATACseq_analysis/pathway_analysis/up/sabah_1_m_up_heatmap.pdf", width = 8, height = 10)
ggplot(sabah_1_paths_plot, aes(y = Term_name, x = 1, fill = Fold_Enrichment)) +
  geom_tile(width = 0.5) +
  scale_fill_gradient(low = "white", high = "#B2182B", limits = c(0, NA))+
  theme_bw()
dev.off()

sabah_1_paths_plotd = sabah_1_paths_d%>%
  rename(`# Term Name` = "Term_name", "Binom Fold Enrichment" = "Fold_Enrichment")

pdf("/working_groups/boylelab/shared/ZulyPava/formating4hpc/ATACseq_analysis/pathway_analysis/down/sabah_1_m_down_heatmap.pdf", width = 8, height = 10)
ggplot(sabah_1_paths_plotd, aes(y = Term_name, x = 1, fill = Fold_Enrichment)) +
  geom_tile(width = 0.5) +
  scale_fill_gradient(low = "white", high = "dodgerblue4", limits = c(0, NA))+
  theme_bw()
dev.off()


##Sabah2

sabah_2_paths_plot = sabah_2_paths%>%
  rename(`# Term Name` = "Term_name", "Binom Fold Enrichment" = "Fold_Enrichment")

pdf("/working_groups/boylelab/shared/ZulyPava/formating4hpc/ATACseq_analysis/pathway_analysis/up/sabah_2_m_up_heatmap.pdf", width = 8, height = 8)
ggplot(sabah_2_paths_plot, aes(y = Term_name, x = 1, fill = Fold_Enrichment)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "#B2182B", limits = c(0, NA))+
  theme_bw()
dev.off()

sabah_2_paths_plotd = sabah_2_pathsD%>%
  rename(`# Term Name` = "Term_name", "Binom Fold Enrichment" = "Fold_Enrichment")

pdf("/working_groups/boylelab/shared/ZulyPava/formating4hpc/ATACseq_analysis/pathway_analysis/down/sabah_2_m_down_heatmap.pdf", width =10, height = 8)
ggplot(sabah_2_paths_plotd, aes(y = Term_name, x = 1, fill = Fold_Enrichment)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "dodgerblue4", limits = c(0, NA))+
  theme_bw()
dev.off()
 
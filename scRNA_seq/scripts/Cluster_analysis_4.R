##Author ZulyP
#Date: 02 10 2024


library(gmodels)
library(here)
library(tidyverse)

##Analysing clusters markers and differences in proportion among days.

##We need all the genes for this, sig and not sig.

psbDE <- read.csv(here(file="output_for_paper/edger/psbulk_allgenes_DayClusters.csv"))


colnames(psbDE)
#Renaming columns
colnames(psbDE)[2] ="celltype"
colnames(psbDE)[4] ="genes"
colnames(psbDE)[5] ="Log2FC"
colnames(psbDE)[9] ="P_adj_val"

head(psbDE)
dim(psbDE)
##removing first column
psbDE <- psbDE[,-1]
## checking clusters name
unique(psbDE$celltype)
##fixing names
psbDE$celltype <- ifelse(psbDE$celltype == "cMonocyte-IL", "cMonocyte-IL1B", psbDE$celltype)

psbDE$day <- ifelse(psbDE$day == "day8", "Day0 vs Day8", psbDE$day)
psbDE$day <- ifelse(psbDE$day == "day36", "Day0 vs Day36", psbDE$day)
psbDE$day <- ifelse(psbDE$day == "day836", "Day8 vs Day36", psbDE$day)
##creating day factors
psbDE$day <- factor(psbDE$day, levels = c("Day0 vs Day8", "Day0 vs Day36", "Day8 vs Day36" ))

##Creating celltype factors
celltype_fact<- c("cMonocyte-IL1B","cMonocyte-S100A8", "cMonocyte-COQ7", "cMonocyte-ISG15", "intMonocyte", "ncMonocyte-CD16", "DCs1", "DCs2", "pDCs" )
##Creating day factors
psbDE$celltype <- factor(psbDE$celltype, levels = celltype_fact)

#Creating dataset with direction, sig, and total not_deg 

psbDE1 <- psbDE%>%
  mutate(
    Dir =case_when(
      Log2FC>0 ~ "Up",
      Log2FC<0 ~ "Down",
      TRUE ~ NA_character_ 
    ),
    Sig = case_when(
      P_adj_val< 0.05 & abs(Log2FC) > 1 ~ "yes",
      TRUE ~ "no"
    ),
  cell_day = paste0(celltype,"_", day)
  )%>%
  filter(day != "Day8 vs Day36")

table(psbDE1$cell_day,psbDE1$Sig, useNA = "ifany")

##Checking that the number of sig genes is the same in both datasets used for plotting:
##read csv prefiltered
psbDESig <- read.csv(here(file="output_for_paper/edger/psbulk_sigGens__all_DayClusters.csv"))
#creating variable per samples
psbDESig$clustday <- paste0(psbDESig$cluster_name,"_",psbDESig$day)
#getting tables for comparison
sig<- as.data.frame(table(psbDESig$clustday))
all <- as.data.frame(table(psbDE$cell_day, psbDE$Sig))
#numbers matched!

##Creating our crosstable
##setting up factors
psbDE1$day <- factor(psbDE1$day, levels = c("Day0 vs Day8", "Day0 vs Day36"))
psbDE1$celltype <- factor(psbDE1$celltype, levels = celltype_fact)
psbDE1 


##Table Setup for Chi-square/Fisher’s Test:
celltype_list <- unique(psbDE1$celltype)

for (name in celltype_list){
  filtered_psbDE1 <- psbDE1%>%
  filter(celltype== name)

    # Run CrossTable for the current cell type
    output <- capture.output(
    gmodels::CrossTable(x=filtered_psbDE1$cell_day, y=filtered_psbDE1$Sig,
                       expected = TRUE,
                       prop.r = FALSE, prop.c = FALSE, prop.t = FALSE, prop.chisq = FALSE,
                       chisq = TRUE, fisher = TRUE,
                       format = "SPSS")
)
    # Create a valid file name by replacing non-alphanumeric characters in the cell type
    valid_ct <- gsub("[^[:alnum:]]", "_", name)
    
    # Construct the full file path using the 'here' package
    file_path <- here(paste0("output_for_paper/ChiS_NoDEG_btw_days_02102025_", valid_ct, ".txt"))
    
    # Save the output to the specified path
    writeLines(output, con = file_path)
    
    cat("Results for Celltype:", name, "saved to", file_path, "\n")  # Indicate completion for each cell type
}

##Individual code
psbDE1 %>%
  filter(celltype=="cMonocyte-COQ7")%>%
  {gmodels::CrossTable(x=.$cell_day, y=.$Sig,
                       expected = TRUE,
                       prop.r = FALSE, prop.c = FALSE, prop.t = FALSE, prop.chisq = FALSE,
                       chisq = TRUE, fisher = TRUE,
                       format = "SPSS")} 
#### cMonocyte-COQ7####
# Cell Contents
# |-------------------------|
#   |                   Count |
#   |         Expected Values |
#   |-------------------------|
#   
#   Total Observations in Table:  24453 
# 
# | .$Sig 
# .$cell_day |       no  |      yes  | Row Total | 
#   -----------------------------|-----------|-----------|-----------|
#   cMonocyte-COQ7_Day0 vs Day36 |    10233  |      266  |    10499  | 
#   | 10364.612  |  134.388  |           | 
#   -----------------------------|-----------|-----------|-----------|
#   cMonocyte-COQ7_Day0 vs Day8 |    13907  |       47  |    13954  | 
#   | 13775.388  |  178.612  |           | 
#   -----------------------------|-----------|-----------|-----------|
#   Column Total |    24140  |      313  |    24453  | 
#   -----------------------------|-----------|-----------|-----------|
#   
#   
#   Statistics for All Table Factors
# 
# 
# Pearson's Chi-squared test 
# ------------------------------------------------------------
# Chi^2 =  228.802     d.f. =  1     p =  1.087924e-51 
# 
# Pearson's Chi-squared test with Yates' continuity correction 
# ------------------------------------------------------------
# Chi^2 =  227.0669     d.f. =  1     p =  2.600289e-51 
# 
#  
# Fisher's Exact Test for Count Data
# ------------------------------------------------------------
#   Sample estimate odds ratio:  0.130034 
# 
# Alternative hypothesis: true odds ratio is not equal to 1
# p =  7.693989e-54 
# 95% confidence interval:  0.09314371 0.178044 
# 
# Alternative hypothesis: true odds ratio is less than 1
# p =  4.006508e-54 
# 95% confidence interval:  0 0.1699016 
# 
# Alternative hypothesis: true odds ratio is greater than 1
# p =  1 
# 95% confidence interval:  0.09828827 Inf 
# 
# 
# 
# Minimum expected frequency: 134.3879 





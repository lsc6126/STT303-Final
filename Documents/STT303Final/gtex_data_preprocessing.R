# Code for cleaning up and merging data
# GTEx data + subject ID (SUBJID) + AGE + specific tissue source (SMTSD)

#Preprocessing data
library(dplyr)
library(tibble)

# Transposing gene expression matrix
gene_expr <- gtex_data |>
  column_to_rownames("Name") |>
  select(-Description) |>
  t() |>
  as.data.frame()

# Adding sample IDs (SAMPID) as a column
gene_expr <- gene_expr |>
  rownames_to_column("SAMPID")

# Merging with Sample metadata to get tissue type
merged_data <-  merge(gene_expr,gtex_metadata[,c("SAMPID","SMTSD","SUBJID")],by= "SAMPID")

# Merging with subject phenotype data to get age ranges
merged_data <- merge(merged_data, gtex_phenotypes[,c("SUBJID","AGE")],by = "SUBJID")

# Reducing age bins to midpoint values (i.e. 30-39->35)
age_map <- c("20-29"=25, "30-39"=35, "40-49"=45, "50-59"=55, "60-69"=65, "70-79"=75)
merged_data$AGE_NUM <- age_map[merged_data$AGE]

#list of genes
gene_columns <- setdiff(colnames(merged_data), c("SAMPID", "SUBJID", "SMTSD", "AGE", "AGE_NUM"))
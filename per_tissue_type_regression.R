# This file is self contained compared to the others and only conducts
# our analysis on the tissue-type correlation with age through linear regression

# Libraries
library(ggplot2)
library(readxl)
library(dplyr)
library(tidyr)
library(broom)

# --- read gct data
read_gct <- function(file) {
  lines <- readLines(file)
  
  # First line is version
  version <- lines[1]
  
  # Second line: number of rows and columns
  dims <- strsplit(lines[2], "\t")[[1]]
  n_rows <- as.integer(dims[1])
  n_cols <- as.integer(dims[2])
  
  # Read the actual data
  df <- read.delim(file, skip = 2, header = TRUE, stringsAsFactors = FALSE)
  
  # Extract expression matrix
  expr <- as.matrix(df[, -(1:2)])  # Remove Name and Description
  rownames(expr) <- df$Name
  
  return(list(
    version = version,
    data = df,
    expression = expr
  ))
}

# sets the working directory (might not work on RStudio in which case change to ~/Downloads
# directory through the session working directory)
setwd("~/Downloads")

# installs the various data used in this analysis
gtex_data <- read_gct("GTEx_Analysis_2022-06-06_v10_RNASeQCv2.4.2_gene_reads.gct")
gtex_attr <- read_xlsx("GTEx_Analysis_v10_Annotations_SampleAttributesDD.xlsx")
gtex_metadata <- read.delim("GTEx_Analysis_v10_Annotations_SampleAttributesDS.txt", sep = "\t", header = TRUE)
gtex_age <- read.delim("GTEx_Analysis_v10_Annotations_SubjectPhenotypesDS.txt", sep = "\t", header = TRUE)

# makes the SUBJID more readable and easier to work with
gtex_age$SUBJID <- gsub("-", ".", gtex_age$SUBJID)

# Gets all tissue types
gtex_tissues <- unique(gtex_metadata$SMTS)

# For the tissue type passed, grabs ID's that share that tissue type, then creates
# and assigns a variable name (such as Blood_IDs instead of Blood) to the IDs
tissue_IDs<- function(tissue) {
  IDs <- gsub("-", ".", gtex_metadata$SAMPID[gtex_metadata$SMTS==tissue])
  tissue <- gsub(" ", "_", tissue)
  var_name <- paste(tissue, "IDs", sep="_")
  assign(var_name, IDs, envir=.GlobalEnv)
  return(var_name)
}

# Makes a list of all tissue types and stores the tissues IDs in them
tissue_list <- NULL
for (i in gtex_tissues) {
  tissue_type <- tissue_IDs(i)
  tissue_list <- c(tissue_list, tissue_type)
  iterable_tissue <- get(tissue_type)
  assign(tissue_type, iterable_tissue, envir=.GlobalEnv)
}

# All tissue types ID's
all_tissues <- c()
for (i in tissue_list) {
  all_tissues <- c(all_tissues, get(i))
}

# All tissues thrown into a daata frame
tissue_df <- data.frame()
for (i in tissue_list) {
  tissue_type <- gsub("_IDs", "", i)
  tissue_type <- gsub("_", " ", tissue_type)
  tissue_df <- rbind(tissue_df, data.frame(tissue=tissue_type, sampleID = get(i)))
}
# selects the part of the ID before the third period then merge with ages on that ID
tissue_df <- tissue_df |> mutate(SUBJID = sub("^([^.]+\\.[^.]+)\\..*", "\\1", sampleID))
ages <- gtex_age |> select(SUBJID, AGE)
age_tissues <- merge(tissue_df, ages, by = "SUBJID")

# selects all genes that are in our ID list that has corresponding ages and selects
# genes with more than 20 expressions
age_sample_list <- as.list(age_tissues$sampleID)
age_expr <- gtex_data$expression[,colnames(gtex_data$expression) %in% age_sample_list]
filtered_age_expr <- age_expr[rowSums(age_expr > 20) >= 20, ]

# assigns probabilities to all genes that have ages attached then grabs the top 500 of them,
# then it grabs the top genes and transposes them.
gene_vars_age <- apply(filtered_age_expr, 1, var)
top_genes_age <- names(sort(gene_vars_age, decreasing = TRUE)[1:500])
filtered_expr_age <- filtered_age_expr[top_genes_age, ]
expr_age <- filtered_expr_age
expr_t_age <- t(expr_age)

# takes the transposed genes and turns it into a data frame and conducts a
# probability comparison on the data.
expr_t_age_df <- as.data.frame(expr_t_age)
expr_t_age_df$sampleID <- rownames(expr_t_age_df)
expr_age_pca_result <- prcomp(expr_t_age, center = TRUE, scale. = TRUE)

# selects only the SUBJIDs and ages from the data frame of ages and tissues
# then preforms a left join to merge the dataframes on their ID's
age_merger <- age_tissues |> select(SUBJID, AGE)
age_expr <- left_join(expr_t_age_df, age_tissues, by="sampleID")

# AGE column then gets everything after the dash removed and turned numeric
# ("20-29" becomes 20)
age_expr$AGE <- gsub("-.*", "", age_expr$AGE)
age_expr$AGE <- as.numeric(age_expr$AGE)

# Removes SUBJID column and modifies the data frame to allow Genes and
# Expressions to be represented by columns Gene and Expressions respectively
age_expr_long_prep <- age_expr |> select(-SUBJID)
age_expr_long <- pivot_longer(age_expr_long_prep, cols = -c(sampleID, AGE, tissue), 
                              names_to = "Gene", values_to = "Expression")

# Runs a linear regression on every tissue type
tissue_age_model <- age_expr_long |> 
  group_by(tissue) |> 
  do(tidy(lm(Expression ~ AGE, data = .)))

# grabs the R² values from every tissue type
tissue_r2 <- age_expr_long |> 
  group_by(tissue) |> 
  do(R2 = summary(lm(Expression ~ AGE, data = .))$r.squared)

# unlists the R² values for every tissue type
tissue_r2_clean <- tissue_r2 |> 
  mutate(R2 = unlist(R2))

# Plots the R² values of every tissue type's regression
ggplot(tissue_r2_clean, aes(x = R2, y = tissue)) +
  geom_point(size = 3, color = "blue") +
  theme_minimal() +
  labs(title = "Tissue R² Values", x = "R²", y = "Tissue")

# Significant tissue types based off of the linear regression and drops non age
# related rows
sig_tissue <- tissue_age_model |> filter(p.value < 0.05)
sig_tissue_df <- data.frame(sig_tissue)
sig_tissue_no_intercept <- sig_tissue_df |> filter(term != "(Intercept)")
sig_tissue_no_intercept
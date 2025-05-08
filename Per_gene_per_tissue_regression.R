# Regression: Per Gene, Per Tissue
# Regression run on modified data frame without any form of dimension reduction
# Genes ~ Age
# Gives a list of genes with strongest R^2 and P-value per tissue
# Note: Extremely small R^2 values with large P-values

# Necessary Libraries
library(ggplot2)
library(dplyr)

# Initialize results list
result_list <- list()

# Loop over each tissue
for (tissue in unique(merged_data$SMTSD)) {
  tissue_data <- merged_data %>% filter(SMTSD == tissue)
  
  # Skip tissues with too few samples
  if (nrow(tissue_data) < 10) next
  
  # Loop over each gene and run regression
  for (gene in gene_columns) {
    expr <- tissue_data[[gene]]
    age <- tissue_data$AGE_NUM
    
    if (all(is.na(expr))) next  # Skip if all NA
    
    model <- try(lm(expr ~ age), silent = TRUE)
    
    if (class(model) != "try-error") {
      coef_summary <- summary(model)$coefficients
      if ("age" %in% rownames(coef_summary)) {
        estimate <- coef_summary["age", "Estimate"]
        pval <- coef_summary["age", "Pr(>|t|)"]
        r_squared <- summary(model)$r.squared
        
        result_list[[length(result_list) + 1]] <- data.frame(
          Tissue = tissue,
          Gene = gene,
          Estimate = estimate,
          P_Value = pval,
          R_squared = r_squared,
          stringsAsFactors = FALSE
        )
      }
    }
  }
}

# Combine all results into one data frame
all_results <- do.call(rbind, result_list)

# Optional: rank genes by p-value within each tissue
ranked_results <- all_results %>%
  group_by(Tissue) %>%
  arrange(P_Value, .by_group = TRUE)

r2_ranked_results <- all_results %>%
  group_by(Tissue) %>%
  arrange(R_squared, .by_group = TRUE)

# gene expression vs age
ggplot(all_results, aes(x = R_squared, y = -log10(P_Value), color = Tissue)) +
  geom_point(alpha = 0.6) +
  labs(title = "Gene Expression vs Age by Tissue",
       x = "R² (Strength of Association)",
       y = "-log10(P-value) (Significance)",
       color = "Tissue") +
  theme_minimal()

# top genes by tissue
top_genes <- all_results %>%
  group_by(Tissue) %>%
  slice_min(order_by = P_Value, n = 5) %>%
  ungroup()

# Calling results
top_genes
r2_ranked_results
ranked_results


# For conducting regression analysis on dimension reduced expression data
# Grouping samples based on tissue of origin, conducting PCA, then regression
# KNN to attempt to predict subject age

# Necessary Libraries
library(caret)
library(FNN)
library(dplyr)
library(ggplot2)

# Filter and group data by tissue of origin
filtered_data <- merged_data %>%
  group_by(SMTSD) %>%
  group_split(.keep = TRUE)
# Create a named list where names match tissue types
tissue_data_list <- setNames(filtered_data, map_chr(filtered_data, ~ unique(.x$SMTSD)))

# Function to run PCA, regression, and KNN
# rejects tissues with too few samples
analyze_tissue_expression <- function(tissue_data_list, n_pc = 10, k = 5,
                                      meta_cols = 2, tail_cols = 3,
                                      min_samples = 15) {
  
  # Extract expression matrix
  get_expr_matrix <- function(df) {
    expr <- df[, (meta_cols + 1):(ncol(df) - tail_cols)]
    expr <- expr[, apply(expr, 2, sd) != 0]  # remove constant columns
    expr <- expr[, colMeans(expr != 0) > 0.1]  # keep features expressed in >10% of samples
    as.matrix(expr)
  }
  
  # Run PCA
  run_pca <- function(expr_matrix) {
    scaled <- scale(expr_matrix)
    pca_result <- prcomp(scaled, center = FALSE, scale. = FALSE)
    available_pc <- ncol(pca_result$x)
    if (available_pc < 1) return(NULL)
    n_extract <- min(n_pc, available_pc)
    pca_result$x[, 1:n_extract, drop = FALSE]
  }
  
  # Results containers
  summary_list <- list()
  predictions_list <- list()
  
  # Loop over tissues
  for (tissue_name in names(tissue_data_list)) {
    df <- tissue_data_list[[tissue_name]]
    
    # Skip tissues with missing or uniform AGE_NUM or too few samples
    if (!"AGE_NUM" %in% names(df) || length(unique(df$AGE_NUM)) < 2 || nrow(df) < min_samples) {
      summary_list[[tissue_name]] <- data.frame(
        Tissue = tissue_name, LinearRegression_R2 = NA, LinearRegression_pval = NA, KNN_MSE = NA
      )
      next
    }
    
    # Extract expression matrix
    expr <- tryCatch(get_expr_matrix(df), error = function(e) NULL)
    if (is.null(expr) || ncol(expr) < 1) {
      summary_list[[tissue_name]] <- data.frame(
        Tissue = tissue_name, LinearRegression_R2 = NA, LinearRegression_pval = NA, KNN_MSE = NA
      )
      next
    }
    
    # Run PCA
    pcs <- run_pca(expr)
    if (is.null(pcs) || ncol(pcs) < 1) {
      summary_list[[tissue_name]] <- data.frame(
        Tissue = tissue_name, LinearRegression_R2 = NA, LinearRegression_pval = NA, KNN_MSE = NA
      )
      next
    }
    
    # Keep AGE_NUM separate
    age_vector <- df$AGE_NUM
    available_pc <- ncol(pcs)
    
    # Linear regression
    pcs_df <- as.data.frame(pcs)
    pcs_df$AGE_NUM <- age_vector
    lm_model <- lm(AGE_NUM ~ ., data = pcs_df)
    lm_summary <- summary(lm_model)
    r2 <- lm_summary$r.squared
    
    # Compute p-value from F-statistic
    f_stat <- lm_summary$fstatistic
    p_val <- pf(f_stat[1], f_stat[2], f_stat[3], lower.tail = FALSE)
    
    # KNN regression with train/test split
    set.seed(123)
    train_index <- createDataPartition(age_vector, p = 0.8, list = FALSE)
    
    train_pcs <- pcs[train_index, , drop = FALSE]
    test_pcs <- pcs[-train_index, , drop = FALSE]
    train_age <- age_vector[train_index]
    test_age <- age_vector[-train_index]
    
    n_train <- nrow(train_pcs)
    k_effective <- min(k, n_train - 1)
    #Collect KNN results
    knn_result <- knn.reg(train = train_pcs, test = test_pcs,
                          y = train_age, k = k_effective)
    
    predicted <- knn_result$pred
    actual <- test_age
    mse <- mean((predicted - actual)^2)
    
    # Save results
    summary_list[[tissue_name]] <- data.frame(
      Tissue = tissue_name,
      LinearRegression_R2 = round(r2, 4),
      LinearRegression_pval = signif(p_val, 4),
      KNN_MSE = round(mse, 2)
    )
    
    predictions_list[[tissue_name]] <- data.frame(
      Actual_Age = actual,
      Predicted_Age = predicted
    )
  }
  
  # Combine and return results
  return(list(
    summary = do.call(rbind, summary_list),
    predictions = predictions_list
  ))
}

# Running the function on GTEx data
results <- analyze_tissue_expression(tissue_data_list, n_pc = 10, k = 5)

# Viewing results
results$summary
view(results$summary)
# Combine all tissue predictions into one data frame with tissue labels
predictions_df <- bind_rows(
  lapply(names(results$predictions), function(tissue) {
    df <- results$predictions[[tissue]]
    df$Tissue <- tissue
    return(df)
  })
)

# Plot: Actual vs Predicted Age, faceting by all tissues
ggplot(predictions_df, aes(x = Actual_Age, y = Predicted_Age)) +
  geom_point(alpha = 0.6, size = 1.5) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  facet_wrap(~ Tissue, scales = "free") +
  theme_minimal() +
  labs(
    title = "Predicted vs Actual Age by Tissue",
    x = "Actual Age",
    y = "Predicted Age"
  )

# Plot of actual vs predicted for a particular tissue
tissue_to_plot <- "Muscle - Skeletal"  # change to any valid tissue name
# Filter the data for specified tissue
df_tissue <- predictions_df %>%
  filter(Tissue == tissue_to_plot)
# Generating plot
ggplot(df_tissue, aes(x = Actual_Age, y = Predicted_Age)) +
  geom_jitter(width = 0.5, height = 0.5, alpha = 0.7, size = 2, color = "steelblue") +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  theme_minimal() +
  labs(
    title = paste("Jittered Actual vs Predicted Age for", tissue_to_plot),
    x = "Actual Age",
    y = "Predicted Age"
  )


# Filter out NAs
# NAs due to either small number of samples or small number of PCs, or lack of variabilty in subject age
r2_df <- results$summary %>%
  filter(!is.na(LinearRegression_R2))

# Plot: R^2 values by tissue
ggplot(r2_df, aes(x = Tissue, y = LinearRegression_R2)) +
  geom_jitter(width = 0.2, height = 0, size = 2, alpha = 0.7, color = "steelblue") +
  theme_minimal() +
  labs(
    title = "R² Values by Tissue",
    x = "Tissue",
    y = expression(R^2)
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Get the row with the highest R² (excluding NAs)
best_r2 <- results$summary %>%
  filter(!is.na(LinearRegression_R2)) %>%
  arrange(desc(LinearRegression_R2)) %>%
  slice(1)

# View it
best_r2

# Visualizing the distribution of R^2 values
ggplot(results$summary, aes(x = LinearRegression_R2)) +
  geom_histogram(binwidth = 0.05, fill = "steelblue", color = "black") +
  theme_minimal() +
  labs(
    title = "Distribution of R² Values Across Tissues",
    x = expression(R^2),
    y = "Number of Tissues"
  )

# Summary statistics for R^2 and KNN MSE values
summary(results$summary$LinearRegression_R2)
summary(results$summary$KNN_MSE)


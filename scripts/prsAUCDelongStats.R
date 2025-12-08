#!/usr/bin/env Rscript

# Load required libraries
library(pROC)
library(boot)
library(dplyr)
library(tidyr)
library(argparse)
# ============================================================================
# HELPER: Safely build formulas with special characters
# ============================================================================
build_safe_formula <- function(outcome, predictors) {
  # Wrap predictor names with special characters in backticks
  safe_predictors <- sapply(predictors, function(x) {
    if (grepl("[^a-zA-Z0-9_\\.]", x)) {
      return(paste0("`", x, "`"))
    }
    return(x)
  })
  
  formula_str <- paste(outcome, "~", paste(safe_predictors, collapse = " + "))
  return(as.formula(formula_str))
}


# ============================================================================
# FUNCTION: Calculate Tjur's R²
# ============================================================================
calculate_tjur_r2 <- function(y_true, y_pred) {
  #
  #Tjur's R² = mean(predicted prob | Y=1) - mean(predicted prob | Y=0)
  
  #Args:
  #  y_true: Binary outcome (0/1 or 1/2)
  #  y_pred: Predicted probabilities
  #
  # Convert to 0/1 if needed
  if (all(y_true %in% c(1, 2))) {
    y_true <- y_true - 1
  }
  
  mean_cases <- mean(y_pred[y_true == 1])
  mean_controls <- mean(y_pred[y_true == 0])
  
  tjur_r2 <- mean_cases - mean_controls
  return(tjur_r2)
}

# ============================================================================
# FUNCTION: Bootstrap Standard Error for Tjur's R²
# ============================================================================
tjur_bootstrap <- function(data, indices, formula) {
  # """
  # Bootstrap function for Tjur's R² standard error
  # """
  d <- data[indices, ]
  model <- glm(formula, data = d, family = binomial(link = "logit"))
  y_pred <- predict(model, type = "response")
  tjur <- calculate_tjur_r2(d$PHENOTYPE, y_pred)
  return(tjur)
}

calculate_tjur_se <- function(data, formula, n_boot = 1000) {
  # """
  # Calculate Tjur's R² with bootstrap standard error
  # """
  # Point estimate
  model <- glm(formula, data = data, family = binomial(link = "logit"))
  y_pred <- predict(model, type = "response")
  tjur_r2 <- calculate_tjur_r2(data$PHENOTYPE, y_pred)
  
  # Bootstrap for SE
  set.seed(42)
  boot_results <- boot(data = data, 
                       statistic = tjur_bootstrap, 
                       R = n_boot,
                       formula = formula)
  
  tjur_se <- sd(boot_results$t)
  
  # 95% CI
  ci <- boot.ci(boot_results, type = "perc")
  ci_lower <- ci$percent[4]
  ci_upper <- ci$percent[5]
  
  return(list(
    tjur_r2 = tjur_r2,
    se = tjur_se,
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    z_score = tjur_r2 / tjur_se,
    p_value = 2 * pnorm(-abs(tjur_r2 / tjur_se))
  ))
}

# ============================================================================
# FUNCTION: Calculate AUROC with DeLong Test
# ============================================================================
calculate_auroc_with_delong <- function(data, formula, comparison_roc = NULL) {
  # """
  # Calculate AUROC with optional DeLong test comparison
  # 
  # Args:
  #   data: Data frame with outcomes and predictors
  #   formula: GLM formula
  #   comparison_roc: Optional ROC object to compare against
  # 
  # Returns:
  #   List with AUROC, CI, and DeLong test results
  # """
  # Fit model
  model <- glm(formula, data = data, family = binomial(link = "logit"))
  y_pred <- predict(model, type = "response")
  
  # Convert phenotype to 0/1 if needed
  y_true <- data$PHENOTYPE
  if (all(y_true %in% c(1, 2))) {
    y_true <- y_true - 1
  }
  
  # Calculate ROC
  roc_obj <- roc(y_true, y_pred, quiet = TRUE)
  auc_value <- as.numeric(auc(roc_obj))
  ci_obj <- ci.auc(roc_obj)
  
  result <- list(
    auc = auc_value,
    ci_lower = ci_obj[1],
    ci_upper = ci_obj[3],
    roc_object = roc_obj,
    predictions = y_pred
  )
  
  # DeLong test if comparison provided
  if (!is.null(comparison_roc)) {
    delong_test <- roc.test(roc_obj, comparison_roc, method = "delong")
    result$delong_p_value <- delong_test$p.value
    result$delong_statistic <- delong_test$statistic
  }
  
  return(result)
}

# ============================================================================
# FUNCTION: Pairwise DeLong Comparisons
# ============================================================================
perform_pairwise_comparisons <- function(data, comparisons_list) {
  # """
  # Perform pairwise DeLong tests between specified models
  # 
  # Args:
  #   data: Data frame with PRS scores and phenotype
  #   comparisons_list: List of comparisons, e.g., 
  #                    list(c("model1", "model2"), c("model1", "model3"))
  # 
  # Returns:
  #   Data frame with pairwise comparison results
  # """
  
  pairwise_results <- list()
  
  for (i in seq_along(comparisons_list)) {
    comparison <- comparisons_list[[i]]
    model1_name <- comparison[1]
    model2_name <- comparison[2]
    
    cat(sprintf("\n=== Comparing %s vs %s ===\n", model1_name, model2_name))
    
    # Fit models
    formula1 <- build_safe_formula("PHENOTYPE", model1_name)
    formula2 <- build_safe_formula("PHENOTYPE", model2_name)
    
    model1 <- glm(formula1, data = data, family = binomial(link = "logit"))
    model2 <- glm(formula2, data = data, family = binomial(link = "logit"))
    
    pred1 <- predict(model1, type = "response")
    pred2 <- predict(model2, type = "response")
    
    # Convert phenotype to 0/1 if needed
    y_true <- data$PHENOTYPE
    if (all(y_true %in% c(1, 2))) {
      y_true <- y_true - 1
    }
    
    # Calculate ROCs
    roc1 <- roc(y_true, pred1, quiet = TRUE)
    roc2 <- roc(y_true, pred2, quiet = TRUE)
    
    auc1 <- as.numeric(auc(roc1))
    auc2 <- as.numeric(auc(roc2))
    
    # DeLong test
    delong_test <- roc.test(roc1, roc2, method = "delong")
    
    # Tjur's R² difference
    tjur1 <- calculate_tjur_r2(data$PHENOTYPE, pred1)
    tjur2 <- calculate_tjur_r2(data$PHENOTYPE, pred2)
    tjur_diff <- tjur2 - tjur1
    
    cat(sprintf("AUC %s: %.4f\n", model1_name, auc1))
    cat(sprintf("AUC %s: %.4f\n", model2_name, auc2))
    cat(sprintf("AUC Difference: %.4f\n", auc2 - auc1))
    cat(sprintf("DeLong p-value: %.4e\n", delong_test$p.value))
    cat(sprintf("Tjur R² difference: %.4f\n", tjur_diff))
    
    pairwise_results[[i]] <- data.frame(
      model1 = model1_name,
      model2 = model2_name,
      auc_model1 = auc1,
      auc_model2 = auc2,
      auc_difference = auc2 - auc1,
      delong_z = delong_test$statistic,
      delong_p_value = delong_test$p.value,
      tjur_r2_model1 = tjur1,
      tjur_r2_model2 = tjur2,
      tjur_r2_difference = tjur_diff
    )
  }
  
  return(bind_rows(pairwise_results))
}


# ============================================================================
# MAIN FUNCTION: Compare Multiple PRS Models
# ============================================================================
compare_prs_models <- function(data,
                               prs_columns, 
                               covariate_columns = NULL,
                               phenotype_col = "PHENOTYPE",
                               output_path = NULL,
                               pairwise_comparisons = NULL,
                               pairwise_output = NULL) {
  
  # Ensure phenotype column is named PHENOTYPE for functions
  if (phenotype_col != "PHENOTYPE") {
    data$PHENOTYPE <- data[[phenotype_col]]
  }
  
  # Convert PHENOTYPE from 1/2 to 0/1 if needed
  unique_vals <- unique(data$PHENOTYPE)
  if (setequal(unique_vals, c(1, 2))) {
    data$PHENOTYPE <- data$PHENOTYPE - 1
    cat("Converted PHENOTYPE from 1/2 coding to 0/1 coding (1=control, 2=case → 0=control, 1=case)\n")
  } else if (setequal(unique_vals, c(0, 1))) {
    cat("PHENOTYPE already in 0/1 coding\n")
  } else {
    warning(paste("Unexpected PHENOTYPE values:", paste(unique_vals, collapse=", ")))
  }
  
  results_list <- list()
  roc_objects <- list()  # ADD THIS - was missing
  
  # Store baseline (covariate-only) ROC for comparisons
  baseline_roc <- NULL
  
  # ========================================================================
  # 1. COVARIATE-ONLY MODEL (if covariates provided)
  # ========================================================================
  if (!is.null(covariate_columns) && length(covariate_columns) > 0) {
    cat("\n=== Analyzing Covariate-Only Model ===\n")
    
    formula_cov <- as.formula(paste("PHENOTYPE ~", 
                                    paste(covariate_columns, collapse = " + ")))
    
    # AUROC
    auroc_cov <- calculate_auroc_with_delong(data, formula_cov)
    baseline_roc <- auroc_cov$roc_object
    roc_objects[["covariate"]] <- auroc_cov$roc_object  # Store ROC
    
    # Tjur's R²
    tjur_cov <- calculate_tjur_se(data, formula_cov)
    
    results_list[[1]] <- data.frame(
      model = "covariates_only",
      model_type = "covariate",
      auc = auroc_cov$auc,
      auc_ci_lower = auroc_cov$ci_lower,
      auc_ci_upper = auroc_cov$ci_upper,
      auc_p_value = NA,
      tjur_r2 = tjur_cov$tjur_r2,
      tjur_se = tjur_cov$se,
      tjur_ci_lower = tjur_cov$ci_lower,
      tjur_ci_upper = tjur_cov$ci_upper,
      tjur_z_score = tjur_cov$z_score,
      tjur_p_value = tjur_cov$p_value
    )
    
    cat(sprintf("AUC: %.4f [%.4f - %.4f]\n", 
                auroc_cov$auc, auroc_cov$ci_lower, auroc_cov$ci_upper))
    cat(sprintf("Tjur R²: %.4f (SE: %.4f)\n", 
                tjur_cov$tjur_r2, tjur_cov$se))
  }
  
  # ========================================================================
  # 2. PRS-ONLY MODELS
  # ========================================================================
  cat("\n=== Analyzing PRS-Only Models ===\n")
  
  for (prs_col in prs_columns) {
    cat(sprintf("\n--- %s ---\n", prs_col))
    
    formula_prs <- build_safe_formula("PHENOTYPE", prs_col)
    
    # formula_prs <- as.formula(paste("PHENOTYPE ~", prs_col))
    
    # AUROC with DeLong test vs baseline
    auroc_prs <- calculate_auroc_with_delong(data, formula_prs, baseline_roc)
    roc_objects[[prs_col]] <- auroc_prs$roc_object  # Store ROC
    
    # Tjur's R²
    tjur_prs <- calculate_tjur_se(data, formula_prs)
    
    results_list[[length(results_list) + 1]] <- data.frame(
      model = prs_col,
      model_type = "prs_only",
      auc = auroc_prs$auc,
      auc_ci_lower = auroc_prs$ci_lower,
      auc_ci_upper = auroc_prs$ci_upper,
      auc_p_value = ifelse(is.null(auroc_prs$delong_p_value), 
                           NA, auroc_prs$delong_p_value),
      tjur_r2 = tjur_prs$tjur_r2,
      tjur_se = tjur_prs$se,
      tjur_ci_lower = tjur_prs$ci_lower,
      tjur_ci_upper = tjur_prs$ci_upper,
      tjur_z_score = tjur_prs$z_score,
      tjur_p_value = tjur_prs$p_value
    )
    
    cat(sprintf("AUC: %.4f [%.4f - %.4f]", 
                auroc_prs$auc, auroc_prs$ci_lower, auroc_prs$ci_upper))
    if (!is.null(auroc_prs$delong_p_value)) {
      cat(sprintf(" (p vs baseline: %.4e)", auroc_prs$delong_p_value))
    }
    cat("\n")
    cat(sprintf("Tjur R²: %.4f (SE: %.4f, p: %.4e)\n", 
                tjur_prs$tjur_r2, tjur_prs$se, tjur_prs$p_value))
  }
  
  # ========================================================================
  # 3. COMBINED MODELS (PRS + Covariates)
  # ========================================================================
  if (!is.null(covariate_columns) && length(covariate_columns) > 0) {
    cat("\n=== Analyzing Combined Models (PRS + Covariates) ===\n")
    
    for (prs_col in prs_columns) {
      cat(sprintf("\n--- %s + covariates ---\n", prs_col))
      
      all_predictors <- c(prs_col, covariate_columns)
      formula_combined <- build_safe_formula("PHENOTYPE", all_predictors)
      
      # AUROC with DeLong test vs baseline
      auroc_combined <- calculate_auroc_with_delong(data, formula_combined, baseline_roc)
      roc_objects[[paste0(prs_col, "_plus_cov")]] <- auroc_combined$roc_object  # Store ROC
      
      # Tjur's R²
      tjur_combined <- calculate_tjur_se(data, formula_combined)
      
      results_list[[length(results_list) + 1]] <- data.frame(
        model = paste0(prs_col, "_plus_covariates"),
        model_type = "prs_plus_covariates",
        auc = auroc_combined$auc,
        auc_ci_lower = auroc_combined$ci_lower,
        auc_ci_upper = auroc_combined$ci_upper,
        auc_p_value = ifelse(is.null(auroc_combined$delong_p_value), 
                             NA, auroc_combined$delong_p_value),
        tjur_r2 = tjur_combined$tjur_r2,
        tjur_se = tjur_combined$se,
        tjur_ci_lower = tjur_combined$ci_lower,
        tjur_ci_upper = tjur_combined$ci_upper,
        tjur_z_score = tjur_combined$z_score,
        tjur_p_value = tjur_combined$p_value
      )
      
      cat(sprintf("AUC: %.4f [%.4f - %.4f]", 
                  auroc_combined$auc, auroc_combined$ci_lower, auroc_combined$ci_upper))
      if (!is.null(auroc_combined$delong_p_value)) {
        cat(sprintf(" (p vs baseline: %.4e)", auroc_combined$delong_p_value))
      }
      cat("\n")
      cat(sprintf("Tjur R²: %.4f (SE: %.4f, p: %.4e)\n", 
                  tjur_combined$tjur_r2, tjur_combined$se, tjur_combined$p_value))
    }
  }
  
  # Combine results
  results_df <- bind_rows(results_list)
  
  # ========================================================================
  # 4. PAIRWISE COMPARISONS
  # ========================================================================
  pairwise_df <- NULL
  
  if (!is.null(pairwise_comparisons)) {
    cat("\n\n========================================\n")
    cat("PAIRWISE COMPARISONS\n")
    cat("========================================\n")
    
    pairwise_df <- perform_pairwise_comparisons(data, pairwise_comparisons)
    
    if (!is.null(pairwise_output)) {
      write.csv(pairwise_df, pairwise_output, row.names = FALSE)
      cat(sprintf("\nPairwise comparisons saved to: %s\n", pairwise_output))
    }
  }
  
  # Save main results
  if (!is.null(output_path)) {
    write.csv(results_df, output_path, row.names = FALSE)
    cat(sprintf("\nMain results saved to: %s\n", output_path))
  }
  
  return(list(
    main_results = results_df,
    pairwise_results = pairwise_df,
    roc_objects = roc_objects
  ))
}  


# ============================================================================
# EXAMPLE USAGE
# ============================================================================
# if (interactive() || !exists("data")) {
#   # Example with simulated data
#   set.seed(42)
#   n <- 1000
#   
#   data <- data.frame(
#     IID = 1:n,
#     PHENOTYPE = rbinom(n, 1, 0.3) + 1,  # 1=control, 2=case
#     scaled_prs_main = rnorm(n),
#     scaled_prs_epi = rnorm(n),
#     scaled_prs_cardio = rnorm(n),
#     age = rnorm(n, 50, 10),
#     sex = rbinom(n, 1, 0.5),
#     PC1 = rnorm(n),
#     PC2 = rnorm(n)
#   )
#   
#   # Run comparison
#   results <- compare_prs_models(
#     data = data,
#     prs_columns = c("scaled_prs_main", "scaled_prs_epi", "scaled_prs_cardio"),
#     covariate_columns = c("age", "sex", "PC1", "PC2"),
#     phenotype_col = "PHENOTYPE",
#     output_path = "prs_model_comparison.csv"
#   )
#   
#   print(results)
# }

# ============================================================================
# RUN SCRIPT WITHE INPUT
# ============================================================================
################################ GLOBAL VARIABLES  ########################
# 
parser <- ArgumentParser()
parser$add_argument("--pheno_data", required = TRUE)

#
args <- parser$parse_args()

pheno_data <- args$pheno_data

scores_path = paste0(pheno_data,'/scores')

# Define filepaths
data_file = paste0(scores_path,'/combinedPRSGroups.csv')
cov_file = paste0(scores_path,'/covariate.12.mixed.prs.csv')

#output files
outputAUCFile = paste0(scores_path,"/prs_delong_auc_statistics.csv")
outputPairwiseFile = paste0(scores_path,"/prs_pairwise_comparisons.csv")

# Load data
data <- read.csv(data_file,check.names = FALSE)
cov_data <- read.csv(cov_file,check.names = FALSE)

# Merge dataframes
data <- left_join(data, 
                  cov_data %>% select(IID, scaled_prs_covariate = scaled_prs), 
                  by = "IID")

# Define PRS columns
prs_cols <- c("scaled_prs_main", "scaled_prs_epi", "scaled_prs_cardio",
              "scaled_prs_epi+main","scaled_prs_all")

# Define covariates
covariate_cols <- c("scaled_prs_covariate")

# Define pairwise comparisons
comparisons <- list(
  c("scaled_prs_main", "scaled_prs_epi"),           # main vs epi
  c("scaled_prs_main", "scaled_prs_epi+main"),      # main vs epi+main
  c("scaled_prs_epi", "scaled_prs_epi+main"),       # epi vs epi+main
  c("scaled_prs_main", "scaled_prs_cardio"),        # main vs cardio
  c("scaled_prs_epi", "scaled_prs_cardio"),         # epi vs cardio
  c("scaled_prs_main", "scaled_prs_all"),           # main vs all
  c("scaled_prs_epi", "scaled_prs_all"),            # epi vs all
  c("scaled_prs_cardio", "scaled_prs_all"),         # cardio vs all
  c("scaled_prs_epi+main", "scaled_prs_all")        # epi+main vs all
)

# Run analysis
results <- compare_prs_models(
  data = data,
  prs_columns = prs_cols,
  covariate_columns = covariate_cols,
  phenotype_col = "PHENOTYPE",
  output_path = outputAUCFile,
  pairwise_comparisons = comparisons,
  pairwise_output = outputPairwiseFile
  
)

#!/usr/bin/env Rscript

# Keri Multerer
# Nagelkerke R² with incremental R² for each statistically distinct PRS model
# vs covariate-only baseline.
#
# Input:  combinedPRS.filtered.csv  (one row per individual, scaled_prs_* columns)
#         covariate.12.mixed.prs.csv (scaled_prs covariate column)
# Output: prs_nagelkerke_incremental.csv

library(DescTools)
library(boot)
library(dplyr)
library(argparse)


# ============================================================================
# HELPER: safely wrap predictor names containing special characters
# ============================================================================
build_safe_formula <- function(outcome, predictors) {
  safe_predictors <- sapply(predictors, function(x) {
    if (grepl("[^a-zA-Z0-9_\\.]", x)) paste0("`", x, "`") else x
  })
  as.formula(paste(outcome, "~", paste(safe_predictors, collapse = " + ")))
}


# ============================================================================
# FUNCTION: Nagelkerke R² point estimate for a fitted GLM
# ============================================================================
nagelkerke_r2 <- function(model) {
  as.numeric(PseudoR2(model, which = "Nagelkerke"))
}


# ============================================================================
# FUNCTION: Incremental Nagelkerke R² with bootstrapped CI
#
# Bootstraps the paired difference  R²(cov + PRS) - R²(cov only)
# in the same resample so the CI correctly accounts for correlation
# between the two models sharing the same individuals and covariates.
#
# Parameters:
#   data           - data.frame containing PHENOTYPE (0/1), PRS col, cov cols
#   prs_col        - single scaled_prs_* column name
#   covariate_cols - character vector of covariate column names
#   R              - number of bootstrap replicates (default 2000)
#   conf.level     - confidence level (default 0.95)
#   ci.type        - boot.ci type: "bca" (default, bias-corrected) or "perc"
#   seed           - random seed for reproducibility
#
# Returns a named list:
#   nagelkerke_cov        - baseline (covariates-only) R2
#   nagelkerke_full       - R2 for covariates + this PRS
#   incremental_r2        - point estimate of the difference
#   ci_lower / ci_upper   - bootstrapped CI for the incremental R2
#   conf_level, ci_type, n_bootstrap
# ============================================================================
calculate_incremental_nagelkerke <- function(data, prs_col, covariate_cols,
                                             R = 2000, conf.level = 0.95,
                                             ci.type = "bca", seed = 123) {

  formula_cov  <- build_safe_formula("PHENOTYPE", covariate_cols)
  formula_full <- build_safe_formula("PHENOTYPE", c(prs_col, covariate_cols))

  # Point estimates on the full dataset
  model_cov  <- glm(formula_cov,  data = data, family = binomial)
  model_full <- glm(formula_full, data = data, family = binomial)

  r2_cov  <- nagelkerke_r2(model_cov)
  r2_full <- nagelkerke_r2(model_full)
  delta   <- r2_full - r2_cov

  cat(sprintf(
    "  %-40s  R2_cov=%.4f  R2_full=%.4f  delta=%.4f\n",
    prs_col, r2_cov, r2_full, delta
  ))

  # Bootstrap: fit BOTH models on the same resample and return the difference.
  # Paired bootstrapping preserves the covariance structure and gives correctly
  # centered CIs; bootstrapping each R2 independently would inflate the CI width.
  boot_fn <- function(data, indices) {
    d      <- data[indices, ]
    m_cov  <- glm(formula_cov,  data = d, family = binomial)
    m_full <- glm(formula_full, data = d, family = binomial)
    nagelkerke_r2(m_full) - nagelkerke_r2(m_cov)
  }

  set.seed(seed)
  boot_res <- boot(data = data, statistic = boot_fn, R = R)

  # BCA is preferred (bias-corrected + accelerated); fall back to percentile
  # if BCA fails (can happen with small n or near-zero incremental R2)
  used_ci_type <- ci.type
  ci <- tryCatch(
    boot.ci(boot_res, conf = conf.level, type = ci.type),
    error = function(e) {
      message(sprintf(
        "  [%s] BCA CI failed for %s: %s -- falling back to percentile",
        ci.type, prs_col, conditionMessage(e)
      ))
      used_ci_type <<- "perc"
      boot.ci(boot_res, conf = conf.level, type = "perc")
    }
  )

  bounds <- switch(used_ci_type,
    bca   = ci$bca[4:5],
    perc  = ci$percent[4:5],
    norm  = ci$normal[2:3],
    basic = ci$basic[4:5]
  )

  list(
    nagelkerke_cov  = r2_cov,
    nagelkerke_full = r2_full,
    incremental_r2  = delta,
    ci_lower        = bounds[1],
    ci_upper        = bounds[2],
    conf_level      = conf.level,
    ci_type         = used_ci_type,
    n_bootstrap     = R
  )
}


# ============================================================================
# MAIN FUNCTION: run incremental Nagelkerke for all discovered PRS models
#
# PRS columns are discovered automatically: every scaled_prs_* column that is
# NOT listed in covariate_cols is treated as a distinct model to evaluate.
# ============================================================================
run_nagelkerke_incremental <- function(data,
                                       covariate_cols,
                                       output_path = NULL,
                                       R = 2000,
                                       conf.level = 0.95,
                                       ci.type = "bca",
                                       seed = 123) {

  # Phenotype coding guard: accept 1/2 or 0/1, standardise to 0/1
  unique_vals <- sort(unique(data$PHENOTYPE))
  if (setequal(unique_vals, c(1, 2))) {
    data$PHENOTYPE <- data$PHENOTYPE - 1
    cat("Converted PHENOTYPE: 1/2 -> 0/1\n")
  } else if (setequal(unique_vals, c(0, 1))) {
    cat("PHENOTYPE already 0/1\n")
  } else {
    stop(sprintf(
      "Unexpected PHENOTYPE values: %s. Expected {1,2} or {0,1}.",
      paste(unique_vals, collapse = ", ")
    ))
  }

  # Discover all scaled_prs_* columns, excluding the covariate column(s)
  all_prs_cols <- grep("^scaled_prs_", names(data), value = TRUE)
  prs_cols     <- setdiff(all_prs_cols, covariate_cols)

  if (length(prs_cols) == 0) {
    stop("No scaled_prs_* columns found after excluding covariates.")
  }

  cat(sprintf("\nDiscovered %d PRS models:\n  %s\n",
              length(prs_cols), paste(prs_cols, collapse = "\n  ")))
  cat(sprintf("Covariate baseline: %s\n\n",
              paste(covariate_cols, collapse = ", ")))

  results_list <- vector("list", length(prs_cols))

  for (i in seq_along(prs_cols)) {
    prs_col <- prs_cols[i]
    cat(sprintf("[%d/%d] %s\n", i, length(prs_cols), prs_col))

    res <- tryCatch(
      calculate_incremental_nagelkerke(
        data           = data,
        prs_col        = prs_col,
        covariate_cols = covariate_cols,
        R              = R,
        conf.level     = conf.level,
        ci.type        = ci.type,
        seed           = seed
      ),
      error = function(e) {
        message(sprintf("  ERROR for %s: %s", prs_col, conditionMessage(e)))
        NULL
      }
    )

    if (!is.null(res)) {
      results_list[[i]] <- data.frame(
        model           = prs_col,
        nagelkerke_cov  = res$nagelkerke_cov,
        nagelkerke_full = res$nagelkerke_full,
        incremental_r2  = res$incremental_r2,
        ci_lower        = res$ci_lower,
        ci_upper        = res$ci_upper,
        conf_level      = res$conf_level,
        ci_type         = res$ci_type,
        n_bootstrap     = res$n_bootstrap,
        stringsAsFactors = FALSE
      )
    }
  }

  results_df <- bind_rows(results_list)

  # Sort by descending incremental R2 so the best-performing model is first
  results_df <- results_df[order(-results_df$incremental_r2), ]
  rownames(results_df) <- NULL

  if (!is.null(output_path)) {
    write.csv(results_df, output_path, row.names = FALSE)
    cat(sprintf("\nResults saved to: %s\n", output_path))
  }

  cat("\n=== Incremental Nagelkerke R2 Summary ===\n")
  print(results_df)

  return(results_df)
}


# ============================================================================
# RUN
# ============================================================================
parser <- ArgumentParser()
parser$add_argument("--pheno_data", required = TRUE,
                    help = "Path to phenotype results directory")
args <- parser$parse_args()

pheno_data  <- args$pheno_data
scores_path <- file.path(pheno_data, "scores")

data_file   <- file.path(scores_path, "combinedPRSGroups.filtered.csv")
cov_file    <- file.path(scores_path, "prsScores", "covariate.12.mixed.prs.csv")
output_file <- file.path(scores_path, "prs_nagelkerke_incremental.csv")

# Load and merge
data     <- read.csv(data_file, check.names = FALSE)
cov_data <- read.csv(cov_file,  check.names = FALSE)

data <- left_join(
  data,
  cov_data %>% select(IID, scaled_prs_covariate = scaled_prs),
  by = "IID"
)

covariate_cols <- "scaled_prs_covariate"

results <- run_nagelkerke_incremental(
  data           = data,
  covariate_cols = covariate_cols,
  output_path    = output_file,
  R              = 2000,
  conf.level     = 0.95,
  ci.type        = "bca",
  seed           = 123
)

# Keri Multerer 2026
# Re-run the AUC / DeLong / Nagelkerke statistics on the HOLDOUT cohort using
# the current five ePRS encodings.
#
# WHY A SEPARATE DRIVER
# prsAUCDelongStats.R is kept unchanged so its outputs remain reproducible.
# Its driver differs from the manuscript text in three ways that together
# explain why the published AUCs cannot be reproduced from current files:
#
#   1. It reads combinedPRSGroups.csv, the TRAINING/validation cohort, while
#      the results text reports these numbers "in the independent holdout
#      cohort". Its scaled_prs_main AUC of 0.5911 matches the training-set
#      value exactly and not the holdout value of 0.5823.
#   2. It builds its own composite at runtime by summing raw PRS columns. That
#      was a side experiment -- comparing a model trained on all features with
#      a sum of separately trained models -- and is not part of the manuscript,
#      so it is not scored here.
#   3. It uses the pre-refactor encoding names (scaled_prs_epi,
#      scaled_prs_cardio, scaled_prs_epi+main, scaled_prs_all), which no longer
#      correspond to the five retained models.
#
# This driver reuses the statistical functions from that script unchanged --
# they are sourced, not copied -- and changes only the cohort, the model set
# and the composite definition.
#
# COMPOSITES COMPARED
#   scaled_prs_ensemble  the fitted combination now reported as ePRS_Comp,
#                        with weights estimated on the validation set and
#                        applied once to the holdout. Passed in from
#                        fitted_prs_ensemble.py so both languages score the
#                        identical vector.
#   scaled_prs_combined  the rule-based composite, retained for comparison.
#
# Output goes to holdout-specific filenames so nothing existing is overwritten.

suppressPackageStartupMessages({
  library(pROC); library(boot); library(dplyr); library(tidyr)
  library(argparse); library(DescTools)
})

parser <- ArgumentParser()
parser$add_argument("--pheno_data", required = TRUE,
                    help = "results/<pheno>/combinedAnalysis")
parser$add_argument("--ensemble_file", default = NULL,
                    help = "CSV with IID,ensemble_prob from fitted_prs_ensemble.py")
parser$add_argument("--source_script", default = NULL,
                    help = "prsAUCDelongStats.R (defaults to alongside this file)")
args <- parser$parse_args()

script_dir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
src <- if (!is.null(args$source_script)) args$source_script else
  file.path(script_dir, "prsAUCDelongStats.R")

# Source only the function definitions, stopping before that script's own
# driver so it does not run its training-set analysis as a side effect.
lines <- readLines(src, warn = FALSE)
cut_at <- grep("^parser <- ArgumentParser\\(\\)", lines)[1] - 1
eval(parse(text = paste(lines[1:cut_at], collapse = "\n")))
cat(sprintf("sourced %d function-definition lines from %s\n", cut_at, basename(src)))

scores_path <- file.path(args$pheno_data, "scores")
data <- read.csv(file.path(scores_path, "combinedPRSGroups.holdout.filtered.csv"),
                 check.names = FALSE)
cov_data <- read.csv(file.path(scores_path, "prsScores",
                               "covariate.holdout.mixed.prs.csv"),
                     check.names = FALSE)

# The holdout covariate file stores the score as `prs`; the training file used
# `scaled_prs`. Accept either so the driver works for both phenotypes.
cov_col <- if ("scaled_prs" %in% names(cov_data)) "scaled_prs" else "prs"
cov_data$scaled_prs_covariate <- as.numeric(scale(cov_data[[cov_col]]))
data <- left_join(data, cov_data %>% select(IID, scaled_prs_covariate), by = "IID")

# Current encodings, intersected with what this phenotype retained as distinct.
candidates <- c("scaled_prs_main", "scaled_prs_epi_summed", "scaled_prs_epi_product",
                "scaled_prs_cardio_summed", "scaled_prs_cardio_product")
prs_cols <- candidates[candidates %in% names(data)]
cat(sprintf("retained encodings: %s\n", paste(prs_cols, collapse = ", ")))

if (!is.null(args$ensemble_file)) {
  ens <- read.csv(args$ensemble_file, check.names = FALSE)
  data <- left_join(data, ens %>% select(IID, scaled_prs_ensemble), by = "IID")
  cat("merged fitted ensemble scores\n")
}

composites <- c("scaled_prs_ensemble", "scaled_prs_combined")
composites <- composites[composites %in% names(data)]
all_cols <- c(prs_cols, composites)

# Each composite against the additive baseline, plus the two composites
# against each other where both are present.
comparisons <- lapply(composites, function(c) c("scaled_prs_main", c))
if (all(c("scaled_prs_ensemble", "scaled_prs_combined") %in% names(data)))
  comparisons <- c(comparisons,
                   list(c("scaled_prs_combined", "scaled_prs_ensemble")))

results <- compare_prs_models(
  data = data,
  prs_columns = all_cols,
  covariate_columns = c("scaled_prs_covariate"),
  phenotype_col = "PHENOTYPE",
  output_path = file.path(scores_path, "prs_nagelkerke_auc_statistics.holdout.csv"),
  pairwise_comparisons = comparisons,
  pairwise_output = file.path(scores_path, "prs_pairwise_comparisons.holdout.csv")
)

cat("\nwrote prs_nagelkerke_auc_statistics.holdout.csv and",
    "prs_pairwise_comparisons.holdout.csv\n")

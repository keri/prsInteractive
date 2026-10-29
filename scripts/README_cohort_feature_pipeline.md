# PRS Cohort Important Feature Analysis

Identifies and characterises genomic and clinical features that are statistically
discriminating for each PRS-derived T2D risk cohort, maps them to T2D functional
clusters and disease subtypes, and produces circos plots and donut pie charts
summarising subtype enrichment.

---

## Overview

Five statistically distinct PRS cohorts are compared:

| Cohort | Description |
|---|---|
| `epi_product` | Epistatic PRS (product aggregation) |
| `epi_summed` | Epistatic PRS (sum aggregation) |
| `cardio_product` | Cardio-metabolic PRS (product aggregation) |
| `cardio_summed` | Cardio-metabolic PRS (sum aggregation) |
| `main` | Main additive PRS |

For each cohort, individuals are stratified into three training groups based on
their PRS risk bin:

| Group | Definition | Phenotype |
|---|---|---|
| **HighCases** | Top 20 % of PRS bin | T2D cases (PHENOTYPE == 2) |
| **LowControls** | Bottom 20 % of PRS bin | Non-cases (PHENOTYPE == 1) |
| **ExclusiveHighCases** | Top 20 % in exactly one cohort only | T2D cases (PHENOTYPE == 2) |

Features (genomic SHAP values and raw clinical measures) that discriminate a
cohort's high-risk or exclusive-high-risk group from its within-cohort controls
are then mapped to published T2D functional clusters (Udler 2019, Suzuki 2024,
Kim 2026) and assigned to SAID / SIDD / SIRD / MOD / MARD subtypes via a
concordance-weighted scoring scheme.

---

## Prerequisites

### Software
- Python ≥ 3.9 with packages: `pandas`, `numpy`, `scipy`, `statsmodels`,
  `scikit-learn`, `matplotlib`, `seaborn`
- All scripts in `scripts/` on `PYTHONPATH` (or run from repo root)

### Input files

| File | Description |
|---|---|
| `results/type2Diabetes/combinedAnalysis/scores/combinedPRSGroups.filtered.csv` | PRS bin assignments (training set) |
| `results/type2Diabetes/combinedAnalysis/scores/combinedPRSGroups.holdout.filtered.csv` | PRS bin assignments (holdout set) |
| `results/type2Diabetes/combinedAnalysis/scores/featureScoresReducedFinalModel.combined.filtered.csv` | LASSO/ElasticNet feature coefficients per cohort model |
| `results/type2Diabetes/combinedAnalysis/scores/prsScores/*.mixed.prs.csv` | Per-cohort PRS score files (training set) |
| `results/participant_environment.csv` | Raw (untransformed) clinical and environmental features; used for Mann-Whitney U comparisons |
| `results/type2Diabetes/subtypeResearch/cluster_subtype_map.csv` | Mapping of functional cluster keys to T2D subtypes |
| `results/type2Diabetes/subtypeResearch/*.csv` | Functional cluster loci / clinical trait maps (Udler, Suzuki, Kim, …) |

---

## Scripts

### 1. `calculate_top_features_in_cohort.py`

Trains a multi-output Random Forest on SHAP-derived feature scores for
**HighCases** and **LowControls**, computes feature importance z-scores across
cohorts, and filters for cohort-discriminating features using one of four
cross-cohort strategies.

**Cross-cohort filtering strategies**

| Strategy | Rule |
|---|---|
| `strict` | z > threshold in matching cohort AND all others < 0 |
| `majority` | z > threshold in matching cohort AND < 0 in ≥ `majority_threshold` fraction of others |
| `differential` | Matching cohort z is the highest AND Δz vs mean of others ≥ `delta_threshold` |
| `tiered` *(default)* | Runs all three; labels each feature with the strictest tier it qualifies for |

**Key outputs** (`scores/importantCohortScores/`)

| File | Contents |
|---|---|
| `importantFeaturesAcrossCohortsAndTrainingData.csv` | Raw per-cohort feature importance before filtering |
| `importantFeaturesAcrossCohortsAndTrainingData.Filtered.<strategy>.csv` | Filtered features with specificity tier, z-score, OR/coef |
| `performanceMetricsAcrossCohortsAndTrainingData.csv` | AUC and classification metrics per cohort |
| `probabilitiesAcrossCohortsAndFeaturesUsedInTraining.HighCases.csv` | Predicted probabilities for HighCases |
| `probabilitiesAcrossCohortsAndFeaturesUsedInTraining.LowControls.csv` | Predicted probabilities for LowControls |

---

### 2. `enrich_exclusive_cohort_features.py`

The primary enrichment script. Loads the tiered feature CSV from Step 1,
applies direction and exclusivity filters, tests clinical features with
Mann-Whitney U, and maps all features to functional clusters.

**Feature filters applied**

| Group | Filters |
|---|---|
| HighCases | OR > 1 (coefs > 0) |
| LowControls | OR < 1 (coefs < 0); removes risk-direction contamination |
| ExclusiveHighCases | `n_exclusive_cohorts == 1` AND coefs > 0; OR derived as exp(coefs) |

**Clinical feature comparison**  
For ExclusiveHighCases: exclusive top-20% cases vs within-cohort bottom-20%
controls (Mann-Whitney U with BH FDR correction and rank-biserial r effect
size). Unbalanced (all controls) or balanced (bootstrapped) modes available.

**Key outputs** (`scores/functionalEnrichment/`)

| File | Contents |
|---|---|
| `functional_cluster_enrichment_<strategy>.csv` | One row per feature × cohort with cluster binary columns, OR, effect size, feature source (`genomic_shap` / `clinical`), match method |

---

### 3. `compare_clinical_features_across_cohorts.py`

Broader cross-cohort clinical comparison across all three risk tiers.

| Tier | Comparison design |
|---|---|
| **HighRisk** | Each top-20% cohort vs all other top-20% cohorts combined |
| **LowRisk** | Each bottom-20% cohort vs all other bottom-20% cohorts combined |
| **ExclusiveHighRisk** | Exclusive top-20% cases (one cohort only) vs within-cohort bottom-20% controls |

ExclusiveHighRisk uses the same case-vs-control design as `identify_exclusive_groups()`
in `enrich_exclusive_cohort_features.py`. When a feature is significant in both
HighRisk and ExclusiveHighRisk for the same cohort, the ExclusiveHighRisk result
supersedes the HighRisk result in the priority output.

**Key outputs** (`scores/cohortClinicalComparison/`)

| File | Contents |
|---|---|
| `clinical_comparison_HighRisk_<set>.csv` | HighRisk tier results |
| `clinical_comparison_LowRisk_<set>.csv` | LowRisk tier results |
| `clinical_comparison_ExclusiveHighRisk_<set>.csv` | ExclusiveHighRisk tier results |
| `clinical_comparison_priority_<set>.csv` | Combined priority output (ExclusiveHighRisk supersedes HighRisk) |
| `epi_product_clinical_summary_<set>.csv` | Concise summary for epi_product cohort |

---

### 4. `clusterSubtypeAssignment.py`

Reads the functional enrichment CSV, joins the cluster-subtype map, and scores
each feature against the five T2D subtypes (SAID, SIDD, SIRD, MOD, MARD) using
concordance weighting (+1 concordant, 0 neutral, −1 discordant).

**Key outputs** (`subtypeResearch/enrichment/`)

| File | Contents |
|---|---|
| `subtype_feature_assignments.csv` | One row per feature × cohort with subtype, concordance, and cluster_key |
| `subtype_cohort_enrichment.csv` | Concordance-weighted feature count per subtype × cohort |
| `subtype_cluster_cohort_enrichment.csv` | Enrichment broken down by cluster |
| `subtype_direct_evidence.csv` | Direct-evidence subtype assignments (not via cluster route) |

---

### 5. `circos_subtype_cohort.py`

Generates a circos-style bar chart per cohort showing each feature as an outward
(risk) or inward (protective) bar coloured by feature source and annotated with
T2D subtype concordance dots.

Four output variants per cohort:

| Suffix | Labels | Protective bars |
|---|---|---|
| `_labelled` | Yes | Yes |
| `_clean` | No | Yes |
| `_risk_labelled` | Yes | No |
| `_risk_clean` | No | No |

Bar colour scheme: red = genomic (outward solid / inward hashed `//`);
yellow-green = clinical (outward solid / inward hashed `//`).

**Output:** `subtypeResearch/enrichment/plots/circos_<cohort>_*.png` (20 files)

---

### 6. `subtype_cohort_pies.py`

Donut pie charts showing the concordance-weighted cohort contribution to each
T2D subtype's feature pool.

**Output:** `subtypeResearch/enrichment/plots/subtype_cohort_pies_labelled.png`,
`_clean.png`, plus 10 individual per-subtype PNGs.

---

## Run Order and Commands

> **Quick start:** use the pipeline script (see next section) instead of running
> each command manually.

### Step 1 — SHAP feature importance
```bash
python scripts/calculate_top_features_in_cohort.py \
    --pheno_data                   results/type2Diabetes/combinedAnalysis \
    --feature_scores_file_filtered results/type2Diabetes/combinedAnalysis/scores/featureScoresReducedFinalModel.combined.filtered.csv \
    --filter_strategy              tiered
```

### Step 2 — Functional enrichment
```bash
python scripts/enrich_exclusive_cohort_features.py \
    --pheno_data         results/type2Diabetes/combinedAnalysis \
    --raw_feature_scores results/participant_environment.csv \
    --use_set            holdout \
    --balance_strategy   unbalanced
```

### Step 3 — Clinical comparison *(can run in parallel with Step 2)*
```bash
python scripts/compare_clinical_features_across_cohorts.py \
    --pheno_data        results/type2Diabetes/combinedAnalysis \
    --raw_features_file results/participant_environment.csv \
    --use_set           holdout
```

### Step 4 — Subtype assignment
```bash
python scripts/clusterSubtypeAssignment.py \
    --tiered_csv          results/type2Diabetes/combinedAnalysis/scores/functionalEnrichment/functional_cluster_enrichment_tiered.csv \
    --cluster_subtype_map results/type2Diabetes/subtypeResearch/cluster_subtype_map.csv \
    --out_dir             results/type2Diabetes/subtypeResearch/enrichment
```

### Steps 5 & 6 — Plots *(can run in parallel)*
```bash
python scripts/circos_subtype_cohort.py \
    --assignments results/type2Diabetes/subtypeResearch/enrichment/subtype_feature_assignments.csv \
    --out_dir     results/type2Diabetes/subtypeResearch/enrichment/plots

python scripts/subtype_cohort_pies.py \
    --assignments results/type2Diabetes/subtypeResearch/enrichment/subtype_feature_assignments.csv \
    --out_dir     results/type2Diabetes/subtypeResearch/enrichment/plots
```

---

## Pipeline Script

`run_cohort_feature_pipeline.sh` runs all six steps in order with consistent
path configuration.

```
Usage:
  bash scripts/run_cohort_feature_pipeline.sh [OPTIONS]

Options:
  --pheno-data        PATH   results/type2Diabetes/combinedAnalysis  (default)
  --subtype-dir       PATH   results/type2Diabetes/subtypeResearch    (default)
  --feature-scores    PATH   <pheno-data>/scores/featureScoresReducedFinalModel.combined.filtered.csv
  --raw-features      PATH   results/participant_environment.csv      (default)
  --use-set           STR    holdout | validation | all               (default: holdout)
  --filter-strategy   STR    tiered | strict | majority | differential (default: tiered)
  --balance-strategy  STR    unbalanced | balanced                    (default: unbalanced)
  --skip-calculate           Skip Step 1; use the existing filtered feature CSV
  --only-plots               Skip Steps 1–4; regenerate plots only
  --help
```

**Full run (all steps):**
```bash
bash scripts/run_cohort_feature_pipeline.sh
```

**Skip Step 1** when feature scores are already up to date:
```bash
bash scripts/run_cohort_feature_pipeline.sh --skip-calculate
```

**Regenerate plots only** (Steps 5 & 6):
```bash
bash scripts/run_cohort_feature_pipeline.sh --only-plots
```

**Custom paths / validation set:**
```bash
bash scripts/run_cohort_feature_pipeline.sh \
    --pheno-data results/type2Diabetes/combinedAnalysis \
    --raw-features results/participant_environment.csv \
    --use-set validation \
    --skip-calculate
```

All path defaults can also be set via environment variables:
`PHENO_DATA`, `SUBTYPE_DIR`, `FEATURE_SCORES_FILE_FILTERED`,
`RAW_FEATURE_SCORES`, `USE_SET`, `FILTER_STRATEGY`, `BALANCE_STRATEGY`.

---

## Output Directory Structure

```
results/type2Diabetes/
├── combinedAnalysis/scores/
│   ├── importantCohortScores/
│   │   ├── importantFeaturesAcrossCohortsAndTrainingData.csv
│   │   ├── importantFeaturesAcrossCohortsAndTrainingData.Filtered.tiered.csv
│   │   ├── performanceMetricsAcrossCohortsAndTrainingData.csv
│   │   └── probabilitiesAcrossCohortsAndFeaturesUsedInTraining.{HighCases,LowControls}.csv
│   ├── functionalEnrichment/
│   │   └── functional_cluster_enrichment_tiered.csv
│   └── cohortClinicalComparison/
│       ├── clinical_comparison_{HighRisk,LowRisk,ExclusiveHighRisk}_holdout.csv
│       ├── clinical_comparison_priority_holdout.csv
│       └── epi_product_clinical_summary_holdout.csv
└── subtypeResearch/enrichment/
    ├── subtype_feature_assignments.csv
    ├── subtype_cohort_enrichment.csv
    ├── subtype_cluster_cohort_enrichment.csv
    ├── subtype_direct_evidence.csv
    └── plots/
        ├── circos_{cohort}_{labelled,clean,risk_labelled,risk_clean}.png  (×5 cohorts = 20 files)
        ├── subtype_cohort_pies_labelled.png
        ├── subtype_cohort_pies_clean.png
        └── subtype_pie_{SAID,SIDD,SIRD,MOD,MARD}_{labelled,clean}.png    (10 files)
```

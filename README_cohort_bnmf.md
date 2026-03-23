# Cohort bNMF Subtype Analysis

Unsupervised subtype discovery within high-risk PRS strata across multiple cohorts, followed by biological enrichment using GTEx tissue data and literature mining.

---

## Overview

The pipeline has three stages:

```
Stage 1 — bNMF subtype analysis
  ↓  per-cohort cluster_assignments.csv + feature_loadings.csv

Stage 2 — Cluster enrichment (run per cohort)
  ├── GTEx: tissue expression + eQTL evidence for top genomic features
  └── Literature: DisGeNET · Open Targets · PubTator3 · PubMed · ClinVar

Stage 3 (optional) — Combined bNMF
  ↓  cross-cohort cluster_assignments.csv with shared feature space
```

---

## Prerequisites

### Software

```bash
pip install numpy pandas scipy scikit-learn matplotlib seaborn requests
```

The script imports `prs_bnmf_subtype_analysis.py` for shared NMF routines — both files must be in the same directory.

### Upstream pipeline outputs

Run the cohort analysis pipeline first. The following files must exist before running bNMF:

| File | Source |
|------|--------|
| `{pheno_data}/scores/importantCohortScores/importantFeaturesAcrossCohortsAndTrainingData.Filtered.{strategy}.csv` | `calculate_top_features_in_cohort` |
| `{pheno_data}/scores/cohortClinicalComparison/clinical_comparison_all_{use_set}.csv` | `compare_clinical_features_across_cohorts` |
| `{pheno_data}/scores/combinedPRSGroups.holdout.filtered.csv` | PRS scoring pipeline |
| `{pheno_data}/scores/prsScores/*.mixed.prs.csv` | PRS model scoring |

---

## Stage 1 — bNMF Subtype Analysis

### Quick start

```bash
PHENO_DATA=/path/to/results \
RAW_FEATURES_FILE=/path/to/participant_environment.csv \
./workflows/run_cohort_bnmf_subtype_analysis.sh
```

This runs **both** per-cohort and combined bNMF with default settings.

### Environment variables

#### Required

| Variable | Description |
|----------|-------------|
| `PHENO_DATA` | Path to phenotype results directory |
| `RAW_FEATURES_FILE` | Path to raw (untransformed) clinical features CSV. Must contain one row per IID with raw feature values. |

#### Optional (with defaults)

| Variable | Default | Description |
|----------|---------|-------------|
| `MODE` | `both` | `per_cohort` · `combined` · `both` |
| `FILTER_STRATEGY` | `tiered` | Genomic feature filter: `tiered` · `strict` · `majority` · `differential` |
| `USE_SET` | `holdout` | Individual set: `holdout` · `validation` · `all` |
| `K_MIN` | `2` | Minimum number of clusters to test |
| `K_MAX` | `6` | Maximum number of clusters to test |
| `N_RUNS` | `30` | NMF restarts per k (higher = more stable consensus) |
| `K_SELECT_METHOD` | `elbow` | Optimal-k selection: `elbow` · `max` |
| `L1_RATIO` | `0.1` | L1 sparsity regularisation (0 = L2 only, 1 = L1 only) |
| `ALPHA_W` | `0.1` | Regularisation strength on W (basis) matrix |
| `MIN_EFFECT_SIZE` | `0.0` | Minimum |rank-biserial r| for clinical feature inclusion |
| `SPECIFICITY_TIER_MAX` | `3` | Maximum genomic specificity tier (1 = most cohort-specific) |
| `HIGH_RISK_THRESHOLD` | `8` | bin > this → high-risk cases (decile 9–10, top 20%) |
| `LOW_CONTROL_THRESHOLD` | `3` | bin < this → low-control individuals (decile 1–2, bottom 20%) |
| `COHORTS` | all | Space-separated cohort names to restrict to |
| `NO_HOLDOUT_PRS` | `0` | Set to `1` to exclude holdout PRS files |

> **Bin scale note:** If `combinedPRSGroups.holdout.filtered.csv` uses millile bins (1–1000) instead of decile bins (1–10), the script detects this automatically at load time and normalises. Thresholds are always interpreted on the 1–10 decile scale.

### Individual selection

For each cohort:
- **High-risk cases**: `bin_{cohort} > HIGH_RISK_THRESHOLD` AND `PHENOTYPE == 2`
- **Low controls**: `bin_{cohort} < LOW_CONTROL_THRESHOLD` AND `PHENOTYPE == 1`

`LOW_CONTROL_THRESHOLD = 3` matches the `controls_lt3` criterion used in the upstream cohort clinical/genomic feature analysis.

### Feature types

| Prefix | Source | Description |
|--------|--------|-------------|
| `clin_` | Raw clinical CSV | FDR-significant clinical features from cross-cohort comparison |
| `gen_` | `*.mixed.prs.csv` | Per-individual weighted SHAP contributions for top SHAP-important genomic features |

In combined mode, clinical features are de-duplicated to one column per feature weighted by max effect size across cohorts. PRS anchors are suppressed for cohorts that already have `gen_` features (avoids near-collinearity since PRS ≈ Σ weighted genomic contributions).

### Per-cohort outputs

```
{pheno_data}/scores/cohortBnmf/{cohort}/
  cluster_assignments.csv     — IID, cluster_label, membership_cluster_* (soft), is_high_risk
  feature_loadings.csv        — H matrix: feature weights per cluster
  basis_matrix.csv            — W matrix: per-individual cluster memberships
  cophenetic_curve.csv        — cophenetic correlation at each k tested
  cluster_profile.csv         — per-cluster median / IQR for all features

{pheno_data}/scores/cohortBnmf/
  cohort_bnmf_overview.csv    — one-row summary per cohort (k selected, n_individuals, etc.)

{pheno_data}/figures/cohortBnmf/{cohort}/
  cophenetic_curve.png        — used to verify k selection
  consensus_map_k{k}.png      — consensus matrix heatmap at selected k
  cluster_profile.png         — feature distributions per cluster
  feature_loadings.png        — top features per cluster
```

### Combined mode additional outputs

```
{pheno_data}/scores/combinedCohortBnmf/
  cluster_assignments.csv         — membership + is_high_risk_{cohort} flags per individual
  feature_loadings.csv            — H matrix with shared clin_/gen_ columns
  cohort_cluster_affinity.csv     — proportion of each cohort's high-risk individuals per cluster
  weighted_feature_list.csv       — per-cluster features ranked by loading × cohort weight
  cluster_dominant_cohort.csv     — which cohort drives each cluster

{pheno_data}/figures/combinedCohortBnmf/
  cohort_cluster_affinity.png     — heatmap of cohort–cluster affinities
```

---

## Stage 2 — Cluster Enrichment

Run enrichment for each cohort directory produced in Stage 1. No API keys are required for Open Targets or PubTator3.

### Quick start — single cohort

```bash
CLUSTER_DIR=/path/to/results/scores/cohortBnmf/cardio \
PHENOTYPES="type 2 diabetes,cardiometabolic disease" \
./workflows/run_cluster_enrichment.sh
```

### Quick start — combined

```bash
CLUSTER_DIR=/path/to/results/scores/combinedCohortBnmf \
PHENOTYPES="type 2 diabetes,epilepsy,cardiometabolic disease" \
./workflows/run_cluster_enrichment.sh
```

### Run enrichment for all cohorts

```bash
PHENO_DATA=/path/to/results
PHENOTYPES="type 2 diabetes,epilepsy,cardiometabolic disease"

for cohort_dir in "${PHENO_DATA}/scores/cohortBnmf"/*/; do
    cohort=$(basename "${cohort_dir}")
    echo "=== Enriching ${cohort} ==="
    CLUSTER_DIR="${cohort_dir}" \
    PHENOTYPES="${PHENOTYPES}" \
    ./workflows/run_cluster_enrichment.sh
done
```

### Environment variables

#### Required

| Variable | Description |
|----------|-------------|
| `CLUSTER_DIR` | Path to a bNMF output directory containing `feature_loadings.csv` |
| `PHENOTYPES` | Comma-separated phenotype terms for literature queries (e.g. `"type 2 diabetes,epilepsy"`) |

#### Optional

| Variable | Default | Description |
|----------|---------|-------------|
| `TOP_N` | `20` | Top features per cluster to query |
| `TISSUES` | cardiometabolic panel | Comma-separated GTEx tissue IDs to restrict |
| `NCBI_API_KEY` | none | NCBI Entrez API key (3 req/s without, 10 req/s with) |
| `DISGENET_KEY` | none | DisGeNET API key (unauthenticated by default) |
| `SKIP_GTEX` | `0` | Set to `1` to skip GTEx enrichment |
| `SKIP_LITERATURE` | `0` | Set to `1` to skip literature enrichment |

### Evidence sources

#### GTEx (`enrich_clusters_gtex.py`)

Extracts top `gen_` features per cluster → maps SNP rs-IDs to gene symbols via NCBI dbSNP → queries GTEx REST API for:
- Tissue-specific median TPM per gene
- eQTL effect sizes per SNP × gene × tissue

Composite score: `gtex_score = |eqtl_effect_size| × median_TPM`

#### Literature (`enrich_clusters_literature.py`)

Queries four sources per gene/SNP × phenotype pair:

| Source | What it provides | Key required? |
|--------|-----------------|---------------|
| **DisGeNET** | Gene-disease association (GDA) score | Optional (free tier available) |
| **Open Targets** | Integrated GWAS + literature + functional genomics evidence | No |
| **PubTator3** | Entity-annotated publication co-occurrence (gene/SNP/disease as biological entities) | No |
| **PubMed** | Keyword co-occurrence count via NCBI E-utilities | No (optional key for higher rate limit) |
| **ClinVar** | Variant clinical significance | No |

Combined score weights: DisGeNET × 0.25 + Open Targets × 0.30 + PubTator3 × 0.25 + ClinVar × 0.20

All API responses are cached in `{cluster_dir}/enrichment/.cache/` to avoid repeat queries.

### Enrichment outputs

```
{cluster_dir}/enrichment/gtex/
  gtex_gene_expression.csv          — gene × tissue median TPM
  gtex_eqtl_evidence.csv            — snp, gene, tissue, effect_size, p_value
  gtex_cluster_tissue_summary.csv   — cluster × tissue enrichment scores
  gtex_cluster_tissue_heatmap.png   — heatmap of top tissues per cluster

{cluster_dir}/enrichment/literature/
  gene_evidence_table.csv           — gene × phenotype evidence from all sources
    columns: disgenet_score, open_targets_score, pubtator3_count,
             pubmed_count, clinvar_sig, combined_score
  cluster_evidence_ranked.csv       — per-cluster features ranked by combined evidence
  cluster_literature_summary.csv    — dominant phenotype and top gene per cluster
```

---

## Stage 3 — Optional: Combined bNMF only

To skip per-cohort bNMF and only run the cross-cohort combined analysis:

```bash
PHENO_DATA=/path/to/results \
RAW_FEATURES_FILE=/path/to/participant_environment.csv \
MODE=combined \
./workflows/run_cohort_bnmf_subtype_analysis.sh
```

Then enrich:

```bash
CLUSTER_DIR=/path/to/results/scores/combinedCohortBnmf \
PHENOTYPES="type 2 diabetes,epilepsy" \
./workflows/run_cluster_enrichment.sh
```

---

## Tuning guide

| Parameter | When to adjust |
|-----------|---------------|
| `K_MAX` | Increase if cohort is large (>500 individuals) and you expect more subtypes |
| `N_RUNS` | Increase to 50–100 for final results; use 10–15 for exploratory runs |
| `MIN_EFFECT_SIZE` | Set to 0.2–0.3 to restrict to strongly significant clinical features |
| `SPECIFICITY_TIER_MAX` | Set to 1 or 2 to keep only cohort-specific genomic features |
| `TOP_N` | Increase to 30–50 for enrichment if clusters have many features with similar loadings |
| `TISSUES` | Restrict to relevant tissues to reduce API calls and focus enrichment |

### Verifying cluster quality

```python
import pandas as pd

# Check cluster distribution — should not be collapsed to one cluster
df = pd.read_csv('/path/to/cohortBnmf/cardio/cluster_assignments.csv')
print(df['cluster_label'].value_counts())

# Check membership separation — well-separated clusters have high max membership
mem_cols = [c for c in df.columns if c.startswith('membership_')]
df['max_membership'] = df[mem_cols].max(axis=1)
print(df['max_membership'].describe())
# mean > 0.7 indicates good cluster separation
```

---

## Full example

```bash
# Set paths
export PHENO_DATA=/data/epilepsy_cardio_prs
export RAW_FEATURES_FILE=/data/participant_environment.csv
export PHENOTYPES="type 2 diabetes,epilepsy,cardiometabolic disease"

# Stage 1: bNMF (per-cohort + combined)
PHENO_DATA=${PHENO_DATA} \
RAW_FEATURES_FILE=${RAW_FEATURES_FILE} \
K_MIN=2 \
K_MAX=5 \
N_RUNS=30 \
./workflows/run_cohort_bnmf_subtype_analysis.sh

# Stage 2: Enrich each cohort
for cohort_dir in "${PHENO_DATA}/scores/cohortBnmf"/*/; do
    CLUSTER_DIR="${cohort_dir}" \
    PHENOTYPES="${PHENOTYPES}" \
    TOP_N=20 \
    ./workflows/run_cluster_enrichment.sh
done

# Stage 2: Enrich combined
CLUSTER_DIR="${PHENO_DATA}/scores/combinedCohortBnmf" \
PHENOTYPES="${PHENOTYPES}" \
./workflows/run_cluster_enrichment.sh
```

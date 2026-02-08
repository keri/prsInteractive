#!/usr/bin/env python3

#!/usr/bin/env python3

import pandas as pd
import numpy as np
from statsmodels.stats.contingency_tables import mcnemar
from scipy.stats import pearsonr, spearmanr, hypergeom, norm, ttest_rel, chi2_contingency
from sklearn.metrics import cohen_kappa_score, roc_auc_score, roc_curve
from statsmodels.stats.proportion import proportion_confint
from typing import Tuple, Dict, List
import warnings
import argparse


# ============================================================================
# UTILITY FUNCTIONS
# ============================================================================

def pearson_confidence_interval(corr: float, n: int, alpha: float = 0.05) -> Tuple[float, float]:
        """Calculate confidence interval for Pearson correlation using Fisher Z-transformation."""
        fisher_z = np.arctanh(corr)
        se = 1 / np.sqrt(n - 3)
        z_crit = norm.ppf(1 - alpha / 2)
        z_lower = fisher_z - z_crit * se
        z_upper = fisher_z + z_crit * se
        return np.tanh(z_lower), np.tanh(z_upper)


def create_contingency_table(prs1_high: pd.Series, prs2_high: pd.Series, 
                                                            cases: pd.DataFrame) -> Tuple[np.ndarray, Dict]:
        """Create 2x2 contingency table for McNemar's test."""
        a = cases[prs1_high & prs2_high].shape[0]  # Both high
        b = cases[prs1_high & ~prs2_high].shape[0]  # Only PRS1 high
        c = cases[~prs1_high & prs2_high].shape[0]  # Only PRS2 high
        d = cases[~prs1_high & ~prs2_high].shape[0]  # Both low
    
        table = np.array([[a, b], [c, d]])
        counts = {'a': a, 'b': b, 'c': c, 'd': d}
    
        # Validate
        assert a + b + c + d == len(cases), f"Table sum error: {a+b+c+d} != {len(cases)}"
    
        return table, counts


def print_contingency_table(counts: Dict, prs1_name: str, prs2_name: str, 
                                                        threshold: int, total: int):
        """Print formatted contingency table."""
        a, b, c, d = counts['a'], counts['b'], counts['c'], counts['d']
    
        print(f"\n{'='*70}")
        print(f"McNemar's Test: {prs1_name} vs {prs2_name}")
        print(f"Threshold: bin > {threshold}")
        print(f"{'='*70}")
        print(f"\n2x2 Contingency Table:")
        print(f"                      {prs2_name} High  {prs2_name} Low     Total")
        print(f"{prs1_name} High       {a:6d}         {b:6d}       {a+b:6d}")
        print(f"{prs1_name} Low        {c:6d}         {d:6d}       {c+d:6d}")
        print(f"Total                  {a+c:6d}       {b+d:6d}     {total:6d}")
    
        print(f"\nInterpretation:")
        print(f"  Concordant pairs:   {a+d:6d} ({100*(a+d)/total:.1f}%)")
        print(f"    - Both high:      {a:6d}")
        print(f"    - Both low:       {d:6d}")
        print(f"  Discordant pairs:   {b+c:6d} ({100*(b+c)/total:.1f}%)")
        print(f"    - Only {prs1_name}: {b:6d}")
        print(f"    - Only {prs2_name}: {c:6d}")
    
    
# ============================================================================
# STATISTICAL TEST FUNCTIONS
# ============================================================================
    
def mcnemar_test(prs1_high: pd.Series, prs2_high: pd.Series, 
                                cases: pd.DataFrame, prs1_name: str, prs2_name: str, 
                                threshold: int) -> Dict:
        """Perform McNemar's test for paired binary outcomes."""
        table, counts = create_contingency_table(prs1_high, prs2_high, cases)
        print_contingency_table(counts, prs1_name, prs2_name, threshold, len(cases))
    
        result = mcnemar(table, exact=False)
    
        print(f"\nMcNemar test statistic: {result.statistic:.4f}")
        print(f"P-value: {result.pvalue:.4e}")
    
        if result.pvalue >= 0.05:
                print(f"   → Methods identify similar NUMBERS of cases")
        else:
                winner = prs1_name if counts['b'] > counts['c'] else prs2_name
                print(f"   → {winner} identifies significantly more cases")
            
        return {
                'test': 'McNemar',
                'statistic': result.statistic,
                'pvalue': result.pvalue,
                'conf_int': 'NA'
        }
        

def overlap_uniqueness_test(counts: Dict, n_total: int) -> Dict:
        """Test if overlap is less than expected by chance (hypergeometric test)."""
        print(f"\n{'='*70}")
        print("Overlap Uniqueness Test (Hypergeometric)")
        print(f"{'='*70}")
    
        a, b, c = counts['a'], counts['b'], counts['c']
    
        # Expected overlap under independence
        p_prs1 = (a + b) / n_total
        p_prs2 = (a + c) / n_total
        expected_overlap = n_total * p_prs1 * p_prs2
    
        print(f"Expected overlap (by chance): {expected_overlap:.1f}")
        print(f"Observed overlap: {a}")
    
        # Hypergeometric test: is overlap less than expected?
        p_value = hypergeom.cdf(a, n_total, a+c, a+b)
    
        print(f"P-value: {p_value:.4e}")
    
        if p_value < 0.05:
                print(f"   → YES: Overlap is LESS than chance")
                print(f"   → Methods identify MORE unique cases than expected")
                print(f"   → Methods are COMPLEMENTARY")
        else:
                print(f"   → Overlap is not significantly less than chance")
            
        return {
                'test': 'overlap_uniqueness',
                'statistic': a / expected_overlap if expected_overlap > 0 else np.nan,
                'pvalue': p_value,
                'conf_int': 'NA'
        }
        

def cohen_kappa_test(prs1_high: pd.Series, prs2_high: pd.Series) -> Dict:
        """Calculate Cohen's Kappa for agreement between classifications."""
        print(f"\n{'='*70}")
        print("Cohen's Kappa (Agreement)")
        print(f"{'='*70}")
    
        kappa = cohen_kappa_score(prs1_high, prs2_high)
    
        print(f"Kappa: {kappa:.3f}")
        if kappa < 0.20:
                print(f"   → Slight agreement - methods largely identify different cases")
        elif kappa < 0.40:
                print(f"   → Fair agreement")
        elif kappa < 0.60:
                print(f"   → Moderate agreement")
        else:
                print(f"   → Substantial agreement")
            
        return {
                'test': 'Cohen_Kappa',
                'statistic': kappa,
                'pvalue': 'NA',
                'conf_int': 'NA'
        }
        

def jaccard_index(counts: Dict) -> Dict:
        """Calculate Jaccard index for set overlap."""
        print(f"\n{'='*70}")
        print("Jaccard Index (Set Overlap)")
        print(f"{'='*70}")
    
        a, b, c = counts['a'], counts['b'], counts['c']
    
        intersection = a
        union = a + b + c
        jaccard = intersection / union if union > 0 else 0
    
        print(f"Jaccard: {jaccard:.3f}")
        print(f"   → {jaccard*100:.1f}% of high-risk cases identified by BOTH")
        print(f"   → {(1-jaccard)*100:.1f}% unique to one method")
    
        return {
                'test': 'Jaccard_Overlap',
                'statistic': jaccard,
                'pvalue': 'NA',
                'conf_int': 'NA'
        }
        

def discordance_proportion(counts: Dict, n_total: int) -> Dict:
        """Calculate proportion of discordant pairs with confidence interval."""
        print(f"\n{'='*70}")
        print("Discordance Rate")
        print(f"{'='*70}")
    
        b, c = counts['b'], counts['c']
        n_discordant = b + c
    
        prop_discordant = n_discordant / n_total
        ci_low, ci_high = proportion_confint(n_discordant, n_total, method='wilson')
    
        print(f"Discordant: {n_discordant}/{n_total} ({100*prop_discordant:.1f}%)")
        print(f"95% CI: [{100*ci_low:.1f}%, {100*ci_high:.1f}%]")
        print(f"Breakdown:")
        print(f"  Only PRS1: {b} ({100*b/n_discordant:.1f}% of discordant)")
        print(f"  Only PRS2: {c} ({100*c/n_discordant:.1f}% of discordant)")
    
        return {
                'test': 'discordance_proportion',
                'statistic': prop_discordant,
                'pvalue': 'NA',
                'conf_int': (ci_low, ci_high)
        }
        

def chi_square_independence(table: np.ndarray) -> Dict:
        """Test independence using Chi-square."""
        print(f"\n{'='*70}")
        print("Chi-Square Test (Independence)")
        print(f"{'='*70}")
    
        chi2, p_value, dof, expected = chi2_contingency(table)
    
        print(f"Chi-square: {chi2:.4f}")
        print(f"P-value: {p_value:.4e}")
    
        if p_value < 0.05:
                print(f"   → Classifications are NOT independent")
        else:
                print(f"   → Classifications appear independent")
            
        return {
                'test': 'chi_squared',
                'statistic': chi2,
                'pvalue': p_value,
                'conf_int': f'dof {dof}: expected {expected}'
        }
        

def pearson_correlation_test(data: pd.DataFrame, prs1_name: str, 
                                                            prs2_name: str, test_threshold: int) -> Dict:
        """Calculate Pearson correlation for high-risk cases."""
        print(f"\n{'='*70}")
        print("Pearson Correlation (High-Risk Cases)")
        print(f"{'='*70}")
    
        # Get high risk cases (either PRS is high)
        high_risk = data[
                (data[f'bin_{prs1_name}'] > test_threshold) | 
                (data[f'bin_{prs2_name}'] > test_threshold)
        ]
    
        corr, p_value = pearsonr(
                high_risk[f'scaled_prs_{prs1_name}'],
                high_risk[f'scaled_prs_{prs2_name}']
        )
    
        r_squared = corr ** 2
        n = len(high_risk)
        ci_low, ci_high = pearson_confidence_interval(corr, n)
    
        print(f"N high-risk cases: {n}")
        print(f"Pearson r: {corr:.3f}")
        print(f"R²: {r_squared:.3f}")
        print(f"P-value: {p_value:.4e}")
        print(f"95% CI: [{ci_low:.3f}, {ci_high:.3f}]")
    
        return {
                'test': 'pearson_r',
                'statistic': corr,
                'pvalue': p_value,
                'conf_int': (ci_low, ci_high),
                'n_cases': n
        }
        

def spearman_correlation_test(data: pd.DataFrame, prs1_name: str, 
                                                            prs2_name: str, test_threshold: int) -> Dict:
        """Calculate Spearman correlation for high-risk cases."""
        print(f"\n{'='*70}")
        print("Spearman Correlation (High-Risk Cases)")
        print(f"{'='*70}")
    
        high_risk = data[
                (data[f'bin_{prs1_name}'] > test_threshold) | 
                (data[f'bin_{prs2_name}'] > test_threshold)
        ]
    
        corr, p_value = spearmanr(
                high_risk[f'bin_{prs1_name}'],
                high_risk[f'bin_{prs2_name}']
        )
    
        print(f"Spearman ρ: {corr:.3f}")
        print(f"P-value: {p_value:.4e}")
    
        return {
                'test': 'spearman_r',
                'statistic': corr,
                'pvalue': p_value,
                'conf_int': 'NA'
        }
        

def paired_ttest(data: pd.DataFrame, prs1_name: str, prs2_name: str, 
                                test_threshold: int) -> Dict:
        """Perform paired t-test on high-risk cases."""
        print(f"\n{'='*70}")
        print("Paired T-Test (High-Risk Cases)")
        print(f"{'='*70}")
    
        high_risk = data[
                (data[f'bin_{prs1_name}'] > test_threshold) | 
                (data[f'bin_{prs2_name}'] > test_threshold)
        ]
    
        t_stat, p_value = ttest_rel(
                high_risk[f'scaled_prs_{prs1_name}'],
                high_risk[f'scaled_prs_{prs2_name}']
        )
    
        print(f"t-statistic: {t_stat:.3f}")
        print(f"P-value: {p_value:.4e}")
    
        return {
                'test': 'paired_ttest',
                'statistic': t_stat,
                'pvalue': p_value,
                'conf_int': 'NA'
        }
        

def delong_test(y_true: np.ndarray, scores1: np.ndarray, 
                                scores2: np.ndarray) -> Dict:
        """
        DeLong's test for comparing AUCs of two correlated ROC curves.
        
        Based on:
        DeLong et al. (1988) "Comparing the areas under two or more correlated 
        receiver operating characteristic curves: a nonparametric approach"
        """
        print(f"\n{'='*70}")
        print("DeLong Test (AUC Comparison)")
        print(f"{'='*70}")
    
        # Calculate AUCs
        auc1 = roc_auc_score(y_true, scores1)
        auc2 = roc_auc_score(y_true, scores2)
    
        print(f"AUC 1: {auc1:.4f}")
        print(f"AUC 2: {auc2:.4f}")
        print(f"Δ AUC: {abs(auc1 - auc2):.4f}")
    
        # Get number of positive and negative samples
        n_pos = np.sum(y_true == 1)
        n_neg = np.sum(y_true == 0)
    
        # Get indices
        pos_idx = np.where(y_true == 1)[0]
        neg_idx = np.where(y_true == 0)[0]
    
        # Calculate structural components (V matrices)
        def structural_component(scores, pos_idx, neg_idx):
                """Calculate V10 (structural component for AUC variance)."""
                V = np.zeros((len(pos_idx), len(neg_idx)))
                for i, p in enumerate(pos_idx):
                        for j, n in enumerate(neg_idx):
                                if scores[p] > scores[n]:
                                        V[i, j] = 1
                                elif scores[p] == scores[n]:
                                        V[i, j] = 0.5
                return V
    
        V1 = structural_component(scores1, pos_idx, neg_idx)
        V2 = structural_component(scores2, pos_idx, neg_idx)
    
        # Calculate AUC from V (should match sklearn)
        auc1_check = np.mean(V1)
        auc2_check = np.mean(V2)
    
        # Calculate variances and covariance
        S_01_1 = np.var(np.mean(V1, axis=1))  # Variance across positives
        S_10_1 = np.var(np.mean(V1, axis=0))  # Variance across negatives
    
        S_01_2 = np.var(np.mean(V2, axis=1))
        S_10_2 = np.var(np.mean(V2, axis=0))
    
        # Covariances
        S_01_12 = np.cov(np.mean(V1, axis=1), np.mean(V2, axis=1))[0, 1]
        S_10_12 = np.cov(np.mean(V1, axis=0), np.mean(V2, axis=0))[0, 1]
    
        # Variance of AUC difference
        var_auc1 = (S_01_1 / n_pos) + (S_10_1 / n_neg)
        var_auc2 = (S_01_2 / n_pos) + (S_10_2 / n_neg)
        cov_aucs = (S_01_12 / n_pos) + (S_10_12 / n_neg)
    
        var_diff = var_auc1 + var_auc2 - 2 * cov_aucs
    
        # Z-statistic
        if var_diff > 0:
                z_stat = (auc1 - auc2) / np.sqrt(var_diff)
                p_value = 2 * (1 - norm.cdf(abs(z_stat)))
        else:
                z_stat = np.nan
                p_value = np.nan
                warnings.warn("Variance of AUC difference is non-positive")
            
        print(f"z-statistic: {z_stat:.3f}")
        print(f"P-value: {p_value:.4e}")
    
        if p_value < 0.05:
                winner = "Model 1" if auc1 > auc2 else "Model 2"
                print(f"   → {winner} has significantly better discrimination")
        else:
                print(f"   → No significant difference in discrimination")
            
        return {
                'test': 'DeLong',
                'statistic': z_stat,
                'pvalue': p_value,
                'conf_int': f'AUC1={auc1:.4f}, AUC2={auc2:.4f}',
                'auc1': auc1,
                'auc2': auc2
        }
        

def net_reclassification_index(data: pd.DataFrame, prs1_name: str, 
                                                                prs2_name: str, threshold: int, 
                                                                phenotype_col: str = 'PHENOTYPE') -> Dict:
        """
        Calculate Net Reclassification Index (NRI) and Integrated Discrimination 
        Improvement (IDI).
        
        NRI measures how well the new model (prs2) reclassifies individuals compared 
        to the old model (prs1). Positive NRI means improvement.
        """
        print(f"\n{'='*70}")
        print("Net Reclassification Index (NRI)")
        print(f"{'='*70}")
    
        # Separate cases and controls (assuming PHENOTYPE: 2=case, 1=control)
        cases = data[data[phenotype_col] == 2]
        controls = data[data[phenotype_col] == 1]
    
        # Reclassification for cases (events)
        prs1_high_cases = cases[f'bin_{prs1_name}'] > threshold
        prs2_high_cases = cases[f'bin_{prs2_name}'] > threshold
    
        # Count reclassifications
        up_cases = np.sum(~prs1_high_cases & prs2_high_cases)  # Low→High (good)
        down_cases = np.sum(prs1_high_cases & ~prs2_high_cases)  # High→Low (bad)
    
        # Reclassification for controls (non-events)
        prs1_high_controls = controls[f'bin_{prs1_name}'] > threshold
        prs2_high_controls = controls[f'bin_{prs2_name}'] > threshold
    
        up_controls = np.sum(~prs1_high_controls & prs2_high_controls)  # Low→High (bad)
        down_controls = np.sum(prs1_high_controls & ~prs2_high_controls)  # High→Low (good)
    
        # NRI calculation
        n_cases = len(cases)
        n_controls = len(controls)
    
        nri_events = (up_cases - down_cases) / n_cases
        nri_non_events = (down_controls - up_controls) / n_controls
        nri = nri_events + nri_non_events
    
        # Standard error and p-value
        se_events = np.sqrt((up_cases + down_cases) / n_cases) / n_cases
        se_non_events = np.sqrt((up_controls + down_controls) / n_controls) / n_controls
        se_nri = np.sqrt(se_events**2 + se_non_events**2)
    
        z_stat = nri / se_nri if se_nri > 0 else np.nan
        p_value = 2 * (1 - norm.cdf(abs(z_stat))) if not np.isnan(z_stat) else np.nan
    
        # 95% CI
        ci_low = nri - 1.96 * se_nri
        ci_high = nri + 1.96 * se_nri
    
        print(f"NRI (Events): {nri_events:.3f}")
        print(f"  Cases: {up_cases} moved up, {down_cases} moved down")
        print(f"NRI (Non-events): {nri_non_events:.3f}")
        print(f"  Controls: {down_controls} moved down, {up_controls} moved up")
        print(f"\nTotal NRI: {nri:.3f}")
        print(f"95% CI: [{ci_low:.3f}, {ci_high:.3f}]")
        print(f"z-statistic: {z_stat:.3f}")
        print(f"P-value: {p_value:.4e}")
    
        if p_value < 0.05:
                if nri > 0:
                        print(f"   → {prs2_name} provides significant improvement in reclassification")
                else:
                        print(f"   → {prs1_name} provides better reclassification")
        else:
                print(f"   → No significant difference in reclassification")
            
        # Calculate IDI (Integrated Discrimination Improvement)
        # IDI = improvement in discrimination slope
        mean_prs1_cases = cases[f'scaled_prs_{prs1_name}'].mean()
        mean_prs1_controls = controls[f'scaled_prs_{prs1_name}'].mean()
        mean_prs2_cases = cases[f'scaled_prs_{prs2_name}'].mean()
        mean_prs2_controls = controls[f'scaled_prs_{prs2_name}'].mean()
    
        idi = (mean_prs2_cases - mean_prs2_controls) - (mean_prs1_cases - mean_prs1_controls)
    
        print(f"\nIDI (Integrated Discrimination Improvement): {idi:.4f}")
    
        return {
                'test': 'NRI',
                'statistic': nri,
                'pvalue': p_value,
                'conf_int': (ci_low, ci_high),
                'nri_events': nri_events,
                'nri_non_events': nri_non_events,
                'idi': idi
        }
        

# ============================================================================
# MAIN COMPARISON FUNCTION
# ============================================================================

def run_prs_comparison(cases: pd.DataFrame, prs1_name: str, prs2_name: str,
                                            test_threshold: int, validation_data: pd.DataFrame = None,
                                            phenotype_col: str = 'PHENOTYPE') -> pd.DataFrame:
        """
        Run all statistical comparisons for a pair of PRS scores.
        
        Parameters:
        -----------
        cases : pd.DataFrame
                DataFrame with case data (phenotype == 2)
        prs1_name : str
                Name of first PRS
        prs2_name : str
                Name of second PRS
        test_threshold : int
                Threshold for high-risk classification
        validation_data : pd.DataFrame, optional
                Full validation dataset (for NRI and DeLong)
        phenotype_col : str
                Column name for phenotype
        
        Returns:
        --------
        pd.DataFrame with test results
        """
        results = []
    
        print(f"\n{'#'*70}")
        print(f"# Comparing: {prs1_name} vs {prs2_name}")
        print(f"{'#'*70}\n")
    
        # Create high-risk indicators
        prs1_high = cases[f'bin_{prs1_name}'] > test_threshold
        prs2_high = cases[f'bin_{prs2_name}'] > test_threshold
    
        # Create contingency table (used by multiple tests)
        table, counts = create_contingency_table(prs1_high, prs2_high, cases)
    
        # Run all tests
        tests_to_run = [
                ('mcnemar', lambda: mcnemar_test(prs1_high, prs2_high, cases, 
                                                                                    prs1_name, prs2_name, test_threshold)),
                ('overlap_uniqueness', lambda: overlap_uniqueness_test(counts, len(cases))),
                ('cohen_kappa', lambda: cohen_kappa_test(prs1_high, prs2_high)),
                ('jaccard', lambda: jaccard_index(counts)),
                ('discordance', lambda: discordance_proportion(counts, len(cases))),
                ('chi_square', lambda: chi_square_independence(table)),
                ('pearson', lambda: pearson_correlation_test(cases, prs1_name, 
                                                                                                            prs2_name, test_threshold)),
                ('spearman', lambda: spearman_correlation_test(cases, prs1_name, 
                                                                                                                prs2_name, test_threshold)),
                ('ttest', lambda: paired_ttest(cases, prs1_name, prs2_name, test_threshold))
        ]
    
        # Add DeLong and NRI if full validation data provided
        if validation_data is not None:
                # DeLong test
                try:
                        y_true = (validation_data[phenotype_col] == 2).astype(int)
                        scores1 = validation_data[f'scaled_prs_{prs1_name}'].values
                        scores2 = validation_data[f'scaled_prs_{prs2_name}'].values
                        tests_to_run.append(('delong', lambda: delong_test(y_true, scores1, scores2)))
                except Exception as e:
                        print(f"Warning: Could not run DeLong test: {e}")
                    
                # NRI
                try:
                        tests_to_run.append(('nri', lambda: net_reclassification_index(
                                validation_data, prs1_name, prs2_name, test_threshold, phenotype_col)))
                except Exception as e:
                        print(f"Warning: Could not run NRI: {e}")
                    
        # Execute all tests
        for test_name, test_func in tests_to_run:
                try:
                        result = test_func()
                        result['comparison'] = f'{prs1_name} v {prs2_name}'
                        result['threshold'] = test_threshold
                        results.append(result)
                except Exception as e:
                        print(f"Error running {test_name}: {e}")
                        continue
            
        return pd.DataFrame(results)


def calculate_mcnemar_test(validation_prs_file: str, scores_path: str):
        """
        Main function to run all PRS comparisons.
        
        Parameters:
        -----------
        validation_prs_file : str
                Path to CSV with PRS scores and phenotypes
        scores_path : str
                Path to save output results
        """
        # Load data
        validation_prs = pd.read_csv(validation_prs_file)
        cases = validation_prs[validation_prs['PHENOTYPE'] == 2]
    
        # Define test parameters
        test_threshold = 8
    
        # Define comparisons
        prs_comparisons = [
                ('main', 'epi'),
                ('main', 'epi+main'),
                ('epi', 'cardio'),
                ('epi', 'epi+main'),
                ('main', 'cardio'),
                ('cardio', 'epi+main'),
                ('all', 'combined')
        ]
    
        # Run comparisons
        all_results = []
    
        for prs1, prs2 in prs_comparisons:
                # Handle special case for 'all' vs 'combined'
                if 'all' in [prs1, prs2] or 'combined' in [prs1, prs2]:
                        # Need special handling for combined PRS
                        # (create combined high-risk indicator)
                        if prs2 == 'combined':
                                # Create combined indicator
                                cases_copy = cases.copy()
                                combined_high = (
                                        (cases['bin_main'] > test_threshold) | 
                                        (cases['bin_epi'] > test_threshold) |
                                        (cases['bin_cardio'] > test_threshold) | 
                                        (cases['bin_epi+main'] > test_threshold)
                                )
                                cases_copy['bin_combined'] = combined_high.astype(int) * 10
                            
                                # Use validation_prs for full dataset tests
                                validation_prs_copy = validation_prs.copy()
                                combined_high_full = (
                                        (validation_prs['bin_main'] > test_threshold) | 
                                        (validation_prs['bin_epi'] > test_threshold) |
                                        (validation_prs['bin_cardio'] > test_threshold) | 
                                        (validation_prs['bin_epi+main'] > test_threshold)
                                )
                                validation_prs_copy['bin_combined'] = combined_high_full.astype(int) * 10
                                validation_prs_copy['scaled_prs_combined'] = validation_prs_copy['bin_combined']
                            
                                results = run_prs_comparison(
                                        cases_copy, prs1, prs2, test_threshold, 
                                        validation_data=validation_prs_copy
                                )
                        else:
                                results = run_prs_comparison(
                                        cases, prs1, prs2, test_threshold, 
                                        validation_data=validation_prs
                                )
                else:
                        results = run_prs_comparison(
                                cases, prs1, prs2, test_threshold, 
                                validation_data=validation_prs
                        )
                    
                all_results.append(results)
            
        # Combine all results
        final_results = pd.concat(all_results, ignore_index=True)
    
        # Reorder columns
        column_order = ['comparison', 'test', 'threshold', 'statistic', 
                                        'pvalue', 'conf_int']
        # Add any extra columns that might exist
        extra_cols = [c for c in final_results.columns if c not in column_order]
        final_results = final_results[column_order + extra_cols]
    
        # Save results
        output_file = f'{scores_path}/McNemarStatsTestsAcrossPRSCalculations_refactored.csv'
        final_results.to_csv(output_file, index=False)
    
        print(f"\n{'='*70}")
        print(f"Analysis complete!")
        print(f"Results saved to: {output_file}")
        print(f"Comparisons analyzed: {final_results['comparison'].unique()}")
        print(f"{'='*70}\n")
    
        return final_results


# ============================================================================
# EXAMPLE USAGE
# ============================================================================

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="calculating performance metrics across PRS models....")
    parser.add_argument("--prs_file", help="Path to the input combined PRS file")
    parser.add_argument("--pheno_data", help="path to output directory")
    
    
    # Prefer command-line input if provided; fallback to env var
    pheno_data = args.pheno_data or os.environ.get("PHENO_DATA")
    print(f"[PYTHON] Reading from: {pheno_data}")
    
    scores_path = f'{pheno_data}/scores'
    print(f"[PYTHON] Reading from scores path : {scores_path}")
    
    prs_file = args.prs_file or os.environ.get("COMBINED_PRS_FILE")
    print(f"combined prs file : {prs_file}")
    
#   prs_file = '/path/to/your/combinedPRSGroups.csv'
#   scores_path = '/path/to/output/directory'
    
    results = calculate_mcnemar_test(prs_file, scores_path)
    
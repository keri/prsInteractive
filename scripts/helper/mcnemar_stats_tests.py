#!/usr/prs/env python3


import pandas as pd
import numpy as np
from statsmodels.stats.contingency_tables import mcnemar
from scipy.stats import pearsonr,spearmanr,binom, binomtest, hypergeom, norm, ttest_rel,chi2_contingency
from sklearn.metrics import cohen_kappa_score
from statsmodels.stats.proportion import proportion_confint
#from MLstatkit import Delong_test

#def calculate_nri_pvalue(delta_p_event, delta_p_non_event, n_event, n_non_event):
#   # Standard error of NRI
#   se_nri = np.sqrt((delta_p_event * (1 - delta_p_event)) / n_event + (delta_p_non_event * (1 - delta_p_non_event)) / n_non_event)
#   
#   # Calculate NRI
#   nri = delta_p_event + delta_p_non_event
#   
#   # Z-statistic
#   z_stat = nri / se_nri
#   
#   # P-value (two-tailed)
#   p_value = 2 * (1 - norm.cdf(abs(z_stat)))
#   return nri, z_stat, p_value
    

def pearson_confidence_interval(corr, n, alpha=0.05):
    # Fisher Z-transformation
    fisher_z = np.arctanh(corr)
    # Standard error
    se = 1 / np.sqrt(n - 3)
    # Z critical value
    z_crit = norm.ppf(1 - alpha / 2)
    # Confidence interval in Fisher Z-space
    z_lower = fisher_z - z_crit * se
    z_upper = fisher_z + z_crit * se
    # Transform back to correlation space
    lower = np.tanh(z_lower)
    upper = np.tanh(z_upper)
    return lower, upper


def calculate_mcnemar_test(validationPRSFile,scoresPath):

    validationPRS = pd.read_csv(validationPRSFile)

    cases = validationPRS[validationPRS['PHENOTYPE'] == 2]
    controls = validationPRS[validationPRS['PHENOTYPE'] == 1]
    
    yTrue = validationPRS[['IID','PHENOTYPE']]

    results = pd.DataFrame()

    test_threshold = 8

    for risk_tuple in [('main','epi'),('epi','cardio'),('main','cardio'),('all','combined')]:
        resultsTemp = pd.DataFrame()
        
        test = []
        threshold= []
        stat = []
        n_list = []
        pvalue = []
        conf_int = []
        comparison = []
        
        
        prs1 = risk_tuple[0]
        prs2 = risk_tuple[1]
        
        print(f'########## TEST {prs1} V {prs2} ###################################')

        
    #   mcnemar_pvalue = .06

        
        
        if 'all' in [prs1,prs2]:
            print('')
            print('########## TEST ALL FEATURES MODELLED TOGETHER V SEPARATE PRS COMBINED ###################################')
            print('')
            print('')
#           print('######################  DeLong Test ################################')
#           print('')
#           test.append('DeLong Test')
#           threshold.append('None')
#           
#           prs1 = 'all'
#           prs2 = 'combined'
#           
##           z_score, p_value = calculate_delong(scoresPath, yTrue, prs1, prs2)
#           
#           comparison.append('all v prs_combined')
#           stat.append(z_score)
#           pvalue.append(p_value)
#           conf_int.append('NA')
#           print("P-value:", p_value)
            
            print('')
            print('######################  McNemar Test ################################')
            print('')
            test.append(f'McNemar')
            threshold.append(test_threshold)
                        
            
            all_high = cases[f'bin_all'] > test_threshold
            combined_high = (cases[f'bin_main'] > test_threshold) | (cases[f'bin_epi'] > test_threshold)  | (cases[f'bin_cardio'] > test_threshold) | (cases[f'bin_epi+main'] > test_threshold)
            prs1_high = all_high
            prs2_high = combined_high
            # McNemar's 2x2 table
            # a: Both methods identify as high risk
            a = cases[all_high & combined_high].shape[0]
            
            # b: Main OR epi identifies, but combined does NOT (discordant)
            b = cases[all_high & ~combined_high].shape[0]
            
            # c: Combined identifies, but main OR epi does NOT (discordant)  
            c = cases[~all_high & combined_high].shape[0]
            
            # d: Neither method identifies as high risk
            d = cases[~all_high & ~combined_high].shape[0]
            
            # Create contingency table for McNemar's test
            mcnemar_table = np.array([[a, b],
                                    [c, d]])
            
            # Validate table sums to total cases
            total = a + b + c + d
            assert total == len(cases)
            print(f"McNemar table error: {total} != {len(cases)}")
            
            n_total = len(cases)
            n_concordant = a + d
            n_discordant = b + c
            
            
            # Print formatted table
            print(f"\n{'='*60}")
            print(f"McNemar's Test: Separate PRS vs combined")
            print(f"Threshold: bin > {test_threshold}")
            print(f"{'='*60}")
            print(f"\n2x2 Contingency Table:")
            print(f"                         combined High  combined Low     Total")
            print(f"separate PRS High           {a:6d}         {b:6d}       {a+b:6d}")
            print(f"separate PRS Low            {c:6d}         {d:6d}       {c+d:6d}")
            print(f"Total                       {a+c:6d}       {b+d:6d}     {total:6d}")
            
            print(f"\nInterpretation:")
            print(f"  Concordant pairs:   {a+d:6d} ({100*(a+d)/total:.1f}%)")
            print(f"    - Both high:             {a:6d}")
            print(f"    - Both low:              {d:6d}")
            print(f"  Discordant pairs:   {b+c:6d} ({100*(b+c)/total:.1f}%)")
            print(f"    - Only separate PRS:     {b:6d}")
            print(f"    - Only combined PRS:     {c:6d}")
            
            # Perform McNemar test
            result = mcnemar(mcnemar_table, exact=False)  # Use exact test for small samples; set exact=False for large samples
            comparison.append('all v prs_combined')
            stat.append(result.statistic)
            pvalue.append(result.pvalue)
            conf_int.append('NA')
            print("McNemar test statistic :", result.statistic)
            print("P-value:", result.pvalue)
            
            if result.pvalue >= 0.05:
                print(f"   → Methods identify similar NUMBERS of cases")
            else:
                winner = 'all' if b > c else 'combined'
                print(f"   → {winner} identifies significantly more cases")
                
                
            print('################### P-value for uniqueness based on overlap of high risk cases ####################')
            print('')
            
            test.append('overlap uniqueness')
            threshold.append(test_threshold)
            # Summary
            total_unique = b + c
            total_high = a + b + c
            
            
            p_main = (a + b) / n_total
            p_epi = (a + c) / n_total
            expected_overlap = n_total * p_main * p_epi
            
            print(f"\n2. Is overlap LESS than expected by chance?")
            print(f"   Expected overlap: {expected_overlap:.1f}")
            print(f"   Observed overlap: {a}")
            
            p_hyper_lower = hypergeom.cdf(a, n_total, a+c, a+b)
            print(f"   P-value: {p_hyper_lower:.4e}")
            
            if p_hyper_lower < 0.05:
                print(f"   → YES: Overlap is LESS than chance")
                print(f"   → Methods identify MORE unique cases than expected")
                print(f"   → Methods are COMPLEMENTARY")
            else:
                print(f"   → NO: Overlap is not less than chance")
                
            # Test 3: Is proportion of unique cases significant?
            prop_unique = total_unique / total_high
            p_binom = binomtest(total_unique, total_high, 0.5)
            
            print(f"\n3. Is proportion of UNIQUE cases significantly high?")
            print(f"   Proportion unique: {prop_unique:.3f}")
            print(f"   Testing against 50%")
            print(f"   P-value: {p_binom.pvalue:.3f}")
            
            if prop_unique > 0.5 and p_binom.pvalue < 0.05:
                print(f"   → YES: More unique than shared")
                print(f"   → Methods provide COMPLEMENTARY information")
            elif prop_unique < 0.5 and p_binom.pvalue < 0.05:
                print(f"   → NO: More shared than unique")
                print(f"   → Methods are largely REDUNDANT")
            else:
                print(f"   → Unique vs shared is balanced ≈50/50")
                
            print(f"\n{'='*70}")
            print(f"\nOverall Conclusion:")
            
            if p_hyper_lower < 0.05 or (prop_unique > 0.6 and p_binom.pvalue < 0.05):
                print(f"  ✓ Methods identify substantially different cases")
                print(f"  ✓ They provide COMPLEMENTARY risk information")
                print(f"  ✓ Using both captures {total_unique} additional cases")
            elif result.pvalue < 0.05:
                print(f"  ⚠ One method identifies more unique cases than the other")
            else:
                print(f"  → Methods have moderate overlap")
                print(f"  → Both identify some unique and some shared cases")
                
            print(f"{'='*70}\n")
            
            comparison.append('all v prs_combined')
            stat.append(prop_unique)
            pvalue.append(f'overlap (#of non-unique cases) less than expected? {p_hyper_lower} | proportion of overlap signficant? {p_binom.pvalue}')
            conf_int.append('NA')
            
            
            
            print('')
            print('################# Cohen Kappa test ################')
            print('')
            
            test.append(f'Cohen Kappa')
            threshold.append(test_threshold)
            kappa = cohen_kappa_score(prs1_high, prs2_high)
            
            print(f"\n2. Cohen's Kappa (How much do they AGREE?):")
            print(f"   Kappa: {kappa:.3f}")
            if kappa < 0.20:
                print(f"   → Slight agreement - methods largely identify different cases")
            elif kappa < 0.40:
                print(f"   → Fair agreement")
            elif kappa < 0.60:
                print(f"   → Moderate agreement")
            else:
                print(f"   → Substantial agreement")
            
            comparison.append('all v prs_combined')
            stat.append(kappa)
            pvalue.append('NA')
            conf_int.append('NA')
            print("Kappa test statistic :", kappa)
            print("P-value: None calculated")
            
                
            # 4. Jaccard Index (Overlap)
            print('########## Jaccard overlap index ##############')
            print("")
            test.append('Jaccard Overlap')
            threshold.append(test_threshold)
            
            intersection = a
            union = a + b + c
            jaccard = intersection / union if union > 0 else 0
            
            print(f"\n3. Jaccard Index (Set OVERLAP):")
            print(f"   Jaccard: {jaccard:.3f}")
            print(f"   → {jaccard*100:.1f}% of high-risk cases identified by BOTH")
            print(f"   → {(1-jaccard)*100:.1f}% unique to one method")
            
            comparison.append('all v prs_combined')
            stat.append(jaccard)
            pvalue.append('NA')
            conf_int.append('NA')
            print("Jaccard test statistic :", jaccard)
            print("P-value: None calculated")
            
            print('########## Calculating Discordance #############')
            print('')
            test.append('discordance proportion')
            threshold.append(test_threshold)
            
            # 5. Discordance
            prop_discordant = n_discordant / n_total
            ci_low, ci_high = proportion_confint(n_discordant, n_total, method='wilson')
            
            print(f"\n4. Discordance Rate:")
            print(f"   Discordant: {n_discordant}/{n_total} ({100*prop_discordant:.1f}%)")
            print(f"   95% CI: [{100*ci_low:.1f}%, {100*ci_high:.1f}%]")
            print(f"   Breakdown:")
            print(f"     Only {prs1_high}: {b} ({100*b/n_discordant:.1f}% of discordant)")
            print(f"     Only {prs2_high}: {c} ({100*c/n_discordant:.1f}% of discordant)")
            
            comparison.append('all v prs_combined')
            stat.append(prop_discordant)
            pvalue.append('NA')
            conf_int.append((ci_low,ci_high))
            print("Discordant test statistic :", prop_discordant)
            print("P-value: None calculated")
            
            print('############# chi-squared test ###################')
            print('')
            test.append('chi squared')
            threshold.append(test_threshold)
            # 6. Chi-square independence
            chi2, p_chi2, dof, expected = chi2_contingency(mcnemar_table)
            
            print(f"\n5. Chi-Square Test (Are they INDEPENDENT?):")
            print(f"   Chi-square: {chi2:.4f}")
            print(f"   P-value: {p_chi2:.4e}")
            if p_chi2 < 0.05:
                print(f"   → Classifications are NOT independent")
            else:
                print(f"   → Classifications appear independent")
                
            print(f"\n{'='*70}\n")
            
            comparison.append('all v prs_combined')
            stat.append(chi2)
            pvalue.append(p_chi2)
            conf_int.append(f'dof {dof}: expected {expected}')
            print("chi-squared test statistic :", chi2)
            print(f"P-value: {p_chi2}")
            
            print('')
            print('######################  pearson correlation coefficient ################################')
            print('')
            
            test.append('pearson_r')
            threshold.append(test_threshold)
            # Identify high-risk cases
            highRiskCases = cases[
                (cases['bin_main'] > test_threshold) | 
                (cases['bin_epi'] > test_threshold) | 
                (cases['bin_cardio'] > test_threshold) | 
                (cases['bin_epi+main'] > test_threshold) | 
                (cases['bin_all'] > test_threshold)
            ].copy()
            
            # Combine separate models by averaging
            separate_models = ['main', 'epi', 'cardio', 'epi+main']
            highRiskCases['scaled_prs_combined'] = highRiskCases[
                [f'scaled_prs_{model}' for model in separate_models]
            ].mean(axis=1)
            
            # Calculate Pearson correlation
            prs1_label = 'all'
            prs2_label = 'combined'
            
            pearson_corr, pearson_pvalue = pearsonr(
                highRiskCases[f'scaled_prs_{prs1_label}'],
                highRiskCases[f'scaled_prs_{prs2_label}']
            )
            
            # R-squared
            r_squared = pearson_corr ** 2
            
            print(f'n high risk cases for {prs1_label} OR {prs2_label} ',highRiskCases.shape[0])
            n_list.append(f'n high risk cases for {prs1_label} OR {prs2_label} :'+str(highRiskCases.shape[0]))
            comparison.append(f'{prs1_label} or {prs2_label}:high risk')
            
            stat.append(pearson_corr)
            pvalue.append(pearson_pvalue)
            
            print(f"Pearson r-squared  of high risk cases for {prs1_label} v {prs2_label} : ", r_squared)
            print(f"Pearson Correlation  of high risk cases for {prs1_label} v {prs2_label} : ", pearson_corr)
            print(f"Pearson Correlation pvalue  of {prs1_label} v {prs2_label} : ", pearson_pvalue)
            
            n = highRiskCases[f'scaled_prs_{prs1_label}'].shape[0]
            conf_int_test = pearson_confidence_interval(pearson_corr, n)
            conf_int.append(conf_int_test)
            
            print("95% Confidence Interval:", conf_int_test)
            print('')
            print('######################  spearmanr correlation coefficient ################################')
            print('')
            test.append('spearman_r')
            threshold.append(test_threshold)
            
            highRiskCases['bin_combined'] = highRiskCases[
                [f'bin_{model}' for model in separate_models]
            ].mean(axis=1)
            
            
            
            # Spearman correlation
            spearmanr_corr_test, spearman_pvalue_test = spearmanr(highRiskCases[f'bin_{prs1_label}'], highRiskCases[f'bin_{prs2_label}'])
            print("Spearman Correlation:", spearmanr_corr_test)
            print("Spearman pvalue:", spearman_pvalue_test)
            print('')
            
            n_list.append(f'n high risk cases for {prs1_label} OR {prs2_label} :'+str(highRiskCases.shape[0]))
            comparison.append(f'{prs1_label} or {prs2_label}:high risk')
            
            stat.append(spearmanr_corr_test)
            pvalue.append(spearman_pvalue_test)
            conf_int.append('NA')
            
            print('')
            print('######## t test for high risk cases ##############')
            print('')
            test.append('paired_ttest')
            threshold.append(test_threshold)
            
            t_stat, p_value = ttest_rel(highRiskCases[f'scaled_prs_{prs1_label}'], highRiskCases[f'scaled_prs_{prs2_label}'])
            print("Paired t-test P-value for high risk cases :", p_value)
            print('t statistic = ',t_stat)
            print('')
            
            n_list.append(f'n high risk cases for {prs1_label} OR {prs2_label} :'+str(highRiskCases.shape[0]))
            comparison.append(f'{prs1_label} or {prs2_label}:high risk')
            
            stat.append(t_stat)
            pvalue.append(p_value)
            conf_int.append('NA')
            

                
        elif prs1 in ['epi', 'main'] and prs2 in ['epi', 'main']:
            
            print('')
            print('########## TEST EPI+MAIN FEATURES MODELLED TOGETHER V SEPARATE EPI, MAIN PRS  ###################################')
            print('')
#           print('######################  DeLong Test ################################')
#           print('')
#           test.append('DeLong Test')
#           threshold.append('None')
#           
#           prs1Delong = 'epi+main'
#           prs2Delong = 'main|epi'
#           
#           z_score, p_value = calculate_delong(scoresPath, yTrue, prs1Delong, prs2Delong)
#           
#           comparison.append(f'{prs1}or{prs2} v G+(GxG)')
#           stat.append(z_score)
#           pvalue.append(p_value)
#           conf_int.append('NA')
#           print("P-value:", p_value)
#           
#           print('')
            
            print('')

            print('######################  McNemar Test ################################')
            print('')
            # Create boolean indicators
            main_or_epi_high = (cases[f'bin_{prs1}'] > test_threshold) | (cases[f'bin_{prs2}'] > test_threshold)
            combined_high = cases[f'bin_epi+main'] > test_threshold
            
            prs1_high = main_or_epi_high
            prs2_high = combined_high
            
            # McNemar's 2x2 table
            # a: Both methods identify as high risk
            a = cases[main_or_epi_high & combined_high].shape[0]
            
            # b: Main OR epi identifies, but combined does NOT (discordant)
            b = cases[main_or_epi_high & ~combined_high].shape[0]
            
            # c: Combined identifies, but main OR epi does NOT (discordant)  
            c = cases[~main_or_epi_high & combined_high].shape[0]
            
            # d: Neither method identifies as high risk
            d = cases[~main_or_epi_high & ~combined_high].shape[0]
            
            # Create contingency table for McNemar's test
            mcnemar_table = np.array([[a, b],
                                    [c, d]])
                                
            # Validate table sums to total cases
            total = a + b + c + d
            assert total == len(cases), f"McNemar table error: {total} != {len(cases)}"
            
            n_total = len(cases)
            n_concordant = a + d
            n_discordant = b + c
            
            # Print formatted table
            print(f"\n{'='*60}")
            print(f"McNemar's Test: {prs1}|{prs2} vs combined")
            print(f"Threshold: bin > {test_threshold}")
            print(f"{'='*60}")
            print(f"\n2x2 Contingency Table:")
            print(f"                        Combined High  Combined Low      Total")
            print(f"{prs1}|{prs2} High          {a:6d}         {b:6d}       {a+b:6d}")
            print(f"{prs1}|{prs2} Low           {c:6d}         {d:6d}       {c+d:6d}")
            print(f"Total                       {a+c:6d}       {b+d:6d}     {total:6d}")
            
            print(f"\nInterpretation:")
            print(f"  Concordant pairs:   {a+d:6d} ({100*(a+d)/total:.1f}%)")
            print(f"    - Both high:              {a:6d}")
            print(f"    - Both low:               {d:6d}")
            print(f"  Discordant pairs:   {b+c:6d} ({100*(b+c)/total:.1f}%)")
            print(f"    - Only {prs1}|{prs2}:     {b:6d}")
            print(f"    - Only epi+main:          {c:6d}")
            
            
            # Perform McNemar test
            result = mcnemar(mcnemar_table, exact=False)  # Use exact test for small samples; set exact=False for large samples
            test.append('McNemar')
            comparison.append(f'{prs1}or{prs2} v G+(GxG)')
            stat.append(result.statistic)
            pvalue.append(result.pvalue)
            conf_int.append('NA')
            threshold.append(test_threshold)
            print("McNemar test statistic :", result.statistic)
            print("P-value:", result.pvalue)
            
            print('################### P-value for uniqueness based on overlap of high risk cases ####################')
            print('')
            
            test.append('overlap uniqueness')
            threshold.append(test_threshold)
            # Summary
            total_unique = b + c
            total_high = a + b + c
            
            
            p_main = (a + b) / n_total
            p_epi = (a + c) / n_total
            expected_overlap = n_total * p_main * p_epi
            
            print(f"\n2. Is overlap LESS than expected by chance?")
            print(f"   Expected overlap: {expected_overlap:.1f}")
            print(f"   Observed overlap: {a}")
            
            p_hyper_lower = hypergeom.cdf(a, n_total, a+c, a+b)
            print(f"   P-value: {p_hyper_lower:.4e}")
            
            if p_hyper_lower < 0.05:
                print(f"   → YES: Overlap is LESS than chance")
                print(f"   → Methods identify MORE unique cases than expected")
                print(f"   → Methods are COMPLEMENTARY")
            else:
                print(f"   → NO: Overlap is not less than chance")
                
            # Test 3: Is proportion of unique cases significant?
            prop_unique = total_unique / total_high
            p_binom = binomtest(total_unique, total_high, 0.5)
            
            print(f"\n3. Is proportion of UNIQUE cases significantly high?")
            print(f"   Proportion unique: {prop_unique:.3f}")
            print(f"   Testing against 50%")
            print(f"   P-value: {p_binom.pvalue:.3f}")
            
            if prop_unique > 0.5 and p_binom.pvalue < 0.05:
                print(f"   → YES: More unique than shared")
                print(f"   → Methods provide COMPLEMENTARY information")
            elif prop_unique < 0.5 and p_binom.pvalue < 0.05:
                print(f"   → NO: More shared than unique")
                print(f"   → Methods are largely REDUNDANT")
            else:
                print(f"   → Unique vs shared is balanced (≈50/50)")
                
            print(f"\n{'='*70}")
            print(f"\nOverall Conclusion:")
            
            if p_hyper_lower < 0.05 or (prop_unique > 0.6 and p_binom.pvalue < 0.05):
                print(f"  ✓ Methods identify substantially different cases")
                print(f"  ✓ They provide COMPLEMENTARY risk information")
                print(f"  ✓ Using both captures {total_unique} additional cases")
            elif result.pvalue < 0.05:
                print(f"  ⚠ One method identifies more unique cases than the other")
            else:
                print(f"  → Methods have moderate overlap")
                print(f"  → Both identify some unique and some shared cases")
                
            print(f"{'='*70}\n")
            
            comparison.append(f'{prs1}or{prs2} v G+(GxG)')
            stat.append(prop_unique)
            pvalue.append(f'overlap (#of non-unique cases) less than expected? {p_hyper_lower} | proportion of overlap signficant? {p_binom.pvalue}')
            conf_int.append('NA')
            
            
            
            print('################# Cohen Kappa test ################')
            print('')
            
            test.append(f'Cohen Kappa')
            threshold.append(test_threshold)
            kappa = cohen_kappa_score(prs1_high, prs2_high)
            
            print(f"\n2. Cohen's Kappa (How much do they AGREE?):")
            print(f"   Kappa: {kappa:.3f}")
            if kappa < 0.20:
                print(f"   → Slight agreement - methods largely identify different cases")
            elif kappa < 0.40:
                print(f"   → Fair agreement")
            elif kappa < 0.60:
                print(f"   → Moderate agreement")
            else:
                print(f"   → Substantial agreement")
                
            comparison.append(f'{prs1}or{prs2} v G+(GxG)')
            stat.append(kappa)
            pvalue.append('NA')
            conf_int.append('NA')
            print("Kappa test statistic :", kappa)
            print("P-value: None calculated")
            
            
            # 4. Jaccard Index (Overlap)
            print('########## Jaccard overlap index ##############')
            print('')
            test.append('Jaccard Overlap')
            threshold.append(test_threshold)
            
            intersection = a
            union = a + b + c
            jaccard = intersection / union if union > 0 else 0
            
            print(f"\n3. Jaccard Index (Set OVERLAP):")
            print(f"   Jaccard: {jaccard:.3f}")
            print(f"   → {jaccard*100:.1f}% of high-risk cases identified by BOTH")
            print(f"   → {(1-jaccard)*100:.1f}% unique to one method")
            
            comparison.append(f'{prs1}or{prs2} v G+(GxG)')
            stat.append(jaccard)
            pvalue.append('NA')
            conf_int.append('NA')
            print("Jaccard test statistic :", jaccard)
            print("P-value: None calculated")
            
            print('########## Calculating Discordance #############')
            print('')
            test.append('discordance proportion')
            threshold.append(test_threshold)
            
            # 5. Discordance
            prop_discordant = n_discordant / n_total
            ci_low, ci_high = proportion_confint(n_discordant, n_total, method='wilson')
            
            print(f"\n4. Discordance Rate:")
            print(f"   Discordant: {n_discordant}/{n_total} ({100*prop_discordant:.1f}%)")
            print(f"   95% CI: [{100*ci_low:.1f}%, {100*ci_high:.1f}%]")
            print(f"   Breakdown:")
            print(f"     Only {prs1_high}: {b} ({100*b/n_discordant:.1f}% of discordant)")
            print(f"     Only {prs2_high}: {c} ({100*c/n_discordant:.1f}% of discordant)")
            
            comparison.append(f'{prs1}or{prs2} v G+(GxG)')
            stat.append(prop_discordant)
            pvalue.append('NA')
            conf_int.append((ci_low,ci_high))
            print("Discordant test statistic :", prop_discordant)
            print("P-value: None calculated")
            
            print('############# chi-squared test ###################')
            print('')
            test.append('chi squared')
            threshold.append(test_threshold)
            # 6. Chi-square independence
            chi2, p_chi2, dof, expected = chi2_contingency(mcnemar_table)
            
            print(f"\n5. Chi-Square Test (Are they INDEPENDENT?):")
            print(f"   Chi-square: {chi2:.4f}")
            print(f"   P-value: {p_chi2:.4e}")
            if p_chi2 < 0.05:
                print(f"   → Classifications are NOT independent")
            else:
                print(f"   → Classifications appear independent")
                
            print(f"\n{'='*70}\n")
            
            comparison.append(f'{prs1}or{prs2} v G+(GxG)')
            stat.append(chi2)
            pvalue.append(p_chi2)
            conf_int.append(f'dof {dof}: expected {expected}')
            print("chi-squared test statistic :", chi2)
            print(f"P-value: {p_chi2}")
            
            

            print('')
            print('######################  pearson correlation coefficient ################################')
            print('')
            
            test.append('pearson_r')
            threshold.append(test_threshold)
            
            highRiskCases = cases[(cases[f'bin_{prs1}'] > test_threshold) | (cases[f'bin_{prs2}'] > test_threshold) | (cases[f'bin_epi+main'] > test_threshold)]
            
            separate_models = ['main', 'epi']
            highRiskCases['scaled_prs_combined'] = highRiskCases[
                [f'scaled_prs_{model}' for model in separate_models]
            ].mean(axis=1)
            # Calculate Pearson correlation
            prs1_label = 'epi+main'
            prs2_label = 'epi|main'
            
            pearson_corr, pearson_pvalue = pearsonr(
                highRiskCases[f'scaled_prs_{prs1_label}'],
                highRiskCases['scaled_prs_combined']
                
            )
            
            # R-squared
            r_squared = pearson_corr ** 2
            
            print(f'n high risk cases for {prs1_label} OR {prs2_label} ',highRiskCases.shape[0])
            n_list.append(f'n high risk cases for {prs1_label} OR {prs2_label} :'+str(highRiskCases.shape[0]))
            comparison.append(f'{prs1_label} or {prs2_label}:high risk')
            
            stat.append(pearson_corr)
            pvalue.append(pearson_pvalue)
            
            print(f"Pearson r-squared  of high risk cases for {prs1_label} v {prs2_label} : ", r_squared)
            print(f"Pearson Correlation  of high risk cases for {prs1_label} v {prs2_label} : ", pearson_corr)
            print(f"Pearson Correlation pvalue  of {prs1_label} v {prs2_label} : ", pearson_pvalue)
            
            n = highRiskCases[f'scaled_prs_{prs1_label}'].shape[0]
            conf_int_test = pearson_confidence_interval(pearson_corr, n)
            conf_int.append(conf_int_test)
            
            print("95% Confidence Interval:", conf_int_test)
            
            print('')
            print('######################  spearmanr correlation coefficient ################################')
            print('')
            test.append('spearman_r')
            threshold.append(test_threshold)
            
            highRiskCases['bin_combined'] = highRiskCases[
                [f'bin_{model}' for model in separate_models]
            ].mean(axis=1)
            
            
            # Spearman correlation
            spearmanr_corr_test, spearman_pvalue_test = spearmanr(highRiskCases[f'bin_{prs1_label}'], highRiskCases['bin_combined'])
            print("Spearman Correlation:", spearmanr_corr_test)
            print("Spearman pvalue:", spearman_pvalue_test)
            print('')
            
            n_list.append(f'n high risk cases for {prs1_label} OR {prs2_label} :'+str(highRiskCases.shape[0]))
            comparison.append(f'{prs1_label} or {prs2_label}:high risk')
            
            stat.append(spearmanr_corr_test)
            pvalue.append(spearman_pvalue_test)
            conf_int.append('NA')
            print('')
            print('######## t test for high risk cases dataset ##############')
            print('')
            
            test.append('paired_ttest')
            threshold.append(test_threshold)
            
            t_stat, p_value = ttest_rel(highRiskCases[f'scaled_prs_{prs1_label}'], highRiskCases['scaled_prs_combined'])
            print("Paired t-test P-value for high risk cases :", p_value)
            print('t statistic = ',t_stat)
            print('')
            
            n_list.append(f'n high risk cases for {prs1_label} OR {prs2_label} :'+str(highRiskCases.shape[0]))
            comparison.append(f'{prs1_label} or {prs2_label}:high risk')
            
            stat.append(t_stat)
            pvalue.append(p_value)
            conf_int.append('NA')
            
            ##### COMPARE MAIN AND EPI SEPARATELY  ######
            print('')
            print(f'#################### TEST FOR : {prs1} v {prs2}  #####################################')
            print('')
            print('')
#           print('######################  DeLong Test ################################')
#           print('')
#           test.append('DeLong Test')
#           threshold.append('None')
#           
#           prs1 = 'epi'
#           prs2 = 'main'
#           
#           z_score, p_value = calculate_delong(scoresPath, yTrue, prs1, prs2)
#           
#           comparison.append(f'{prs1} v {prs2}')
#           stat.append(z_score)
#           pvalue.append(p_value)
#           conf_int.append('NA')
#           print("P-value:", p_value)
#           
#           print('')
            print('')

            print(f'######################  McNemar Test ################################')
            print('')
            
            test.append('McNemar')
            threshold.append(test_threshold)
            
            prs1_high = cases[f'bin_{prs1}'] > test_threshold
            prs2_high = cases[f'bin_{prs2}'] > test_threshold
            
            a = cases[prs1_high & prs2_high].shape[0]
            b = cases[prs1_high & ~prs2_high].shape[0]
            c = cases[~prs1_high & prs2_high].shape[0]
            d = cases[~prs1_high & ~prs2_high].shape[0]
            
            # Create contingency table for McNemar's test
            mcnemar_table = np.array([[a, b],
                                    [c, d]])
                                
            # Validate table sums to total cases
            total = a + b + c + d
            assert total == len(cases), f"McNemar table error: {total} != {len(cases)}"
            
            n_total = len(cases)
            n_concordant = a + d
            n_discordant = b + c
            
            # Print formatted table
            print(f"\n{'='*60}")
            print(f"McNemar's Test: {prs1} vs {prs2}")
            print(f"Threshold: bin > {test_threshold}")
            print(f"{'='*60}")
            print(f"\n2x2 Contingency Table:")
            print(f"                    {prs2} High   {prs2} Low      Total")
            print(f"{prs1} High           {a:6d}        {b:6d}       {a+b:6d}")
            print(f"{prs1} Low            {c:6d}        {d:6d}       {c+d:6d}")
            print(f"Total                 {a+c:6d}      {b+d:6d}     {total:6d}")
            
            print(f"\nInterpretation:")
            print(f"  Concordant pairs:   {a+d:6d} ({100*(a+d)/total:.1f}%)")
            print(f"    - Both high:      {a:6d}")
            print(f"    - Both low:       {d:6d}")
            print(f"  Discordant pairs:   {b+c:6d} ({100*(b+c)/total:.1f}%)")
            print(f"    - Only {prs1}:    {b:6d}")
            print(f"    - Only {prs2}:    {c:6d}")
            
            
            comparison.append(f'{prs1} v {prs2}')
            
            # Perform McNemar test
            result = mcnemar(mcnemar_table, exact=False)  # Use exact test for small samples; set exact=False for large samples
            stat.append(result.statistic)
            pvalue.append(result.pvalue)
            conf_int.append('NA')
            print("McNemar test statistic :", result.statistic)
            print("P-value:", result.pvalue)
            
            print('################### P-value for uniqueness based on overlap of high risk cases ####################')
            print('')
            
            test.append('overlap uniqueness')
            threshold.append(test_threshold)
            # Summary
            total_unique = b + c
            total_high = a + b + c
            
            
            p_main = (a + b) / n_total
            p_epi = (a + c) / n_total
            expected_overlap = n_total * p_main * p_epi
            
            print(f"\n2. Is overlap LESS than expected by chance?")
            print(f"   Expected overlap: {expected_overlap:.1f}")
            print(f"   Observed overlap: {a}")
            
            p_hyper_lower = hypergeom.cdf(a, n_total, a+c, a+b)
            print(f"   P-value: {p_hyper_lower:.4e}")
            
            if p_hyper_lower < 0.05:
                print(f"   → YES: Overlap is LESS than chance")
                print(f"   → Methods identify MORE unique cases than expected")
                print(f"   → Methods are COMPLEMENTARY")
            else:
                print(f"   → NO: Overlap is not less than chance")
                
            # Test 3: Is proportion of unique cases significant?
            prop_unique = total_unique / total_high
            p_binom = binomtest(total_unique, total_high, 0.5)
            
            print(f"\n3. Is proportion of UNIQUE cases significantly high?")
            print(f"   Proportion unique: {prop_unique:.3f}")
            print(f"   Testing against 50%")
            print(f"   P-value: {p_binom.pvalue:.3f}")
            
            if prop_unique > 0.5 and p_binom.pvalue < 0.05:
                print(f"   → YES: More unique than shared")
                print(f"   → Methods provide COMPLEMENTARY information")
            elif prop_unique < 0.5 and p_binom.pvalue < 0.05:
                print(f"   → NO: More shared than unique")
                print(f"   → Methods are largely REDUNDANT")
            else:
                print(f"   → Unique vs shared is balanced (≈50/50)")
                
            print(f"\n{'='*70}")
            print(f"\nOverall Conclusion:")
            
            if p_hyper_lower < 0.05 or (prop_unique > 0.6 and p_binom.pvalue < 0.05):
                print(f"  ✓ Methods identify substantially different cases")
                print(f"  ✓ They provide COMPLEMENTARY risk information")
                print(f"  ✓ Using both captures {total_unique} additional cases")
            elif result.pvalue < 0.05:
                print(f"  ⚠ One method identifies more unique cases than the other")
            else:
                print(f"  → Methods have moderate overlap")
                print(f"  → Both identify some unique and some shared cases")
                
            print(f"{'='*70}\n")
            
            comparison.append(f'{prs1} v {prs2}')
            stat.append(prop_unique)
            pvalue.append(f'overlap (#of non-unique cases) less than expected? {p_hyper_lower} | proportion of overlap signficant? {p_binom.pvalue}')
            conf_int.append('NA')
            
            print('################# Cohen Kappa test ################')
            print('')
            
            test.append(f'Cohen Kappa')
            threshold.append(test_threshold)
            kappa = cohen_kappa_score(prs1_high, prs2_high)
            
            print(f"\n2. Cohen's Kappa (How much do they AGREE?):")
            print(f"   Kappa: {kappa:.3f}")
            if kappa < 0.20:
                print(f"   → Slight agreement - methods largely identify different cases")
            elif kappa < 0.40:
                print(f"   → Fair agreement")
            elif kappa < 0.60:
                print(f"   → Moderate agreement")
            else:
                print(f"   → Substantial agreement")
                
            comparison.append(f'{prs1} v {prs2}')
            stat.append(kappa)
            pvalue.append('NA')
            conf_int.append('NA')
            print("Kappa test statistic :", kappa)
            print("P-value: None calculated")
            
            
            # 4. Jaccard Index (Overlap)
            print('########## Jaccard overlap index ##############')
            print('')
            test.append('Jaccard Overlap')
            threshold.append(test_threshold)
            
            intersection = a
            union = a + b + c
            jaccard = intersection / union if union > 0 else 0
            
            print(f"\n3. Jaccard Index (Set OVERLAP):")
            print(f"   Jaccard: {jaccard:.3f}")
            print(f"   → {jaccard*100:.1f}% of high-risk cases identified by BOTH")
            print(f"   → {(1-jaccard)*100:.1f}% unique to one method")
            
            comparison.append(f'{prs1} v {prs2}')
            stat.append(jaccard)
            pvalue.append('NA')
            conf_int.append('NA')
            print("Jaccard test statistic :", jaccard)
            print("P-value: None calculated")
            
            print('########## Calculating Discordance #############')
            print('')
            test.append('discordance proportion')
            threshold.append(test_threshold)
            
            # 5. Discordance
            prop_discordant = n_discordant / n_total
            ci_low, ci_high = proportion_confint(n_discordant, n_total, method='wilson')
            
            print(f"\n4. Discordance Rate:")
            print(f"   Discordant: {n_discordant}/{n_total} ({100*prop_discordant:.1f}%)")
            print(f"   95% CI: [{100*ci_low:.1f}%, {100*ci_high:.1f}%]")
            print(f"   Breakdown:")
            print(f"     Only {prs1_high}: {b} ({100*b/n_discordant:.1f}% of discordant)")
            print(f"     Only {prs2_high}: {c} ({100*c/n_discordant:.1f}% of discordant)")
            
            comparison.append(f'{prs1} v {prs2}')
            stat.append(prop_discordant)
            pvalue.append('NA')
            conf_int.append((ci_low,ci_high))
            print("Discordant test statistic :", prop_discordant)
            print("P-value: None calculated")
            
            print('############# chi-squared test ###################')
            print('')
            test.append('chi squared')
            threshold.append(test_threshold)
            # 6. Chi-square independence
            chi2, p_chi2, dof, expected = chi2_contingency(mcnemar_table)
            
            print(f"\n5. Chi-Square Test (Are they INDEPENDENT?):")
            print(f"   Chi-square: {chi2:.4f}")
            print(f"   P-value: {p_chi2:.4e}")
            if p_chi2 < 0.05:
                print(f"   → Classifications are NOT independent")
            else:
                print(f"   → Classifications appear independent")
                
            print(f"\n{'='*70}\n")
            
            comparison.append(f'{prs1} v {prs2}')
            stat.append(chi2)
            pvalue.append(p_chi2)
            conf_int.append(f'dof {dof}: expected {expected}')
            print("chi-squared test statistic :", chi2)
            print(f"P-value: {p_chi2}")
            
            print('')
            print('')
            print('######################  pearson correlation coefficient ################################')
            print('')
            
            test.append('pearson_r')
            threshold.append(test_threshold)
            highRiskCases = cases[(cases[f'bin_{prs1}'] > test_threshold) | (cases[f'bin_{prs2}'] > test_threshold)]
            
            # Calculate Pearson correlation
            prs1_label = prs1
            prs2_label = prs2
            
            pearson_corr, pearson_pvalue = pearsonr(
                highRiskCases[f'scaled_prs_{prs1}'],
                highRiskCases[f'scaled_prs_{prs2}']
                
            )
            
            # R-squared
            r_squared = pearson_corr ** 2
            
            print(f'n high risk cases for {prs1_label} OR {prs2_label} ',highRiskCases.shape[0])
            n_list.append(f'n high risk cases for {prs1_label} OR {prs2_label} :'+str(highRiskCases.shape[0]))
            comparison.append(f'{prs1_label} or {prs2_label}:high risk')
            
            stat.append(pearson_corr)
            pvalue.append(pearson_pvalue)
            
            print(f"Pearson r-squared  of high risk cases for {prs1_label} v {prs2_label} : ", r_squared)
            print(f"Pearson Correlation  of high risk cases for {prs1_label} v {prs2_label} : ", pearson_corr)
            print(f"Pearson Correlation pvalue  of {prs1_label} v {prs2_label} : ", pearson_pvalue)
            
            n = highRiskCases[f'scaled_prs_{prs1}'].shape[0]
            conf_int_test = pearson_confidence_interval(pearson_corr, n)
            conf_int.append(conf_int_test)
            
            print("95% Confidence Interval:", conf_int_test)
            print('')
            print('######################  spearmanr correlation coefficient ################################')
            print('')
            test.append('spearman_r')
            threshold.append(test_threshold)
            
            # Spearman correlation
            spearmanr_corr_test, spearman_pvalue_test = spearmanr(highRiskCases[f'bin_{prs1_label}'], highRiskCases[f'bin_{prs2_label}'])
            print("Spearman Correlation:", spearmanr_corr_test)
            print("Spearman pvalue:", spearman_pvalue_test)
            print('')
            
            n_list.append(f'n high risk cases for {prs1_label} OR {prs2} :'+str(highRiskCases.shape[0]))
            comparison.append(f'{prs1_label} or {prs2_label}:high risk')
            
            stat.append(spearmanr_corr_test)
            pvalue.append(spearman_pvalue_test)
            conf_int.append('NA')
            print('')
            print('######## t test for high risk cases dataset ##############')
            print('')
        
            test.append('paired_ttest')
            threshold.append(test_threshold)
            
            t_stat, p_value = ttest_rel(highRiskCases[f'scaled_prs_{prs1_label}'], highRiskCases[f'scaled_prs_{prs2}'])
            print("Paired t-test P-value for high risk cases :", p_value)
            print('t statistic = ',t_stat)
            print('')
            
            n_list.append(f'n high risk cases for {prs1_label} OR {prs2_label} :'+str(highRiskCases.shape[0]))
            comparison.append(f'{prs1_label} or {prs2_label}:high risk')
            
            stat.append(t_stat)
            pvalue.append(p_value)
            conf_int.append('NA')
            
        else:
            print('')
            print(f'#################### TEST FOR : {prs1} v {prs2}  #####################################')
            print('')
            print('')
#           print('######################  DeLong Test ################################')
#           print('')
#           test.append('DeLong Test')
#           threshold.append('None')
#           
#           z_score, p_value = calculate_delong(scoresPath, yTrue, prs1, prs2)
#           
#           comparison.append(f'{prs1} v {prs2}')
#           stat.append(z_score)
#           pvalue.append(p_value)
#           conf_int.append('NA')
#           print("P-value:", p_value)
#           
#           print('')
            print('')
            print(f'######################  McNemar Test ################################')
            print('')
            
            prs1_high = cases[f'bin_{prs1}'] > test_threshold
            prs2_high = cases[f'bin_{prs2}'] > test_threshold
            
            a = cases[prs1_high & prs2_high].shape[0]
            b = cases[prs1_high & ~prs2_high].shape[0]
            c = cases[~prs1_high & prs2_high].shape[0]
            d = cases[~prs1_high & ~prs2_high].shape[0]
            
            # Create contingency table for McNemar's test
            mcnemar_table = np.array([[a, b],
                                    [c, d]])
            
            # Validate table sums to total cases
            total = a + b + c + d
            assert total == len(cases), f"McNemar table error: {total} != {len(cases)}"
            
            
            n_total = len(cases)
            n_concordant = a + d
            n_discordant = b + c
            
            # Print formatted table
            print(f"\n{'='*60}")
            print(f"McNemar's Test: {prs1} vs {prs2}")
            print(f"Threshold: bin > {test_threshold}")
            print(f"{'='*60}")
            print(f"\n2x2 Contingency Table:")
            print(f"                    {prs2} High    {prs2} Low     Total")
            print(f"{prs1} High          {a:6d}        {b:6d}       {a+b:6d}")
            print(f"{prs1} Low           {c:6d}        {d:6d}       {c+d:6d}")
            print(f"Total              {a+c:6d}      {b+d:6d}       {total:6d}")
            
            print(f"\nInterpretation:")
            print(f"  Concordant pairs:   {a+d:6d} ({100*(a+d)/total:.1f}%)")
            print(f"    - Both high:      {a:6d}")
            print(f"    - Both low:       {d:6d}")
            print(f"  Discordant pairs:   {b+c:6d} ({100*(b+c)/total:.1f}%)")
            print(f"    - Only {prs1}:     {b:6d}")
            print(f"    - Only {prs2}:     {c:6d}")
            
        
            result = mcnemar(mcnemar_table, exact=False)  # Use exact test for small samples; set exact=False for large samples
            test.append('McNemar')
            comparison.append(f'{prs1} v {prs2}')
            stat.append(result.statistic)
            pvalue.append(result.pvalue)
            conf_int.append('NA')
            threshold.append(test_threshold)
            print("McNemar test statistic :", result.statistic)
            print("P-value:", result.pvalue)
            
            
            
            
            print('################### P-value for uniqueness based on overlap of high risk cases ####################')
            print('')
            
            test.append('overlap uniqueness')
            threshold.append(test_threshold)
            # Summary
            total_unique = b + c
            total_high = a + b + c
            
            
            p_main = (a + b) / n_total
            p_epi = (a + c) / n_total
            expected_overlap = n_total * p_main * p_epi
            
            print(f"\n2. Is overlap LESS than expected by chance?")
            print(f"   Expected overlap: {expected_overlap:.1f}")
            print(f"   Observed overlap: {a}")
            
            p_hyper_lower = hypergeom.cdf(a, n_total, a+c, a+b)
            print(f"   P-value: {p_hyper_lower:.4e}")
            
            if p_hyper_lower < 0.05:
                print(f"   → YES: Overlap is LESS than chance")
                print(f"   → Methods identify MORE unique cases than expected")
                print(f"   → Methods are COMPLEMENTARY")
            else:
                print(f"   → NO: Overlap is not less than chance")
                
            # Test 3: Is proportion of unique cases significant?
            prop_unique = total_unique / total_high
            p_binom = binomtest(total_unique, total_high, 0.5)
            
            print(f"\n3. Is proportion of UNIQUE cases significantly high?")
            print(f"   Proportion unique: {prop_unique:.3f}")
            print(f"   Testing against 50%")
            print(f"   P-value: {p_binom.pvalue:.3f}")
            
            if prop_unique > 0.5 and p_binom.pvalue < 0.05:
                print(f"   → YES: More unique than shared")
                print(f"   → Methods provide COMPLEMENTARY information")
            elif prop_unique < 0.5 and p_binom.pvalue < 0.05:
                print(f"   → NO: More shared than unique")
                print(f"   → Methods are largely REDUNDANT")
            else:
                print(f"   → Unique vs shared is balanced (≈50/50)")
                
            print(f"\n{'='*70}")
            print(f"\nOverall Conclusion:")
            
            if p_hyper_lower < 0.05 or (prop_unique > 0.6 and p_binom.pvalue < 0.05):
                print(f"  ✓ Methods identify substantially different cases")
                print(f"  ✓ They provide COMPLEMENTARY risk information")
                print(f"  ✓ Using both captures {total_unique} additional cases")
            elif result.pvalue < 0.05:
                print(f"  ⚠ One method identifies more unique cases than the other")
            else:
                print(f"  → Methods have moderate overlap")
                print(f"  → Both identify some unique and some shared cases")
                
            print(f"{'='*70}\n_total")
            
            comparison.append(f'{prs1} v {prs2}')
            stat.append(prop_unique)
            pvalue.append(f'overlap (#of non-unique cases) less than expected? {p_hyper_lower} | proportion of overlap signficant? {p_binom.pvalue}')
            conf_int.append('NA')

            
            
            print('################# Cohen Kappa test ################')
            print('')
            
            test.append(f'Cohen Kappa')
            threshold.append(test_threshold)
            kappa = cohen_kappa_score(prs1_high, prs2_high)
            
            print(f"\n2. Cohen's Kappa (How much do they AGREE?):")
            print(f"   Kappa: {kappa:.3f}")
            if kappa < 0.20:
                print(f"   → Slight agreement - methods largely identify different cases")
            elif kappa < 0.40:
                print(f"   → Fair agreement")
            elif kappa < 0.60:
                print(f"   → Moderate agreement")
            else:
                print(f"   → Substantial agreement")
                
            comparison.append(f'{prs1} v {prs2}')
            stat.append(kappa)
            pvalue.append('NA')
            conf_int.append('NA')
            print("Kappa test statistic :", kappa)
            print("P-value: None calculated")
            
            
            # 4. Jaccard Index (Overlap)
            print('########## Jaccard overlap index ##############')
            print('')
            test.append('Jaccard Overlap')
            threshold.append(test_threshold)
            
            intersection = a
            union = a + b + c
            jaccard = intersection / union if union > 0 else 0
            
            print(f"\n3. Jaccard Index (Set OVERLAP):")
            print(f"   Jaccard: {jaccard:.3f}")
            print(f"   → {jaccard*100:.1f}% of high-risk cases identified by BOTH")
            print(f"   → {(1-jaccard)*100:.1f}% unique to one method")
            
            comparison.append(f'{prs1} v {prs2}')
            stat.append(jaccard)
            pvalue.append('NA')
            conf_int.append('NA')
            print("Jaccard test statistic :", jaccard)
            print("P-value: None calculated")
            
            print('########## Calculating Discordance #############')
            print('')
            test.append('discordance proportion')
            threshold.append(test_threshold)
            
            # 5. Discordance
            prop_discordant = n_discordant / n_total
            ci_low, ci_high = proportion_confint(n_discordant, n_total, method='wilson')
            
            print(f"\n4. Discordance Rate:")
            print(f"   Discordant: {n_discordant}/{n_total} ({100*prop_discordant:.1f}%)")
            print(f"   95% CI: [{100*ci_low:.1f}%, {100*ci_high:.1f}%]")
            print(f"   Breakdown:")
            print(f"     Only {prs1_high}: {b} ({100*b/n_discordant:.1f}% of discordant)")
            print(f"     Only {prs2_high}: {c} ({100*c/n_discordant:.1f}% of discordant)")
            
            comparison.append(f'{prs1} v {prs2}')
            stat.append(prop_discordant)
            pvalue.append('NA')
            conf_int.append((ci_low,ci_high))
            print("Discordant test statistic :", prop_discordant)
            print("P-value: None calculated")
            
            print('############# chi-squared test ###################')
            print('')
            test.append('chi squared')
            threshold.append(test_threshold)
            # 6. Chi-square independence
            chi2, p_chi2, dof, expected = chi2_contingency(mcnemar_table)
            
            print(f"\n5. Chi-Square Test (Are they INDEPENDENT?):")
            print(f"   Chi-square: {chi2:.4f}")
            print(f"   P-value: {p_chi2:.4e}")
            if p_chi2 < 0.05:
                print(f"   → Classifications are NOT independent")
            else:
                print(f"   → Classifications appear independent")
                
            print(f"\n{'='*70}\n")
            
            comparison.append(f'{prs1} v {prs2}')
            stat.append(chi2)
            pvalue.append(p_chi2)
            conf_int.append(f'dof {dof}: expected {expected}')
            print("chi-squared test statistic :", chi2)
            print(f"P-value: {p_chi2}")
            
            
            print('')
            print('')
            print('######################  pearson correlation coefficient ################################')
            print('')
            
            test.append('pearson_r')
            threshold.append(test_threshold)
            highRiskCases = cases[(cases[f'bin_{prs1}'] > test_threshold) | (cases[f'bin_{prs2}'] > test_threshold)]
            
            # Calculate Pearson correlation
            prs1_label = prs1
            prs2_label = prs2
            
            pearson_corr, pearson_pvalue = pearsonr(
                highRiskCases[f'scaled_prs_{prs1}'],
                highRiskCases[f'scaled_prs_{prs2}']
                
            )
            
            # R-squared
            r_squared = pearson_corr ** 2
            
            print(f'n high risk cases for {prs1_label} OR {prs2_label} ',highRiskCases.shape[0])
            n_list.append(f'n high risk cases for {prs1_label} OR {prs2_label} :'+str(highRiskCases.shape[0]))
            comparison.append(f'{prs1_label} or {prs2_label}:high risk')
            
            stat.append(pearson_corr)
            pvalue.append(pearson_pvalue)
            
            print(f"Pearson r-squared  of high risk cases for {prs1_label} v {prs2_label} : ", r_squared)
            print(f"Pearson Correlation  of high risk cases for {prs1_label} v {prs2_label} : ", pearson_corr)
            print(f"Pearson Correlation pvalue  of {prs1_label} v {prs2_label} : ", pearson_pvalue)
            
            n = highRiskCases[f'scaled_prs_{prs1}'].shape[0]
            conf_int_test = pearson_confidence_interval(pearson_corr, n)
            conf_int.append(conf_int_test)
            
            print("95% Confidence Interval:", conf_int_test)
            print('')
            print('######################  spearmanr correlation coefficient ################################')
            print('')
            test.append('spearman_r')
            threshold.append(test_threshold)
            
            # Spearman correlation
            spearmanr_corr_test, spearman_pvalue_test = spearmanr(highRiskCases[f'bin_{prs1_label}'], highRiskCases[f'bin_{prs2_label}'])
            print("Spearman Correlation:", spearmanr_corr_test)
            print("Spearman pvalue:", spearman_pvalue_test)
            print('')
            
            n_list.append(f'n high risk cases for {prs1_label} OR {prs2} :'+str(highRiskCases.shape[0]))
            comparison.append(f'{prs1_label} or {prs2_label}:high risk')
            
            stat.append(spearmanr_corr_test)
            pvalue.append(spearman_pvalue_test)
            conf_int.append('NA')
            print('')
            print('######## t test for high risk cases dataset ##############')
            print('')
            
            test.append('paired_ttest')
            threshold.append(test_threshold)
            
            t_stat, p_value = ttest_rel(highRiskCases[f'scaled_prs_{prs1_label}'], highRiskCases[f'scaled_prs_{prs2}'])
            print("Paired t-test P-value for high risk cases :", p_value)
            print('t statistic = ',t_stat)
            print('')
            
            n_list.append(f'n high risk cases for {prs1_label} OR {prs2_label} :'+str(highRiskCases.shape[0]))
            comparison.append(f'{prs1_label} or {prs2_label}:high risk')
            
            stat.append(t_stat)
            pvalue.append(p_value)
            conf_int.append('NA')
            
            
        ###### populate results temp with statistics  #########
        
        resultsTemp['comparison'] = comparison
        resultsTemp['test'] = test
        resultsTemp['threshold'] = threshold
        resultsTemp['statistic'] = stat
#       resultsTemp['n_in_stat'] = n_list
        resultsTemp['pvalue'] = pvalue
        resultsTemp['95_CI'] = conf_int
        
        results = pd.concat([results,resultsTemp],ignore_index=True)
        
    print('McNemar, ttest, and pearson correlation analysis is completed for combinations:')
    print(results['comparison'].unique())
    
    results.to_csv(f'{scoresPath}/McNemarStatsTestsAcrossPRSCalculations.csv',index=False)
        
    
    
#   
#prsFile = '/Users/kerimulterer/prsInteractive/results/type2Diabetes/summedEpi/scores/combinedPRSGroups.csv'
#scoresPath = '/Users/kerimulterer/prsInteractive/results/type2Diabetes/summedEpi/scores'
#
#
#calculate_mcnemar_test(prsFile,scoresPath)
    
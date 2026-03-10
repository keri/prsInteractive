#helper functions for downloading and processing data used in model training/testing and creating PRS 

import pandas as pd
import numpy as np
import time
import plotly.graph_objects as go
import plotly.io as pio
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib_venn import venn3, venn2
import seaborn as sns
from sklearn import metrics
from itertools import combinations
from typing import List, Dict, Tuple, Optional
import os
import scipy as sp
from scipy.stats import pearsonr
import sys

# =============================================================================
# COLOR AND MARKER DEFINITIONS
# =============================================================================

# ============================================================================
# NEW FUNCTIONS FOR COMBINED ANALYSIS - Added for product/summed comparisons
# ============================================================================

# Extended color scheme for product/summed variants
COHORT_COLORS_EXTENDED = {
    'main': '#E69F00',      # Orange
    'epi': '#56B4E9',       # Sky blue
    'epi+main': '#CC79A7',  # Pinkish purple
    'cardio': '#009E73',    # Bluish green
    'all': '#F0E442',       # Yellow
    'combined': '#D55E00',  # Vermillion
    # Product variants (same as base)
    'epi_product': '#56B4E9',
    'cardio_product': '#009E73',
    'epi+main_product': '#CC79A7',
    'all_product': '#F0E442',
    'main_product': '#E69F00',
    # Summed variants (darker shades)
    'epi_summed': '#0073A8',      # Darker blue
    'cardio_summed': '#006B52',   # Darker green
    'epi+main_summed': '#9F5580', # Darker purple
    'all_summed': '#C7B800',      # Darker yellow
    'main_summed': '#C77D00',     # Darker orange
}

COHORT_MARKERS_EXTENDED = {
    'main': 'o',       # Circle
    'epi': 's',        # Square
    'epi+main': 'p',   # Pentagon
    'cardio': '^',     # Triangle
    'all': 'x',        # X
    'combined': 'D',   # Diamond
    # Product/summed use same markers as base
    'epi_product': 's',
    'epi_summed': 's',
    'cardio_product': '^',
    'cardio_summed': '^',
    'epi+main_product': 'p',
    'epi+main_summed': 'p',
    'all_product': 'x',
    'all_summed': 'x',
    'main_product': 'o',
    'main_summed': 'o',
}


# Color-blind friendly palette (Wong 2011)
COHORT_COLORS = {
    'main': '#E69F00',      # Orange
    'epi': '#56B4E9',       # Sky blue
    'epi+main': '#CC79A7',  # Pinkish purple
    'cardio': '#009E73',    # Bluish green
    'all': '#F0E442',       # Yellow
    'combined': '#D55E00'   # Vermillion
}

COHORT_MARKERS = {
    'main': 'o',       # Circle
    'epi': 's',        # Square
    'epi+main': 'p',   # Pentagon
    'cardio': '^',     # Triangle up
    'all': 'x'         # X
}

# Marker sizes for different contexts
MARKER_SIZES = {
    'control': 40,
    'case_low_risk': 40,
    'case_high_risk_single': 50,
    'case_high_risk_multi': 50
}

# Control colors
CONTROL_COLOR = '#BEBEBE'  # Gray
LOW_RISK_COLOR = '#333333' #black

# ── Draw order: background → foreground ───────────────────────────────
# main always foreground; then epi variants; then cardio variants behind
DRAW_PRIORITY = [
    'cardio_product',
    'cardio_summed',
    'epi_product',
    'epi_summed',
    'main',
]

# ── Pre-assign display colour by model priority (main = highest) ───────
COLOUR_PRIORITY = ['main', 'epi_summed', 'epi_product', 'cardio_summed', 'cardio_product']


# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

def get_model_color(model: str) -> str:
    """Get color for a model.
    
    Checks COHORT_COLORS first (base models), then COHORT_COLORS_EXTENDED
    (product/summed variants), then falls back to black.
    """
    return COHORT_COLORS.get(model) or COHORT_COLORS_EXTENDED.get(model, '#000000')


def get_model_marker(model: str) -> str:
    """Get marker for a model.
    
    Checks COHORT_MARKERS first (base models), then COHORT_MARKERS_EXTENDED
    (product/summed variants), then falls back to circle.
    """
    return COHORT_MARKERS.get(model) or COHORT_MARKERS_EXTENDED.get(model, 'o')

def get_model_color_extended(model: str) -> str:
    """Get color for a model, with fallback."""
    return COHORT_COLORS_EXTENDED.get(model,'#888888')

def get_model_marker_extended(model: str) -> str:
    """Get marker for a model, with fallback."""
    return COHORT_MARKERS_EXTENDED.get(model, 'o')

def _draw_rank(model: str) -> int:
    """Lower rank = drawn first (background). Higher rank = drawn last (foreground)."""
    for i, pattern in enumerate(DRAW_PRIORITY):
        if model == pattern:
            return i
    # Any model not explicitly listed goes behind cardio (rank -1)
    return -1

def _assign_display_color(row):
    """Return the colour of the highest-priority model that flags this case."""
    for model in COLOUR_PRIORITY:
        if model in row.index and row[model]:
            return get_model_color(model)
    return LOW_RISK_COLOR



def create_case_control_histogram(df,pheno_col,continous_col,figPath,figsize=(12,6)):
    """
    Plot the distribution of a cases and controls
    
    Parameters:
    -----------
    df : pandas.DataFrame
        DataFrame containing prs/bin calculations
    pheno_col : str
        Column name of the binarized version of the phenotype
    continous_col : str
        column name of the measure to be plotted
    figsize : tuple
        Figure size
        
    Returns:
    --------
    matplotlib.figure.Figure
        Distribution plot
    """
    plt.figure(figsize=figsize)
    
    if not set(df[pheno_col].unique()).issubset({0, 1}):
        df['phenotype'] = df['PHENOTYPE'] - 1
        pheno_col = 'phenotype'
        
        
    # Split data
    cases = df[df[pheno_col] == 1][[continous_col]]
    controls = df[df[pheno_col] == 0][[continous_col]]
    
    # Plot histograms
    plt.hist(controls, bins=30, alpha=0.6, label='Controls', color='skyblue', edgecolor='black')
    plt.hist(cases, bins=30, alpha=0.6, label='Cases', color='salmon', edgecolor='black')
    
    case_mean = cases[continous_col].mean()
    # Plot threshold line
    plt.axvline(x=case_mean, color='black', linestyle='--', linewidth=2)
    
    
    plt.title(f'Distribution of cases/controls across: {continous_col}')
    plt.xlabel(continous_col)
    plt.ylabel('Frequency')
    plt.grid(alpha=0.3)
    
    plt.savefig(f'{figPath}/histogramCaseControl.{continous_col}.png')
    plt.close()
    
    
def create_sankey_plot_clinical_data(df,figPath,use_epi_main=False):
    
    if use_epi_main:
        prs_mathods = ['cardio','epi','main','epi+main']
    else:
        prs_methods = ['cardio', 'epi', 'main'] 
        
    print(f"Total holdout samples: {len(df)}")
    print(f"Cases in holdout: {df['status'].sum()}")
    
    # check which binary columns exist:
    binary_cols = [col for col in df.columns if col.endswith('_binary')]
    print(f"\nAvailable binary columns: {binary_cols}")
    
    for binar_col in binary_cols:
        # Filter for cases with low clinical risk (binary = 1)
        df_filtered = df[(df['PHENOTYPE'] == 1) & (df[clinical_binary_col] == 1)].copy()
        
        # Count flows from each PRS method to combined
        flows = []
        for method in prs_methods:
            high_risk_col = f'{method}_high_risk'
            
            if high_risk_col not in df_filtered.columns:
                print(f"Warning: {high_risk_col} not found in data")
                continue
            
            # Count individuals high-risk in this method -> combined high-risk
            high_to_high = len(df_filtered[(df_filtered[high_risk_col] == 1) & 
                                            (df_filtered['combined_high_risk'] == 1)])
            
            # Count individuals high-risk in this method -> combined low-risk
            high_to_low = len(df_filtered[(df_filtered[high_risk_col] == 1) & 
                                            (df_filtered['combined_high_risk'] == 0)])
            
            flows.append({
                'source': f'{method.upper()} High Risk',
                'target': 'Combined High Risk',
                'value': high_to_high,
                'color': COHORT_COLORS[method]
            })
            
            if high_to_low > 0:
                flows.append({
                    'source': f'{method.upper()} High Risk',
                    'target': 'Combined Low Risk',
                    'value': high_to_low,
                    'color': COHORT_COLORS[method]
                })
                
            print(f"{method.upper()}: {df_filtered[high_risk_col].sum()} high-risk individuals")
            print(f"  -> Combined High: {high_to_high}")
            print(f"  -> Combined Low: {high_to_low}")
            
            print(f"\nCombined high-risk: {df_filtered['combined_high_risk'].sum()}")
            
            # Create Sankey diagram
            # Build unique labels
            all_labels = []
            label_dict = {}
            
            for flow in flows:
                if flow['source'] not in label_dict:
                    label_dict[flow['source']] = len(all_labels)
                    all_labels.append(flow['source'])
                if flow['target'] not in label_dict:
                    label_dict[flow['target']] = len(all_labels)
                    all_labels.append(flow['target'])
                    
            # Build source, target, value, and color lists
            sources = [label_dict[flow['source']] for flow in flows]
            targets = [label_dict[flow['target']] for flow in flows]
            values = [flow['value'] for flow in flows]
            link_colors = [flow['color'] for flow in flows]
            
            # Create node colors
            node_colors = []
            for label in all_labels:
                if 'CARDIO' in label:
                    node_colors.append(COHORT_COLORS['cardio'])
                elif 'EPI' in label and 'main' not in label.lower():
                    node_colors.append(COHORT_COLORS['epi'])
                elif 'MAIN' in label:
                    node_colors.append(COHORT_COLORS['main'])
                elif 'Combined High' in label:
                    node_colors.append('#D55E00')  # Red-orange for combined high risk
                elif 'Combined Low' in label:
                    node_colors.append('#999999')  # Gray for combined low risk
                else:
                    node_colors.append('#CCCCCC')
                    
            # Create Sankey diagram
            fig = go.Figure(data=[go.Sankey(
                node=dict(
                    pad=15,
                    thickness=20,
                    line=dict(color='black', width=0.5),
                    label=all_labels,
                    color=node_colors
                ),
                link=dict(
                    source=sources,
                    target=targets,
                    value=values,
                    color=link_colors
                )
            )])
            
            fig.update_layout(
                title=dict(
                    text=f"PRS High-Risk Classification Flow to Combined PRS<br>" + 
                        f"<sub>Cases with Low Clinical Risk (n={len(df_filtered)})</sub>",
                    x=0.5,
                    xanchor='center'
                ),
                font=dict(size=12, family='Arial'),
                width=1000,
                height=600,
                margin=dict(l=20, r=20, t=80, b=20)
            )
            
            # Save figure
            fig.write_html(f'{figPath}/Sankey{binary_col}.PRStoCombinedHighRisk.html')
            print("\nSankey plot saved as 'sankey_combined_prs_flow.html'")
            
            pio.write_image(fig,f'{figPath}/Sankey{binary_col}.PRStoCombinedHighRisk.png')
            
            # Also display if in interactive environment
            fig.show()
            
            # Print summary statistics
            print("\n" + "="*60)
            print("SUMMARY STATISTICS")
            print("="*60)
            
            for method in prs_methods:
                high_risk_col = f'{method}_high_risk'
                if high_risk_col in df_filtered.columns:
                    n_high = df_filtered[high_risk_col].sum()
                    overlap_combined = len(df_filtered[(df_filtered[high_risk_col] == 1) & 
                                                        (df_filtered['combined_high_risk'] == 1)])
                    if n_high > 0:
                        overlap_pct = (overlap_combined / n_high) * 100
                        print(f"{method.upper()}: {overlap_pct:.1f}% of high-risk -> Combined high-risk")
                        
            combined_high = df_filtered['combined_high_risk'].sum()
            print(f"\nTotal combined high-risk: {combined_high} ({combined_high/len(df_filtered)*100:.1f}%)")
            
            
            
            
            
def create_venn_diagram(df,figPath,image_str):
    
    epiORmainCases = set(df[((df['epiMain'] == 1) & (df['group'] == 'Cases'))]['feature'].tolist())
    epiCases = set(df[((df['iprs'] == 1) & (df['group'] == 'Cases'))]['feature'].tolist())
    mainCases = set(df[((df['prs'] == 1) & (df['group'] == 'Cases'))]['feature'].tolist())
    
    epiORmainBoth = set(df[((df['epiMain'] == 1) & (df['group'] == 'Everyone'))]['feature'].tolist())
    epiBoth = set(df[((df['iprs'] == 1) & (df['group'] == 'Everyone'))]['feature'].tolist())
    mainBoth = set(df[((df['prs'] == 1) & (df['group'] == 'Everyone'))]['feature'].tolist())
    
#   epiORmainControls = set(df[((df['epiMain'] == 1) & (df['group'] == 'Controls'))]['feature'].tolist())
#   epiControls = set(df[((df['iprs'] == 1) & (df['group'] == 'Controls'))]['feature'].tolist())
#   mainControls = set(df[((df['prs'] == 1) & (df['group'] == 'Controls'))]['feature'].tolist())
    
#   additiveControls = set(df[((df['additive'] == 1) & (df['group'] == 'Controls'))]['feature'].tolist())
#   epiControls = set(df[((df['iprs'] == 1) & (df['group'] == 'Controls'))]['feature'].tolist())
#   mainControls = set(df[((df['prs'] == 1) & (df['group'] == 'Controls'))]['feature'].tolist())
    
#   additive = set(df[df['additive'] == 1]['feature'].tolist())
#   epi = set(df[df['iprs'] == 1]['feature'].tolist())
#   main = set(df[df['prs'] == 1]['feature'].tolist())
    
    plt.figure(figsize=(10,10))
    
    venn = venn3([epiORmainCases, epiCases, mainCases],('epi+main iPRS', 'epi iPRS','main PRS'),alpha = .5)
#   venn = venn([epiCases, mainCases],('epi iPRS','main PRS'),alpha = .5)
#   plt.title("Important Features For High Risk Cases",fontsize=20)
    # Set colors and scales
    venn.get_patch_by_id('100').set_color('purple')  # additive
    venn.get_patch_by_id('010').set_color('blue')  # epi
    venn.get_patch_by_id('001').set_color('red')  # main
#   # Set the offset for labels in subsets
    #venn.subset_labels[0].set_x((venn.subset_labels[0].get_position()[0] - .02))
#   venn.subset_labels[1].set_x((venn.subset_labels[1].get_position()[0] + .05))
#   venn.subset_labels[3].set_x((venn.subset_labels[3].get_position()[0] - .03))
    
    # Set the offset for labels in subsets
#   venn.set_labels[0].set_y((venn.set_labels[1].get_position()[1]))
#   venn.set_labels[1].set_y((venn.set_labels[1].get_position()[1]))
#   venn.set_labels[2].set_y((venn.set_labels[1].get_position()[1]))
    
    plt.savefig(f'{figPath}/groupedImportantFeaturesCases.{image_str}.Venn.png')
    
    
    ###########################  PLOT CASES AND EVERYONE ON THE SAME PLOT BUT DIFFERENT AXIS  ###################
    plt.figure(figsize=(10,10))
    
    plt.subplot(2,1,1)
    
    venn = venn3([epiORmainCases, epiCases, mainCases],('iPRS epi+main', 'iPRS epi','PRS main'), alpha = .5)
#   venn = venn3([epiCases, mainCases],('epi iPRS','main PRS'), alpha = .5)
    plt.title("Cases",fontsize=20)
    # Set colors and scales
    venn.get_patch_by_id('100').set_color('purple')  # additive
    venn.get_patch_by_id('010').set_color('blue')  # epi
    venn.get_patch_by_id('001').set_color('red')  # main
    
    plt.subplot(2,1,2)
    
    venn = venn3([epiORmainBoth, epiBoth, mainBoth],('iPRS epi+main', 'iPRS epi','PRS main'),alpha = .5)
#   venn = venn3([epiBoth, mainBoth],('epi iPRS','main PRS'),alpha = .5)
    plt.title("Everyone",fontsize=20)
    # Set colors and scales
    venn.get_patch_by_id('100').set_color('purple')  # additive
    venn.get_patch_by_id('010').set_color('blue')  # epi
    venn.get_patch_by_id('001').set_color('red')  # main
    
    plt.savefig(f'{figPath}/groupedImportantFeaturesEveryoneCasesTogether.{image_str}.Venn.png')
    
    #####################  PLOT EVERYONE HIGH V LOW PER GROUP ####################################################
    
    plt.figure(figsize=(10,10))
    
    venn = venn3([epiORmainBoth, epiBoth, mainBoth],('iPRS epi+main', 'iPRS epi','PRS main'),alpha = .5)
#   venn = venn([epiCases, mainCases],('epi iPRS','main PRS'),alpha = .5)
#   plt.title("Important Features For High Risk Cases",fontsize=20)
    # Set colors and scales
    venn.get_patch_by_id('100').set_color('purple')  # additive
    venn.get_patch_by_id('010').set_color('blue')  # epi
    venn.get_patch_by_id('001').set_color('red')  # main
#   # Set the offset for labels in subsets
    #venn.subset_labels[0].set_x((venn.subset_labels[0].get_position()[0] - .02))
#   venn.subset_labels[1].set_x((venn.subset_labels[1].get_position()[0] + .05))
#   venn.subset_labels[3].set_x((venn.subset_labels[3].get_position()[0] - .03))
    
    # Set the offset for labels in subsets
#   venn.set_labels[0].set_y((venn.set_labels[1].get_position()[1]))
#   venn.set_labels[1].set_y((venn.set_labels[1].get_position()[1]))
#   venn.set_labels[2].set_y((venn.set_labels[1].get_position()[1]))
    
    plt.savefig(f'{figPath}/groupedImportantFeaturesCasesAndControls.{image_str}.Venn.png')
    
    
    
    
def create_violin_plot(df,model_type,figurePath):
    
    '''df has scaled prs for epi, main, and epi+main created from separate models combined into 1 dataframe
       columns = scaled_prs, model,phenotype'''
    title = 'KDE with boxplot of cases/controls for different model types'
    dfCopy = df.copy()
    dfCopy['phenotype'] = 99
    
    
    dfCopy.loc[dfCopy['PHENOTYPE'] == 1, 'phenotype'] = 'Without T2D'
    dfCopy.loc[dfCopy['PHENOTYPE'] == 2, 'phenotype'] = 'With T2D'
    
    fig, ax = plt.subplots(figsize=(10,10))
    sns.violinplot(data=dfCopy,x='model',y='scaled_prs',hue='phenotype',ax=ax,alpha=0.3,palette={"Without T2D": "b", "With T2D": "r"})
    sns.despine(left=True)
    for violin in ax.collections:
        violin.set_alpha(0.6)
    fig.suptitle(title,fontsize=18,fontweight='bold')
    ax.set_xlabel("Model Type",size = 16,alpha=0.7)
    ax.set_ylabel("Scaled PRS/iPRS",size = 16,alpha=0.7)
    leg = plt.legend()
    ax.get_legend().set_visible(False)
    
    fig.savefig(f'{figurePath}/{model_type}.violinplot.png')
    
def create_prevalence_plot(df,model_type,figurePath,value_type='scaled_prs'):
    dfCopy = df.copy()
    #sort scaled prs for entire dataset and break up into centiles
    dfCopy.sort_values([value_type],inplace=True)
    
    #get min max values to create centiles
#   prsMin = dfCopy['scaled_prs'].min()-.0001
#   print(prsMin)
#   prsMax = dfCopy['scaled_prs'].max()+.0001
#   print(prsMax)
    
#   bins = np.linspace(prsMin,prsMax,num=11)
#   len(bins)
    
    try:
        dfCopy['centile'] = pd.qcut(dfCopy[value_type], 10, labels=list(range(1,11)),duplicates='drop')
        
        dfCount = dfCopy.groupby(['centile','PHENOTYPE']).count()
        dfCount.reset_index(inplace=True)
        dfCountCases = dfCount[dfCount['PHENOTYPE'] == 2]
        dfCountCases['case_count'] = dfCountCases[value_type]
        
        dfTotal = dfCopy.groupby(['centile']).count()
        dfTotal.reset_index(inplace=True)
        dfTotal['total'] = dfTotal[value_type]
        
        dfPrevalence = dfTotal[['centile','total']].merge(dfCountCases[['centile','case_count']],on=['centile'])
        
        dfPrevalence['prevalence'] = dfPrevalence['case_count'] / dfPrevalence['total']
        
        title = f'prevalence over deciles {model_type} '
        ax = dfPrevalence.sort_values(['centile']).plot(x='centile',y='prevalence',kind='scatter',s=100,figsize=(10,10))
        ax.plot()
        ax.set_xlabel('PRS Decile Bins')
        ax.set_title(title)
        xticks = list(range(1,11))
        xlabels = [str(i) for i in xticks]
        ax.set_xticks(xticks,labels=xlabels)
#       ax.set_xticklabels(dfPrevalence['centile'])
        plt.savefig(f'{figurePath}/{model_type}.prevalencePlot.png')
        plt.close()
    except ValueError:
        pass
        
def calculate_tpr_fpr(cm):
    '''cm = np.array([[TP,FP],[FN,TN]])
       TPR = #TP / (#TP + #FN)
       FPR = #FP / (#FP + #TN)'''
    tp = cm[0][0]
    fn = cm[1][0]
    fp = cm[0][1]
    tn = cm[1][1]
    TPR = tp / (tp + fn)
    FPR = fp / (fp + tn)
    return(TPR,FPR)

    
def standard_confusion_matrix(y_true, y_predict):
    '''returns TN, FP, FN, TN based on predictions'''
    TP = np.sum((y_true == 1) & (y_predict == 1))
    TN = np.sum((y_true == 0) & (y_predict == 0))
    FP = np.sum((y_true == 0) & (y_predict == 1))
    FN = np.sum((y_true == 1) & (y_predict == 0))
    
    return (np.array([[TP,FP],[FN,TN]]),TP,FP,TN,FN)

def cb_matrix(cost_TP, cost_FP, cost_FN, cost_TN):
    return (np.array([[cost_TP, cost_FP],[cost_FN, cost_TN]]))

def y_predict_threshold(y_probs, threshold):
    '''returns matrix with predictions based on threshold'''
    changed = np.array(y_probs > threshold).astype(int)
    return(changed)


def calculate_auc_metrics(y_true, y_prs,thresholds):
    '''creates an array of predicted y's for every every threshold from prs range.
        Returns tpr,fpr for every threshold'''
    tprList = []
    fprList = []
    
    tpList = []
    tnList = []
    fpList = []
    fnList = []
    
    
    for threshold in thresholds:
        new_predictions = y_predict_threshold(y_prs, threshold)
        confusion_matrix,tp,fp,tn,fn = standard_confusion_matrix(y_true,new_predictions)
        tpr,fpr = calculate_tpr_fpr(confusion_matrix)
        tprList.append(tpr)
        fprList.append(fpr)
        tpList.append(tp)
        tnList.append(tn)
        fpList.append(fp)
        fnList.append(fn)
        
        
    return(thresholds,tprList,fprList,tpList,tnList,fpList,fnList)

    
def true_positive(df):
    y_true = np.array(df['PHENOTYPE'])
    #change the phenotypes from 1,2 to 0,1
    y_true[y_true == 1] = 0
    y_true[y_true == 2] = 1
    return(y_true)

def get_thresholds(df):
    #get min max values to create centiles
    prsMin = df['scaled_prs'].min()-.0001
    prsMax = df['scaled_prs'].max()+.0001
    
    thresholds = np.linspace(prsMin,prsMax,num=1000)
    return(thresholds)

def get_best_threshold(dfAUC):
    dfAUC['diff_TPR_FPR'] = dfAUC['TPR'] - dfAUC['FPR']
    best_threshold = dfAUC['diff_TPR_FPR'].max()
    
    dfAUC['color'] = 'orange'
    #dfAUC.loc[dfAUC['diff_TPR_FPR'] == best_threshold,'color'] = 'darkblue'
    dfAUC['size'] = 30
    #dfAUC.loc[dfAUC['diff_TPR_FPR'] == best_threshold,'size'] = 200
    return(dfAUC,best_threshold)

#def draw_auc_plots_prs_iprs(filePRS,fileIPRS,figurePath,model_type):
def draw_auc_plots_prs_iprs(prsDf,figurePath,threshold):
    '''filepathPRS/IPRS : str path/to/AUC_metricTable generated in create_auc_graph'''
#   df = pd.read_csv(filepathPRS)
#   df1 = pd.read_csv(filepathIPRS)
    
    model_type = f'mainVmain+epi_iPRS1_threshold{threshold}'
    #get threshold values for both iPRS and prs
    thresholdMax = prsDf['scaled_prs_main'].max()
    thresholdMin = prsDf['scaled_prs_main'].min()
    if prsDf['scaled_prs'].max() > thresholdMax:
        thresholdMax = prsDf['scaled_prs'].max()
    if prsDf['scaled_prs'].min() < thresholdMin:
        thresholdMin = prsDf['scaled_prs'].min()
        
    prsMin = round(thresholdMin -.0001,2)
    prsMax = round(thresholdMax +.0001,2)
    
    thresholds = np.linspace(prsMin,prsMax,num=1000)
    
    df = prsDf[['scaled_prs_main','PHENOTYPE']]
    df.rename(columns={'scaled_prs_main':'scaled_prs'},inplace=True)
    df1 = prsDf[['scaled_prs','PHENOTYPE']]
    
    #AUC dataframes returned have columns : [prs (thresholds),
    dfAUC2Main,aucMain = create_auc_dataframe(df,thresholds=thresholds)
    dfAUC2Main['model'] = 'main'
    dfAUC2MainEpi,aucMainEpi = create_auc_dataframe(df1,thresholds=thresholds)
    dfAUC2MainEpi['model'] = 'epi+main'
    
    dfAUC3 = pd.concat([dfAUC2Main,dfAUC2MainEpi],ignore_index=True)
    
    dfAUC3.to_csv(f'{figurePath}/AUC_metrics_table_{model_type}.csv',index=False)
    
#   fprListMain = dfAUC2Main['FPR']
#   tprListMain = dfAUC2Main['TPR']
    best_threshold_main = dfAUC2Main['diff_TPR_FPR'].max()
#   aucMain = metrics.auc(np.array(fprListMain),np.array(tprListMain))
    
#   fprListMainEpi = dfAUC2MainEpi['FPR']
#   tprListMainEpi = dfAUC2MainEpi['TPR']
    best_threshold_main_epi = dfAUC2MainEpi['diff_TPR_FPR'].max()
    aucMainEpi = metrics.auc(np.array(fprListMainEpi),np.array(tprListMainEpi))
    aucMainEpi = round(aucMainEpi,2)
    
    
#   dfAUC2Main['color'] = 'orange'
#   dfAUC2MainEpi['color']
#   dfAUC2MainEpi.loc[dfAUC['diff_TPR_FPR'] == best_threshold,'color'] = 'darkblue'
#   dfAUC2Main['size'] = 10
    
    fig = plt.figure(figsize=(10,10))
    ax1 = fig.add_subplot(111)
    ax2 = ax1.twinx()
    
    dfAUC2Main.plot(x='FPR',y='TPR',c='orange',alpha=.6,ax=ax1,label='Main')
    dfAUC2MainEpi.plot(x='FPR',y='TPR',c='blue',alpha=.6,ax=ax2,label='Main+Epi')
    
    ax1.set_xlabel('False Positive Rate : #FP/(#FP + #TN)')
    ax1.set_ylabel('True Positive Rate : #TP / (#TP + #FN)')
    ax1.set_title(f'Main AUC = {aucMain}  Main+Epi AUC = {aucMainEpi} threshold {threshold}')
    lines = ax1.get_lines() + ax2.get_lines()
    
    fig.legend(lines,[l.get_label() for l in lines],loc='upper center',bbox_to_anchor=(0,0,1,1))
    
    plt.savefig(f'{figurePath}/{model_type}.AUC.png')
    
# =============================================================================
# PAIRWISE CORRELATION PLOTS
# =============================================================================
    
def draw_pairwise_correlation_plot(
    cases: pd.DataFrame,
    controls: pd.DataFrame,
    model1: str,
    model2: str,
    output_path: str,
    plot_title: str = None,
    r_squared: float = None,
    pvalue: float = None,
    figsize: Tuple[int, int] = (10, 10),
    save: bool = True
) -> plt.Figure:
    """
    Create a correlation scatter plot for two PRS models.
    
    Parameters:
    -----------
    cases : pd.DataFrame
        Cases data with columns: scaled_prs_{model}, color, size, alpha, marker
    controls : pd.DataFrame
        Controls data with columns: scaled_prs_{model}, color, size, alpha
    model1 : str
        Name of first model (x-axis)
    model2 : str
        Name of second model (y-axis)
    output_path : str
        Directory to save the plot
    plot_title : str, optional
        Custom plot title
    r_squared : float, optional
        Pearson r-squared value to display
    pvalue : float, optional
        P-value to display
    figsize : Tuple[int, int], default=(10, 10)
        Figure size
    save : bool, default=True
        Whether to save the plot
        
    Returns:
    --------
    plt.Figure
        The matplotlib figure object
    """
    fig, ax = plt.subplots(figsize=figsize)
    
    # Column names for PRS scores
    model1_col = f'scaled_prs_{model1}'
    model2_col = f'scaled_prs_{model2}'
    
    # Scatter plot for controls (uniform marker)
    ax.scatter(
        x=controls[model1_col].values,
        y=controls[model2_col].values,
        marker=".",
        c=controls['color'],
        s=controls['size'],
        alpha=controls['alpha'],
        label='Controls',
        zorder=1
    )
    
    # Scatter plot for cases (colored by high-risk group, model-specific markers)
    # Group by marker type to create separate scatter plots
    for marker_type in cases['marker'].unique():
        marker_cases = cases[cases['marker'] == marker_type]
        if len(marker_cases) == 0:
            continue
        
        # Further group by color for legend entries
        for color in marker_cases['color'].unique():
            color_marker_cases = marker_cases[marker_cases['color'] == color]
            if len(color_marker_cases) == 0:
                continue
            
            # Determine label
            if color == '#666666':
                label = 'Cases (low risk)'
            else:
                # Find which model this color represents
                model_name = [k for k, v in COHORT_COLORS.items() if v == color]
                if model_name:
                    label = f'High on {model_name[0]}'
                else:
                    label = 'High risk'
                    
            ax.scatter(
                x=color_marker_cases[model1_col].values,
                y=color_marker_cases[model2_col].values,
                marker=marker_type,
                c=color_marker_cases['color'],
                s=color_marker_cases['size'],
                alpha=color_marker_cases['alpha'],
                label=label,
                edgecolors='black',
                linewidths=0.5,
                zorder=3 if color != 'black' else 2
            )
            
    # Labels and styling
    ax.set_xlabel(f'{model1.upper()} PRS (standardized)', fontsize=12, fontweight='bold')
    ax.set_ylabel(f'{model2.upper()} PRS (standardized)', fontsize=12, fontweight='bold')
    ax.grid(True, alpha=0.3, linestyle='--')
    
    # Title
    if plot_title is None:
        plot_title = f'PRS Correlation: {model2.upper()} vs {model1.upper()}'
    ax.set_title(plot_title, fontsize=14, fontweight='bold', pad=15)
    
    # Legend
    handles, labels = ax.get_legend_handles_labels()
    # Remove duplicate labels
    by_label = dict(zip(labels, handles))
    ax.legend(by_label.values(), by_label.keys(), 
             loc='upper left', fontsize=10, framealpha=0.95)
    
    # Add statistics if provided
    if r_squared is not None:
        r_squared_round = round(r_squared, 3)
        r = np.sqrt(r_squared) if r_squared >= 0 else -np.sqrt(abs(r_squared))
        plt.text(
            0.98, 0.02,
            f'r = {r:.3f}\nr² = {r_squared_round}',
            transform=ax.transAxes,
            fontsize=11,
            verticalalignment='bottom',
            horizontalalignment='right',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.9, edgecolor='black'),
            family='monospace'
        )
        
    if pvalue is not None:
        pval_text = 'p < 0.001' if pvalue < 0.001 else f'p = {pvalue:.4f}'
        y_pos = 0.10 if r_squared is not None else 0.02
        plt.text(
            0.98, y_pos,
            pval_text,
            transform=ax.transAxes,
            fontsize=11,
            verticalalignment='bottom',
            horizontalalignment='right',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.9, edgecolor='black'),
            family='monospace'
        )
        
    plt.tight_layout()
    
    # Save plot
    if save:
        os.makedirs(output_path, exist_ok=True)
        filename = f'{model2}_vs_{model1}.correlation_plot.png'
        filepath = os.path.join(output_path, filename)
        plt.savefig(filepath, dpi=300, bbox_inches='tight')
        print(f"  ✓ Saved: {filename}")
        
    return fig



def create_all_pairwise_plots(
    combined_prs: pd.DataFrame,
    models_to_keep: List[str],
    output_path: str,
    threshold: int = 8,
#   priority_order: List[str] = None,
    include_all_model: bool = False,
    plot_prefix: str = '',
    verbose: bool = True
) -> Dict[str, Dict]:
    """
    Create correlation plots for all pairwise combinations of PRS models.
    
    Parameters:
    -----------
    combined_prs : pd.DataFrame
        DataFrame with columns: IID, PHENOTYPE, scaled_prs_{model}, bin_{model}
    models_to_keep : List[str]
        List of model names to include
    output_path : str
        Directory to save plots
    threshold : int, default=8
        Threshold for high-risk classification
    priority_order : List[str], optional
        Order of priority for color assignment
    include_all_model : bool, default=False
        Whether to include 'all' model in coloring
    plot_prefix : str, default=''
        Prefix for plot filenames
    verbose : bool, default=True
        Print progress information
        
    Returns:
    --------
    Dict[str, Dict]
        Dictionary mapping model pairs to their statistics
    """
    # Add 'all' model if requested
    models_for_colors = models_to_keep.copy()
#   if include_all_model and 'all' in combined_prs.columns:
#       models_for_colors.append('all')
    
    # Separate cases and controls
    cases = combined_prs[combined_prs['PHENOTYPE'] == 2].copy()
    controls = combined_prs[combined_prs['PHENOTYPE'] == 1].copy()
    
    
    # Set controls to default color
    controls['color'] = CONTROL_COLOR
    controls['size'] = MARKER_SIZES['control']
    controls['alpha'] = 0.6
    controls['marker'] = '.'
    
    cases['color'] = LOW_RISK_COLOR
    cases['size'] = MARKER_SIZES['case_low_risk']
    cases['alpha'] = 0.6
    cases['marker'] = '+'
    
    
    # Generate all pairwise combinations
    model_pairs = list(combinations(models_to_keep, 2))
    
    if verbose:
        print(f"\n{'='*70}")
        print(f"Creating {len(model_pairs)} pairwise correlation plots")
        print(f"Models: {', '.join(models_to_keep)}")
        print(f"Output: {output_path}")
        print(f"{'='*70}\n")
        
    # Store statistics for each pair
    pair_statistics = {}
    
    # Create plots for each pair
    for i, (model1, model2) in enumerate(model_pairs, 1):
        # Check if columns exist
        model1_col = f'scaled_prs_{model1}'
        model2_col = f'scaled_prs_{model2}'
        
        if model1_col not in combined_prs.columns or model2_col not in combined_prs.columns:
            if verbose:
                print(f"  ⚠️  Skipping {model1} vs {model2}: columns not found")
            continue
        
        # Calculate correlation on full dataset
        valid_mask = (
            combined_prs[model1_col].notna() & 
            combined_prs[model2_col].notna()
        )
        
        corr, pval = pearsonr(
            combined_prs.loc[valid_mask, model1_col],
            combined_prs.loc[valid_mask, model2_col]
        )
        
        #assign cases color based on high_risk cases for each prs model
        cases.loc[valid_mask, model1_col]['color'] = COHORT_COLORS['case_high_risk_single']
        cases.loc[valid_mask, model1_col]['size'] = MARKER_SIZES['case_high_risk_single']
        cases.loc[valid_mask, model1_col]['alpha'] = 0.85
        
        #assign cases color based on high_risk cases for each prs model
        cases.loc[valid_mask, model2_col]['color'] = COHORT_COLORS['case_high_risk_single']
        cases.loc[valid_mask, model2_col]['size'] = MARKER_SIZES['case_high_risk_single']
        cases.loc[valid_mask, model2_col]['alpha'] = 0.85
        
        r_squared = corr ** 2
        n_samples = valid_mask.sum()
        
        # Store statistics
        pair_key = f'{model1}_vs_{model2}'
        pair_statistics[pair_key] = {
            'correlation': corr,
            'r_squared': r_squared,
            'pvalue': pval,
            'n_samples': n_samples
        }
        
        # Create plot title
        plot_title = f'{plot_prefix}{model2.upper()} vs {model1.upper()} PRS Correlation'
        if include_all_model:
            plot_title += ' (with all model)'
            
        # Create plot
        if verbose:
            print(f"  [{i}/{len(model_pairs)}] {model2} vs {model1}: r={corr:.3f}, p={pval:.2e}")
            
        #assign color to cases
            
            
        fig = draw_pairwise_correlation_plot(
            cases=cases,
            controls=controls,
            model1=model1,
            model2=model2,
            output_path=output_path,
            plot_title=plot_title,
            r_squared=r_squared,
            pvalue=pval
        )
        
        plt.close(fig)
        
    if verbose:
        print(f"\n{'='*70}")
        print(f"✅ Created {len(pair_statistics)} correlation plots")
        print(f"{'='*70}\n")
        
    return pair_statistics


# =============================================================================
# COMPREHENSIVE HIGH-RISK PLOT
# =============================================================================

def plot_all_high_risk_cases(
        combined_prs: pd.DataFrame,
        models_to_keep: List[str],
        output_path: str,
        x_model: str = 'epi',
        y_model: str = 'main',
        threshold: int = 8,
        model_order: List[str] = None,
        include_low_risk: bool = True,
        include_all_model: bool = False,
        figsize: Tuple[int, int] = (14, 10),
        show_model_counts: bool = True,
        verbose: bool = True
) -> Tuple[plt.Figure, pd.DataFrame]:
        """
        Plot high-risk cases per model, layered background → foreground.

        Each model is drawn as its own scatter layer using its color.
        Alpha distinguishes variant type:
            - product / base : 0.85 (solid)
            - summed         : 0.45 (semi-transparent, so product shows through)
            - combined       : 0.90

        Marker shape distinguishes exclusivity:
            - * star  : case is high-risk on this model ONLY
            - P plus  : case is high-risk on ≥2 models (shared)

        Legend is compact: one entry per model (color swatch + name).
        A separate mini-legend explains star vs plus.
        """
        # ── Alpha per variant type ─────────────────────────────────────────────
        def _alpha(model):
                if model == 'main':
                        return 0.50
                if model.endswith('_summed'):
                        return 0.60
                return 0.85   # product, base, or main
    
        if verbose:
                print(f"\n{'='*70}")
                print(f"Creating high-risk plot: {y_model.upper()} vs {x_model.upper()}")
                print(f"{'='*70}\n")
            
        # ── Build model_order ─────────────────────────────────────────────────
        if model_order is None:
                model_order = models_to_keep.copy()
#       else:
        draw_order = sorted(
            [m for m in model_order if m in models_to_keep],
            key=_draw_rank   # ascending: background first, foreground last
        )
    
                            
        if include_all_model and 'all' not in draw_order:
                draw_order = ['all'] + draw_order
            
        if verbose:
                print(f"  Plot order (background→foreground): {' → '.join(draw_order)}")
            
        # ── Cases / controls ──────────────────────────────────────────────────
        cases    = combined_prs[combined_prs['PHENOTYPE'] == 2].copy()
        controls = combined_prs[combined_prs['PHENOTYPE'] == 1].copy()
    
        # Verify axes exist
        for col in [f'scaled_prs_{x_model}', f'scaled_prs_{y_model}']:
                if col not in combined_prs.columns:
                        raise ValueError(f"Column {col!r} not found in combined_prs")
                    
        # ── High-risk membership per model ────────────────────────────────────
        high_risk_by_model = {}
        for model in draw_order:
                bin_col = f'bin_{model}'
                if bin_col in cases.columns:
                        high_risk_by_model[model] = cases[bin_col] > threshold
                else:
                        if verbose:
                                print(f"  WARNING: {bin_col} not found — skipping {model}")
                            
        if not high_risk_by_model:
                raise ValueError("No bin_ columns found for any model in model_order")
            
        # For each case, count how many models flag it as high-risk
        hr_df = pd.DataFrame(high_risk_by_model, index=cases.index).fillna(False)
        n_flagged = hr_df.sum(axis=1)           # 0 = low risk, 1 = exclusive, ≥2 = shared
    
        # ── Summary assignments DataFrame (for CSV export) ────────────────────
        assignments_df = pd.DataFrame({
                'IID':            cases.get('IID', pd.Series(cases.index, index=cases.index)),
                'n_models':       n_flagged,
                'high_on_models': hr_df.apply(lambda r: ', '.join(r.index[r]), axis=1),
                f'prs_{x_model}': cases[f'scaled_prs_{x_model}'].values,
                f'prs_{y_model}': cases[f'scaled_prs_{y_model}'].values,
        }, index=cases.index)
    
        # hr_df rows = cases, cols = model bool flags
        display_colors = hr_df.apply(_assign_display_color, axis=1)
        
    
        # ── Figure ────────────────────────────────────────────────────────────
        fig, ax = plt.subplots(figsize=figsize)
    
        # Controls (background)
        ax.scatter(
                controls[f'scaled_prs_{x_model}'],
                controls[f'scaled_prs_{y_model}'],
                c=CONTROL_COLOR, s=MARKER_SIZES['control'],
                alpha=0.35, marker='.', zorder=0, rasterized=True
        )
    
        # Low-risk cases
        if include_low_risk:
                low_risk_idx = n_flagged[n_flagged == 0].index
                if len(low_risk_idx):
                        ax.scatter(
                                cases.loc[low_risk_idx, f'scaled_prs_{x_model}'],
                                cases.loc[low_risk_idx, f'scaled_prs_{y_model}'],
                                c=LOW_RISK_COLOR, s=30, alpha=.85,
                                marker='P', zorder=1, label=f'Low risk (n={len(low_risk_idx)})',
                                rasterized=True
                        )
                    
                    
        # ── Single scatter pass for all high-risk cases ────────────────────────
        high_risk_any_idx = n_flagged[n_flagged >= 1].index
        exclusive_idx     = n_flagged[n_flagged == 1].index
        shared_idx        = n_flagged[n_flagged >= 2].index
        
        # ── Priority color assignment (main > epi > cardio) ───────────────────
        # Iterate draw_order low→high priority so the highest (main) overwrites last.
        point_color = pd.Series('', index=high_risk_any_idx, dtype=object)
        for model in draw_order:          # cardio → epi → main
            if model not in high_risk_by_model:
                continue
            is_high  = high_risk_by_model[model]
            flagged  = is_high[is_high].index.intersection(high_risk_any_idx)
            point_color.loc[flagged] = get_model_color(model)
    
        # ── Two scatter passes: shape encodes exclusivity, color encodes priority ──
        if len(exclusive_idx):
            ax.scatter(
                cases.loc[exclusive_idx, f'scaled_prs_{x_model}'],
                cases.loc[exclusive_idx, f'scaled_prs_{y_model}'],
                c=point_color.loc[exclusive_idx].tolist(),
                s=55, alpha=0.85, marker='*',
                zorder=5, edgecolors='none', rasterized=True
            )
            
        if len(shared_idx):
            ax.scatter(
                cases.loc[shared_idx, f'scaled_prs_{x_model}'],
                cases.loc[shared_idx, f'scaled_prs_{y_model}'],
                c=point_color.loc[shared_idx].tolist(),
                s=45, alpha=0.85, marker='P',
                zorder=5, edgecolors='none', rasterized=True
            )
            
        # ── Legend: highest priority first, counts = total flagged by each model ──
        legend_model_handles = []
        for model in reversed(draw_order):      # main → epi → cardio in legend
            if model not in high_risk_by_model:
                continue
            is_high  = high_risk_by_model[model]
            n_excl   = int(is_high[is_high & (n_flagged == 1)].sum())
            n_shared = int(is_high[is_high & (n_flagged >= 2)].sum())
            n_total  = n_excl + n_shared
            n_shown  = int((point_color == get_model_color(model)).sum())
            if verbose:
                print(f"  {model}: {n_total} high-risk  "
                        f"(exclusive={n_excl}, shared={n_shared}, shown in colour={n_shown})")
            if n_total:
                color = get_model_color(model)
                legend_model_handles.append(
                    plt.scatter([], [], c=color, alpha=0.85,
                                s=50, marker='s',
                                label=f'{model.upper()}  (n={n_total})')
                )
        

                    
        # ── Compact legend ────────────────────────────────────────────────────
        # Model colors
        model_legend = ax.legend(
                handles=legend_model_handles,
                loc='upper left',
                fontsize=8,
                framealpha=0.9,
                edgecolor='#CCCCCC',
                title='Models',
                title_fontsize=9,
                markerscale=1.1,
        )
        ax.add_artist(model_legend)
    
        # Marker shape legend (star = exclusive, plus = shared)
        from matplotlib.lines import Line2D
        shape_handles = [
                Line2D([0], [0], marker='*', color='#666', linestyle='none',
                                markersize=9, label='Exclusive (1 model)'),
                Line2D([0], [0], marker='P', color='#666', linestyle='none',
                                markersize=8, label='Shared (≥2 models)'),
        ]
        ax.legend(handles=shape_handles, loc='upper right',
                            fontsize=8, framealpha=0.9, edgecolor='#CCCCCC',
                            title='Marker', title_fontsize=9)
    
        # ── Stats box ─────────────────────────────────────────────────────────
        if show_model_counts:
                n_any   = (n_flagged > 0).sum()
                n_multi = (n_flagged > 1).sum()
                lines = ['High-risk counts:']
                for model in model_order:
                        if model in high_risk_by_model:
                                n = high_risk_by_model[model].sum()
                                lines.append(f'  {model}: {n}')
                lines += [f'', f'Total: {n_any}', f'Multi-model: {n_multi}']
                ax.text(
                        0.98, 0.02, '\n'.join(lines),
                        transform=ax.transAxes, fontsize=9,
                        va='bottom', ha='right',
                        bbox=dict(boxstyle='round', facecolor='white',
                                            alpha=0.95, edgecolor='#CCCCCC'),
                        family='monospace'
                )
            
        # ── Axes / title ──────────────────────────────────────────────────────
        ax.set_xlabel(f'{x_model.upper()} PRS (standardized)', fontsize=13, fontweight='bold')
        ax.set_ylabel(f'{y_model.upper()} PRS (standardized)', fontsize=13, fontweight='bold')
        ax.set_title(
                f'High-Risk Cases: {y_model.upper()} vs {x_model.upper()} PRS\n'
                f'Threshold: bin > {threshold}',
                fontsize=15, fontweight='bold', pad=15
        )
        ax.grid(True, alpha=0.25, linestyle='--')
        ax.axhline(0, color='k', lw=0.4, alpha=0.4)
        ax.axvline(0, color='k', lw=0.4, alpha=0.4)
        plt.tight_layout()
    
        # ── Save ──────────────────────────────────────────────────────────────
        os.makedirs(output_path, exist_ok=True)
        suffix = '_with_all' if include_all_model else ''
        filename  = f'all_high_risk_cases_{y_model}_vs_{x_model}{suffix}.png'
        csv_name  = f'high_risk_case_assignments_{y_model}_vs_{x_model}{suffix}.csv'
        plt.savefig(os.path.join(output_path, filename),  dpi=300, bbox_inches='tight')
        assignments_df.to_csv(os.path.join(output_path, csv_name), index=False)
    
        if verbose:
                print(f"  ✓ {filename}")
                print(f"  ✓ {csv_name}\n")
            
        return fig, assignments_df

# =============================================================================
# CORRELATION MATRIX
# =============================================================================

def create_correlation_matrix_plot(
    combined_prs: pd.DataFrame,
    models_to_keep: List[str],
    output_path: str,
    figsize: Tuple[int, int] = None,
    verbose: bool = True
) -> plt.Figure:
    """
    Create a correlation matrix heatmap for all models.
    
    Parameters:
    -----------
    combined_prs : pd.DataFrame
        DataFrame with scaled_prs_{model} columns
    models_to_keep : List[str]
        List of model names
    output_path : str
        Directory to save plot
    figsize : Tuple[int, int], optional
        Figure size (auto-calculated if None)
    verbose : bool, default=True
        Print progress
        
    Returns:
    --------
    plt.Figure
        The matplotlib figure object
    """
    if verbose:
        print(f"\n{'='*70}")
        print(f"Creating correlation matrix heatmap")
        print(f"{'='*70}\n")
        
    # Prepare data for correlation matrix
    prs_columns = [f'scaled_prs_{model}' for model in models_to_keep]
    
    # Filter to available columns
    available_cols = [col for col in prs_columns if col in combined_prs.columns]
    
    if not available_cols:
        raise ValueError("No PRS columns found in data")
        
    # Calculate correlation matrix
    corr_matrix = combined_prs[available_cols].corr()
    
    # Rename columns/index to model names only
    rename_map = {f'scaled_prs_{m}': m.upper() for m in models_to_keep}
    corr_matrix.rename(columns=rename_map, index=rename_map, inplace=True)
    
    # Auto-calculate figure size
    if figsize is None:
        n = len(corr_matrix)
        figsize = (max(8, n * 1.5), max(6, n * 1.2))
        
    # Create figure
    fig, ax = plt.subplots(figsize=figsize)
    
    
    
    # Create custom colormap (blue for negative, white for 0, red for positive)
    colors = ['#0571b0', 'white', '#ca0020']
    n_bins = 100
    cmap = mcolors.LinearSegmentedColormap.from_list('custom', colors, N=n_bins)
    
    im = ax.imshow(corr_matrix, cmap=cmap, vmin=-1, vmax=1, aspect='auto')
    
    # Add colorbar
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label('Pearson Correlation (r)', rotation=270, labelpad=20, fontweight='bold')
    
    # Add correlation values as text
    for i in range(len(corr_matrix)):
        for j in range(len(corr_matrix)):
            value = corr_matrix.iloc[i, j]
            color = 'white' if abs(value) > 0.5 else 'black'
            ax.text(j, i, f'{value:.2f}', ha='center', va='center', 
                   color=color, fontsize=11, fontweight='bold')
            
    # Set ticks and labels
    ax.set_xticks(range(len(corr_matrix)))
    ax.set_yticks(range(len(corr_matrix)))
    ax.set_xticklabels(corr_matrix.columns, rotation=45, ha='right', fontsize=11, fontweight='bold')
    ax.set_yticklabels(corr_matrix.index, fontsize=11, fontweight='bold')
    
    ax.set_title('PRS Model Correlation Matrix', fontsize=14, fontweight='bold', pad=20)
    
    plt.tight_layout()
    
    # Save
    os.makedirs(output_path, exist_ok=True)
    filename = 'correlation_matrix_heatmap_across_models.png'
    filepath = os.path.join(output_path, filename)
    plt.savefig(filepath, dpi=300, bbox_inches='tight')
    
    if verbose:
        print(f"  ✓ Saved: {filename}\n")
        print(f"{'='*70}\n")
        
    return fig


# =============================================================================
# MAIN PLOTTING FUNCTION
# =============================================================================

def create_all_prs_plots(
    combined_prs: pd.DataFrame,
    models_to_keep: List[str],
    output_path: str,
    threshold: int = 8,
    priority_order: List[str] = ['cardio','epi','epi+main','main'],
    create_pairwise: bool = False,
    create_comprehensive: bool = True,
    create_matrix: bool = True,
    include_all: bool = False,
    comprehensive_axes: List[Tuple[str, str]] = None,
    verbose: bool = True
) -> Dict:
    """
    Create all PRS visualization plots in one function call.
    
    This is the main entry point for creating a complete visualization suite.
    
    Parameters:
    -----------
    combined_prs : pd.DataFrame
        DataFrame with PRS scores and bins
    models_to_keep : List[str]
        List of model names to include
    output_path : str
        Directory to save all plots
    threshold : int, default=8
        Threshold for high-risk classification
    priority_order : List[str], optional
        Priority order for color assignment
    create_pairwise : bool, default=True
        Whether to create pairwise correlation plots
    create_comprehensive : bool, default=True
        Whether to create comprehensive high-risk plots
    create_matrix : bool, default=True
        Whether to create correlation matrix
    comprehensive_axes : List[Tuple[str, str]], optional
        List of (x_model, y_model) pairs for comprehensive plots
        Default: [('epi', 'main')]
    verbose : bool, default=True
        Print progress information
        
    Returns:
    --------
    Dict
        Dictionary with keys:
        - 'pairwise_stats': Statistics from pairwise plots
        - 'comprehensive_assignments': DataFrames from comprehensive plots
        - 'output_path': Where plots were saved
    """
    results = {
        'pairwise_stats': {},
        'comprehensive_assignments': {},
        'output_path': output_path
    }
    
    print("\n" + "="*70)
    print("PRS VISUALIZATION SUITE")
    print("="*70)
    print(f"Models: {', '.join(models_to_keep)}")
    print(f"Threshold: bin > {threshold}")
    print(f"Output: {output_path}")
    print("="*70)
    
    # Create output directory
    os.makedirs(output_path, exist_ok=True)
    
    if include_all and 'all' in combined_prs.columns:
        models_to_keep.append('all')
        priority_order.append('all')
        
    # Build priority order from models_to_keep using base order as a guide.
    # models_to_keep may contain suffixed names like 'cardio_product', 'epi_summed'
    # which won't match the base priority_order list directly.
    # Instead: expand each base name to cover its product/summed variants in order.
    priority_order_filtered = []
    for base in priority_order:
        for model in models_to_keep:
            # Match exact base name OR suffixed variants (e.g. cardio → cardio_product, cardio_summed)
            if model == base or model.startswith(base + '_'):
                if model not in priority_order_filtered:
                    priority_order_filtered.append(model)
                    
                
    # Append any remaining models_to_keep not covered by priority_order
    for model in models_to_keep:
        if model not in priority_order_filtered:
            priority_order_filtered.append(model)
            
    print(f'priority order for combined plot: {priority_order_filtered}')
    
    # 1. Pairwise correlation plots
    if create_pairwise:
        pair_stats = create_all_pairwise_plots(
            combined_prs=combined_prs,
            models_to_keep=models_to_keep,
            output_path=output_path,
            threshold=threshold,
            include_all_model=include_all,
#           priority_order=None,
            verbose=verbose
        )
        results['pairwise_stats'] = pair_stats
        
    # 2. Comprehensive high-risk plots
    if create_comprehensive:
        if comprehensive_axes is None:
            # Default: create main vs epi if both are in models
            if 'main' in models_to_keep and 'epi_product' in models_to_keep:
                comprehensive_axes = [('epi_product', 'main')]
            elif 'main' in models_to_keep and 'epi_summed' in models_to_keep:
                comprehensive_axes = [('epi_summed', 'main')]
            elif len(models_to_keep) >= 2:
                comprehensive_axes = [(models_to_keep[-2], models_to_keep[-1])]
            else:
                comprehensive_axes = []
                    
        for x_model, y_model in comprehensive_axes:
            if x_model in models_to_keep and y_model in models_to_keep:
                fig, assignments = plot_all_high_risk_cases(
                    combined_prs=combined_prs,
                    models_to_keep=models_to_keep,
                    output_path=output_path,
                    x_model=x_model,
                    y_model=y_model,
                    threshold=threshold,
                    model_order=priority_order_filtered,
                    include_all_model = False,
                    show_model_counts = True,
                    figsize = (14, 10),
                    verbose=verbose
                )
                
                
                plt.close(fig)
                results['comprehensive_assignments'][f'{y_model}_vs_{x_model}'] = assignments
                
    # 3. Correlation matrix
    if create_matrix:
        fig = create_correlation_matrix_plot(
            combined_prs=combined_prs,
            models_to_keep=models_to_keep,
            output_path=output_path,
            verbose=verbose
        )
        plt.close(fig)
        
    print("="*70)
    print("✅ ALL PLOTS CREATED SUCCESSFULLY")
    print("="*70)
    print(f"\nOutput directory: {output_path}")
    if create_pairwise:
        print(f"  - {len(results.get('pairwise_stats', {}))} pairwise correlation plots")
    if create_comprehensive:
        print(f"  - {len(results.get('comprehensive_assignments', {}))} comprehensive high-risk plots")
    if create_matrix:
        print(f"  - 1 correlation matrix heatmap")
    print("\n")
    
    return results

def create_auc_dataframe(df,thresholds=[]):
    y_true = true_positive(df)
    
    #prs will be used as the yHat based on threshold
    yPRS = np.array(df['scaled_prs'])
    
    if len(thresholds) > 0:
        pass
    else:
        #the thresholds are scaled prs broken up into 1K increments
        thresholds = get_thresholds(df)
        
    thresholds,tprList,fprList,tpList,tnList,fpList,fnList = calculate_auc_metrics(y_true,yPRS,thresholds)
    dfAUC = pd.DataFrame({'scaled_prs':thresholds,'FPR':fprList,'TPR':tprList,'TP':tpList,'FP':fpList,
                         'TN':tnList,'FN':fnList})
    dfAUC2,best_threshold = get_best_threshold(dfAUC)
    auc = metrics.auc(np.array(fprList),np.array(tprList))
    
    return(dfAUC2,auc)


def create_auc_graph(df,model_type,figurePath):
    
    dfAUC2,auc = create_auc_dataframe(df)
    auc = round(auc,2)
    fprList = dfAUC2['FPR']
    tprList = dfAUC2['TPR']
    
    dfAUC2.to_csv(f'{figurePath}/AUC_metrics_table_{model_type}.csv',index=False)
    
    fig,ax = plt.subplots(figsize=(10,10))
    ax.scatter(x=fprList,y=tprList,marker='.',c=dfAUC2['color'].values,s=dfAUC2['size'].values)
    ax.set_xlabel('False Positive Rate : #FP/(#FP + #TN)')
    ax.set_ylabel('True Positive Rate : #TP / (#TP + #FN)')
    ax.set_title(f'AUC = {auc}')
    plt.savefig(f'{figurePath}/{model_type}.AUC.png')
    plt.close()
#   return(best_threshold)
    
def create_density_plot(df,model_type,figurePath,prs_col='scaled_prs'):
# def create_density_plot(df,model_type,figurePath,threshold):
    try:
        meanDiff = round(df.groupby(['PHENOTYPE']).mean()[prs_col].diff().loc[2],2)
        print('mean diff between cases and controls is : ',meanDiff)
        meanControl =  df.groupby(['PHENOTYPE']).mean()[prs_col].loc[1]
        meanCase =  df.groupby(['PHENOTYPE']).mean()[prs_col].loc[2]
        
        title = f'Combined PRS for {model_type} - mean diff = {meanDiff}'
        ###split the data into groups based on types
        g = df.groupby('PHENOTYPE')
        
        ## From here things change as I make use of the seaborn library
        yes = g.get_group(2)
        no = g.get_group(1)
        
        fig, ax = plt.subplots(figsize=(10,10))
        
        #
        ax = sns.kdeplot(data=yes[prs_col], label='cases', ax=ax,color='red',shade=True)
        ax = sns.kdeplot(data=no[prs_col], label='controls', ax=ax,color='blue',shade=True)
        
        ax.legend(['cases','controls'])
        
        ax.set_title(title)
        fig.savefig(f'{figurePath}/{model_type}.densityPlot.png')
        plt.close()
    except KeyError:
        pass
        
        
        
def create_box_plot(prsCopy,model_type,figurePath,box_value='scaled_prs'):
    print('creating box plots ...')
    try:
        meanDiff = prsCopy.groupby(['PHENOTYPE'])[box_value].mean().diff().loc[2]
        print(meanDiff)
        ax = prsCopy.boxplot(column=box_value,by='PHENOTYPE',figsize=(10,10))
        ax.plot()
        title = f'prs mean diff for {model_type} = {meanDiff}'
        plt.title(title)
        plt.savefig(f'{figurePath}/{model_type}.boxplot.{box_value}.png')
        plt.close()
    except KeyError: #means there are no cases to plot
        pass
    
    
def create_saturation_graph(df,model_type,figurePath):
    '''df dataframe with columns : number_features:int, risk_direction:str, mean_diff:float'''
    
    dfCopy = df.copy()
    dfCopy('number_features', inplace=True)
    #group data by product and display sales as line chart
    dfCopy.groupby('risk_direction')['mean_diff'].plot(figsize=(10,10),legend=True)
    image_str = f'{model_type}.saturationPlot.png'
    plt.savefig(f'{figurePath}/{image_str}')
    
def create_density_plot_with_individuals(df,dfInd,model_type,figurePath,threshold):
    
    title = f'High Low Group with individual risk change when epi is added'
    ###split the data into groups based on types
    g = df.groupby('PHENOTYPE')
    
    ## From here things change as I make use of the seaborn library
    yes = g.get_group(2)
    no = g.get_group(1)
    model_type_temp = f'{model_type}.HighLowGroupCases'
    
    fig, ax = plt.subplots(figsize=(10,10))
    model_type = f'{model_type}.HighLowGroupCases'
    ax = sns.kdeplot(data=yes['scaled_prs'], label='with type 2 diabetes', ax=ax,color='red',shade=True)
    ax.scatter(dfInd['scaled_prs'].values,dfInd['y'].values,color='black',marker="^",s=40)
    ax.scatter(dfInd['scaled_prs_main'].values,dfInd['y'].values,color='black',marker="o",s=40)
    ax.axvline(x=threshold,linestyle='--')
    #plot the line graph connecting the two iprs,prs markers
    for i in range(dfInd.shape[0]):
        ax.plot([dfInd.loc[i]['scaled_prs'],dfInd.loc[i]['scaled_prs_main']],[dfInd.loc[i]['y'],dfInd.loc[i]['y']],c='black',alpha=.8)
        
        
    ax.legend(['With T2D','Main + Epi','Main Only'])
    
    ax.set_title(title)
    fig.savefig(f'{figurePath}/{model_type_temp}.diffWithScaledPRS{threshold}.densityPlot.png')
    plt.close()
    
    model_type_temp = f'{model_type}.HighLowGroup'
    fig, ax = plt.subplots(figsize=(10,10))
    
    ax = sns.kdeplot(data=yes['scaled_prs'], label='with type 2 diabetes', ax=ax,color='red',shade=True)
    ax = sns.kdeplot(data=no['scaled_prs'], label='without type 2 diabetes', ax=ax,color='blue',shade=True)
    ax.scatter(dfInd['scaled_prs'].values,dfInd['y'].values,c='black',marker="^",s=40)
    ax.scatter(dfInd['scaled_prs_main'].values,dfInd['y'].values,c='black',marker="o",s=40)
    ax.axvline(x=threshold,linestyle='--')
    #plot the line graph connecting the two iprs,prs markers
    for i in range(dfInd.shape[0]):
        ax.plot([dfInd.loc[i]['scaled_prs'],dfInd.loc[i]['scaled_prs_main']],[dfInd.loc[i]['y'],dfInd.loc[i]['y']],c='black',alpha=.8)
        
        
    ax.legend(['Cases','Controls','Main + Epi','Main Only'])
    
    ax.set_title(title)
    fig.savefig(f'{figurePath}/{model_type_temp}.diffWithScaledPRS{threshold}.densityPlot.png')
    plt.close()
    
    model_type_temp = f'{model_type}.HighLowGroupControlsTogether'
    fig, ax = plt.subplots(figsize=(10,10))
    
    ax = sns.kdeplot(data=df['scaled_prs'], label='entire High Low population', ax=ax,color='red',shade=True)
    ax.scatter(dfInd['scaled_prs'].values,dfInd['y'].values,c='black',marker="^",s=40)
    ax.scatter(dfInd['scaled_prs_main'].values,dfInd['y'].values,c='black',marker="o",s=40)
    ax.axvline(x=threshold,linestyle='--')
    #plot the line graph connecting the two iprs,prs markers
    for i in range(dfInd.shape[0]):
        ax.plot([dfInd.loc[i]['scaled_prs'],dfInd.loc[i]['scaled_prs_main']],[dfInd.loc[i]['y'],dfInd.loc[i]['y']],c='black',alpha=.8)
        
        
    ax.legend(['High Low Risk Group','Main + Epi','Main Only'])
    
    ax.set_title(title)
    fig.savefig(f'{figurePath}/{model_type_temp}.diffWithPRSThreshold{threshold}.densityPlot.png')
    plt.close()
    
    
def create_density_volcano_plot(df,dfInd,model_type,threshold,figurePath):
    dfInd2 = dfInd.sort_values(by=['scaled_prs_main'],ascending=False)
#	dfInd2['y'] = dfInd2['normalized_diff'].values
    meanDiff = df.groupby(['PHENOTYPE']).mean()['scaled_prs_main'].diff().loc[2]
    meanControl =  df.groupby(['PHENOTYPE']).mean()['scaled_prs_main'].loc[1]
    meanCase =  df.groupby(['PHENOTYPE']).mean()['scaled_prs_main'].loc[2]
    
#   title = f'Combined PRS for {model_type} : median diff = {meanDiff}'
    ###split the data into groups based on types
    g = df.groupby('PHENOTYPE')
    
    ## From here things change as I make use of the seaborn library
    yes = g.get_group(2)
    no = g.get_group(1)
    
#   fig, ax = plt.subplots(figsize=(10,10))
#   ax = sns.kdeplot(data=yes['scaled_prs'], label='with type 2 diabetes', ax=ax,color='white',shade=True)
#   ax = sns.kdeplot(data=no['scaled_prs'], label='without type 2 diabetes', ax=ax,color='white',shade=True)
#   ax = sns.kdeplot(data=yes['scaled_prs'], label='with type 2 diabetes', ax=ax,color='red',shade=True)
#   ax = sns.kdeplot(data=no['scaled_prs'], label='without type 2 diabetes', ax=ax,color='blue',shade=True)
    
    #uncomment this for high cases > threshold and controls < threshold
#   dfInd2CasesControls = dfInd2[(dfInd2['PHENOTYPE'] == 2)&(dfInd2['scaled_prs'] >1)|(dfInd2['PHENOTYPE'] == 1)&(dfInd2['scaled_prs'] < -1)]
    #use this for cases only
#   dfInd2CasesControls = dfInd2[(dfInd2['PHENOTYPE'] == 2)&(dfInd2['scaled_prs'] >threshold)|(dfInd2['PHENOTYPE'] == 2)&(dfInd2['scaled_prs'] < -1)]
    
    for c,cText in zip([99,1,2,-99],['CaseAndControl','Controls','Cases','LowControlsHighCases']):
        fig, ax = plt.subplots(figsize=(10,10))
        ax = sns.kdeplot(data=yes['scaled_prs_main'], label='with type 2 diabetes', ax=ax,color='red',shade=True)
        ax = sns.kdeplot(data=no['scaled_prs_main'], label='without type 2 diabetes', ax=ax,color='blue',shade=True)
        
        if c == 99: #cases AND controls 
            dfInd3 = dfInd2.copy()
        elif c in [1,2]: #cases OR controls graphs
            dfInd3 = dfInd2[((dfInd2['PHENOTYPE'] == c) & (dfInd2['scaled_prs'] < -threshold)) | ((dfInd2['PHENOTYPE'] == c) & (dfInd2['scaled_prs'] > threshold))]
        
        else: #low controls AND high cases graph
            dfInd3 = dfInd2[((dfInd2['PHENOTYPE'] == 1) & (dfInd2['scaled_prs'] < -threshold)) | ((dfInd2['PHENOTYPE'] == 2) & (dfInd2['scaled_prs'] > threshold))]
        
        
        
        #uncomment this to include all cases AND controls in plot 
    #   dfInd2CasesControls = dfInd2[(dfInd2['scaled_prs'] >1)|(dfInd2['scaled_prs'] < -1)]
#       xticks = dfInd2CasesControls['scaled_prs_main'].unique().tolist()
        ax.scatter(dfInd3['scaled_prs'].values,dfInd3['normalized_diff'].values,c=dfInd3['color'].values,alpha=.2)
#       ax.axvline(x=threshold,linestyle='--',c='black')
#       ax.axvline(x=-(threshold),linestyle='--',c='black')
    #   ax.legend(['With T2D','Without T2D','Epi decreases Risk','Epi increases risk'])
        #ax.legend(['Epi decreases Risk','Epi increases risk'])
        title = f'{cText} Scatter WRT iPRS Main Only Overlay Main Density:  mean diff = {meanDiff}'
        ax.set_title(title)
        fig.savefig(f'{figurePath}/{model_type}.WRTiPRSOverlayMain.{cText}Only.volcanoDensityPlot.png')
        plt.close()
        
def create_forest_plot(dfInd,model_type,figurePath,threshold):
    
    title = f'Cases with respect to PRS treshold = {threshold}'
    
    fig, ax = plt.subplots(figsize=(10,20))
    nPeople = dfInd.shape[0]
    ax.scatter(dfInd['scaled_prs'].values,dfInd['y'].values,c='purple',marker="^",s=30)
    ax.scatter(dfInd['scaled_prs_main'].values,dfInd['y'].values,c='cyan',marker="o",s=30)
#   if 'Regardless' in model_type:
#       ax.axvline(x=1,linestyle='--',c='purple')
#       ax.axvline(x=-1,linestyle='--',c='purple')
#   else:
#       ax.axvline(x=threshold,linestyle='--',c='purple')
    
    #plot the line graph connecting the two iprs,prs markers
    for i in range(dfInd.shape[0]):
        ax.plot([dfInd.loc[i]['scaled_prs'],dfInd.loc[i]['scaled_prs_main']],[dfInd.loc[i]['y'],dfInd.loc[i]['y']],c='black',alpha=.8)
        
    ax.set_ylim(0, nPeople)
#   ax.legend(['Main + Epi','Main Only'],prop = {"size": 30},markerscale=3)
    ax.set_xlabel('Difference scaled iPRS - scaled PRS',fontsize='40')
    ax.set_ylabel('Individuals',fontsize='40')
    #ax.set_title(title,fontsize='20')
    fig.savefig(f'{figurePath}/{model_type}.diffWithPRS.forestPlot.png')
    plt.close()
    
def create_important_feature_forest_plot(df,cohortFeatures,figurePath):
    '''graph the features with the that are important within each cohort compared with shapley values for other cohorts in a forest plot
    
    input : df : dataframe with Shapley values for all cohorts to graph columns = [feature:str  model:str(shapley_values) cohort:str(epi,main,epi+main,cardio)  feature_importance:float (Shapley value)  split:float  ranking:float  support:bool]
            cohortFeatures : dataframe with features within each cohort to plot columns = [model:str, feature:str]
    '''
    
    
    
    fig, ax = plt.subplots(figsize=(30,30))
    
    #plot the line graph connecting the two iprs,prs markers
    y=0
    labels = []
    cohorts = ['epi+main','epi','main','cardio']
    
    for cohort in cohorts:
        #get feature for that cohort
#       cohortPlot = df[df['cohort'] == cohort].sort_values([cohort]).set_index('feature')
        
        #get the features important in that cohort (met the importance threshold when passed into function)
        cohortFeatureToPlot = cohortFeatures[cohortFeatures['cohort'] == cohort].sort_values(['feature_importance'],ascending=False)['feature'].tolist()
        
        #get all of the Shapley values for all of the cohorts
        cohortDf = df[df['feature'].isin(cohortFeatureToPlot)]
        
        #sort the feature by feature importance and take the top value, this will have cohort as well
        cohortDf.sort_values(['feature','feature_importance'],ascending=False,inplace=True)
        filteredCohortFeatures = cohortDf.drop_duplicates(subset=['feature'],keep='first')
        
        #filter the features that are highest only in that cohort
        finalCohortFeatures = filteredCohortFeatures[filteredCohortFeatures['cohort'] == cohort].sort_values(['feature_importance'],ascending=False)['feature'].tolist()
        
        #set index for the dataframe with Shapley values for cohorts for features specific to that cohort
        cohortDf.set_index('feature',inplace=True)
        
#       features = cohortPlot.index.tolist()
        labels += finalCohortFeatures
        
        for f in finalCohortFeatures:
            featureToPlot = df[df['feature'] == f]
            
            y+=1
            mainDot = featureToPlot[featureToPlot['cohort'] == 'main']['feature_importance']
            epiDot = featureToPlot[featureToPlot['cohort'] == 'epi']['feature_importance']
            epiMainDot = featureToPlot[featureToPlot['cohort'] == 'epi+main']['feature_importance']
            cardioDot = featureToPlot[featureToPlot['cohort'] == 'cardio']['feature_importance']
            
            ax.scatter(mainDot,y,c='red',marker="D",s=100)
            ax.scatter(epiDot,y,c='blue',marker="o",s=100)
            ax.scatter(epiMainDot,y,c='purple',marker="X",s=100)
            ax.scatter(cardioDot,y,c='#f5a142',marker="v",s=100)
            
            ax.plot([mainDot,epiDot,epiMainDot,cardioDot],[y,y,y,y],c='black',alpha=.8)
            
    dfLimits = df[df['feature'].isin(labels)]
    ax.set_ylim(0, len(labels)+1)
    ax.tick_params(axis='y',labelsize=16,direction="out")
    ax.set_yticks(range(1,len(labels)+1),labels=labels)
    ax.set_xlim(dfLimits['feature_importance'].min()-1,dfLimits['feature_importance'].max()+1)
    plt.tight_layout()
    
    fig.savefig(f'{figurePath}/importantCohortFeatures.forestPlot.png')
    
    
def create_optimized_prevalence_plot(df,figurePath,image_str):
    '''prevalence plot across deciles for iPRS,PRS, and optimized
        ]


        figurePath = fig path for test or holdout set
        '''
    
    
    
    ####################################################################
    #                            Prevalence Plot                       #
    ####################################################################
    
    fig,ax = plt.subplots(figsize=(10,10))
    
    
    
    ###########################  PREVALENCE PLOT ################################
    title = 'Centile Prevalence Across PRS models '
    
    for model in df['model'].unique():
        plotDf = df[df['model'] == model]
        marker = plotDf['marker'].tolist()[0]
        color = plotDf['color'].tolist()[0]
        ax.scatter(x=plotDf['prevalence'],y=np.array(plotDf['bin'])-1,c=color,s=200,marker=marker,alpha=.7)
        ax.plot(plotDf['prevalence'],np.array(plotDf['bin'])-1,color=color,alpha=.4,linestyle='dashed')
        
    ax.set_ylabel('deciles')
    ax.set_xlabel('prevalence')
    ax.set_title(title)
    #	xlabels = [f"{i}\n{j}" for i,j in enumerate(bin_limits,start=1)]
    #	xlabels = df['bin_limits'].unique().tolist()
    yticks = list(range(0,10))
    ylabels = [str(i) for i in range(1,11)]
    ax.set_yticks(yticks,labels=ylabels)
    
    plt.savefig(f'{figurePath}/prevalenceAcrossModels{image_str}.scatter.png')
    
def plot_holdout_with_combined(
        holdout_df: pd.DataFrame,
        models_to_keep: List[str],
        output_path: str,
        threshold: int = 8,
        figsize: Tuple[int, int] = (14, 10),
        show_model_counts: bool = True,
        verbose: bool = True,
) -> Tuple[plt.Figure, pd.DataFrame]:
        """
        Plot holdout cases with x=scaled_prs_main, y=scaled_prs_combined.

        scaled_prs_combined provides the y-axis (the meta-score space); it is NOT
        used as a coloring category.  Cases are colored solely by which individual
        model(s) flagged them, using exactly the same per-model layer approach,
        priority order, alpha, and marker conventions as plot_all_high_risk_cases:

            Priority: cardio_product → cardio_summed → epi_product → epi_summed → … → main
            Alpha:    product/base = 0.85  |  summed = 0.45
            Marker:   * star = exclusive to one model  |  P plus = shared ≥2 models

        Cases not flagged by any individual model are shown as low-risk (gray P).

        Parameters
        ----------
        holdout_df : pd.DataFrame
                Holdout data.  Must contain scaled_prs_combined, scaled_prs_main,
                and bin_{model} columns for each model in models_to_keep.
        models_to_keep : List[str]
                Individual models from the validation filtering step.
        output_path : str
                Directory to save figure and CSV.
        threshold : int
                Bin threshold for high-risk classification (default 8).
        figsize : Tuple
                Figure size.
        show_model_counts : bool
                Show per-model count stats box.
        verbose : bool
                Print progress.

        Returns
        -------
        Tuple[plt.Figure, pd.DataFrame]
                Figure and per-case assignments DataFrame.
        """
        # Alpha gradient: first model in order (background) = most opaque,
        # last model (foreground) = most transparent so top layers don't
        # completely obscure what's beneath them.
        ALPHA_MAX  = 0.85   # first model (background)
        ALPHA_MIN  = 0.20   # last model  (foreground)
    
        def _positional_alpha(position: int, n_models: int) -> float:
                """Linear gradient: position 0 → ALPHA_MAX, position n-1 → ALPHA_MIN."""
                if n_models <= 1:
                        return ALPHA_MAX
                return ALPHA_MAX - (position / (n_models - 1)) * (ALPHA_MAX - ALPHA_MIN)
    
        # ── Validate required columns ─────────────────────────────────────────
        for col in ('scaled_prs_combined', 'scaled_prs_main'):
                if col not in holdout_df.columns:
                        print(f"  ERROR: {col} not found in holdout_df")
                        return None, None
            
        x_col = 'scaled_prs_main'
        y_col = 'scaled_prs_combined'
    
        # ── Build priority order — no 'combined' entry ────────────────────────
        base_order = ['cardio', 'epi', 'epi+main', 'main']
        model_order = []
        for base in base_order:
                for model in models_to_keep:
                        if (model == base or model.startswith(base + '_')) and model not in model_order:
                                model_order.append(model)
        for model in models_to_keep:
                if model not in model_order:
                        model_order.append(model)
        # 'combined' intentionally excluded — it provides the y-axis, not a color
                    
        if verbose:
                print(f"\nCreating holdout scatter: scaled_prs_combined vs scaled_prs_main")
                print(f"  Plot order (background→foreground): {' → '.join(model_order)}")
            
        # ── Cases / controls ──────────────────────────────────────────────────
        cases    = holdout_df[holdout_df['PHENOTYPE'] == 2].copy()
        controls = holdout_df[holdout_df['PHENOTYPE'] == 1].copy()
    
        # ── High-risk membership — use pre-computed {model}_high_risk bool cols ─
        # These are authoritative flags already set by assign_holdout_risk_bin;
        # do NOT re-derive from bin_{model} > threshold.
        high_risk_by_model = {}
        for model in model_order:
                hr_col = f'{model}_high_risk'
                if hr_col in cases.columns:
                        high_risk_by_model[model] = cases[hr_col].astype(bool)
                        if verbose:
                                print(f"  {model}: {high_risk_by_model[model].sum()} high-risk")
                else:
                        if verbose:
                                print(f"  WARNING: {hr_col} not found — skipping {model}")
                            
        # Membership matrix — sum excludes combined_high_risk per spec
        hr_df     = pd.DataFrame(high_risk_by_model, index=cases.index).fillna(False)
        n_flagged = hr_df.sum(axis=1)   # 0 = low risk, 1 = exclusive, ≥2 = shared
    
        # model_assigned = highest-priority model that flags each case (the color
        # that will actually be painted on the point in the plot)
        def _assigned_model(row):
                # model_order is background→foreground; last True wins (highest priority)
                flagged = [m for m in model_order if row.get(m, False)]
                return flagged[-1] if flagged else ''
        model_assigned = hr_df.apply(_assigned_model, axis=1)
    
        # ── Assignments DataFrame (for CSV) ───────────────────────────────────
        assignments_df = pd.DataFrame({
                'IID':            cases['IID'].values if 'IID' in cases.columns else cases.index,
                'n_models':       n_flagged.values,
                'high_on_models': hr_df.apply(lambda r: ', '.join(r.index[r]), axis=1).values,
                'model_assigned': model_assigned.values,
                'prs_combined':   cases[y_col].values,
        }, index=cases.index)
    
        # ── Figure ────────────────────────────────────────────────────────────
        fig, ax = plt.subplots(figsize=figsize)
    
        # Controls
        ax.scatter(
                controls[x_col], controls[y_col],
                c=CONTROL_COLOR, s=MARKER_SIZES['control'],
                alpha=0.35, marker='.', zorder=0, rasterized=True,
        )
    
        # Low-risk cases (no individual model flags them)
        low_risk_idx = n_flagged[n_flagged == 0].index
        if len(low_risk_idx):
                ax.scatter(
                        cases.loc[low_risk_idx, x_col],
                        cases.loc[low_risk_idx, y_col],
                        c='black', s=30, alpha=0.85,
                        marker='P', zorder=1, rasterized=True,
                        label=f'Low risk (n={len(low_risk_idx):,})'
                )
            
        # ── One scatter layer per model — identical logic to plot_all_high_risk_cases ─
        legend_model_handles = []
    
        # Only count models that are actually present in the data
        active_models = [m for m in model_order if m in high_risk_by_model]
        n_active = len(active_models)
    
        if verbose:
                print(f"  Alpha gradient: {ALPHA_MAX:.2f} (background) → {ALPHA_MIN:.2f} (foreground)")
            
        for position, model in enumerate(active_models):
                z_order = position + 2
            
                color = get_model_color(model)
                alpha = _positional_alpha(position, n_active)
            
                is_high      = high_risk_by_model[model]
                excl_idx     = is_high[is_high & (n_flagged == 1)].index
                shared_idx   = is_high[is_high & (n_flagged >= 2)].index
                n_excl       = len(excl_idx)
                n_shared     = len(shared_idx)
                n_total      = n_excl + n_shared
            
                if verbose:
                        print(f"  {model} (alpha={alpha:.2f}): excl(*)={n_excl}  shared(P)={n_shared}")
                    
                if n_excl:
                        ax.scatter(
                                cases.loc[excl_idx, x_col],
                                cases.loc[excl_idx, y_col],
                                c=color, s=55, alpha=alpha,
                                marker='*', zorder=z_order,
                                edgecolors='none', rasterized=True
                        )
                    
                if n_shared:
                        ax.scatter(
                                cases.loc[shared_idx, x_col],
                                cases.loc[shared_idx, y_col],
                                c=color, s=45, alpha=alpha,
                                marker='P', zorder=z_order,
                                edgecolors='none', rasterized=True
                        )
                    
                if n_total:
                        legend_model_handles.append(
                                plt.scatter([], [], c=color, alpha=max(alpha, 0.7),
                                                        s=50, marker='s',
                                                        label=f'{model.upper()}  (n={n_total})')
                        )
                    
        # ── Compact two-part legend (mirrors plot_all_high_risk_cases) ────────
        model_legend = ax.legend(
                handles=legend_model_handles,
                loc='upper left', fontsize=8, framealpha=0.9,
                edgecolor='#CCCCCC', title='Models', title_fontsize=9,
        )
        ax.add_artist(model_legend)
    
        from matplotlib.lines import Line2D
        shape_handles = [
                Line2D([0], [0], marker='*', color='#555', linestyle='none',
                                markersize=9, label='Exclusive (1 model)'),
                Line2D([0], [0], marker='P', color='#555', linestyle='none',
                                markersize=8, label='Shared (≥2 models)'),
        ]
        ax.legend(handles=shape_handles, loc='upper right',
                            fontsize=8, framealpha=0.9, edgecolor='#CCCCCC',
                            title='Marker', title_fontsize=9)
    
        # ── Stats box ─────────────────────────────────────────────────────────
        if show_model_counts:
                n_any   = (n_flagged > 0).sum()
                n_multi = (n_flagged > 1).sum()
                lines = ['High-risk counts:']
                for model in model_order:
                        if model in high_risk_by_model:
                                lines.append(f'  {model}: {high_risk_by_model[model].sum()}')
                lines += ['', f'Total: {n_any}', f'Multi-model: {n_multi}']
                ax.text(
                        0.98, 0.02, '\n'.join(lines),
                        transform=ax.transAxes, fontsize=9,
                        va='bottom', ha='right',
                        bbox=dict(boxstyle='round', facecolor='white',
                                            alpha=0.95, edgecolor='#CCCCCC'),
                        family='monospace'
                )
            
        # ── Axes / title ──────────────────────────────────────────────────────
        ax.set_xlabel('MAIN PRS (standardized)', fontsize=13, fontweight='bold')
        ax.set_ylabel('Combined PRS (standardized)', fontsize=13, fontweight='bold')
        ax.set_title(
                f'Holdout: High-Risk Cases by Model in Combined PRS Space\n'
                f'Threshold: bin > {threshold}  |  N = {len(cases):,} cases',
                fontsize=15, fontweight='bold', pad=15
        )
        ax.grid(True, alpha=0.25, linestyle='--')
        ax.axhline(0, color='k', lw=0.4, alpha=0.4)
        ax.axvline(0, color='k', lw=0.4, alpha=0.4)
        plt.tight_layout()
    
        # ── Save ──────────────────────────────────────────────────────────────
        os.makedirs(output_path, exist_ok=True)
        filename = 'holdout_high_risk_main_vs_combined.png'
        csv_name = 'holdout_high_risk_case_assignments_main_vs_combined.csv'
        plt.savefig(os.path.join(output_path, filename), dpi=300, bbox_inches='tight')
        assignments_df.to_csv(os.path.join(output_path, csv_name), index=False)
    
        if verbose:
                print(f"  ✓ {filename}")
                print(f"  ✓ {csv_name}\n")
            
        return fig, assignments_df    
    
if __name__ == '__main__':
    
#   pheno = sys.argv[1]
#   pheno = 'type2Diabetes'
#   scoresPath = f'/your/results/{pheno}/productEpi/scores' 
#   figPath = f'/your/results/{pheno}/productEpi/figures'
#   combinedPRSFile = f'{scoresPath}/combinedPRSGroups.csv'
#   combinedPRS = pd.read_csv(combinedPRSFile)
    create_qq_plot_groups(combinedPRS,figPath,use_epi_main=True)
    
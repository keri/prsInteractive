#Keri Multerer Jun 2023
#helper functions for downloading and processing data used in model training/testing and creating PRS 


import pandas as pd
import numpy as np
import time
import matplotlib.pyplot as plt
from matplotlib_venn import venn3, venn2
import seaborn as sns
from sklearn import metrics
import scipy as sp
import sys

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

def create_prevalence_plot(df,model_type,figurePath):
    dfCopy = df.copy()
    #sort scaled prs for entire dataset and break up into centiles
    dfCopy.sort_values(['scaled_prs'],inplace=True)

    #get min max values to create centiles
#   prsMin = dfCopy['scaled_prs'].min()-.0001
#   print(prsMin)
#   prsMax = dfCopy['scaled_prs'].max()+.0001
#   print(prsMax)

#   bins = np.linspace(prsMin,prsMax,num=11)
#   len(bins)
    
    try:
        dfCopy['centile'] = pd.qcut(dfCopy['scaled_prs'], 10, labels=list(range(1,11)),duplicates='drop')
        
        dfCount = dfCopy.groupby(['centile','PHENOTYPE']).count()
        dfCount.reset_index(inplace=True)
        dfCountCases = dfCount[dfCount['PHENOTYPE'] == 2]
        dfCountCases['case_count'] = dfCountCases['scaled_prs']
    
        dfTotal = dfCopy.groupby(['centile']).count()
        dfTotal.reset_index(inplace=True)
        dfTotal['total'] = dfTotal['scaled_prs']
    
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
    
    
#def create_qq_plot_groups(combinedPRS,figurePath,suffix):
#   ''' use input to color cases in different colors based on group:
#   input : dataframe columns = ['PHENOTYPE', 'prs', 'scaled_prs', 'prs_main', 'scaled_prs_main',
#   'epi_helps', 'diff', 'color', 'normalized_diff', 'bin_EpiMain',
#   'bin_Main', 'prs_Final', 'scaled_prs_Final', 'bin_Final', 'bin_limits',
#   'caseANDControlGroup', 'caseORControlGroup'] 
#   output : correlation plot with groups colors'''
#   
#   cases = combinedPRS[combinedPRS['PHENOTYPE'] == 2]
#   cases['color'] = 'black'
#   controls = combinedPRS[combinedPRS['PHENOTYPE'] == 1]
#   controls['color'] = '#776cb1'
#   
#   #assign colors based on categories
#   #category 1 : high risk cases > 8 main decile AND epi helps == 0
#   #category 2 : high risk cases > 8 epi decile AND epi helps == 1
#   #category 3 : cases for which model doesn't work well main decile < 9 AND epi helps == 0
#   #category 4 : cases for which model doesn't work well epi decils < 9 AND epi helps == 1
#   
#   for group in ['high risk cases PRS','high risk cases iPRS','cases predicted low risk PRS','cases predicted low risk iPRS','all','high risk','low risk']:
#       if group == 'high risk cases PRS':
#           fig,ax = plt.subplots(figsize=(10,10))
#           cases.loc[((cases['bin_Main'] == 9) & (cases['epi_helps'] == 0)),'color'] == 'rgb(211,0,63,0)'
#           cases.loc[((cases['bin_Main'] == 10) & (cases['epi_helps'] == 0)),'color'] == 'rgb(201,0,22,0)'
#           
#           model_type = f'HighRiskPRS.{suffix}'
#       #   ax.scatter(x=controls['scaled_prs'].values,y=controls['scaled_prs_main'].values,marker=".",c='#776cb1',s=40,alpha=.6,label='Controls')
#           ax.scatter(x=controls['scaled_prs'].values,y=controls['scaled_prs_main'].values,marker=".",c=controls['color'].values,s=40,alpha=.6,label='Controls')
#       #   ax.scatter(x=epiCases['scaled_prs'].values,y=mainCases['scaled_prs'].values,marker="+",c='black',s=40,alpha=.9,label='Cases')
#           ax.scatter(x=cases['scaled_prs'].values,y=cases['scaled_prs'].values,marker="+",c=cases['color'].values,s=40,alpha=.9,label='Cases')
#       #   ax.plot([main['scaled_prs'].min(),main['scaled_prs'].max()],[linreg.intercept+linreg.slope*epi['scaled_prs'].min(),linreg.intercept+linreg.slope*epi['scaled_prs'].max()])
#           
#           ax.set_xlabel('epi')
#           ax.set_ylabel('main only')
#           ax.grid(True) 
#           ax.set_title(f'Standardized High Risk PRS Cases : Main V Epi')
#           handles,labels = ax.get_legend_handles_labels()
#           ax.legend(handles,labels,loc='upper left',fontsize=15)
#       #   rsquared = round((linreg.rvalue*linreg.rvalue),2)
#           ax.text(-3.8, 3, f'{group}',fontsize=15)
#           plt.savefig(f'{figurePath}/PRSHigh.QQColorPlot.png')
#           ax.clear()
#           cases['color'] = 'black'
#           
#       if group == 'high risk cases iPRS':
#           cases.loc[((cases['bin_EpiMain'] == 9) & (cases['epi_helps'] == 1)),'color'] == 'rgb(0, 0, 205,.6)'
#           cases.loc[((cases['bin_EpiMain'] == 10) & (cases['epi_helps'] == 1)),'color'] == 'rgb(0, 0, 128,.6)'
#           
#           model_type = f'HighRiskiPRS.{suffix}'
#       #   ax.scatter(x=controls['scaled_prs'].values,y=controls['scaled_prs_main'].values,marker=".",c='#776cb1',s=40,alpha=.6,label='Controls')
#           ax.scatter(x=controls['scaled_prs'].values,y=controls['scaled_prs_main'].values,marker=".",c=controls['color'].values,s=40,alpha=.6,label='Controls')
#       #   ax.scatter(x=epiCases['scaled_prs'].values,y=mainCases['scaled_prs'].values,marker="+",c='black',s=40,alpha=.9,label='Cases')
#           ax.scatter(x=cases['scaled_prs'].values,y=cases['scaled_prs'].values,marker="+",c=cases['color'].values,s=40,alpha=.9,label='Cases')
#       #   ax.plot([main['scaled_prs'].min(),main['scaled_prs'].max()],[linreg.intercept+linreg.slope*epi['scaled_prs'].min(),linreg.intercept+linreg.slope*epi['scaled_prs'].max()])
#           
#           ax.set_xlabel('epi')
#           ax.set_ylabel('main only')
#           ax.grid(True) 
#           ax.set_title(f'Standardized High Risk iPRS Cases : Main V Epi')
#           handles,labels = ax.get_legend_handles_labels()
#           ax.legend(handles,labels,loc='upper left',fontsize=15)
#       #   rsquared = round((linreg.rvalue*linreg.rvalue),2)
#           ax.text(-3.8, 3, f'{group}',fontsize=15)
#           plt.savefig(f'{figurePath}/iPRSHigh.QQColorPlot.png')
#           ax.clear()
#           cases['color'] = 'black'
#           
#       if group == 'high risk':
#           fig,ax = plt.subplots(figsize=(10,10))
#           cases.loc[((cases['bin_Main'] == 9) & (cases['epi_helps'] == 0)),'color'] == 'rgb(211,0,63,0)'
#           cases.loc[((cases['bin_Main'] == 10) & (cases['epi_helps'] == 0)),'color'] == 'rgb(201,0,22,0)'
#           cases.loc[((cases['bin_EpiMain'] == 9) & (cases['epi_helps'] == 1)),'color'] == 'rgb(0, 0, 205,.6)'
#           cases.loc[((cases['bin_EpiMain'] == 10) & (cases['epi_helps'] == 1)),'color'] == 'rgb(0, 0, 128,.6)'
#           model_type = f'HighRiskPRSANDiPRS.{suffix}'
#       #   ax.scatter(x=controls['scaled_prs'].values,y=controls['scaled_prs_main'].values,marker=".",c='#776cb1',s=40,alpha=.6,label='Controls')
#           ax.scatter(x=controls['scaled_prs'].values,y=controls['scaled_prs_main'].values,marker=".",c=controls['color'].values,s=40,alpha=.6,label='Controls')
#       #   ax.scatter(x=epiCases['scaled_prs'].values,y=mainCases['scaled_prs'].values,marker="+",c='black',s=40,alpha=.9,label='Cases')
#           ax.scatter(x=cases['scaled_prs'].values,y=cases['scaled_prs'].values,marker="+",c=cases['color'].values,s=40,alpha=.9,label='Cases')
#       #   ax.plot([main['scaled_prs'].min(),main['scaled_prs'].max()],[linreg.intercept+linreg.slope*epi['scaled_prs'].min(),linreg.intercept+linreg.slope*epi['scaled_prs'].max()])
#           
#           ax.set_xlabel('epi')
#           ax.set_ylabel('main only')
#           ax.grid(True) 
#           ax.set_title(f'Standardized High Risk PRS & iPRS Cases : Main V Epi')
#           handles,labels = ax.get_legend_handles_labels()
#           ax.legend(handles,labels,loc='upper left',fontsize=15)
#       #   rsquared = round((linreg.rvalue*linreg.rvalue),2)
#           ax.text(-3.8, 3, f'{group}',fontsize=15)
#           plt.savefig(f'{figurePath}/iPRS_PRSHigh.QQColorPlot.png')
#           ax.clear()
#           cases['color'] = 'black'
#           
#       if group == 'all':
#           fig,ax = plt.subplots(figsize=(10,10))
#           cases.loc[((cases['bin_Main'] == 9) & (cases['epi_helps'] == 0)),'color'] == 'rgb(211,0,63,0)'
#           cases.loc[((cases['bin_Main'] == 10) & (cases['epi_helps'] == 0)),'color'] == 'rgb(201,0,22,0)'
#           cases.loc[((cases['bin_EpiMain'] == 9) & (cases['epi_helps'] == 1)),'color'] == 'rgb(0, 0, 205,.6)'
#           cases.loc[((cases['bin_EpiMain'] == 10) & (cases['epi_helps'] == 1)),'color'] == 'rgb(0, 0, 128,.6)'
#           cases.loc[((cases['bin_Main'] < 9) & (cases['epi_helps'] == 0)),'color'] == 'rgb(240,128,128)'
#           cases.loc[((cases['bin_EpiMain'] < 9) & (cases['epi_helps'] == 1)),'color'] == 'rgb(176, 224, 230,.6)'
#           model_type = f'HighANDLowRiskPRSANDiPRS.{suffix}'
#       #   ax.scatter(x=controls['scaled_prs'].values,y=controls['scaled_prs_main'].values,marker=".",c='#776cb1',s=40,alpha=.6,label='Controls')
#           ax.scatter(x=controls['scaled_prs'].values,y=controls['scaled_prs_main'].values,marker=".",c=controls['color'].values,s=40,alpha=.6,label='Controls')
#       #   ax.scatter(x=epiCases['scaled_prs'].values,y=mainCases['scaled_prs'].values,marker="+",c='black',s=40,alpha=.9,label='Cases')
#           ax.scatter(x=cases['scaled_prs'].values,y=cases['scaled_prs'].values,marker="+",c=cases['color'].values,s=40,alpha=.9,label='Cases')
#       #   ax.plot([main['scaled_prs'].min(),main['scaled_prs'].max()],[linreg.intercept+linreg.slope*epi['scaled_prs'].min(),linreg.intercept+linreg.slope*epi['scaled_prs'].max()])
#           
#           ax.set_xlabel('epi')
#           ax.set_ylabel('main only')
#           ax.grid(True) 
#           ax.set_title(f'Standardized High & Low Risk PRS & iPRS Cases : Main V Epi')
#           handles,labels = ax.get_legend_handles_labels()
#           ax.legend(handles,labels,loc='upper left',fontsize=15)
#       #   rsquared = round((linreg.rvalue*linreg.rvalue),2)
#           ax.text(-3.8, 3, f'{group}',fontsize=15)
#           plt.savefig(f'{figurePath}/iPRS_PRSHighANDLow.QQColorPlot.png')
#           ax.clear()
#           cases['color'] = 'black'
#           
#       if group == 'cases predicted low risk PRS':
#           fig,ax = plt.subplots(figsize=(10,10))
#           cases.loc[((cases['bin_Main'] < 9) & (cases['epi_helps'] == 0)),'color'] == 'rgb(240,128,128)'
#           model_type = f'LowRiskPRS.{suffix}'
#       #   ax.scatter(x=controls['scaled_prs'].values,y=controls['scaled_prs_main'].values,marker=".",c='#776cb1',s=40,alpha=.6,label='Controls')
#           ax.scatter(x=controls['scaled_prs'].values,y=controls['scaled_prs_main'].values,marker=".",c=controls['color'].values,s=40,alpha=.6,label='Controls')
#       #   ax.scatter(x=epiCases['scaled_prs'].values,y=mainCases['scaled_prs'].values,marker="+",c='black',s=40,alpha=.9,label='Cases')
#           ax.scatter(x=cases['scaled_prs'].values,y=cases['scaled_prs'].values,marker="+",c=cases['color'].values,s=40,alpha=.9,label='Cases')
#       #   ax.plot([main['scaled_prs'].min(),main['scaled_prs'].max()],[linreg.intercept+linreg.slope*epi['scaled_prs'].min(),linreg.intercept+linreg.slope*epi['scaled_prs'].max()])
#           
#           ax.set_xlabel('epi')
#           ax.set_ylabel('main only')
#           ax.grid(True) 
#           ax.set_title(f'Standardized Low Risk PRS Cases : Main V Epi')
#           handles,labels = ax.get_legend_handles_labels()
#           ax.legend(handles,labels,loc='upper left',fontsize=15)
#       #   rsquared = round((linreg.rvalue*linreg.rvalue),2)
#           ax.text(-3.8, 3, f'{group}',fontsize=15)
#           plt.savefig(f'{figurePath}/PRSLow.QQColorPlot.png')
#           ax.clear()
#           cases['color'] = 'black'
#           
#       if group == 'cases predicted low risk iPRS':
#           fig,ax = plt.subplots(figsize=(10,10))
#           cases.loc[((cases['bin_EpiMain'] < 9) & (cases['epi_helps'] == 1)),'color'] == 'rgb(176, 224, 230,.6)'
#           model_type = f'LowRiskiPRS.{suffix}'
#       #   ax.scatter(x=controls['scaled_prs'].values,y=controls['scaled_prs_main'].values,marker=".",c='#776cb1',s=40,alpha=.6,label='Controls')
#           ax.scatter(x=controls['scaled_prs'].values,y=controls['scaled_prs_main'].values,marker=".",c=controls['color'].values,s=40,alpha=.6,label='Controls')
#       #   ax.scatter(x=epiCases['scaled_prs'].values,y=mainCases['scaled_prs'].values,marker="+",c='black',s=40,alpha=.9,label='Cases')
#           ax.scatter(x=cases['scaled_prs'].values,y=cases['scaled_prs'].values,marker="+",c=cases['color'].values,s=40,alpha=.9,label='Cases')
#       #   ax.plot([main['scaled_prs'].min(),main['scaled_prs'].max()],[linreg.intercept+linreg.slope*epi['scaled_prs'].min(),linreg.intercept+linreg.slope*epi['scaled_prs'].max()])
#           
#           ax.set_xlabel('epi')
#           ax.set_ylabel('main only')
#           ax.grid(True) 
#           ax.set_title(f'Standardized Low Risk iPRS Cases : Main V Epi')
#           handles,labels = ax.get_legend_handles_labels()
#           ax.legend(handles,labels,loc='upper left',fontsize=15)
#       #   rsquared = round((linreg.rvalue*linreg.rvalue),2)
#           ax.text(-3.8, 3, f'{group}',fontsize=15)
#           plt.savefig(f'{figurePath}/iPRSLow.QQColorPlot.png')
#           ax.clear()
#           cases['color'] = 'black'
#           
#       if group == 'low risk':
#           fig,ax = plt.subplots(figsize=(10,10))
#           cases.loc[((cases['bin_Main'] < 9) & (cases['epi_helps'] == 0)),'color'] == 'rgb(240,128,128)'
#           cases.loc[((cases['bin_EpiMain'] < 9) & (cases['epi_helps'] == 1)),'color'] == 'rgb(176, 224, 230,.6)'
#           model_type = f'LowRiskPRSANDiPRS.{suffix}'
#       #   ax.scatter(x=controls['scaled_prs'].values,y=controls['scaled_prs_main'].values,marker=".",c='#776cb1',s=40,alpha=.6,label='Controls')
#           ax.scatter(x=controls['scaled_prs'].values,y=controls['scaled_prs_main'].values,marker=".",c=controls['color'].values,s=40,alpha=.6,label='Controls')
#       #   ax.scatter(x=epiCases['scaled_prs'].values,y=mainCases['scaled_prs'].values,marker="+",c='black',s=40,alpha=.9,label='Cases')
#           ax.scatter(x=cases['scaled_prs'].values,y=cases['scaled_prs'].values,marker="+",c=cases['color'].values,s=40,alpha=.9,label='Cases')
#       #   ax.plot([main['scaled_prs'].min(),main['scaled_prs'].max()],[linreg.intercept+linreg.slope*epi['scaled_prs'].min(),linreg.intercept+linreg.slope*epi['scaled_prs'].max()])
#           
#           ax.set_xlabel('epi')
#           ax.set_ylabel('main only')
#           ax.grid(True) 
#           ax.set_title(f'Standardized Low Risk PRS & iPRS Cases : Main V Epi')
#           handles,labels = ax.get_legend_handles_labels()
#           ax.legend(handles,labels,loc='upper left',fontsize=15)
#       #   rsquared = round((linreg.rvalue*linreg.rvalue),2)
#           ax.text(-3.8, 3, f'{group}',fontsize=15)
#           plt.savefig(f'{figurePath}/iPRS_PRSLow.QQColorPlot.png')
#           ax.clear()
#           cases['color'] = 'black'
    
    
#def create_prs_iprs_qq_plots(combinedPRS,figurePath,suffix):
#   '''filepaths to prs scores using columns = [IID,prs,scaled_prs,PHENOTYPE]'''
#   
#   
#   epiMain = combinedPRS[combinedPRS['model'] == 'epi+main']
#   main = combinedPRS[combinedPRS['model'] == 'main']
#   epi = combinedPRS[combinedPRS['model'] == 'epi']
#   
#   mainCases = combinedPRS[(combinedPRS['model'] == 'main') & (combinedPRS['PHENOTYPE'] == 2)]
#   mainControls = combinedPRS[(combinedPRS['model'] == 'main') & (combinedPRS['PHENOTYPE'] == 1)]
#   epiCases = combinedPRS[(combinedPRS['model'] == 'epi') & (combinedPRS['PHENOTYPE'] == 2)]
#   epiControls = combinedPRS[(combinedPRS['model'] == 'epi') & (combinedPRS['PHENOTYPE'] == 1)]
#   epiMainCases = combinedPRS[(combinedPRS['model'] == 'epi+main') & (combinedPRS['PHENOTYPE'] == 2)]
#   epiMainControls = combinedPRS[(combinedPRS['model'] == 'epi+main') & (combinedPRS['PHENOTYPE'] == 1)]
#   
#   suffix = suffix
#   
#   ############ Main v Epi+Main ##############################
#   
#   fig,ax = plt.subplots(figsize=(10,10))
#   linreg = sp.stats.linregress(epiMain['scaled_prs'].values, main['scaled_prs'].values)
#   model_type = f'mainVmain+epi.{suffix}'
#   ax.scatter(x=epiMainControls['scaled_prs'].values,y=mainControls['scaled_prs'].values,marker=".",c='#776cb1',s=40,alpha=.6,label='Controls')
#   ax.scatter(x=epiMainCases['scaled_prs'].values,y=mainCases['scaled_prs'].values,marker="+",c='black',s=40,alpha=.9,label='Cases')
##   ax.plot([epiMain['scaled_prs'].min(),epiMain['scaled_prs'].max()],[linreg.intercept+linreg.slope*main['scaled_prs'].min(),linreg.intercept+linreg.slope*main['scaled_prs'].max()])
#   
#   ax.set_xlabel('epi + main')
#   ax.set_ylabel('main only')
#   ax.grid(True) 
#   ax.set_title(f'Standardized Risk Scores : Main V Main+Epi')
#   rsquared = round((linreg.rvalue*linreg.rvalue),2)
#   ax.text(-3.8, 3, f'r squared = {rsquared}',fontsize=15)
#   handles,labels = ax.get_legend_handles_labels()
#   ax.legend(handles,labels,loc='upper left',fontsize=15)
#
#   plt.savefig(f'{figurePath}/{model_type}.QQWithR2.png')
#   ax.clear()
#   
#   ############ Main v Epi ##############################
#   
#   fig,ax = plt.subplots(figsize=(10,10))
#   linreg = sp.stats.linregress(main['scaled_prs'].values, epi['scaled_prs'].values)
#   
#   model_type = f'{suffix}'
#   ax.scatter(x=epiControls['scaled_prs'].values,y=mainControls['scaled_prs'].values,marker=".",c='#776cb1',s=40,alpha=.6,label='Controls')
#   ax.scatter(x=epiCases['scaled_prs'].values,y=mainCases['scaled_prs'].values,marker="+",c='black',s=40,alpha=.9,label='Cases')
##   ax.plot([main['scaled_prs'].min(),main['scaled_prs'].max()],[linreg.intercept+linreg.slope*epi['scaled_prs'].min(),linreg.intercept+linreg.slope*epi['scaled_prs'].max()])
#   
#   ax.set_xlabel('epi')
#   ax.set_ylabel('main only')
#   ax.grid(True) 
#   ax.set_title(f'Standardized Risk Scores : Main V Epi')
#   handles,labels = ax.get_legend_handles_labels()
#   ax.legend(handles,labels,loc='upper left',fontsize=15)
#   rsquared = round((linreg.rvalue*linreg.rvalue),2)
#   ax.text(-3.8, 3, f'r squared = {rsquared}',fontsize=15)
#   plt.savefig(f'{figurePath}/{model_type}.QQWithR2.png')
#   
#   ax.clear()
#   
#   ############ Epi v Epi+Main ##############################
#   
#   fig,ax = plt.subplots(figsize=(10,10))
#   
#   model_type = f'epiVMainepi.{suffix}'
#   linreg = sp.stats.linregress(epiMain['scaled_prs'].values, epi['scaled_prs'].values)
#   
#   ax.scatter(x=epiMainControls['scaled_prs'].values,y=epiControls['scaled_prs'].values,marker=".",c='#776cb1',s=40,alpha=.6,label='Controls')
#   ax.scatter(x=epiMainCases['scaled_prs'].values,y=epiCases['scaled_prs'].values,marker="+",c='black',s=40,alpha=.9,label='Cases')
##   ax.plot([epiMain['scaled_prs'].min(),epiMain['scaled_prs'].max()],[linreg.intercept+linreg.slope*epi['scaled_prs'].min(),linreg.intercept+linreg.slope*epi['scaled_prs'].max()])
#   
#   ax.set_xlabel('epi + main')
#   ax.set_ylabel('epi only')
#   ax.grid(True) 
#   ax.set_title(f'Standardized Risk Scores : Epi V Epi + Main')
#   handles,labels = ax.get_legend_handles_labels()
#   labels.append('r squared')
#   ax.legend(handles,labels,loc='upper left',fontsize=15)
#   rsquared = round((linreg.rvalue*linreg.rvalue),2)
#   ax.text(-3.8, 3, f'r squared = {rsquared}',fontsize=15)
#   plt.savefig(f'{figurePath}/{model_type}.QQWithR2.png')
#   
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

def create_density_plot(df,model_type,figurePath):
# def create_density_plot(df,model_type,figurePath,threshold):
    try:
        meanDiff = round(df.groupby(['PHENOTYPE']).mean()['scaled_prs'].diff().loc[2],2)
        print('mean diff between cases and controls is : ',meanDiff)
        meanControl =  df.groupby(['PHENOTYPE']).mean()['scaled_prs'].loc[1]
        meanCase =  df.groupby(['PHENOTYPE']).mean()['scaled_prs'].loc[2]

        title = f'Combined PRS for {model_type} - mean diff = {meanDiff}'
        ###split the data into groups based on types
        g = df.groupby('PHENOTYPE')
    
        ## From here things change as I make use of the seaborn library
        yes = g.get_group(2)
        no = g.get_group(1)
    
        fig, ax = plt.subplots(figsize=(10,10))
    
        #
        ax = sns.kdeplot(data=yes['scaled_prs'], label='cases', ax=ax,color='red',shade=True)
        ax = sns.kdeplot(data=no['scaled_prs'], label='controls', ax=ax,color='blue',shade=True)
    
        ax.legend(['cases','controls'])
    
        ax.set_title(title)
        fig.savefig(f'{figurePath}/{model_type}.densityPlot.png')
        plt.close()
    except KeyError:
        pass

def create_box_plot(prsCopy,model_type,figurePath):
    try:
        meanDiff = prsCopy.groupby(['PHENOTYPE']).mean()['scaled_prs'].diff().loc[2]
        print(meanDiff)
        ax = prsCopy.boxplot(column='scaled_prs',by='PHENOTYPE',figsize=(10,10))
        ax.plot()
        title = f'prs mean diff for {model_type} = {meanDiff}'
        plt.title(title)
        plt.savefig(f'{figurePath}/{model_type}.boxplot.png')
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
    ####################   GET CARDIO AND GENO SEPARATE AND TOGETHER  ###########
    
    cardio_models = df[df['model'].str.contains('cardio')]
    geno_models = df[~df['model'].str.contains('cardio')]
    
    
    ####################################################################
    #                            Prevalence Plot                       #
    ####################################################################

    fig,ax = plt.subplots(figsize=(10,10))

    

    ###########################  PREVALENCE PLOT ################################
    title = 'Prevalence Across Centiles optimized'
    
    for model in cardio_models['model'].unique():
        plotDf = cardio_models[cardio_models['model'] == model]
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
    
    plt.savefig(f'{figurePath}/PRScrCardioFeatures{image_str}.scatter.png')
    
    ##########################  PLOT GENO NO CARDIO  ###############################
    fig,ax = plt.subplots(figsize=(10,10))
    title = 'PRScr geno models across deciles'
    
    for model in geno_models['model'].unique():
        plotDf = geno_models[geno_models['model'] == model]
        marker = plotDf['marker'].tolist()[0]
        color = plotDf['color'].tolist()[0]
        ax.scatter(x=plotDf['prevalence'],y=np.array(plotDf['bin'])-1,color=color,s=200,marker=marker,alpha=.7)
        ax.plot(plotDf['prevalence'],np.array(plotDf['bin'])-1,color=color,alpha=.4,linestyle='dashed')
        
    ax.set_ylabel('deciles')
    ax.set_xlabel('prevalence')
    ax.set_title(title)
    #	xlabels = [f"{i}\n{j}" for i,j in enumerate(bin_limits,start=1)]
    #	xlabels = df['bin_limits'].unique().tolist()
    yticks = list(range(0,10))
    ylabels = [str(i) for i in range(1,11)]
    ax.set_yticks(yticks,labels=ylabels)
    
    plt.savefig(f'{figurePath}/PRScrGenoFeatures{image_str}.scatter.png')
    
    
    ##########################  PLOT BOTH GENO AND CARDIO TOGETHER  ###############################
    
    title = 'PRScr models across deciles with ExGxG'
    fig,ax = plt.subplots(figsize=(10,10))
    #models are 3 PRScr models, cardio_main, cardio
    models = [col for col in df['model'].unique() if 'PRScr' in col] + ['cardio']
    
#   for model in df['model'].unique():
    for model in models:
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
    
    plt.savefig(f'{figurePath}/PRScrAll{image_str}.scatter.png')
    
    
    ##########################  PLOT ALL WITH PRScr geno and PRScr_epi  ###############################
    
    fig,ax = plt.subplots(figsize=(10,10))
    data_types_to_plot = [col for col in df['model'].unique() if '+cardio' not in col]
    for model in data_types_to_plot:
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
    
    plt.savefig(f'{figurePath}/PRScrGenoWithCardio{image_str}.scatter.png')
    
    
    
    #def draw_pie(r1,r2,xpos,ypos,size,ax=None):
    #	x = [0] + np.cos(np.linspace(0, 2 * np.pi * r1, 10)).tolist()
    #	y = [0] + np.sin(np.linspace(0, 2 * np.pi * r1, 10)).tolist()
    #	xy1 = np.column_stack([x, y])
    #	
    #	x = [0] + np.cos(np.linspace(2 * np.pi * r1, 2 * np.pi, 10)).tolist()
    #	y = [0] + np.sin(np.linspace(2 * np.pi * r1, 2 * np.pi, 10)).tolist()
    #	xy2 = np.column_stack([x, y])
    #	ax.scatter([xpos],[ypos],marker=xy1,s=size,facecolor='blue',alpha=.7)
    #	ax.scatter([xpos],[ypos],marker=xy2,s=size,facecolor='red',alpha=.7)
    #	return(ax)
    
def plot_important_features_modelling(pheno):
    filePath = f'/Users/kerimulterer/ukbiobank/{pheno}/tanigawaSet/prs/reducedSHAP'

    # get the prevalence plot for quintiles and cardio metabolic features
    df = pd.read_csv(f'/Users/kerimulterer/ukbiobank/{pheno}/tanigawaSet/prs/reducedSHAP/holdout/CohortAssignedCosine.cardioOnlyQuintiles.prevalenceCombinedPRSGroupsWithcardio_main.holdout.csv')
    figurePath = f'/Users/kerimulterer/ukbiobank/{pheno}/tanigawaSet/figures/validation/sectionModels/paperFigures/filteredMainAllEpi/forestPlots'
    df = pd.read_csv(f'{filePath}/featureImportance/groupedCohortFeatureImportanceAcrossMethods.csv')
    
    #get the features that have at least one cohort with low pvalue and support == True
    importantFeaturesDf = df[(df['model'] == 'one_way_anova') & (df['fdr_bh_rejected'] == True) & (df['support'] == True)]
    importantFeatures = importantFeaturesDf['feature'].tolist()
    importantFeatures = list(set(importantFeatures))
    
    dfPlot = df[df['feature'].isin(importantFeatures)]
    
    dfPlot = dfPlot[dfPlot['model'] == 'shapley_values'][['cohort','feature','feature_importance']]
    
    #cohort with the highest shap value then filter the features that don't work for that cohort
    #i.e. epi interactions in main cohort
    importantCohorts = dfPlot.pivot(index='feature', columns='cohort', values='feature_importance').reset_index()
    
    importantCohorts['cohort'] = importantCohorts[['main','epi','epiMain']].idxmax(axis=1)
    
    filteredFeaturesDf = importantCohorts[(importantCohorts['cohort'] == 'epiMain') | ((importantCohorts['cohort'] == 'main') & (~importantCohorts['feature'].str.contains(','))) | ((importantCohorts['cohort'] == 'epi') & (importantCohorts['feature'].str.contains(',')))]

    
    filteredFeatures = filteredFeaturesDf['feature'].tolist()
    
    dfPlot = dfPlot[dfPlot['feature'].isin(filteredFeatures)]
    finalFeatures = dfPlot.merge(importantCohorts[['feature','cohort']],on=['feature','cohort'],how='right')
    finalFeatures.dropna(inplace=True)
    finalFeatures2 = finalFeatures[['feature','cohort','feature_importance']].merge(df[df['model'] == 'shapley_values'],on=['feature','cohort','feature_importance'],how='inner')

    finalFeatures2.to_csv(f'{filePath}/featureImportance/groupedCohortFeatureImportanceAcrossMethods.Filtered.csv',index=False)
    create_important_feature_forest_plot(dfPlot,filteredFeaturesDf,figurePath)
    
    

    def main(pheno):
        pass
    
    
if __name__ == '__main__':
    
#   pheno = sys.argv[1]
    pheno = 'type2Diabetes'
    
    main(pheno)
    
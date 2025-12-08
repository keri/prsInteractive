#!/usr/bin/env python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from helper.create_migration_plot import *

COHORT_COLORS = {
    'main': '#E69F00',      # Orange
    'epi': '#56B4E9',    # Sky blue
    'epi+main': '#CC79A7',    # Pinkish purple
    'cardio': '#009E73',   # Bluish green
    'all': '#F0E442'   # Bluish green
    
}

def calculate_prevalance(df,prscr_group):
    dfCopy = df.copy()
    dfCopy = dfCopy[dfCopy['PHENOTYPE'] == 2]
    
    #get the total count for cases in high and low group
    dfCases = dfCopy.groupby(['quintile']).count().rename(columns={'IID':'total_cases'})[['total_cases']]
    
    #count the number of cases in high and environmental group
    n_prs_cases_low = dfCopy[(dfCopy['quintile'] == 1) & (dfCopy[f'decile_bin_{prscr_group}'] > 8)].shape[0]
    n_prs_cases_high = dfCopy[(dfCopy['quintile'] == 5) & (dfCopy[f'decile_bin_{prscr_group}'] > 8)].shape[0]
    
    #get the count in each decile
    dfTotal = df.groupby(['quintile']).count().rename(columns={'IID':'total'})[['total']]
#	dfTotal.rename(columns={'color':'total'},inplace=True)[['total']]
    #get the prevalence in each decile
    dfTotal['prevalence'] = round((dfCases['total_cases'] / dfTotal['total']),2)
    dfTotal['n_cases'] = dfCases['total_cases']
    dfTotal.reset_index(inplace=True)
    dfTotal.loc[dfTotal['quintile'] == 1,'n_cases_env'] =  dfTotal['n_cases'] - n_prs_cases_low
    dfTotal.loc[dfTotal['quintile'] == 1,'n_cases_prs'] =  n_prs_cases_low
    
    dfTotal.loc[dfTotal['quintile'] == 5,'n_cases_env'] =  dfTotal['n_cases'] - n_prs_cases_high
    dfTotal.loc[dfTotal['quintile'] == 5,'n_cases_prs'] =  n_prs_cases_high
    

    dfTotal = dfTotal[['prevalence','quintile','n_cases','n_cases_env','n_cases_prs']]
    return(dfTotal)

def create_graph(nodeDf,linkDf,mainSource,quintile,figPath):
#	if mainSource == 'PRStoiPRS':

#   node_position_y = [0,.1,.2,.3,.4,.5,.6,.7,.8,.9,0,.1,.2,.3,.4,.5,.6,.7,.8,.9]
    node_positions_x = [.0001,.0001,.0001,.0001,.0001,.0001,.0001,.0001,.0001,.0001,.99,.99,.99,.99,.99,.99,.99,.99,.99,.99]
    #get the node sizes for each node
    node_sizes = linkDf.groupby(['source']).sum()['IID'].tolist() + linkDf.groupby(['target']).sum()['IID'].tolist()
    
    # Normalize the node sizes for relative positioning (y axis adjustment)
    total_flow = sum(linkDf['IID'].tolist())
    
#   node_heights = [size / total_flow for size in node_sizes]  # Node height based on flow
    
    node_positions_y = [
        0.00001,  # Source A (normalized placement)
        node_sizes[0]/total_flow,  # Source B placed after Source A's height
        sum(node_sizes[:2])/total_flow,  # Source C after Source B's height
        sum(node_sizes[:3])/total_flow,
        sum(node_sizes[:4])/total_flow,  # Target B after Target A's height
        sum(node_sizes[:5])/total_flow,  # Target C after Target B's height
        sum(node_sizes[:6])/total_flow,
        sum(node_sizes[:7])/total_flow,
        sum(node_sizes[:8])/total_flow,
        sum(node_sizes[:9])/total_flow,
        0.00001,  # Target A (normalized placement)
        node_sizes[10]/total_flow,  # Source B placed after Source A's height
        sum(node_sizes[10:12])/total_flow,  # Source C after Source B's height
        sum(node_sizes[10:13])/total_flow,
        sum(node_sizes[10:14])/total_flow,  # Target B after Target A's height
        sum(node_sizes[10:15])/total_flow,  # Target C after Target B's height
        sum(node_sizes[10:16])/total_flow,
        sum(node_sizes[10:17])/total_flow,
        sum(node_sizes[10:18])/total_flow,
        sum(node_sizes[10:19])/total_flow
    ]
    
    
    linkDf.sort_values(['source'],inplace=True)
#	if epi_contribution == 'decreases':
#		nodeDf.sort_values(['model','source'],inplace=True)
    graphData = {}
    binSums = pd.concat([linkDf.groupby(['source']).sum()[['IID']].reset_index(), linkDf.groupby(['target']).sum()[['IID']].reset_index()],ignore_index=True)
    

    node = {
        'pad' : 2,
        'thickness' : 15,
        'line' : dict(color = "white", width = 0),
#		"label" :  [x.split(' ')[-1] for x in nodeDf['node_label']],
#		"label" :  nodeDf['node_label'].tolist(),
        'color' :  nodeDf['node_color'].tolist(),
        'customdata' : np.stack((nodeDf['node_label'], binSums['IID']), axis=-1),
        'y': node_positions_y,
        'x': node_positions_x,
        'hovertemplate' : ' Total # of people in node %{label} = %{customdata[1]}<extra></extra>'
        
    }
    
    graphData['node'] = node
    
    link = {
        "source" :  linkDf['source'].tolist(),
        "target" :  linkDf['target'].tolist(),
        "value" :  linkDf['value'].tolist(),
        'color' :  linkDf['color'].tolist(),
        'customdata' : linkDf['IID'].tolist(),
        'hovertemplate' : '%{customdata} people from node %{source.customdata[0]} to %{target.customdata[0]}<extra></extra>'
    }
        
    graphData['link'] = link

    
    layout = go.Layout(
        autosize=False,
        width=1000,
        height=1000
    )
    
    fig = go.Figure(data=[go.Sankey(
        arrangement = "snap",
        valueformat = ".0f",
        valuesuffix = "",
        # Define nodes
        node = graphData['node'],
        # Add links
        link = graphData['link'])],layout=layout)
    
#   fig = go.Figure(go.Sankey(
#       valueformat = ".0f",
#       # Define nodes
#       node = graphData['node'],
#       # Add links
#       link = graphData['link']))
    
    fig.data[0].node.x = node_positions_x
    # Update layout and display the plot
    fig.update_layout(font_size=10)

    
    pio.write_image(fig,f'{figPath}/SankeyDiagramBinMigrationQuintile{quintile}.{mainSource}.png')
    fig.write_html(f'{figPath}/SankeyDiagramBinMigrationQuintile{quintile}.{mainSource}.html')


def main_clinical_features(df,env_factor,figPath,prscr_group):
    '''create migration plot from individuals stratified by quintiles 1-4 (denoted as quintile 1) , 5 and 0 (combined quintile)
    
    input :

    df : cases that are > 8 decile in PRScr_all
        columns = [IID,bin_main, bin_epi, bin_epi+main,bin_cardio,bin_PRScr_all,quintile]
    
    env_factor : str
    
    '''
    node_dict = {'decile1':'rgba(255, 0, 0, 0)','decile2':'rgba(255, 0, 0, 0)','decile3':'rgba(255, 0, 0, 0)','decile4':'rgba(255, 0, 0, 0)','decile5':'rgba(255, 0, 0, 0)','decile6':'rgba(255, 0, 0, 0)','decile7':'rgba(255, 0, 0, 0)','decile8':'rgba(255, 0, 0, 0)','decile9':'rgb(255, 0, 255)','decile10':'rgb(255, 0, 255)'}
    
    
    
    for q in [1,5]:
        
        individualCases = df[df['quintile']==q]
        
        finalLinkDf,finalNodeDf = create_graph_data(individualCases,prscr_group)
        
        create_graph(finalNodeDf, finalLinkDf, f'{prscr_group}.{env_factor}',q,figPath)
        
    #get migration for all quintiles combined (overall migration)

    finalLinkDf,finalNodeDf = create_graph_data(df,prscr_group)
    
    create_graph(finalNodeDf, finalLinkDf, f'{prscr_group}.Overall',0,figPath)


def main(phenoData,covarFile,envFile,use_epi_main=False):

    figPath = f'{phenoData}/figures'
    scoresPath = f'{phenoData}/scores'
    

    #columns =  ['prevalence', 'bin', 'total', 'n_cases', 'model', 'method', 'quintile','env_factor', 'marker', 'color']
#   prevalenceFile = f'{dataPath}/holdout/CohortAssignedMax.EnvQuintiles.prevalenceCombinedPRSGroupsWithcardio_main1Kbins.holdout.csv'
    
    combinedPRSFile = f'{scoresPath}/combinedPRSGroups.holdout.csv'
#   ORFile = f'{dataPath}/holdout/CohortAssignedMax.{env_factor}Quintiles.ORPRSGroups.holdout.csv'
    #create the PRScr_mix column


    
    color_dict = {'main':'#C90016','epi':'#19f0f2','epi+main':'#39c700','cardio':'#f2de19','PRScr_geno':'#f73bd7','PRScr_all':'#9c2287','cardio_main':'#c4771f','PRScr_mix':'#f542bc','all':'#f69697'}
    marker_dict = {'main':'D','epi':'o','epi+main':'X','cardio':'v','PRScr_geno':'8','PRScr_all':'8','cardio_main':'v','PRScr_mix':'8'}
    
    
    combinedPRS = pd.read_csv(combinedPRSFile)
    
    if use_epi_main:
        cohorts = ['decile_bin_main','decile_bin_epi+main','decile_bin_epi','decile_bin_cardio']
        colors = ['#f2de19','#19f0f2','#C90016','#39c700']
    else:
        cohorts = ['decile_bin_main','decile_bin_epi','decile_bin_cardio']
        colors = ['#f2de19','#19f0f2','#C90016']
    
    combinedPRS['prscr_mix'] = combinedPRS[cohorts].idxmax(axis=1) #cohort = row[columns].idxmax()
    combinedPRS['prscr_mix'] = combinedPRS['PRScr_mix'].apply(lambda x : x.replace('decile_bin_',''))

#   for prscr in ['prscr_all','prscr_mix']:
    for prscr in ['prscr_mix']:
            
        #create a copy to be used for 0 quintile
        #before data is filtered downstream
        combinedNoEnv = combinedPRS.copy()
        prevalenceDf = pd.DataFrame()
        
#       prevalenceDfTemp = prevalenceDf[(prevalenceDf['bin'] > 8) & (prevalenceDf['model'].isin([prscr]))]
#       #filter the PRScr values using the correct method
#       prevalenceDfTemp = prevalenceDfTemp[(prevalenceDfTemp['method'] == 'PRScr_methods_over_1K_bins_using_scaled_prs')]
        if pheno == 'type2Diabetes':
            env_factors = ['BMI','HbA1c','Glucose']
            fig,(ax1,ax2,ax3) = plt.subplots(1,3,figsize=(10,10),squeeze=True,sharey=True)
            axs = (ax1,ax2,ax3)
            ylimBreakdown = 1900
            ylimPrev = .50
        else:
            env_factors = ['metabolic rate','urea','HC']
            fig,(ax1,ax2,ax3) = plt.subplots(1,3,figsize=(10,10),squeeze=True,sharey=True)
            axs = (ax1,ax2,ax3)
            ylimBreakdown = 300
            ylimPrev = .055
    
        ######################  CREATE PREVALENCE BAR CHART ##########################

        columns_to_use = cohorts + ['IID','PRScr_all','scaled_prs_PRScr_all','decile_bin_PRScr_all','PRScr_mix','scaled_prs_PRScr_mix','decile_bin_PRScr_mix','PHENOTYPE']
        
        #################### for PRScr geno and PRScr all get the combined number of people in deciles 9 and 10 within each PRScr  #########
        stats = pd.DataFrame()
        binnedCountFlagged = pd.DataFrame()

        
        
        for i, ax in enumerate(axs):
            print(i)
            env_factor = env_factors[i]
            print(env_factor)
            binnedDf = combinedPRS[columns_to_use+[f'{env_factor}_quintile',env_factor]]
            binnedDf.rename(columns={f'{env_factor}_quintile':'quintile'},inplace=True)
            
            
    #       binnedDf = pd.read_csv(f'{dataPath}/holdout/CohortAssignedMax.{env_factor}Quintiles.combinedPRSGroupsWithcardio_main.holdout.csv',usecols=columns_to_use+[env_factor])
            
            if env_factor == 'BMI':
                binnedDf['quintile'] = binnedDf[env_factor].apply(lambda x: 1 if x < 25  else 5)
            elif env_factor == 'Glucose':
                binnedDf['quintile'] = binnedDf[env_factor].apply(lambda x: 1 if x < 12 else 5)
            elif env_factor == 'HbA1c': #env factor = HbA1c
                binnedDf['quintile'] = binnedDf[env_factor].apply(lambda x: 1 if x < 41 else 5)
            else:    
                #combine quintiles 1-4 and by assigning to quintile 1
                binnedDf['quintile'] = binnedDf['quintile'].apply(lambda x: 1 if x < 5 else x)

            prevalence = calculate_prevalance(binnedDf,prscr)
            prevalence['env_factor'] = env_factor
            prevalenceDf = pd.concat([prevalence,prevalenceDf],ignore_index=True)
            
            env_factor_stats = binnedDf[[env_factor,'quintile']].groupby(['quintile']).mean()
            env_factor_stats.columns = ['mean']
            env_factor_stats['min'] = binnedDf[[env_factor,'quintile']].groupby(['quintile']).min()
            env_factor_stats['max'] = binnedDf[[env_factor,'quintile']].groupby(['quintile']).max()
            env_factor_stats['std'] = binnedDf[[env_factor,'quintile']].groupby(['quintile']).std()
            env_factor_stats['env_factor'] = env_factor
            stats = pd.concat([env_factor_stats.reset_index(),stats],ignore_index=True)
            
            
            #use cases only for count
            binnedDf = binnedDf[binnedDf['PHENOTYPE'] == 2]
            
            binnedCount = pd.DataFrame()
#           binnedCount2 = pd.DataFrame()
    
            
    #       for prscr_tuple in [('PRScr_all',['bin_main','bin_epi','bin_epi+main','bin_cardio']),('PRScr_geno',['bin_main','bin_epi','bin_epi+main'])]:
    #       for prscr_tuple in [('PRScr_all',['bin_main','bin_epi','bin_epi+main','bin_cardio'])]:
#           for prscr_tuple in [('PRScr_all',cohorts),('PRScr_mix',cohorts)]:
#           prscr_group = prscr_tuple[0]
#           cohorts = prscr_tuple[1]
            
            ############  CREATE MIGRATION PLOT HERE  ###############
            main_clinical_features(binnedDf,env_factor,figPath,prscr)
            
            #filter data to high risk decile 9 & 10 and combine
            binnedDfPRScr = binnedDf[binnedDf[f'decile_bin_{prscr}'] > 8]
            binnedDfPRScr.rename(columns={prscr:'PRScr_group'},inplace=True)
                        
            
            # number of people for which all risk calculations flag as high risk
            binnedDfPRScr.set_index('IID',inplace=True)
            binnedDfPRScr['n_cohort_flagged'] = 0
            
            
            #get the people who are flagged using a cohort exclusively
#           for c in cohorts:
#               temp = binnedDfPRScr[cohorts][binnedDfPRScr[c] > 8]
#           temp = temp[temp[cohorts] > 8][cohorts].count(axis=1)
#           binnedDfPRScr.loc[temp.index,'n_cohort_flagged'] = temp    
            binnedDfPRScr['n_cohort_flagged'] = binnedDfPRScr[binnedDfPRScr[cohorts] > 8][cohorts].count(axis=1)

            # calculate number of people within each PRScr all cohort 
            binnedDfPRScr = binnedDfPRScr[binnedDfPRScr[f'decile_bin_{prscr}'] > 8]
            binnedDfCount = binnedDfPRScr.reset_index()[['quintile','PRScr_group','IID']].groupby(['PRScr_group','quintile']).count().reset_index()
            binnedCount = pd.concat([binnedDfCount,binnedCount])
            
            if i == 0:
                #get count of all quintiles together and assign as quintile 0
                binnedDfPRScrAll = binnedDfPRScr.reset_index()[['PRScr_group','IID']].groupby(['PRScr_group']).count().reset_index()
                binnedDfPRScrAll['quintile'] = 0
                binnedCount = pd.concat([binnedCount,binnedDfPRScrAll])

            binnedCount['model'] = prscr
#           binnedCount['env_factor'] = env_factor
            
            
            ############################  NUMBER OF PEOPLE FOR WHICH ALL COHORTS FLAG AS HIGH RISK AND EXCLUSIVELY WITHIN EACH COHORT  ##################
#           
#           # calculate the number of people within each cohort and quintile which would be missed with any other 
#           binnedDfCount2= binnedDfPRScr.reset_index()[['quintile','PRScr_group','n_cohort_flagged','IID']].groupby(['PRScr_group','quintile','n_cohort_flagged']).count()
#           binnedDfCount2.reset_index(inplace=True)
#           binnedDfCount2 = binnedDfCount2[binnedDfCount2['n_cohort_flagged'].isin([1,len(cohorts)])]
#           binnedCount2 = pd.concat([binnedDfCount2,binnedCount2],ignore_index=True)
#           
#           if i == 0:
#               #get count of all quintiles together and assign as quintile 0
#               binnedDfCount2 = binnedDfPRScr.reset_index()[['PRScr_group','n_cohort_flagged','IID']].groupby(['PRScr_group','n_cohort_flagged']).count().reset_index()
#               binnedDfCount2['quintile'] = 0
#               binnedDfCount2 = binnedDfCount2[binnedDfCount2['n_cohort_flagged'].isin([1,len(cohorts)])]
#               binnedCount2 = pd.concat([binnedDfCount2,binnedCount2],ignore_index=True)
#               
#               
#           binnedCount2['model'] = prscr
#           
#               
##               binnedDfCount2 = pd.concat([binnedDfCount2,binnedDfPRScrAll],ignore_index=True)
#           
##           binnedCount = pd.concat([binnedDfCount,binnedCount],ignore_index=True)
##           binnedCount2 = pd.concat([binnedDfCount2,binnedCount2],ignore_index=True)
#           binnedCount2['env_factor'] = env_factor
#           binnedCount2.rename(columns={'IID':'n_cases'},inplace=True)
    
#           binnedCountFlagged = pd.concat([binnedCount2,binnedCountFlagged],ignore_index=True)
            
            ######################  Plot the stacked count for each group in a bar chart  ##########
            dfPlot = binnedCount.pivot_table(index='quintile', columns=['model','PRScr_group'], values='IID', aggfunc='sum')
            
            
            # Number of quintiles and width of bars
            n_quintiles = len(dfPlot)
            bar_width = 0.4  # Width of each bar
            
            # X locations for the groups (quintiles)
            indices = np.arange(n_quintiles)
    
            
            # For Model A (first model)
            dfPlot[prscr].plot(kind='bar', stacked=True, ax=ax, width=bar_width, position=1, color=colors,legend=False)
            
            # For Model B (second model)
    #       dfPlot['PRScr_geno'].plot(kind='bar', stacked=True, ax=ax, width=bar_width, position=0, color=['blue','purple','red'],legend=False)
            
            # Customize the plot
            ax.set_title(env_factor, fontsize=10)
            ax.set_xlabel('Quintile', fontsize=12)
    #       ax.set_title(env_factor, fontsize=10)
            ax.set_xlim(left=-.5)
            ax.set_ylim(0,ylimBreakdown)
            
            # Set x-ticks to represent the quintiles
            ax.set_xticks(indices)
            ax.set_xticklabels(dfPlot.index)
            
            ax1.set_ylabel('# of Cases', fontsize=12)
            
            # Remove x-axis labels except on the last plot
        #   plt.setp((ax0, ax1, ax2, ax3)[:-1], xticks=[])
            
        #   Adjust layout for better spacing
        plt.tight_layout()
    
        plt.xticks(rotation=0)  # Keep the x-axis labels horizontal
    
        plt.savefig(f'{figPath}/{prscr}.barChartAcrossQuintilesForDecile9&10.CohortBreakdown.PerClinicalMeasure.png')
        
        stats.to_csv(f'{dataPath}/holdout/{prscr}.environmentFactorStats.PerClinicalMeasure.csv',index=False)
#       binnedCount2.to_csv(f'{dataPath}/holdout/{prscr}.cohortPerEnvironmentalFactorFlagged.PerClinicalMeasure.csv',index=False)
            
        
        ######### create prevalence plot  ##########

                
        #combine quintiles 1-4 and by assigning to quintile 1
        combinedNoEnv['quintile'] = combinedNoEnv[f'decile_bin_{prscr}'].apply(lambda x: 1 if x < 9 else 5)
        prevalence = calculate_prevalance(combinedNoEnv,prscr)
        prevalence['env_factor'] = 'None'
        prevalenceDf = pd.concat([prevalenceDf,prevalence])
        
        prevalenceDf.to_csv(f'{dataPath}/holdout/{prscr}.prevalenceAndCountAcrossEnvironmentalFactors.csv',index=False)
        
        # Pivot the DataFrame for grouped bar plotting
        pivot_env_count = prevalenceDf.pivot(index="quintile", columns="env_factor", values="n_cases_env")
        pivot_prs_count = prevalenceDf.pivot(index="quintile",columns="env_factor",values="n_cases_prs")
        
        # Create the bar chart
        bar_width = 0.35  # Width of the bars
        x = range(len(prevalenceDf['env_factor'].unique()))  # X-axis positions

        
        # Create figure and axis
        fig, ax = plt.subplots(figsize=(8, 6))
        
        bar_pos1 = [pos - bar_width / 2 for pos in x]
        bar_pos2 = [pos + bar_width / 2 for pos in x]
        
        for g,pos1,pos2 in zip(prevalenceDf['env_factor'].unique(),bar_pos1,bar_pos2):
            #get the prs count and env count dataframe
            non_prs_count = pivot_env_count[g]
            prs_count = pivot_prs_count[g]
            
            
            # Bar positions
            bar1 = ax.bar(pos1, prs_count.loc[1], bar_width, label='Q2', color='lightGrey')
            bar1_prs = ax.bar(pos1, non_prs_count.loc[1], bar_width, bottom=prs_count.loc[1], label='Q1', color='lightBlue')
            
            bar2 = ax.bar(pos2, prs_count.loc[5], bar_width, label='Q2', color='lightGrey')
            bar2_prs = ax.bar(pos2, non_prs_count.loc[5], bar_width, bottom=prs_count.loc[5], label='Q2', color='DarkBlue')
            

        # Add labels and title
        ax.set_xlabel("env group")
#        ax.set_ylabel("Number of Cases")
#        ax.set_title("Stacked Bar Chart: Case Distribution and PRS Contribution")
        ax.set_xticks(x)
        ax.set_xticklabels(prevalenceDf['env_factor'].unique().tolist())
#       ax.legend()
        
        # Show plot
#       plt.show()
        
        
#       
#       # Create the bar chart
#       bar_width = 0.35  # Width of the bars
#       x = range(len(pivot_env_count))  # X-axis positions
#       
#       plt.figure(figsize=(8, 5))
#       
#       # Bar positions
#       plt.bar([pos - bar_width / 2 for pos in x], pivot_df[1], bar_width, label='Q1', color='lightBlue')
#       plt.bar([pos + bar_width / 2 for pos in x], pivot_df[5], bar_width, label='Q5', color='DarkBlue')
#       
#       # Add labels, title, and legend
#       plt.xlabel('Environmental Factor Groups', fontsize=12)
#       plt.ylabel('Prevalence', fontsize=12)
#
#       plt.xticks(ticks=x, labels=pivot_df.index, fontsize=10)
#       plt.yticks(fontsize=10)
#       plt.legend(fontsize=10)
#       #plt.grid(axis='y', linestyle='--', alpha=0.6)
#       
#       # Adjust layout and show the chart
#       plt.tight_layout()
#          
#           
        fig.savefig(f'{figPath}/{prscr}.barChartAcrossQuintiles.PrevalencePlot.PRSCountLow.png')
            
        plt.show() 
    
#   plt.show()
    
    
if __name__ == '__main__':
    
#   pheno = sys.argv[1]
    pheno = 'celiacDisease'
    
    main(pheno,use_epi_main=False)
import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime
from scipy.stats import chisquare
import ast

def basic_stats(df,start_date=datetime(2013,1,1) , end_date =datetime.now(),save = False,save_name='n_unique.xlsx',select_infection=[]):
    if len(select_infection)>0:
        df = df[df['Infection'].isin(select_infection)]
    
    #Select data between two dates
    mask_date = (df['Date'] > start_date) & (df['Date'] <= end_date)
    df = df.loc[mask_date]
    
    #change name of dead and alive for result formatting
    df['survival_original'] = df['survival_original'].replace({1:'Dead',0:'Alive'})
    
    #Group by Infection for the analysis
    group_by_infection = df.groupby(['Infection'])
    
    #RESULT: N_Experiment; N_Mice
    n_unique = group_by_infection.nunique()
    n_unique_infos = n_unique.loc[:,:'Mouse_ID']
    n_unique_infos = n_unique_infos.rename(columns={'ID_Experiment':'N_Experiment','Mouse_ID':'N_Mice'})
    
    #RESULT: Date_first_experiment; Date_last_experiment
    min_date = pd.Series(group_by_infection['Date'].min().dt.strftime('%d-%m-%Y'),name='Date_first_experiment')
    max_date = pd.Series(group_by_infection['Date'].max().dt.strftime('%d-%m-%Y'),name='Date_last_experiment')
    
    #RESULT: Alive; Dead; Alive_in_%
    dead_alive = group_by_infection['survival_original'].value_counts().sort_index(ascending=False).unstack()
    dead_alive['mortality_in_%'] = round(dead_alive['Dead']/(dead_alive['Alive']+dead_alive['Dead'])*100,1)
    
    
    result = pd.concat([n_unique_infos,dead_alive,min_date,max_date],axis=1)
    if save:
        result.to_excel(save_name)
    return result

def survival_per_percentage(df, start_date=datetime(2013,1,1) , end_date =datetime.now(),save=False,plot=False,save_name = 'Count_percentage_death_per_THR_2023.xlsx'):
    mask_date = (df['Date'] > start_date) & (df['Date'] <= end_date)
    df = df.loc[mask_date]
    Infection = ['C. albicans','Staphylococcus aureus','H1N1','Pseudomonas aeruginosas'] + ['S. pneumoniae','E. coli','Listeria'] 
    df = df[df['Infection'].isin(Infection)]
    df = df.loc[:,['Infection','survival_original','survival_0.3','survival_0.25','survival_0.2','survival_0.15','survival_0.1']]
    df = df.replace({0:'survive',1:'dead'})
    count_total = df.groupby(['Infection']).apply(lambda x: x.value_counts().sort_index())
    count_original = df.groupby(['Infection'])['survival_original'].apply(lambda x: x.value_counts().sort_index())
    count_30 = df.groupby(['Infection'])['survival_0.3'].apply(lambda x: x.value_counts().sort_index())
    count_25 = df.groupby(['Infection'])['survival_0.25'].apply(lambda x: x.value_counts().sort_index())
    count_20 = df.groupby(['Infection'])['survival_0.2'].apply(lambda x: x.value_counts().sort_index())
    count_15 = df.groupby(['Infection'])['survival_0.15'].apply(lambda x: x.value_counts().sort_index())
    count_10 = df.groupby(['Infection'])['survival_0.1'].apply(lambda x: x.value_counts().sort_index())
    mix = pd.concat([count_original,count_30,count_25,count_20,count_15,count_10],axis=1)

    if plot:
        absolute_death = mix.apply(lambda x: x.survival_original-x,axis=1)
        idx = pd.IndexSlice
        absolute_death = absolute_death.loc[idx[:,'survive'],'survival_0.3':]
        absolute_death = absolute_death.droplevel(1)
        absolute_death = absolute_death.rename(columns={"survival_0.3": "THR_30", "survival_0.25": "THR_25",'survival_0.2':'THR_20','survival_0.15':'THR_15','survival_0.1':'THR_10'})
        absolute_death.T.plot.bar()
        plt.xticks( rotation='horizontal')
        plt.show()
    if save:
        mix.to_excel(save_name)
    return mix
    #count_total.to_excel('Count_percentage_death_per_THR_StepVisualisation_YEAR2023.xlsx')

def weight_loss_analysis(df,column_thr = 'survival_original'):
    slice = [(2,0.85),(0.85,0.8),(0.8,0.75),(0.75,0.7),(0.7,0)]
    name_column_df = ["<15%", "15%-20%","20%-25%","25%-30%",">30%"]
    result = {}
    
    for index, thr in enumerate(slice):
        result_Dead_Alive={}
        df_slice = df[(df['max_loss_weight_percentage']<thr[0]) & (df['max_loss_weight_percentage']>=thr[1])]
        result_Dead_Alive['Dead'] = len(df_slice[df_slice[column_thr] ==1][column_thr])
        result_Dead_Alive['Alive'] = len(df_slice[df_slice[column_thr] ==0][column_thr])
        result[name_column_df[index]] = result_Dead_Alive
    
    df_result = pd.DataFrame(result)
    return df_result.T


def plot_proportion_of_mice_by_weightLoss(df,title='Title',x_label='X label', y_label = 'Y label',save=False,path_to_save='',show=False):
    print("graph_weght_loss")
    ax = df.plot.bar(stacked=True,rot=0,color=['indianred','royalblue'],figsize=(14,8))
    Alive = df['Alive'].values.tolist()
    Dead = df['Dead'].values.tolist()
    total = sum(Alive+Dead)
    for index,p in enumerate(ax.patches):
        if index < 5:
            tota_in_p = round((Alive[index] + Dead[index])/total*100,1)
            print(Alive[index])
            print(Dead[index])
            ax.annotate(str(round(Alive[index])), ((p.get_x()+0.19) , (Alive[index])/2 + Dead[index]-5))
            ax.annotate(str(tota_in_p)+' %', ((p.get_x()+0.19) , Alive[index]+Dead[index]+2))
        else:
            ax.annotate(str(round(Dead[index-5])), ((p.get_x()+0.19) , Dead[index-5]/3))
    ax.set_title(title)
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    if save:
        plt.savefig(path_to_save)

    if show:
        plt.show()

def find_new_maximums(data):
    max_so_far = float('-inf')
    new_maximums = []
    for value in data:
        if value > max_so_far:
            max_so_far = value
        new_maximums.append(max_so_far)
    return new_maximums

def lambda_weight(x):
    x = pd.to_numeric(x,errors='coerce')
    y = x[:x.last_valid_index()]
    if y[0] == max(y):
        worse_weight_npercentage = min(y/max(y))
    else:
        new_maximums  = find_new_maximums(y)
        worse_weight_npercentage = min(y.values/new_maximums)
    return worse_weight_npercentage

            

def weight_change_hundred_percent(df):
    df_data = df.loc[:,'weight_T_infection':'weight_T32']
    series_weight = df_data.apply(lambda x: lambda_weight(x),axis=1)
    df['continuous_maximum'] = series_weight
    return df
    
def chi_square_computing(x):
    idx = pd.IndexSlice
    result = {}
    for name in x.columns.values[1:]:
        observed = x.loc[idx[:,:],name].values
        print(name)
        print(observed)
        
        expected = x.loc[idx[:,:],'survival_original'].values
        print(expected)
        chisq, p = chisquare(observed, expected)
        result[name] = {'chisq':chisq,'p_value':p}
    df = pd.DataFrame(result)
    print(df)
    return df
def chi_square_analysis(df):
    print(df)
    percentage = df/df.groupby('Infection').sum()
    result = df.groupby('Infection').apply(lambda x: chi_square_computing(x))
    print(result)
    result.to_excel('Chi_squared_analysis_YEAR2018.xlsx')

def chi_square_in_lambda_function(x):
    observed = x[x['threshold']=='survival_original'][['dead','survive']].values.tolist()[0]
    result = []
    for index, rows in x.iterrows():
        data = rows[['dead','survive']]
        chi2, p = chisquare(observed,data)
        result += [p]
    print(pd.Series(result,x.index))
    return pd.Series(result,x.index)

def find_number_of_mice_that_die_due_to_maxlossweight(df,start_date=datetime(2013,1,1) , end_date =datetime.now(),):
    print(df)
    mask_date = (df['Date'] > start_date) & (df['Date'] <= end_date)
    df = df.loc[mask_date]
    group = df.groupby(by="Infection")
    result = {}
    for index, df_group in group:
        print(index)
        if index in ("C. albicans", "H1N1","S. pneumoniae"):
            print(f'30 percent {index}')
            Total = len(df_group)
            Died = len(df_group[df_group['survival_original']==1])
            Under_thr = len(df_group[df_group['max_loss_weight_percentage']<0.7])
            score_6 = len(df_group[(df_group['score_max']>=6) & (df_group['max_loss_weight_percentage']>=0.7)])
            result_infection = {"Numbe_of_mice":Total,"Number_of_death":Died,"Death_due_to_weight":Under_thr,"Death_due_to_score_uniquely":score_6}
            result[index] = result_infection
        else:
            print(f'20 percent {index}')
            Total = len(df_group)
            Died = len(df_group[df_group['survival_original']==1])
            Under_thr = len(df_group[df_group['max_loss_weight_percentage']<0.8])
            score_6 = len(df_group[(df_group['score_max']>=6) & (df_group['max_loss_weight_percentage']>=0.8)])
            result_infection = {"Numbe_of_mice":Total,"Number_of_death":Died,"Death_due_to_weight":Under_thr,"Death_due_to_score_uniquely":score_6}
            result[index] = result_infection
    df_result = pd.DataFrame(result)
    return df_result.T

def string_to_list(s):
        if '-' in s:
            return []
        elif '[' in s:
            s = s.replace("nan", "0.0")
            return ast.literal_eval(s)
        else:
            return [float(x.strip()) for x in s.split(',')]

if __name__ == '__main__':
    df = pd.read_excel('./DF_for_simulation/df_2023_max14days_new_time_calculation_H0.xlsx',index_col=0,parse_dates=['Date'])
    

    #remove data from anesthesia
    df = df[~df.loc[:,"weight_T_infection":"weight_T14"].apply(lambda x: x.astype(str).str.contains('ane')).any(axis=1)]
    

    wl_proportion = weight_loss_analysis(df)
    
    #plot_proportion_of_mice_by_weightLoss(wl_proportion,show=False)
    survival = survival_per_percentage(df=df,start_date=datetime(2010,5,1),save_name='test.xlsx',save=False,plot=False)
    survival = survival.reset_index().rename(columns={"level_0":"infection","level_1":"survival"})
    survival = survival.melt(id_vars=["Infection","survival"],value_name="number_of_mice",var_name="threshold")


    # Pivot the table to have 'Infection' as columns and calculate the ratio
    pivot_table = survival.pivot(index=['Infection',"threshold"], columns='survival', values='number_of_mice')
    pivot_table['Ratio'] = pivot_table['dead'] / (pivot_table['dead'] + pivot_table['survive'])
    pivot_table = pivot_table.reset_index(level=1)
    pivot_table['supplementary_death'] = pivot_table.groupby('Infection').apply(lambda x: x['dead'] - x[x['threshold']=="survival_original"]['dead']).values
    
    #--------------
    #add chi2 tests
    #--------------
    """pivot_table = pivot_table.reset_index()
    chi_result = pivot_table.groupby("Infection").apply(lambda x: chi_square_in_lambda_function(x))
    chi_result = chi_result.reset_index()
    pivot_table["p_values"] = chi_result[0]
    pivot_table_of_interest = pivot_table[pivot_table['Infection'].isin(("C. albicans","Listeria","S. pneumoniae","H1N1"))]
    pivot_table_of_interest = pivot_table_of_interest.sort_values(by=["Infection",'threshold'],ascending=False)
    pivot_table_of_interest = pivot_table_of_interest[~pivot_table_of_interest['threshold'].isin(["survival_original","survival_0.3"])]
    pivot_table_of_interest = pivot_table_of_interest.drop(["dead","survive"],axis=1)
    pivot_table_of_interest = pivot_table_of_interest.pivot(index="Infection",columns='threshold',values=['Ratio',"supplementary_death","p_values"])
    pivot_table_of_interest = pivot_table_of_interest.reorder_levels([1,0],axis=1)
    columns_name = pivot_table_of_interest.columns.tolist()
    columns_name = columns_name[::-1]
    print(columns_name)
    pivot_table_of_interest = pivot_table_of_interest.reindex(columns=columns_name)
    pivot_table_of_interest.to_excel("./result_chi_squared/infection_of_interest_NO_ANE_4Infections.xlsx")
    print(pivot_table_of_interest)
    
    Infection_simulation = ['C. albicans','Staphylococcus aureus','H1N1'] + ['S. pneumoniae','E. coli','Listeria']
    infection_of_interest =  ['C. albicans','S. pneumoniae','Listeria',"H1N1"]
    df_of_interest = df[df['Infection'].isin(infection_of_interest)]
    df_of_interest_NO_LD = df_of_interest[~df_of_interest['Experiment'].str.contains("/LD")]
    bas_stat = basic_stats(df_of_interest,start_date=datetime(2010,5,1),save=False,save_name='')
    bas_stat_NO_LD = basic_stats(df_of_interest_NO_LD,start_date=datetime(2010,5,1),save=False,save_name='')
    print(df_of_interest['Strain'].value_counts()) # count strains
    print(df_of_interest[df_of_interest['Experiment'].str.contains("Training")].groupby("Infection")['Experiment'].nunique())
    print(bas_stat)"""
    #df_of_interest['Strain'].value_counts().to_excel("./result_number_mice_per_infection/STRAIN.xlsx")
    #bas_stat.to_excel("./result_number_mice_per_infection/data_of_interest_summary_NoANE_withLD.xlsx")
    #bas_stat_NO_LD.to_excel("./result_number_mice_per_infection/data_of_interest_summary_NoANE_NoLD.xlsx")

    #-----------------WEIGHT 100% IN CONTINUOUS
    """df_weight_continuous = weight_change_hundred_percent(df_simulation)
    print(df_weight_continuous)
    df_weight_continuous['<15%'] = df_weight_continuous['survival_original']
    df_weight_continuous['<20%'] = df_weight_continuous['survival_original']
    df_weight_continuous.loc[df_weight_continuous['continuous_maximum']<0.85, '<15%'] = 1
    df_weight_continuous.loc[df_weight_continuous['continuous_maximum']<0.80, '<20%'] = 1
    print(df_weight_continuous)


    mask_date = (df_weight_continuous['Date'] > datetime(2018,5,1)) & (df_weight_continuous['Date'] <= datetime.now())
    df_weight_continuous = df_weight_continuous.loc[mask_date]
    df_weight_continuous = df_weight_continuous[df_weight_continuous["Infection"].isin(['H1N1','S. pneumoniae','Listeria','E. coli'])]
    diff = df_weight_continuous.loc[~(df_weight_continuous['survival_0.2'] == df_weight_continuous['<20%'])]
    print(len(df_weight_continuous))
    print(diff[['Mouse_ID','survival_0.2','<20%','max_loss_weight_percentage','continuous_maximum']])
    print(len(diff[['Mouse_ID','survival_0.2','<20%']]))"""


    #FOR CHI2
    #number of mice for the last autorisaton
    """df_table_survival = pd.read_excel('./results_death_per_thr/table_survival_per_thr_YEAR2018.xlsx',index_col=[0,1])
    chi_square_analysis(df_table_survival)"""

    #number of mice that die uniquely due to maximum weight loss
    print(df)
    df_of_interest = df[df['Infection'].isin(("C. albicans","Listeria","S. pneumoniae","H1N1"))]
    print(df_of_interest['Scores'])
    # Custom function to convert string to list
    

# Convert strings to lists
    df_of_interest['Scores'] = df_of_interest['Scores'].apply(string_to_list)
    df_of_interest['score_max'] = df_of_interest['Scores'].apply(lambda x: max(x) if x else 0)
    df_res_with_scores = find_number_of_mice_that_die_due_to_maxlossweight(df_of_interest,start_date=datetime(2012,5,1))
    print(df_res_with_scores)
    #df_res_with_scores.to_excel("./autre/mortality_by_score_or_weight.xlsx")
#--------------OLD CODE-------------

"""    
def formatting_loss_weight_df(df):
    serie_number_mice = df['original_survival'].value_counts()
    serie_percentage_mice = df['original_survival'].value_counts(normalize=True)
    serie_percentage_mice.index = [n+'_in_%' for n in serie_percentage_mice.index.values]
    df_result_percentage = pd.concat([serie_number_mice,serie_percentage_mice])
    return df_result_percentage

def weight_loss_percentage(df,percentage):
    #Name of index for graphic
    index = ["{}% - {}%".format(n*100,percentage[idx+1]*100) if idx < len(percentage)-1 else '>{}%'.format(100*n) for idx, n in enumerate(percentage)]
    index = ["<{}%".format(percentage[0]*100)] + index
    
    dict_weight_loss = {}
    left_born = 2
    for id_percentage, n in enumerate(percentage):
        right_born = 1-n
        df_per = df[(df['max_loss_weight_percentage']<left_born) & (df['max_loss_weight_percentage']>=right_born)]
        df_result_percentage = formatting_loss_weight_df(df_per)
        dict_weight_loss[index[id_percentage]]= df_result_percentage
        left_born = right_born
    
    df_per = df[df['max_loss_weight_percentage']<=left_born]
    df_result_percentage = formatting_loss_weight_df(df_per)
    dict_weight_loss[index[-1]]= df_result_percentage
    df_w_l = pd.DataFrame(dict_weight_loss).T
    return df_w_l
    

def graph_weight_loss_percentage(df,title='Title',x_label='X label', y_label = 'Y label',save=False,path_to_save='',show=False):
    print("graph_weght_loss")
    ax = df.loc[:,['Dead','Alive']].plot.bar(stacked=True,rot=0,color=['indianred','royalblue'],figsize=(14,8))
    Alive = df['Alive'].values.tolist()
    Dead = df['Dead'].values.tolist()
    total = sum(Alive+Dead)
    for index,p in enumerate(ax.patches):
        if index < 5:
            tota_in_p = round((Alive[index] + Dead[index])/total*100,1)
            print(Alive[index])
            print(Dead[index])
            ax.annotate(str(round(Alive[index])), ((p.get_x()+0.19) , (Alive[index])/2 + Dead[index]-5))
            ax.annotate(str(tota_in_p)+' %', ((p.get_x()+0.19) , Alive[index]+Dead[index]+2))
        else:
            ax.annotate(str(round(Dead[index-5])), ((p.get_x()+0.19) , Dead[index-5]/3))
    ax.set_title(title)
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    if save:
        plt.savefig(path_to_save)

    if show:
        plt.show()


    w_l_p_graphique = False

    if w_l_p_graphique:
        
        df_weight_loss_chronique = df[df['Infection'].isin(['C. albicans','Staphylococcus aureus','H1N1'])]
        df_weight_loss_acute = df[df['Infection'].isin(['S. pneumoniae','E. coli','Listeria'])]
        df_weight_loss_Training = df[df['Experiment'].str.contains('Training')]
        
        df_training_wl = weight_loss_percentage(df_weight_loss_Training,[0.15,0.20,0.25,0.30])
        df_chronique_wl = weight_loss_percentage(df_weight_loss_chronique,[0.15,0.20,0.25,0.30])
        df_acute_wl = weight_loss_percentage(df_weight_loss_acute,[0.15,0.20,0.25,0.30])
        graph_weight_loss_percentage(df_wl_p,title="Survivability of mice at end of experiments in function of weight loss\n ALL",
        x_label='weight loss',y_label='Survivability',save=False,path_to_save='./graphiques/basic_stats/bar_aLL_absolute_values',show=False) #ALL EXPERIMENT
        
        #graph_weight_loss_percentage(df_training_wl,title="Survivability of mice at end of experiments in function of weight loss\n TRAINING",
        #x_label='weight loss',y_label='Survivability',save=True,path_to_save='./graphiques/basic_stats/bar_training_absolute_values',show=True)#TRAINING
        
        #graph_weight_loss_percentage(df_acute_wl,title="Survivability of mice at end of experiments in function of weight loss\n ACUTE",
        #x_label='weight loss',y_label='Survivability',save=True,path_to_save='./graphiques/basic_stats/bar_acute_absolute_values',show=True)#ACUTE
        
        #graph_weight_loss_percentage(df_chronique_wl,title="Survivability of mice at end of experiments in function of weight loss\n CHRONIQUE",
        #x_label='weight loss',y_label='Survivability',save=True,path_to_save='./graphiques/basic_stats/bar_chronique_absolute_values',show=True)#CHRONIQUE
    """
import pandas as pd
from lifelines.statistics import multivariate_logrank_test
import math
from datetime import datetime

def format_df_for_p_values(df, end_of_df_column_name = 'weight_T32'):
    """
    Add a column: max_loss_weight_percentage. Compute the maximum loss weight percentage of each mouse during the experiment.
    Invert Group name: Inside the column Group, if the name is a string, it will invert the name. Example: 0A --> A0
    Add column: exp. Take the first value if it's a string, otherwise the same value. Example: A0 --> A; 1 --> 1
    *Arguments*
    -df: dataframe object
    -end_of_df_column_name:
    *Return*
    nothing, but add columns to df: max_loss_weight_percentage and exp
    """

    df['max_loss_weight_percentage'] = df.loc[:,'weight_T1':end_of_df_column_name].apply(lambda x: pd.to_numeric(x, errors='coerce')).min(axis=1)/df['weight_T_infection']
    df['Group'] = df['Group'].apply(lambda x: x[1]+x[0] if isinstance(x,str) else x)
    df['exp'] = df['Group'].apply(lambda x:x[0] if isinstance(x,str) else x)



def true_events(x,max_weight):
    """
    Change mice survival datas to 0 when the event of death was not recorded.
    For example, if the mouse was sacrified at a specific weight loss threshold
    *Arguments*
    -x: serie object (pandas object): survival columns
    -max_weight: serie object (pandas object), max_weight (axis = 0) of all the mice in the df. 
    To find the current max weight -> max_weight[id_of_interest] = max_weight[x.name]

    *Return*
    x: serie object
    """
    if math.isnan(max_weight[x.name]): #When the mouse died the day after infection, there is no max_loss_weight --> return x
        return x
    else:
        idx_max = round((1-max_weight[x.name])*100) #1-percentage max obtained --> directly correspond to the id to keep as each column represent 1%
        if idx_max >= len(x): 
            x[:-1]=0 #everythinf has survive --> 0
        else:
            x[(len(x)-idx_max):-1]=0   #after a specific threshold, we can't say if the mice would be dead later --> censored at 0
        return x

def censored_survival(df,columns_survival,end_of_df_column_name = 'weight_T32'):
    """
    Create a new dataframe (censored) who count differenty how a death event occure
    In the original df, a death event occure when the mouse reach a specific threshold or when the scientist decid to sacrifice the mouse.
    In the censored df, a death event occure only if the mouse has be found dead. If criteria of sacrifice are met (like a weight loss threshold),
    it is not consiered as a true death.

    To resume: If the mouse is found dead --> survival = 1 (dead); if the mouse is sacrified --> survival = 0 (censored)
    *Arguments*
    -df: dataframe object
    -columns_survival: 
    -end_of_df_column_name:
    *Return*
    new_df, datafrme object (pandas object)
    """

    new_df = df.copy()
    new_df.loc[(new_df.loc[:,'weight_T1':end_of_df_column_name].values == 'sacrified').any(1),'original_survival'] = 0 #change the column original survival to 0 if the world "sacrified is found"
    new_df.loc[new_df.loc[:,'original_survival'] == 0,columns_survival] = 0 #
    new_df[columns_survival] = new_df[columns_survival].apply(lambda x: true_events(x,new_df.max_loss_weight_percentage),axis = 1)
    return new_df

def choose_group_and_df_Type(df,censored,training, Infection=[],Infection_excluded=[],start_date = datetime(2000,1,1),end_date = datetime.now(),**kwarg):
    """
    Function to chose the datas to analyze in the dataframe.
    You can chose to censore the data or not. If a mouse is sacrified, it will be censure: envent = 0 and not 1
    You can chose only experiment of training
    You can chose specific infection, for example S. pneuom, Listeria, etc.
    You can also chose to exclude some infections

    To resume: If the mouse is found dead --> survival = 1 (dead); if the mouse is sacrified --> survival = 0 (censored)
    *Arguments*
    -df: dataframe object
    -censored: Bool, chose if we want to censore or not the data
    -training: int 1: only keep training data, -1 remove training data, other number: do nothing
    -Infection: list object, keep only the desired infection
    -Infection_excluded: list object, remove the desired infection
    *Return*
    df, datafrme object (pandas object)
    """
    df = df.copy()
    mask_date = (df['Date'] > start_date) & (df['Date'] <= end_date)
    df = df.loc[mask_date]
    if censored:
        columns_survival = df.columns[df.columns.str.contains('survival')].values.tolist()
        df = censored_survival(df,columns_survival,kwarg['end_of_df_column_name'])
    
    if training == 1:
        df = df[df['Experiment'].str.contains("Training")]
    elif training == -1:
        df = df[~df['Experiment'].str.contains("Training")]
    
    if len(Infection) >0:
        df = df[df['Infection'].isin(Infection)]
    if len(Infection_excluded)>0:
        df = df[~df['Infection'].isin(Infection_excluded)]
        
    return df



def add_p_values(p_value,alpha = 0.05):
    """
    function that return 1 if the p value is lower than a chosen alpha

    *Arguments*
    -p_value: p value to evaluate
    -alpha: alpha chosen
    *Return*
    - 1 if p < alpha; 0 otherwise
    """
    try:
        p_value > 0
    except ValueError:
        print("negative value")
    if p_value < alpha:
        return 1
    else:
        return 0

def log_rank_analysis(df,survival_percentage,time_percentage):
    """
    Compute a log rank test with a dataframe of two groups

    *Arguments*
    -df: Dataframe object with a column containing the groups name Group
    -survival_percentage: Name of column in df that contain the datas of survival (1 dead, 0 survive)
    -time_percentage: Name of column in df that contain the time datas (time until an event)

    *Return*
    - result_multi.p_value: p_value of logr rank test (float)
    """
    E = df[survival_percentage].values
    G = df['Group'].values
    T = df[time_percentage].values
    result_multi = multivariate_logrank_test(T, G, E)
    return result_multi.p_value



def log_rank_group_multi_group_in_binome(df,time_col='time_original',survival_col='survival_original'):
    """
    Perform a log rank analysis for each sub group (A,B,C,...) in the experiment.
    Rreturn a serie with p_value attached to each experiemnt.

    *Arguments*
    -df: Dataframe object with a column containing the different sub experiment name exp
    -time_col: Name of column in df that contain the datas of survival (1 dead, 0 survive)
    -survival_col: Name of column in df that contain the time datas (time until an event)

    *Return*
    - result_sub_grp: serie of p value per experiment
    """

    sub_group = df.groupby('exp')
    result_sub_grp = sub_group.apply(lambda x : log_rank_analysis(x,time_percentage=time_col,survival_percentage=survival_col))
    return result_sub_grp

def log_rank_multi_group_with_control(df,time_col='time_original',survival_col='survival_original'):
    """
    Perform a log rank analysis for each variable with a control. Group = 0 --> always a control groupe
    Example: Unique Group = [0,1,2,3]
        Log rank: 0_1 --> pvalue 1 
        Log rank: 0_2 --> pvalue 2 
        Log rank: 0_3 --> pvalue 3 

    *Arguments*
    -df: Dataframe object with a column containing the different sub experiment name exp
    -time_col: Name of column in df that contain the datas of survival (1 dead, 0 survive)
    -survival_col: Name of column in df that contain the time datas (time until an event)

    *Return*
    - panda Serie: serie of p value per variable
    """
    dict_result_multi = {}
    for n in df['Group'].unique():
        if n ==0:
            continue
        else:
            experiment = df[df['Group'].isin([0,n])]
            p_val = log_rank_analysis(experiment,time_percentage=time_col,survival_percentage=survival_col)
            dict_result_multi['0_{}'.format(n)] = p_val
    return pd.Series(dict_result_multi)
            
            


def count_p_values_stats(df,time_col='time_original',survival_col='survival_original',alpha=0.05):
    """
    Compute the number of p_values statistically significative (alpha = 0.05) for one given ID experiment

    *Arguments*
    -df: Dataframe object of one Experiment ID
    -time_col: Name of column in df that contain the datas of survival (1 dead, 0 survive)
    -survival_col: Name of column in df that contain the time datas (time until an event)
    -aplha: significance level. probablity of rejecting the null hypothesis. Default 5%
    *Return*
    -count_p_values: number of p_values that are below alpha (reject the null hypothesis)
    -number_of_experiment: number of p_values computed = number of log rank tested
    """

    groups = df['Group'].unique().tolist() #number of different variable/group in the experiemnt
    is_numeric = all([isinstance(item, int) for item in groups]) #test if there is non numberic group: example A, B,... --> corresponding to subexperiemnt
    count_p_values = 0 
    number_of_experiment = 0
    
    if (0 not in groups) and (is_numeric) and (len(groups)>2): #test if there is an experiment without a control neither sub group and greater than 2
        print('bad experiment', df.name)                       #Experiment with group: [1,2,3] are removed as we can't decide which test to effectuate
        return count_p_values, number_of_experiment
    else:
        serie_result = pd.Series([],dtype='float64')

        if len(groups)==2: #if only two groups, do a simple log rank test
            p_val = multivariate_logrank_test(df[time_col],df['Group'].astype(str),df[survival_col]).p_value
            count_p_values = add_p_values(p_value=p_val,alpha=alpha)
            number_of_experiment = 1
            return count_p_values, number_of_experiment

        if not is_numeric: #if there is sub group (A,B,C...) perform log rank analysis for each sub group
            serie = log_rank_group_multi_group_in_binome(df,time_col=time_col,survival_col=survival_col)
            serie_result = pd.concat([serie_result,serie])

        if 0 in groups: #if there is more than one variable perform all combinaison of log rank test with the control group
            serie = log_rank_multi_group_with_control(df,time_col=time_col,survival_col=survival_col)
            serie_result = pd.concat([serie_result,serie])

        serie_result = serie_result.dropna()
        bonferonie_corr = len(serie_result)
        count_p_values = serie_result[serie_result < alpha/bonferonie_corr].count() #alpha/len(serie_result) -->  bonferroni correction
        number_of_experiment = len(serie_result)
            
        
        

    return count_p_values, number_of_experiment

def analyse_p_value(df,time_col='time_original',survival_col='survival_original',resume=True):
    """
    Compute the number of p_values statistically significative (alpha = 0.05) for all ID experiment

    *Arguments*
    -df: whole dataframe
    -time_col: Name of column in df that contain the datas of survival (1 dead, 0 survive)
    -survival_col: Name of column in df that contain the time datas (time until an event)
    *Return*
    result.sum() : serie object, sum of all statisticaly significative p_values
    """
    group = df.groupby('ID_Experiment')
    result = group.apply(lambda x: count_p_values_stats(x,time_col=time_col,survival_col=survival_col))
    result = result.apply(pd.Series)
    result.columns=['p_stats','number_p']
    if resume:
        return result.sum()
    else:
        return result


def count_p_values_per_thr(df,save = False,path_to_save = 'test.xlsx',resume=True):
        """
        Compute the number of p_values statistically significative (alpha = 0.05) for all Threshold of sacrifice.

        *Arguments*
        -df: whole dataframe
        -save: bool, if the result serie need to be save in a xlsx file
        -path_to_save: name of the file
        *Return*
        -None
        """
        
        #create list of column name
        survival = df_to_analyse.columns[df.columns.str.contains('survival')].values.tolist()
        pop = survival.pop()
        survival = [pop] + survival
        time = df_to_analyse.columns[df.columns.str.contains('time')].values.tolist()
        pop = time.pop()
        time = [pop] + time
        
        dict_percentage = {}
        zip_column = zip(time,survival)
        
        if resume:
            for index, columns in enumerate(zip_column):
                print(columns)
                sum_p = analyse_p_value(df,time_col=columns[0],survival_col=columns[1])
                dict_percentage[columns[0]] = sum_p['p_stats']
                """if index ==0:
                    break"""
            dict_percentage['number_p'] = sum_p['number_p']
            serie_per_p_value = pd.Series(dict_percentage)
            print(serie_per_p_value)
            if save:
                serie_per_p_value.to_excel(path_to_save)
        else:
            list_df = []
            for index, columns in enumerate(zip_column):
                print(columns)
                result_p_value_per_experiment = analyse_p_value(df,time_col=columns[0],survival_col=columns[1],resume=resume)
                name = columns[0].split('_')[1]
                result_p_value_per_experiment = result_p_value_per_experiment.rename(columns={'p_stats':'p_stats_{}'.format(name),'number_p':'number_p_{}'.format(name)})
                list_df += [result_p_value_per_experiment]
            result = pd.concat(list_df,axis=1)
            if save:
                result.to_excel(path_to_save)
def add_p_values_with_H0(df,p_value,alpha=0.05):
    if any(df['H0']<0) and p_value>alpha:
        return 1
    elif any(df['H0']>0) and p_value<alpha:
        return 1
    else:
        return 0

def compute_p_values(df,time_col='time_original',survival_col='survival_original',alpha=0.05):
    dict_result_multi = {"count_p_values":0,"number_of_experiment":0}

    count_p_values = 0
    number_of_experiment = 0
    groups = df['exp'].unique().tolist()
    if len(groups) > 0 and (0 not in groups):
        return pd.Series(dict_result_multi)
    if len (groups) == 2:
        if any(df['H0']<0):
            return pd.Series(dict_result_multi)
        p_value = multivariate_logrank_test(df[time_col],df['exp'].astype(str),df[survival_col]).p_value
        dict_result_multi["count_p_values"] =  add_p_values_with_H0(df,p_value,alpha=alpha)
        dict_result_multi["number_of_experiment"] = 1
        return pd.Series(dict_result_multi)
    else:
        for n in df['exp'].unique():
            if n ==0:
                continue
            else:
                experiment = df[df['exp'].isin([0,n])]
                if any(experiment['H0']<0):
                    continue
                p_value = log_rank_analysis(experiment,time_percentage=time_col,survival_percentage=survival_col)
                count_p_values += add_p_values_with_H0(experiment,p_value,alpha=alpha)
                number_of_experiment += 1
        dict_result_multi['count_p_values']=count_p_values
        dict_result_multi['number_of_experiment'] = number_of_experiment
        return pd.Series(dict_result_multi)


def stats_with_H0(df,time_col='time_original',survival_col='survival_original',alpha=0.05):
    """
    Compute the number of p_values statistically significative (alpha = 0.05) for one given ID experiment

    *Arguments*
    -df: Dataframe object of one Experiment ID
    -time_col: Name of column in df that contain the datas of survival (1 dead, 0 survive)
    -survival_col: Name of column in df that contain the time datas (time until an event)
    -aplha: significance level. probablity of rejecting the null hypothesis. Default 5%
    *Return*
    -count_p_values: number of p_values that are below alpha (reject the null hypothesis)
    -number_of_experiment: number of p_values computed = number of log rank tested
    """

    group = df.groupby(['ID_Experiment','sub_exp'])
    result = group.apply(lambda x: compute_p_values(x,time_col=time_col,survival_col=survival_col))
    
    return result.sum()


if __name__ == '__main__':

    #df = pd.read_excel('./results_weight_percentage_clean_goodDate_14days.xlsx')
    #format_df_for_p_values(df,end_of_df_column_name='weight_T14')

    #df = pd.read_excel('./Dataframe_for_p_value_analysis.xlsx') #FOR ORIGINAL DATA
    #df = pd.read_excel('./data_14days_for_simulation.xlsx') #14 DAYS DATA
    df = pd.read_excel('./DF_for_simulation/df_2023_max14days_new_time_calculation_H0.xlsx',index_col=0)
    columns_time = df.columns[df.columns.str.contains('time')].values.tolist()
    columns_survival = df.columns[df.columns.str.contains('survival')].values.tolist()
    df = df[~df['Experiment'].str.contains("/LD")]
    df = df[df['ID_Experiment']!="ID_exp_89"]

    chronique = ['C. albicans','Staphylococcus aureus','H1N1'] #chronique
    acute = ['S. pneumoniae','E. coli','Listeria'] #acute
    model_of_interest  = ['C. albicans','S. pneumoniae','Listeria',"H1N1"]
    df_to_analyse = choose_group_and_df_Type(df,censored = False, training = 0,Infection=['H1N1'])
    #Analyse p values
        #count_p_values_per_thr(df_to_analyse,save=False,path_to_save='./result_p_values/test/test.xlsx')
    
    CALCULATE_P_VALUES_WITH_H0 = True
    if CALCULATE_P_VALUES_WITH_H0:
        columns_time = [columns_time[-1]] + columns_time[:-1]
        columns_survival = [columns_survival[-1]] + columns_survival[:-1]
        result = {}
        for col_time, col_survival in zip(columns_time,columns_survival):
            name = col_time.split("_")[1]
            print(name)
            stats = stats_with_H0(df_to_analyse,col_time,col_survival)
            result[name] = stats
        result_df = pd.DataFrame(result).T
        print(result_df)
        result_df.to_excel("./result_p_values/WITH_H0/H1N1.xlsx")
    
    CALCULATE_P_VALUES_PER_MODEL = False
    if CALCULATE_P_VALUES_PER_MODEL:
        for model in ['E. coli', 'S. pneumoniae', 'Listeria','Staphylococcus aureus','H1N1','C. albicans']:
            df_unique_model = choose_group_and_df_Type(df,censored=False,training=0,Infection=[model],start_date=datetime(2010,4,11))
            count_p_values_per_thr(df_unique_model,save=True,path_to_save='./result_p_values/models/{}_2023_new_time_update.xlsx'.format(model),resume=True)

    #OLD CODE
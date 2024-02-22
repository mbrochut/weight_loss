import pandas as pd
import os
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

import researchpy as rp
import statsmodels.api as sm
import scipy.stats as stats
import statsmodels.formula.api as smf

def keep_weight_post_infection(x,weight_end = "weight_T14"):
    dates = x['Dates']
    t_infection = x['Time_infection']
    datas = x['weight_T_infection':weight_end]

    #find date la plus proche de J infection
    new_time_infection = dates[dates <= t_infection][-1]
    location_of_TI = dates.get_loc(new_time_infection)
    if np.isnan(datas[location_of_TI]):
        return datas 
    shifted_series = pd.Series([np.nan] * len(datas), index=datas.index)
    if location_of_TI == 0:
        return datas
    else:
        shifted_series[:-location_of_TI] = datas.values.tolist()[location_of_TI:]
    # Shift the values of the input series by the specified index
    return shifted_series


def Mixed_Effects_Models(df,chosen_infection = 'S. pneumoniae',time_to_exclude = 8):
    df_infection = df[df['Infection'] == chosen_infection]
    df_infection = df_infection[~df_infection['Time'].isin([f"T{n}" for n in range(time_to_exclude,15,1)])]#remove unused data
    
    time_point = df_infection['Time'].unique()
    weight_point_to_integer = dict(zip(time_point,[n for n in range(len(time_point))]))
    df_infection['Time'] = df_infection['Time'].replace(weight_point_to_integer)
    df_infection = df_infection[df_infection['group_death'].isin(['Control_alive','Challenge_alive'])]
    model = smf.mixedlm("weight ~ Time",df_infection,groups=df_infection['group_death'],missing="drop").fit()
    print(model.summary())

def plot_weight_by_infection(df,infection="Listeria",group="group_death",plot = False,save = False,path_to_save = "weight_over_time.png",errorbar = ("sd"),max_time="T14"):
    df = df[df['Infection'] == infection]
    if group == "group_death":
        hue_order = ['Control_alive','Control_dead','Challenge_alive','Challenge_dead']
    else:
        hue_order = None
    plt.figure(figsize=(12,8))
    sns.lineplot(data=df,x = "Time",y = "weight",hue=group,errorbar=errorbar,err_style="band", hue_order=hue_order,legend=False, estimator=np.mean) #, estimator=np.median
    plt.xticks(rotation=45, ha='right', fontsize=15)
    plt.yticks(fontsize=15)
    plt.xlabel('Time [days]', fontsize=18)
    plt.ylabel('Weight loss [%]', fontsize=18)
    plt.title(infection, fontsize=22)
    plt.subplots_adjust(top=0.9,bottom=0.2)
    plt.xlim("Tinfection", max_time)
    plt.ylim(65, 105)
    if plot:
        plt.show()
    if save:
        plt.savefig(path_to_save+infection+".png")

def count_add_death(df):
    print(df)
    for thr in [75,80,85,90]:
        df[f'add_thr_{thr}'] = df.apply(lambda x: 1 if x['survival_original']==0 and x['min_weight']<thr else 0,axis=1)
    counts = df.groupby('Infection')[["add_thr_75", "add_thr_80", "add_thr_85", "add_thr_90"]].sum()
    print(counts)
if __name__ == "__main__":

    #READ DF
    df = pd.read_excel("./DF_for_simulation/df_2023_max14days_new_time_calculation_H0.xlsx")
    #DATA FOR ANALYZE
    #PREPARE DF FOR ANALYSIS, change dates column to datetimindex and transform weight datas to numeric only
    df.loc[:,"weight_T_infection":"weight_T14"] = df.loc[:,"weight_T_infection":"weight_T14"].apply(pd.to_numeric,errors='coerce')
    serie_dates = df['Time_point'].apply(lambda x: pd.to_datetime(x.split(','),dayfirst=True))
    df['Dates'] = serie_dates
    data = df.apply(lambda x: keep_weight_post_infection(x),axis=1)

    
    #normalization by weight at T infection
    normalize = data.div(data['weight_T_infection'],axis=0)*100
    df_normalize = df.copy()
    df_normalize.loc[:,"weight_T_infection":"weight_T14"] = normalize
    
    df_normalize['min_weight'] = df_normalize.loc[:,"weight_T_infection":"weight_T14"].min(axis=1)
    df_normalize['t_origin'] = df_normalize['time_original']
    count_add_death(df_normalize)

    #FIND COLUMN TO ANALYZE
    columns = df_normalize.loc[:,"weight_T_infection":"weight_T14"].columns.tolist()
    columns_index = df_normalize.loc[:, ~df_normalize.columns.isin(columns)]
    column_time = [n for n in df_normalize.columns.tolist() if "time_" in n]
    column_time = [column_time[-1]] + column_time[:-1]
    columns_index_time = ['ID_Experiment','Mouse_ID','Date','Infection','Group','exp','survival_original','t_origin']
    #TRANSFORM TO TIDY DATA
    df_longer = df_normalize.melt(id_vars=columns_index,value_vars=columns,var_name="weight_point",value_name="weight")
    df_longer_time = df_normalize.melt(id_vars=columns_index_time,value_vars=column_time,var_name="Threshold",value_name="Time")
    df_longer_time.to_csv("./DF_for_simulation/time_clean.csv")

    #INFECTION TO ANALYZE
    chosen_infection = ['E. coli',"C. albicans", "Listeria", "S. pneumoniae", "H1N1", "Pseudomonas aeruginosas", "Staphylococcus aureus"]
    df_longer = df_longer[df_longer['Infection'].isin(chosen_infection) ]

    #find group control
    df_longer['control'] = df_longer['Group'].apply(lambda x: x if isinstance(x,int) else int(x[0]))
    df_longer['control'] = df_longer['control'].apply(lambda x: x if x == 0 else 1)


    #RENAME FOR PLOTING AND ADD GROUP DEATH
    rename_group = {"0_0":"Control_alive","0_1":"Control_dead","1_0":"Trained_alive","1_1":"Trained_dead"}
    rename_group2 = {"0_0":"Control_alive","0_1":"Control_dead","1_0":"Challenge_alive","1_1":"Challenge_dead"}
    df_longer['group_death'] = df_longer['control'].astype(str) + "_" +df_longer['survival_original'].astype(str)
    #choose which rename to use
    df_longer = df_longer.replace({"group_death":rename_group2})
    
    #change weight name for better visualisation
    df_longer['weight_point'] = df_longer['weight_point'].apply(lambda x: "".join(x.split("_")[1:]))
    df_longer = df_longer.rename(columns={"weight_point":"Time"})
    columns_to_export = ['ID_Experiment','Mouse_ID',"Infection",'Time_infection','time_original','survival_original','max_loss_weight_percentage','exp','min_weight','Time','weight','control','group_death']
    #df_longer.loc[:,columns_to_export].to_excel("./DF_for_simulation/weigth_clean.xlsx")
    
    print(df_longer)
    df_longer_listeria = df_longer[df_longer['Infection']=='Listeria']
    count = df_longer_listeria[df_longer_listeria['weight']<80].groupby('Mouse_ID')['weight'].count()
    count = count[count >=2]
    count.to_excel("./autre/Listeria_under_20.xlsx")
    #type of data to analyze: infection and training
    df_longer_infection = df_longer[df_longer['Infection'] == 'S. pneumoniae']
    df_longer_training = df_longer[df_longer['Experiment'].str.contains("Training")]
    

    #Mixed_Effects_Models(df_longer,"C. albicans",10)
    #PLOTING

    #plot_weight_by_infection(df_longer,plot=True)
    
    
    #plot by infection
    """sns.relplot(data=df_longer,kind="line",x = "Time",y = "weight",hue="Infection",errorbar="ci",err_style="band") #, estimator=np.median
    plt.xticks(rotation=45, ha='right')
    plt.show()"""


    #plot all Infections

    """infections = df_longer['Infection'].unique()
    print(infections)
    infection_of_interest = zip(["T12","T7","T9"],['C. albicans','Listeria','S. pneumoniae'])
    for max_time, infection in infection_of_interest:
        plot_weight_by_infection(df=df_longer,infection=infection,save=True,path_to_save="./graphiques/NEW_ANALYSIS/weight_over_time/",max_time=max_time)"""
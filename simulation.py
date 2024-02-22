import pandas as pd
import numpy as np
import statistics as stat
from lifelines.statistics import multivariate_logrank_test
import time
import multiprocessing as mp

def take_data_for_log_rank(df,n_mice,col_time='time_original',col_survival="survival_original"):
        df_control = df[df['control']==0].sample(n=n_mice,replace=True)
        df_experiment = df[df['control']==1].sample(n=n_mice,replace=True)
        E_c = df_control[col_survival].values
        G_c = df_control['control'].values
        T_c = df_control[col_time].values
        E_e = df_experiment[col_survival].values
        G_e = df_experiment['control'].values
        T_e = df_experiment[col_time].values
        event = np.concatenate((E_c,E_e))
        group = np.concatenate((G_c,G_e))
        time_of_event = np.concatenate((T_c ,T_e))
        return event,group,time_of_event

def simulation_log_rank(df,length=100,n_mice=10,col_time='time_original',col_survival="survival_original"):
    stat_p_values = 0
    for simu in range(length):
        event, group, time_of_event = take_data_for_log_rank(df,n_mice,col_time,col_survival)
        log_rank = multivariate_logrank_test(time_of_event,group,event)
        p_value = log_rank.p_value
        if p_value < 0.05:
            stat_p_values +=1
        else:
            continue
    return stat_p_values/length


def multiple_simulation_with_results(df,N_simulation,N_repetition,N_mice,cpu=4,col_time='time_original',col_survival="survival_original"):
        pool = mp.Pool(cpu)
        power = pool.starmap(simulation_log_rank,[(df,N_simulation,N_mice,col_time,col_survival) for n in range(N_repetition)])
        pool.close()
        result = {"Power":power,
        "Repetition_No":np.arange(1,N_repetition+1,1),
        "N_simulation":N_repetition*[N_simulation],
        "N_mice":N_repetition*[N_mice],
        "Threshold":N_repetition*[col_time.split("_")[1]]}
        df_result = pd.DataFrame(result)
        return df_result

def automatic_simulation_per_infection(df,N_simulation,N_repetition,mice_range,cpu=4,col_time='time_original',col_survival="survival_original"):
    infections = df['Infection'].unique().tolist()
    for infection in infections:
        df_infection = df[df['Infection']==infection]
        counts = df_infection.groupby(["Infection"])['control'].value_counts()
        print(counts)
        t1 = time.time()
        df_total = []
        for N_mice in mice_range:
            print(N_mice)
            tsim1 = time.time()
            df_result = multiple_simulation_with_results(df_infection,
            N_simulation=N_simulation,
            N_repetition=N_repetition,
            N_mice=N_mice,
            cpu = cpu,
            col_time=col_time,
            col_survival=col_survival)

            df_result['Infection'] = infection
            df_total +=[df_result]
            tsim2 = time.time()
            print("Simulation took: ",tsim2-tsim1,"s")
        df_simulation = pd.concat(df_total,axis=0).reset_index(inplace=False,drop=True)
        print(df_simulation)
        t2 = time.time()
        print("Complete simulation took: ",t2-t1,"s")
        df_simulation.to_excel(f"./simulations/new_simulations/{infection}_Nsimulation{N_simulation}_Nrepetition{N_repetition}_Nmice400.xlsx")



def automatic_simulation_per_infection_fix_mice(df,N_simulation,N_repetition,column_threshold,dict_mice_range,cpu=4):
    infections = df['Infection'].unique().tolist()
    for infection in infections:
        df_infection = df[df['Infection']==infection]
        counts = df_infection.groupby(["Infection"])['control'].value_counts()
        print(counts)
        t1 = time.time()
        df_total = []
        N_mice = dict_mice_range[infection]
        print("N mice 80% = ",N_mice)
        for col_survival,col_time in column_threshold:
            print(col_survival,col_time)
            tsim1 = time.time()

            df_result = multiple_simulation_with_results(df_infection,
            N_simulation=N_simulation,
            N_repetition=N_repetition,
            N_mice=N_mice,
            cpu = cpu,
            col_time=col_time,
            col_survival=col_survival)

            df_result['Infection'] = infection
            df_total +=[df_result]
            tsim2 = time.time()
            print("Simulation took: ",tsim2-tsim1,"s")
        df_simulation = pd.concat(df_total,axis=0).reset_index(inplace=False,drop=True)
        print(df_simulation)
        t2 = time.time()
        print("Complete simulation took: ",t2-t1,"s")
        df_simulation.to_excel(f"./simulations/new_simulations/{infection}_Nsimulation{N_simulation}_Nrepetition{N_repetition}_fixNmice{N_mice}.xlsx")


if __name__ == '__main__':
    print('SIMULATION ANALYSIS')

    df = pd.read_excel('./data/df_for_analysis.xlsx',index_col=0)
    df['Time_infection'] = df['Time_infection'].dt.strftime("%Y-%m-%d")
    df['Time_point'] = df['Time_point'].apply(lambda x: pd.to_datetime(x.split(','),dayfirst=True,exact=False))
   
    #remove lethal dose
    df = df[~df['Experiment'].str.contains("/LD")]
    #Dataset for Article
    infection_of_interest = ["C. albicans","Listeria","S. pneumoniae","H1N1","E. coli","Pseudomonas aeruginosas","Staphylococcus aureus"]

    df = df[df['Infection'].isin(infection_of_interest)]
    
    #add info of group control or not
    df['control'] = df['exp'].apply(lambda x: 0 if x==0 else 1)

    #list of column name for survival and time for each threshold
    columns_time = df.columns[df.columns.str.contains('time')].values.tolist()
    columns_survival = df.columns[df.columns.str.contains('survival')].values.tolist()

    columns_survival = [columns_survival[-1]] + columns_survival[:-1] #PUT ORIGINAL DATA AT FIRST
    columns_time = [columns_time[-1]] + columns_time[:-1] #PUT ORIGINAL DATA AT FIRST
    column_threshold = list(zip(columns_survival,columns_time)) #USEFUL TO LOOP OVER CONDITION: THRESHOLD
  
    #FOR CHECKING NUMBER OF MICE IN EACH SIMULATION
    

    mice_100 = np.logspace(1, 2, num=20).astype(int).tolist()
    mice_200 = np.arange(110,210,10).astype(int).tolist()
    mice_600 = np.arange(210,410,10).astype(int).tolist()
    mice_range = mice_100 + mice_200 + mice_600

    mice_range = mice_range[0:10]
    chose_infection = "C. albicans"
    N_simulation = 10
    N_repetition = 2
    #df = df[df['Infection']==chose_infection]
    counts = df.groupby(["Infection"])['control'].value_counts()
    print(counts)
    counts.to_excel("./result_number_mice_per_infection/count_control_group.xlsx")

    #automatic_simulation_per_infection(df,N_simulation=5000,N_repetition=20,mice_range = mice_range,cpu=12)


    t1 = time.time()

    #"C. albicans","Listeria","S. pneumoniae","H1N1","E. coli","Pseudomonas aeruginosas","Staphylococcus aureus"
    infection_of_interest = ["Listeria","S. pneumoniae","E. coli","Staphylococcus aureus"]
    N_mice_80 = [35,56,140,126]
    dict_mice_range = dict(zip(infection_of_interest,N_mice_80))
    print(dict_mice_range)
    df = df[df['Infection'].isin(infection_of_interest)]
    
    """automatic_simulation_per_infection_fix_mice(df,
    N_simulation=5000,
    N_repetition=20,
    column_threshold=column_threshold,
    dict_mice_range = dict_mice_range,
    cpu=11)"""
    
    """df_total = []
    for N_mice in mice_range:
        print(N_mice)
        tsim1 = time.time()
        df_result = multiple_simulation_with_results(df,N_simulation=N_simulation,N_repetition=N_repetition,N_mice=N_mice,cpu = 6,col_survival='survival_0.15',col_time='time_0.15')
        df_total +=[df_result]
        tsim2 = time.time()
        print("Simulation took: ",tsim2-tsim1,"s")
    df_simulation = pd.concat(df_total,axis=0).reset_index(inplace=False,drop=True)
    print(df_simulation)
    t2 = time.time()
    print("Complete simulation took: ",t2-t1,"s")
    df_simulation.to_excel(f"./simulations/new_simulations/{chose_infection}_Nsimulation{N_simulation}_Nrepetition{N_repetition}_TESTTING_THR.xlsx")"""
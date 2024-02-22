import pandas as pd
import numpy as np
from lifelines.statistics import multivariate_logrank_test
import time
import multiprocessing as mp

def take_data_for_log_rank(df,n_mice,col_time='time_original',col_survival="survival_original"):
    """
    Extracts data for log-rank test.

    *Arguments*
    - df: DataFrame containing the data.
    - n_mice: Number of mice to sample.
    - col_time: Column name containing the time data.
    - col_survival: Column name containing the survival data.

    *Returns*
    - event: Concatenated array of survival data.
    - group: Concatenated array of group labels.
    - time_of_event: Concatenated array of time data.
    """
    # Sample from control and experiment groups
    df_control = df[df['control']==0].sample(n=n_mice,replace=True)
    df_experiment = df[df['control']==1].sample(n=n_mice,replace=True)

    # Extract relevant columns from the DataFrame
    E_c = df_control[col_survival].values
    G_c = df_control['control'].values
    T_c = df_control[col_time].values
    E_e = df_experiment[col_survival].values
    G_e = df_experiment['control'].values
    T_e = df_experiment[col_time].values

    # Concatenate arrays
    event = np.concatenate((E_c,E_e))
    group = np.concatenate((G_c,G_e))
    time_of_event = np.concatenate((T_c ,T_e))
    return event,group,time_of_event

def simulation_log_rank(df,length=100,n_mice=10,col_time='time_original',col_survival="survival_original"):
    """
    Runs a simulation of the log-rank test.

    *Arguments*
    - df: DataFrame containing the data.
    - length: Number of simulations to run.
    - n_mice: Number of mice to sample for each simulation.
    - col_time: Column name containing the time data.
    - col_survival: Column name containing the survival data.

    *Returns*
    - stat_p_values: The proportion of simulations with a p-value < 0.05.
    """
    
    # Initialize counter for significant p-values
    stat_p_values = 0

    # Run simulations
    for simu in range(length):
        # Extract data for log-rank test
        event, group, time_of_event = take_data_for_log_rank(df,n_mice,col_time,col_survival)

        # Compute log-rank test
        log_rank = multivariate_logrank_test(time_of_event,group,event)
        p_value = log_rank.p_value

        # Check if p-value is significant
        if p_value < 0.05:
            stat_p_values +=1
        else:
            continue
    return stat_p_values/length


def multiple_simulation_with_results(df,N_simulation,N_repetition,N_mice,cpu=4,col_time='time_original',col_survival="survival_original"):
    """
    Performs multiple simulations of the log-rank test and returns the results.

    *Arguments*
    - df: DataFrame containing the data.
    - N_simulation: Number of simulations to run for each repetition.
    - N_repetition: Number of repetitions.
    - N_mice: Number of mice to sample for each simulation.
    - cpu: Number of CPU cores to use for multiprocessing.
    - col_time: Column name containing the time data.
    - col_survival: Column name containing the survival data.

    *Returns*
    - df_result: DataFrame containing the simulation results.
    """
    # Initialize multiprocessing pool
    pool = mp.Pool(cpu)

    # Perform simulations
    power = pool.starmap(simulation_log_rank,[(df,N_simulation,N_mice,col_time,col_survival) for n in range(N_repetition)])

    # Close multiprocessing pool
    pool.close()

    # Create result DataFrame
    result = {"Power":power,
    "Repetition_No":np.arange(1,N_repetition+1,1),
    "N_simulation":N_repetition*[N_simulation],
    "N_mice":N_repetition*[N_mice],
    "Threshold":N_repetition*[col_time.split("_")[1]]}
    df_result = pd.DataFrame(result)
    return df_result

def automatic_simulation_per_infection(df,N_simulation,N_repetition,mice_range,cpu=4,col_time='time_original',col_survival="survival_original"):
    """
    Performs automatic simulations for each infection type and for multiple number of mice.

    *Arguments*
    - df: DataFrame containing the data.
    - N_simulation: Number of simulations to run for each repetition.
    - N_repetition: Number of repetitions.
    - mice_range: Range of number of mice to sample for each simulation.
    - cpu: Number of CPU cores to use for multiprocessing.
    - col_time: Column name containing the time data.
    - col_survival: Column name containing the survival data.

    *Returns*
    - None
    """
    # Get unique infection types
    infections = df['Infection'].unique().tolist()

    # Loop through each infection type
    for infection in infections:

        # Filter DataFrame for the current infection
        df_infection = df[df['Infection']==infection]

        # print number of control that exist for this type of infection. Useful to check if there is enough eterogenous data in the simulation
        counts = df_infection.groupby(["Infection"])['control'].value_counts()
        print(counts)

        # Start timer
        t1 = time.time()

        # Initialize list to store DataFrames
        df_total = []

        # Loop through each number of mice in mice_range --> Each number of mice N will have M number of simulation with repetition
        for N_mice in mice_range:
            print(N_mice)

            # Start timer for simulation
            tsim1 = time.time()

            # Perform multiple simulations (M times) for the choosen number of mice
            df_result = multiple_simulation_with_results(df_infection,
            N_simulation=N_simulation,
            N_repetition=N_repetition,
            N_mice=N_mice,
            cpu = cpu,
            col_time=col_time,
            col_survival=col_survival)

            # Add infection column to results DataFrame
            df_result['Infection'] = infection

            # Append results DataFrame to df_total list
            df_total +=[df_result]

            # End timer for simulation
            tsim2 = time.time()

            # Print simulation time
            print("Simulation took: ",tsim2-tsim1,"s")

        # Concatenate DataFrames in df_total list    
        df_simulation = pd.concat(df_total,axis=0).reset_index(inplace=False,drop=True)

        # End timer
        t2 = time.time()

        # Print complete simulation time
        print("Complete simulation took: ",t2-t1,"s")

        # Save simulation results to Excel file
        df_simulation.to_excel(
            f"./results/simulations/N_mice_variable_THR_fix/{infection}_Nsimulation{N_simulation}_Nrepetition{N_repetition}_Nmice400.xlsx"
        )



def automatic_simulation_per_infection_fix_mice(df,N_simulation,N_repetition,column_threshold,dict_mice_range,cpu=4):
    """
    Performs automatic simulations for each infection type with fixed number of mice.

    *Arguments*
    - df: DataFrame containing the data.
    - N_simulation: Number of simulations to run for each repetition.
    - N_repetition: Number of repetitions.
    - column_threshold: List of tuples containing the survival and time column names.
    - dict_mice_range: Dictionary containing the infection and corresponding number of mice.
    - cpu: Number of CPU cores to use for multiprocessing.

    *Returns*
    - None
    """
    # Get unique infection types
    infections = df['Infection'].unique().tolist()

    # Loop through each infection type
    for infection in infections:
        # Filter DataFrame for the current infection
        df_infection = df[df['Infection']==infection]

        # print number of control that exist for this type of infection. Useful to check if there is enough eterogenous data in the simulation
        counts = df_infection.groupby(["Infection"])['control'].value_counts()
        print(counts)

        # Start timer
        t1 = time.time()

        # Initialize list to store DataFrames
        df_total = []

        # Get number of mice at the power of 80% for current infection from dict_mice_range
        N_mice = dict_mice_range[infection]
        print("N mice 80% = ",N_mice)
        
        # Loop through each (col_survival, col_time) tuple in column_threshold
        for col_survival,col_time in column_threshold:
            print(col_survival,col_time)

            # Start timer for simulation
            tsim1 = time.time()

            # Perform multiple simulations (M times) for the choosen Threshold (number of mice N is fixed)
            df_result = multiple_simulation_with_results(df_infection,
            N_simulation=N_simulation,
            N_repetition=N_repetition,
            N_mice=N_mice,
            cpu = cpu,
            col_time=col_time,
            col_survival=col_survival)

            # Add infection column to results DataFrame
            df_result['Infection'] = infection

            # Append results DataFrame to df_total list
            df_total +=[df_result]

            # End timer for simulation
            tsim2 = time.time()

            # Print simulation time
            print("Simulation took: ",tsim2-tsim1,"s")

        # Concatenate DataFrames in df_total list
        df_simulation = pd.concat(df_total,axis=0).reset_index(inplace=False,drop=True)

        # End timer
        t2 = time.time()
        # Print complete simulation time
        print("Complete simulation took: ",t2-t1,"s")

        # Save simulation results to Excel file
        df_simulation.to_excel(
            f"./results/simulations/N_mice_fix_THR_variable/{infection}_Nsimulation{N_simulation}_Nrepetition{N_repetition}_fixNmice{N_mice}.xlsx"
        )

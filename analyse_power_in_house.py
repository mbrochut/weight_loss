import pandas as pd
from lifelines.statistics import multivariate_logrank_test



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




def add_p_values_with_H0(df,p_value,alpha=0.05):
    """
    Add p-values to the count if the null hypothesis is rejected based on the given alpha level.

    *Arguments*
    - df: DataFrame object containing the survival data.
    - p_value: p-value computed from a statistical test.
    - alpha: Significance level for hypothesis testing (default is 0.05).

    *Return*
    - 1 if the null hypothesis is rejected (p-value < alpha), otherwise 0.
    """

    # Check if any value in the 'H0' column is less than 0
    if any(df['H0'] < 0):
        # If the p-value is greater than the significance level alpha, return 1
        if p_value > alpha:
            return 1
    # Check if any value in the 'H0' column is greater than 0
    elif any(df['H0'] > 0):
        # If the p-value is less than the significance level alpha, return 1
        if p_value < alpha:
            return 1
    # If neither condition is met, return 0
    return 0


def compute_p_values(df,time_col='time_original',survival_col='survival_original',alpha=0.05):

    """
    Compute p-values for survival data using the multivariate log-rank test or the log-rank analysis for each group in the 'exp' column.
    
    *Arguments*
    - df: DataFrame object containing the survival data.
    - time_col: Name of the column in df containing the time data (default is 'time_original').
    - survival_col: Name of the column in df containing the survival data (default is 'survival_original').
    - alpha: Significance level for hypothesis testing (default is 0.05).

    *Return*
    - A pandas Series with the count of p-values and the number of experiments.
    """

    # Initialize a dictionary to store the results
    dict_result_multi = {"count_p_values": 0, "number_of_experiment": 0}

    # Initialize variables to keep track of count and number of experiments
    count_p_values = 0
    number_of_experiment = 0

    # Get unique values in the 'exp' column
    groups = df['exp'].unique().tolist()

    # Check if there are more than one unique group and the 'exp' column doesn't contain 0
    if len(groups) > 0 and (0 not in groups):
        # Return a pandas Series with the results
        return pd.Series(dict_result_multi)

    # Check if there are exactly two unique groups
    if len(groups) == 2:
        # Check if any value in the 'H0' column is less than 0
        if any(df['H0'] < 0):
            # Return a pandas Series with the results
            return pd.Series(dict_result_multi)
        
        # Compute the p-value using the multivariate log-rank test
        p_value = multivariate_logrank_test(df[time_col], df['exp'].astype(str), df[survival_col]).p_value
        
        # Add p-values to the count if the null hypothesis is rejected
        dict_result_multi["count_p_values"] = add_p_values_with_H0(df, p_value, alpha=alpha)
        
        # Update the number of experiments
        dict_result_multi["number_of_experiment"] = 1
        
        # Return a pandas Series with the results
        return pd.Series(dict_result_multi)

    # If there are more than two unique groups
    else:
        # Iterate through each unique group
        for n in df['exp'].unique():
            # Skip the group if it is 0
            if n == 0:
                continue
            else:
                # Filter the dataframe for the current group and group 0
                experiment = df[df['exp'].isin([0, n])]
                
                # Check if any value in the 'H0' column is less than 0
                if any(experiment['H0'] < 0):
                    # Continue to the next iteration
                    continue
                
                # Perform the log-rank analysis and compute the p-value
                p_value = log_rank_analysis(experiment, time_percentage=time_col, survival_percentage=survival_col)
                
                # Add p-values to the count if the null hypothesis is rejected
                count_p_values += add_p_values_with_H0(experiment, p_value, alpha=alpha)
                
                # Increment the number of experiments
                number_of_experiment += 1
        
        # Update the dictionary with the count and number of experiments
        dict_result_multi['count_p_values'] = count_p_values
        dict_result_multi['number_of_experiment'] = number_of_experiment
        
        # Return a pandas Series with the results
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
    -sum of compute_p_values: number of p_values that are below alpha (reject the null hypothesis) and the number of experiments
    """

    group = df.groupby(['ID_Experiment','sub_exp'])
    result = group.apply(lambda x: compute_p_values(x,time_col=time_col,survival_col=survival_col))
    
    return result.sum()


def analyse_per_model(df, columns_time, columns_survival, model="Listeria", path_to_export=None):
    """
    Analyze survival statistics for a specific model using time and survival columns.
    
    *Arguments*
    - df: DataFrame object containing the dataset.
    - columns_time: List of column names containing time data.
    - columns_survival: List of column names containing survival data.
    - model: Model name for which statistics are computed (default is "Listeria").
    - path_to_export: Path to export the results as an Excel file (optional).

    *Return*
    - result_df: DataFrame containing the computed statistics.
    """
    # Filter the DataFrame for the specific model
    df_to_analyze = df[df['Infection'] == model]
    
    # Initialize an empty dictionary to store the results
    result = {}
    
    # Iterate over the columns of time and survival
    for col_time, col_survival in zip(columns_time, columns_survival):
        # Extract the name from the column name
        name = col_time.split("_")[1]
        print("Thr applied: ", name)
        
        # Compute statistics with H0 for the specified model, time, and survival columns
        stats = stats_with_H0(df_to_analyze, col_time, col_survival)
        
        # Store the statistics in the result dictionary
        result[name] = stats
    
    # Convert the result dictionary to a DataFrame and transpose it
    result_df = pd.DataFrame(result).T
    
    # Export the DataFrame to an Excel file if a path is provided
    if path_to_export:
        result_df.to_excel(path_to_export)
    
    # Return the resulting DataFrame
    return result_df




if __name__ == '__main__':

    df = pd.read_excel('./data/df_for_analysis.xlsx',index_col=0)

    df = df[~df['Experiment'].str.contains("/LD")]


    columns_time = df.columns[df.columns.str.contains('time')].values.tolist()
    columns_survival = df.columns[df.columns.str.contains('survival')].values.tolist()

    columns_time = [columns_time[-1]] + columns_time[:-1]
    columns_survival = [columns_survival[-1]] + columns_survival[:-1]

    
    # For one model
    # analyse_per_model(df,columns_time,columns_survival,path_to_export="./results/in_house/Listeria.xlsx")

    # for all models
    """for model in ['C. albicans','S. pneumoniae','Listeria',"H1N1"]:
        analyse_per_model(df,columns_time,columns_survival,model=model,path_to_export=f"./results/in_house/{model}.xlsx")"""
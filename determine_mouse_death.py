import pandas as pd
import numpy as np

def convert_time_point_in_str_format(df,column='Time_point'):
    """
    Change datetimeIndex inside the list of column "Time_point" into string separate by comma for a better export
    
    *Arguments*
    -df: dataframe object

    *Return*
    nothing, but change the column "Time_point" in the df (list[string] --> string: date_1, date_2, ..., date_n)
    """
    list_dt = []
    for dates in df[column]:
        dt = dates.strftime("%d-%m-%Y")
        dt = ', '.join(dt.values.tolist())
        list_dt += [dt]
    df[column] = list_dt

def add_time(x,last_time):
    """
    calculate for a given mouse (represented by the serie object x) the survival time in days
    *Arguments*
    -x: serie object (pandas object)
    -last_time: datetime object representing the last time recorded in the experiment (end of experiment)

    *Return*
    time: datetime object in days
    """
    last_index_in_Dates = x.index.get_loc(x.last_valid_index())-x.index.get_loc('weight_T_infection')
    last_index = x.index.get_loc(x.last_valid_index())
    try:
        x.Dates[last_index_in_Dates]
    except IndexError:
        print(x['ID_Experiment'])



    first_date_after_infection = x.Dates[x.Dates >x.Time_infection][0]
    initial_time = first_date_after_infection-x.Time_infection
    if x.loc["weight_T_infection":].astype(str).str.contains("ane").any():#if the mouse die due to anesthesia
        return 0

    if x.last_valid_index_original == x.index.get_loc("weight_T_infection"): #if the mouse died the day after first record
        return initial_time.days
    elif x.last_valid_index_original == last_index: #if last record data correspond to original last record data
        time = x.Dates[last_index_in_Dates] -x.Time_infection
        if x.str.contains('sacrified').any():
            return  time.days + 0.5
        return time.days + 1
    elif x.Dates[last_index_in_Dates] <= x.Time_infection: #if the last record time is before time of infection because of THR
        return initial_time.days +0.5
    else:
        time = x.Dates[last_index_in_Dates] -x.Time_infection #if the threshold has an impact on the survivability
        return time.days + 1.5



def add_survival(x,last_time):
    """
    calculate for a given mouse (represented by the serie object x) if the mouse has survived or not during the experiment
    dead = 1
    survive = 0
    *Arguments*
    -x: serie object (pandas object)
    -last_time: datetime object representing the last time recorded in the experiment (end of experiment)

    *Return*
    survival: int (1 or 0) representing if the mouse is dead or alive
    """
    last_index = x.index.get_loc(x.last_valid_index())-x.index.get_loc('weight_T_infection')

    try:
        x.Dates[last_index]
    except IndexError:
        print(x['ID_Experiment'])
    survival = 1
    if x.loc["weight_T_infection":].astype(str).str.contains("dead|sacrified|ane").any() or x.Dates[last_index] <last_time:
        survival = 1
    elif x.str.contains('survive').any() or x.Dates[last_index] ==last_time:
        survival = 0
    else:
        print(x.ID_Experiment,x.index)
    return survival

def survival_time_percentage_weight(x,percentage):
    """
    keep data untile a specific weight loss on a mouse (threshold of sacrifice)
    If a mouse lose a given percentage of its weight during the experiment, datas are kept until the threshold.
    The percentage of weight loss is calculated reltively to the time of infection
    *Arguments*
    -x: serie object (pandas object)
    -percentage: maximum loss weight that a mouse can take.

    *Return*
    x: serie object. Data untile a threshold given by percentage
    """
    y = x.loc['weight_T_infection':] #take datas of serie x
    y = y[~y.astype(str).str.isalpha()] #only take numerical value (if string, ignore)
    y = y/y.weight_T_infection #change the values of weight in percentage (in function of weight at Time of infection)
    z = y[y<1-percentage] #find index max where the weight is under the threshold given by percentag (weight loss max)
    if len(z) !=0:
        x.loc[z.index[0]:] = np.nan #keep data untile the given index, replace the end of the serie by nan values
        return x
    else:
        return x
    
    

    
def add_survival_and_time(df,percentage=0):
    """
    Apply the function add_time, add_survival in function of ID_Experiment and Time_point.
    It is important to group the function by ID_experiment and Time_point in order to correctly take the last_time of each group!
    The maximum of threshold authorized can be pass to the function. It can change how long a mouse can survive or not
    *Arguments*
    -df: dataframe object (pandas object)
    -percentage: int(optional), maximum loss weight that a mouse can take.

    *Return*
    total_time: serie object. serie of length N for DF of length N, representing the number of days a mouse survive during the experiemnt.
    total_survival: serie object. serie of length N for DF of length N, representing if a mouse has survive or not (1, 0)
    """
    groups = df.groupby(['ID_Experiment','Time_point'])
    total_time = pd.Series()
    total_survie = pd.Series()
    for name, group in groups:
        group = group.dropna(how='all',axis=1)
        last_time = group.reset_index().at[0,'Dates'].values[-1]
        if percentage!=0:
            group = group.apply(lambda x : survival_time_percentage_weight(x,percentage),axis=1)
        
        time = group.apply(lambda x : add_time(x,last_time),axis=1)
        survival = group.apply(lambda x : add_survival(x,last_time),axis=1)
        total_time = pd.concat([total_time, time])
        total_survie = pd.concat([total_survie,survival])
    return total_time, total_survie

def multiple_percentage_weigth(df,list_percentage):
    """
    Function to performe multiple time the function add_survival_and_time and store the result in a new df (df_result).
    *Arguments*
    -df: dataframe object (pandas object)
    -list_percentage:list, list of percentage to evaluate

    *Return*
    df_result: dataframe object, create a new dataframe adding the survival and time of survival for each mouse in the df, for a given list of percentage
    """
    print(list_percentage)
    dict_result = {}
    for percentage in list_percentage:
        print("starting percentage ",percentage)
        total_time, total_survie = add_survival_and_time(df,percentage)
        dict_result["time_{}".format(round(percentage,2))] = total_time
        dict_result["survival_{}".format(round(percentage,2))] = total_survie
    result = pd.DataFrame(dict_result)
    df_result = pd.concat([df,result],axis=1)
    print(result)
    print(df)
    return df_result



def format_df_for_p_values(df, end_of_df_column_name = 'weight_T32'):
    """
    Add a column: max_loss_weight_percentage. Compute the maximum loss weight percentage of each mouse during the experiment.
    Add column: exp. Take the first value if it's a string, otherwise the same value. Example: 1A --> 1 example 2: 4 --> 4
    Add column: sub_exp. Take the second value if it's a string, otherwise put "no". Example: 1A --> A example 2: 4 --> no
    *Arguments*
    -df: dataframe object
    -end_of_df_column_name:
    *Return*
    nothing, but add columns to df: max_loss_weight_percentage and exp
    """

    df['max_loss_weight_percentage'] = df.loc[:,'weight_T1':end_of_df_column_name].apply(lambda x: pd.to_numeric(x, errors='coerce')).min(axis=1)/df['weight_T_infection']
    df['exp'] = df['Group'].apply(lambda x:int(x[0]) if isinstance(x,str) else x)
    df['sub_exp'] = df['Group'].apply(lambda x:x[1] if isinstance(x,str) else "no")

def longer_than_N_days(x,N=14):
    """
    keep for a given mouse (represented by the serie object x) the number of datas recorded until a specific date, by default 14.
    So, an experiment can run only until a given date (N=14) then, all the datas are removed
    *Arguments*
    -x: serie object (pandas object)
    -N: int, number of days max in the experiment

    *Return*
    x: serie object
    """
    
    dt = x.Time_infection +pd.Timedelta("{} days".format(N)) #date maximal of recorded datas (datetime object)
    closest_date_after_Ndays = x.Dates.get_indexer([dt], method='nearest')[0] #found the closest date in the serie x of the maximal date (dt)
    last_valid_index = x.index.get_loc('weight_T_infection') + closest_date_after_Ndays #found the maximal index in the serie
    x.iloc[last_valid_index:] = np.nan #remove unnecessary datas until the last authorized index
    x['Dates'] = x['Dates'][0:closest_date_after_Ndays] #remove unnecessary datas until the last authorized index in Dates column (Datetimeindex)
    return x 



if __name__ == '__main__':
    #LOAD DATA
    df = pd.read_excel("./data/weight_loss_raw_data.xlsx")

    # transform str in DateTimeIndex of Time_point column
    serie_dates = df['Time_point'].apply(lambda x: pd.to_datetime(x.split(', '),dayfirst=True))
    index_to_insert_after = df.columns.get_loc("Time_point") # Get the index of the column to insert after
    df.insert(loc=index_to_insert_after+1,column='Dates',value=serie_dates) # add new column Dates that is a DateTimeIndex of Time_point column
    
    # ADD a column that is the last valid index (i.e not a NaN) for further analysis (to compute survival time)
    last_valid_index_original= df.apply(lambda x: x.index.get_loc(x.last_valid_index())+1,axis=1) #!!!add one as we will add a new column before!!!
    df.insert(loc=index_to_insert_after+2,column='last_valid_index_original',value=last_valid_index_original)

    #remove data longer than N days, by default 14 days after infection
    df_Max14days = df.apply(lambda x: longer_than_N_days(x,14),axis=1)
    df_Max14days = df_Max14days.dropna(axis=1,how='all')
    df_Max14days["last_valid_index_original"]= df_Max14days.apply(lambda x: x.index.get_loc(x.last_valid_index()),axis=1)
    #df_Max14days.to_excel('./data/max_14_days.xlsx')
        
    #ADD TIME AND SURVIVAL FOR EACH PERCENTAGE, FEW MINUTES TO RUN
    time, survie = add_survival_and_time(df_Max14days)
    name_last_col = df_Max14days.columns.values[-1]
    print(name_last_col)
    PERCENTAGES = np.arange(0.3,0.04,-0.01)
    #percentage_testing = [0.3,0.15,0.05]
    df_result = multiple_percentage_weigth(df_Max14days,PERCENTAGES)
    df_result['time_original'] = time
    df_result['survival_original'] = survie
    df_result['Time_point'] = serie_dates

    #CONVERT TIME_POINT IN BETTER FORMAT AND ADD INFORMATION ON SUB EXPERIMENTS
    convert_time_point_in_str_format(df_result)
    format_df_for_p_values(df_result,end_of_df_column_name=name_last_col)
    df_result.to_excel("./data/df_for_analysis.xlsx")
import pandas as pd
import openpyxl
from openpyxl import load_workbook
from datetime import datetime
import ast


class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

def end_of_experiment(df):
    """
    Open the original DF with "-" as experiment separator and return the list of index for each separator 

    *Arguments*
    -df: dataframe object (weight_all_mice_manual_formating _CLEAN.xlsx)

    *Return*
    -index_end_experiment: list of index with "-" as separator
    """
    mice = df['Mouse_ID']
    index_end_experiment = []
    for index, value in enumerate(mice):
        if value =='-':
            index_end_experiment+=[index]
    return index_end_experiment

def create_index_experiment(index_end_exp,length_df=3503):
    """
    Create a list of ID of experiment for each mice in the DF. 
    The purpose is to create a new column in the DF with a unique ID for each experiment and assign it to the mice.

    *Arguments*
    -index_end_exp: index of each "-" separator of index (end of each experiemnt)

    *Return*
    -col_exp: list of ID for each experiment.
    """
    col_exp = index_end_exp[0]*['ID_exp_0'] #indice 0 of list
    for index,end in enumerate(index_end_exp[1:]):
        len_exp = end-index_end_exp[index]-1
        string_experiment = 'ID_exp_{}'.format(index+1)
        col_exp = col_exp + len_exp*[string_experiment]

    length_last_experiment = length_df-index_end_exp[-1] - 1
    col_exp = col_exp + length_last_experiment*['ID_exp_{}'.format(index+2)] #last indice to add
    return col_exp

def controle_Date(df,column = 'Date'):

    """
    Control if there is different dates in a groupby object (id_experiment and group)

    *Arguments*
    -df: dataframe object
    -colum: name of column of the dataframe with date to controle

    *Return*
    nothing
    """

    print(f"{bcolors.HEADER}TEST CONTROL {column}:{bcolors.ENDC}")
    group = df.groupby(['ID_Experiment','Group'])
    check = True
    for experiment in group:
        dates = experiment[1][column].values.tolist()
        result = dates.count(dates[0]) == len(dates)
        if result:
            continue
        else:
            print("problem length of dates at ")
            print(experiment[1].head(1))
            check = False
    if check:
        print(f"{bcolors.OKGREEN}Validate{bcolors.ENDC}")
    else:
        print(f"{bcolors.FAIL}NOT PASS{bcolors.ENDC}")


def controle_unique_name(df,column_name = 'Pre_traitment'):
    """
    Print the list of unique entry in a column (pre_traitment by default). Need to verify by hand if everything is well name

    *Arguments*
    -df: dataframe object
    -colum: name of column of the dataframe to controle unique value

    *Return*
    nothing
    """
    print(f"{bcolors.HEADER}TEST CONTROL UNIQUE NAME {column_name}:{bcolors.ENDC}")
    unique = df[column_name].value_counts()
    print(unique)

def controle_time_pre_traitement(df):
    """
    Control if there is different time in a groupby object (id_experiment and group) for the column Time_pre_traitment

    *Arguments*
    -df: dataframe object

    *Return*
    nothing
    """

    print(f"{bcolors.HEADER}TEST CONTROL TIME PRE TREATMENT:{bcolors.ENDC}")
    group = df.groupby(['ID_Experiment','Group'])
    check = True
    for index, experiment in group:
        unique_date = experiment['Time_pre_traitment'].unique()
        if len(unique_date) >= 2:
            print("PROBLEM TIME IN PRE-TREATMENT")
            print(experiment.head(1))
            print(unique_date)
            check = False
    if check:
        print(f"{bcolors.OKGREEN}Validate{bcolors.ENDC}")
    else:
        print(f"{bcolors.FAIL}NOT PASS{bcolors.ENDC}")


def try_parsing_date(text):
    """
    check multiple date format in string. If found one, return the date found.

    *Arguments*
    -text: string of text containing dates 

    *Return*
    datetime object
    """
    for fmt in ('%Y-%m-%d', '%d.%m.%Y', '%d/%m/%Y','%d-%m-%Y','%d.%m.%y','%d/%m/%y'):
        try:
            return datetime.strptime(text, fmt)
        except ValueError:
            pass
    raise ValueError('no valid date format found',text)


def uniform_time_infection(df):
    """
    Put in datetime format the elements in the column "Time_infection"

    *Arguments*
    -df: dataframe object

    *Return*
    nothing, but change the column "Time_infection" in the df (string --> datetime)
    """
    infecttion = df['Time_infection']
    list_time_infection = []
    for index, value in enumerate(infecttion):
        if isinstance(value,datetime):
            list_time_infection +=[value]
            continue
        else:
            date = try_parsing_date(value)
            list_time_infection +=[date]
    df['Time_infection'] = list_time_infection


def uniform_time_point_in_list(df):
    """
    Change string format of "Time_point" column in list format
    
    *Arguments*
    -df: dataframe object

    *Return*
    nothing, but change the column "Time_point" in the df (string --> list)
    """
    time_point = df['Time_point']
    list_time_point = []
    for index, values in enumerate(time_point):
        if isinstance(values,str) and values[0] =='[':
            x = ast.literal_eval(values)
            list_time_point+=[x]
        elif isinstance(values,str):
            result_split = values.split(", ")
            list_time_point+=[result_split]
        else:
            print("this is not a string",index)
    df['Time_point'] = list_time_point

def change_time_point_in_date_format(df):
    """
    Change string format inside the list of column Time_point in list of datetime
    
    *Arguments*
    -df: dataframe object

    *Return*
    nothing, but change the column "Time_point" in the df (list[string] --> list[datetime])
    """
    time_point = df['Time_point']
    serie_time  = []
    for time in time_point:
        list_time = []
        for point in time:
            date = try_parsing_date(point)
            list_time+=[date]
        serie_time+=[list_time]
    df['Time_point'] = serie_time


def add_time_infection_for_beginning_of_df(df):
    """
    Untile the mouse "TRO-15172", Time_point didn't have the time of infection included in the time point.
    In order to have all the different time point (beginning with the time of infection) we add this time to the list of datetime object
    in the "Time_point" column
    
    *Arguments*
    -df: dataframe object

    *Return*
    nothing, but change the column "Time_point" by adding one datetime to the list
    """
    final_mouse = "TRO-15172"
    print("INDEX : ",df[df['Mouse_ID'].str.contains(final_mouse)].index.values[0])
    index_final_mouse = df[df['Mouse_ID'].str.contains(final_mouse)].index.values[0]
    df.loc[:index_final_mouse,'Time_point'] = df.loc[:index_final_mouse,:].apply(lambda x: [pd.to_datetime(x.Time_infection)]+x.Time_point,axis=1)
    df.loc[:,'Time_point'] =df.loc[:,'Time_point'].apply(lambda x: pd.to_datetime(x))

def check_time(df):
    """
    Control if there is multiple year AND if the time is continuous in Time_point list

    *Arguments*
    -df: dataframe object

    *Return*
    nothing
    """
    print(f"{bcolors.HEADER}TEST CONTROLE CONTINUOUS TIME:{bcolors.ENDC}")
    groups = df.groupby(['ID_Experiment'])
    check = True
    for name, group in groups:
        date_index = group['Time_point'].values
        first_dates = date_index[0]
        year = first_dates.year.tolist()
        dates_ordinal = [int(d.toordinal()) for d in first_dates]
        sort = dates_ordinal[:]
        sort.sort()
        if len(set(year)) >1: #Check if there is multiple year 
            print(name)
            print('year')
            check=False
            continue
        elif dates_ordinal != sort: #Check if time is continuous 
            print(name)
            print('non continuous')
            check=False
    if check:
        print(f"{bcolors.OKGREEN}Validate{bcolors.ENDC}")
    else:
        print(f"{bcolors.FAIL}NOT PASS{bcolors.ENDC}")

def control_length_of_time_point(df,column_name='Time_point'):
    """
    Control if there is enough time point for each weight datas
    *Arguments*
    -df: dataframe object

    *Return*
    nothing
    """
    print(f"{bcolors.HEADER}TEST CONTROLE LENGTH TIME POINT:{bcolors.ENDC}")
    check = True
    groups = df.groupby(['ID_Experiment','Time_point'])
    for name, group in groups:
        datas = group.loc[:,'weight_T_infection':].dropna(axis=1,how='all')
        len_column = len(datas.columns)
        if isinstance(group.reset_index().at[0,column_name],str):
            len_time_point = len(group.reset_index().at[0,column_name].split(', '))
        else:
            len_time_point = len(group.reset_index().at[0,column_name])
            
        if len_column == len_time_point:
            continue
        elif len_time_point > len_column:
            print(f"{bcolors.WARNING}more date point than actual datas at: {bcolors.ENDC}",name)
        else:
            print(name)
            print(datas)
            print("columns length ",len_column)
            print("time length ", len_time_point)
            check = False
    if check:
        print(f"{bcolors.OKGREEN}Validate{bcolors.ENDC}")
    else:
        print(f"{bcolors.FAIL}NOT PASS{bcolors.ENDC}")


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

def string_to_list(s):
        if '-' in s:
            return []
        elif '[' in s:
            s = s.replace("nan", "0.0")
            return ast.literal_eval(s)
        else:
            return [float(x.strip()) for x in s.split(',')]

def convert_scores_in_list_format(df):
    df['Scores'] = df['Scores'].apply(string_to_list)
    return df

if __name__ == '__main__':

    path = "./data/raw_data.xlsx"
    df = pd.read_excel(path) #read an excel file

    index_end = end_of_experiment(df)
    col_experiment = create_index_experiment(index_end,length_df=len(df))
    df = df.dropna(subset=['Strain']) # remove column with no informations (column with "-" separator)
    df = df.reset_index(drop=True) #reset index
    df.insert(loc=0,column='ID_Experiment',value=col_experiment) #Add the Experiment ID for each mouse
    
    #UNIFORMISE NAME IN PRE_TREATMENT COLUMN
    dict_to_replace = {'Propionate':'propionate','zymosan':'training/zymosan',
    'training/zymosan + anakinra':'training/zymosan-anakinra','training':'training/zymosan'}
    df['Pre_traitment'] = df['Pre_traitment'].replace(dict_to_replace)
    
    uniform_time_infection(df) # uniformise time in column time_infection
    uniform_time_point_in_list(df) # uniformise time in column time_point
    change_time_point_in_date_format(df) # change time_point in date format
    add_time_infection_for_beginning_of_df(df) # add time of infection to time point until INDEX :  1089 (then it is already added)
    
    #BEGINNING OF TESTING
    controle_Date(df)
    controle_unique_name(df,column_name = 'Pre_traitment')
    controle_time_pre_traitement(df)
    controle_unique_name(df,column_name = 'Infection')
    controle_Date(df,'Time_infection')

    check_time(df) #test if time point are continuous and not multiple year
    
    #CONVERT TIME_POINT IN BETTER FORMAT
    convert_time_point_in_str_format(df)
    control_length_of_time_point(df)

    df = convert_scores_in_list_format(df)
    df = df[df['Infection'].isin(['C. albicans','S. pneumoniae','Listeria',"H1N1"])]
    name_to_export = 'data/weight_loss_raw_data.xlsx'
    df.to_excel(name_to_export)
    
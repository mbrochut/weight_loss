import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
import os
import seaborn as sns
from ast import literal_eval
import ast
import re

def cubic(x, a, b, c,d):
    return a*x**3 + b*x**2 +c*x + d

def sigmoid(x, x0,k):
    y = 1 / (1 + np.exp(-k*(x-x0)))
    return (y)

def sigmoid2(x, L ,x0, k, b):
    y = L / (1 + np.exp(-k*(x-x0))) + b
    return (y)

def find_roots(x,y):
    s = np.abs(np.diff(np.sign(y))).astype(bool)
    return x[:-1][s] + np.diff(x)[s]/(np.abs(y[1:][s]/y[:-1][s])+1)



def import_dataframe(folder_name,keys=[0,5]):
    All_names = os.listdir(folder_name)
    dict_df = {}
    for name in All_names:
        split_name = name.split('_')
        name_df = "_".join([split_name[key] for key in keys])
        df = pd.read_excel(folder_name+name,index_col=0)
        dict_df[name_df] = df.T
    df = pd.concat(dict_df,axis=1)
    return df




def plot_model(df,name = 'Albicans',save = False,fill_between = False,function_name='sigmoid2',x_func_max=200,show=False,title=False):
    #Selecting datas
    df = df[name]
    idx = pd.IndexSlice
    #complete_name = 'NotCensored_{}'.format(name)
    mean = df.loc[idx[:x_func_max],idx[:,'mean_percentage']] #datas
    SD = df.loc[idx[:x_func_max],idx[:,'SD']] #standard deviation
    
    #reindex for a better visualisation
    
    new_index = ['Thr_15', 'Thr_10', 'Thr_25', 'Thr_20','Thr_original']    
    mean.columns.set_levels(new_index,level=0,inplace=True)
    print(mean)
    SD.columns.set_levels(new_index,level=0,inplace=True)
    #curve fitting
    #choose function
    if function_name == 'cubic':
        func = cubic
    elif function_name == 'sigmoid':
        func = sigmoid      
    elif function_name == 'sigmoid2':
        func = sigmoid2
        
    x_func = np.arange(0,x_func_max,1)
    #plot the graph
    fig, ax = plt.subplots(figsize=(14, 8))
    color = ['blue','red','green','orange','purple']
    for index,experiment in enumerate(['Thr_10', 'Thr_15', 'Thr_20', 'Thr_25','Thr_original']):
        x = mean.index.values
        y = mean[experiment].values.astype(float)[:,0]
        sd = SD[experiment].values.astype(float)[:,0]
        
        if fill_between:
            ax.fill_between(x, y-sd, y+sd, alpha=.5, linewidth=0)

        if function_name == 'sigmoid2':
            p0 = [max(x), np.median(x),1,min(y)] # this is an mandatory initial guess
            popt,cov = curve_fit(func,x, y,p0,method='dogbox',maxfev=10000)
        else:
            popt,cov = curve_fit(func,x, y)

        ax.plot(x, y, linewidth=2,label=experiment)  
        ax.plot(x_func, func(x_func, *popt),color=color[index],alpha=0.4)
    if title:
        ax.set_title('Power, model {}'.format(name))
    ax.set_xlabel("Number of animals used")
    ax.set_ylabel("Power (% p_values < 0.05)")
    ax.legend()
    if save:
        fig.savefig('./graphiques/Simulation Power/N_animal_variableTHR/models/{}_with_{}_fitting.jpg'.format(name,function_name))
    if show:
        plt.show()



def import_models(path,keys = [0,5]):
    model_name = os.listdir(path)
    list_df = []
    for name in model_name:
        df = import_dataframe(path+name+'/',keys=keys)
        list_df +=[df]
    df = pd.concat(list_df,axis=1,keys=model_name)
    return df


def fit_sigmoid(x,y,max_x = 450, bounds = (0,450)):
    p0 = [max(x), np.median(x),1,0] # this is an mandatory initial guess
    popt,cov = curve_fit(sigmoid2,x, y,p0,method='trf',maxfev=10000,bounds=bounds)
    x_func = np.arange(0,max_x,1)
    y_func = sigmoid2(x_func,*popt)
    return x_func, y_func, popt


def find_number_of_mice_at_specific_power(df,model = "Listeria",THR = "THR_original",power = 80,plot=False):

    #!!!BE CAREFUL ALL IS ALL THE SIMULATION COMBINED BUT NOT THE SIMULATION WITH ALL INFECTION AT ONCE!!!
    if model == "all":
        df_model = df
    else:
        df_model = df[df["infection"]==model]
    
    print(df_model)
    mean_thr = df_model.groupby(["N_mice","threshold"]).mean(numeric_only=True)
    mean_thr = mean_thr.reset_index()
    mean_thr = mean_thr[mean_thr['threshold']==THR]

    x = mean_thr['N_mice'].values
    y = mean_thr['percentage'].values

    x_func,y_func,popt = fit_sigmoid(x,y)

    #plot is only used to check if the fitting is correct.
    if plot:
        fig, ax = plt.subplots()
        ax.plot(x, y)
        ax.plot(x_func, sigmoid2(x_func, *popt),alpha=0.4)
        plt.show()


    z = find_roots(x_func,y_func-power)
    return z[0]

def find_number_of_mice_at_specific_power_V2(df,model = "Listeria",power = 80,plot=False):

    #!!!BE CAREFUL ALL IS ALL THE SIMULATION COMBINED BUT NOT THE SIMULATION WITH ALL INFECTION AT ONCE!!!
    if model == "all":
        df_model = df
    elif model == "C. albicans":
        return 0
    else:
        df_model = df[df["Infection"]==model]
    
    print(df_model)
    mean_infection = df_model.groupby(["N_mice"]).mean(numeric_only=True)
    mean_infection = mean_infection.reset_index()

    x = mean_infection['N_mice'].values
    y = mean_infection['Power'].values*100
    print(x,y)
    x_func,y_func,popt = fit_sigmoid(x,y)

    #plot is only used to check if the fitting is correct.
    if plot:
        fig, ax = plt.subplots()
        ax.plot(x, y)
        ax.plot(x_func, sigmoid2(x_func, *popt),alpha=0.4)
        plt.show()


    z = find_roots(x_func,y_func-power)
    
    return z[0]

def find_number_of_mice_all_models(df,columns_infection="infection",V2=False):
    models = df[columns_infection].unique()
    results = {}
    for model in models:
        if V2:
            mice_80 = find_number_of_mice_at_specific_power_V2(df,model)
        else:
            mice_80 = find_number_of_mice_at_specific_power(df,model)
        results[model] = mice_80
    print(results)



def reindex_model_df(df):
    #remove multiindex and put it in tidy format
    df_stack = df.stack([0,1,2]).rename_axis(index=['N_mice','infection','threshold','stats']).reset_index()
    
    #rename threshold for better visualisation
    replace_thr = {"NotCensored_THR0.1":"THR_010",
    "NotCensored_THR0.15":"THR_015",
    "NotCensored_THR0.25":"THR_025",
    "NotCensored_THR0.2":"THR_020",
    "NotCensored_THRoriginal":"THR_original",
    "THR0.1":"THR_010",
    "THR0.15":"THR_015",
    "THR0.25":"THR_025",
    "THR0.2":"THR_020",
    "THRoriginal":"THR_original"}
    df_stack['threshold'] = df_stack['threshold'].replace(replace_thr)
    df_stack = df_stack.rename(columns={0:'values'})
    
    #take only p_values results
    df_p_values = df_stack[df_stack['stats'] == "p_values_good_and_bad"]

    #convert the string containing list of tupple in list of tupple
    df_p_values['values'] = df_p_values['values'].apply(lambda x: ast.literal_eval(x))

    #transform list in 10 columns of values (wide format)
    col_values = ["value_{}".format(n) for n in range(10)]
    df_p_values_2 = pd.DataFrame(df_p_values['values'].to_list(), columns=col_values)
    df_p_values = df_p_values.reset_index(drop=True)
    df_p_values =df_p_values.drop('values',axis=1)
    df_p_values = pd.concat([df_p_values,df_p_values_2],axis=1)
    
    #transform wide format to long format
    df_p_values = df_p_values.melt(id_vars=["N_mice","infection","threshold","stats"],value_vars=col_values,var_name="p_values",value_name="values")
    
    #Add percentage of Pvalues
    df_p_values['percentage'] = df_p_values['values'].apply(lambda x: x[0]/500*100)
    return df_p_values
    
def plot_new_index_model(df,model = "Albicans",N_mice = 120,plot = False,save = None):
    df_model = df[df['infection']==model]
    df_model = df_model[df_model['N_mice']<N_mice]
    plt.figure(figsize=(16,9))
    sns.lineplot(df_model,x="N_mice",y="percentage",hue="threshold",errorbar="sd")
    plt.xticks(rotation=45, ha='right')
    plt.xlabel('Number of Mice', fontsize=16)
    plt.ylabel('Power [%]', fontsize=16)
    plt.title(model, fontsize=20)
    plt.subplots_adjust(top=0.9,bottom=0.15)
    if plot:
        plt.show()
    
    if save:
        print("the file is saved at: ",save)
        plt.savefig(save)


def import_excel_from_simulation(path):
    df = pd.read_excel(path,index_col=0)
    df = df.T
    
    df["p_values_good_and_bad"] = df['p_values_good_and_bad'].apply(lambda x: ast.literal_eval(x))

    col_values = ["value_{}".format(n) for n in range(10)]
    df_extended_from_list = pd.DataFrame(df['p_values_good_and_bad'].to_list(), columns=col_values,index=df.index)
    df_extended_from_list = df_extended_from_list.reset_index().rename(columns={'index':"threshold"})
    df_extended_from_list = df_extended_from_list.melt(id_vars=['threshold'],value_vars=col_values,var_name="p_values",value_name="values")
    df_extended_from_list['percentage'] = df_extended_from_list['values'].apply(lambda x: x[0]/(x[0]+x[1]))

    return df_extended_from_list

def import_simulation_from_folder(path_folder):
    list_file_in_directory= os.listdir(path_folder)
    fix_n_infections = {}
    for file in list_file_in_directory:
        if file[-4:] == "xlsx":
            infection_name = file.split("_")[0]
            df = import_excel_from_simulation(path_folder+file)
            fix_n_infections[infection_name] = df
        else:
            continue
        
    
    df_models_fix =pd.concat(fix_n_infections)
    df_models_fix = df_models_fix.reset_index().rename(columns={"level_0":"infection"})
    return df_models_fix


def graph_fix_animals(df,plot = False,save =  False,path_to_save="./graphiques/NEW_ANALYSIS/Simulation/Fix_N_Mice/Fix_mice.png"):
    plt.figure(figsize=(16,9))
    sns.lineplot(df,x="threshold",y="percentage",errorbar="sd",hue='infection')
    plt.xticks(rotation=45)
    plt.xlabel('Threshold', fontsize=16)
    plt.ylabel('Power [%]', fontsize=16)
    plt.title("Power with fix number of mice", fontsize=20)
    plt.subplots_adjust(top=0.9,bottom=0.15)

    if save:
        plt.savefig(path_to_save)
    if plot:
        plt.show()

def figure_report_variable_N(df):

    sns.lineplot(df,x="N_mice",y="percentage",hue="infection")
    #plt.xticks(rotation=45, ha='right')
    plt.xlabel('Number of Mice', fontsize=16)
    plt.ylabel('Power [%]', fontsize=16)
    plt.title("Simulation", fontsize=20)
    plt.subplots_adjust(top=0.9,bottom=0.15)
    plt.show()

def import_simulations(folder_name):
    All_names = os.listdir(folder_name)
    list_df = []
    for name in All_names:
        if name[-4:] != "xlsx":
            continue
        else:
            print(name)
            split_name = name.split('_')[0]
            df = pd.read_excel(folder_name+name,index_col=0)
            df['Infection'] = split_name
            list_df += [df]
    df = pd.concat(list_df,axis=0).reset_index(drop=True,inplace=False)
    return df

if __name__ == '__main__':

    
    #---------TO GENERATE GRAPH WITH OLD EXCELS-----------
    """path = './simulations/Change_Animals_change_Thr/models/'
    df_models = import_models(path,keys=[0,5])
    new_index_models = reindex_model_df(df_models)
    plot_new_index_model(new_index_models,model='Albicans')
    find_number_of_mice_at_specific_power(new_index_models,model="Ecoli",allModel=False)"""
    
    #---------NEW GRAPHS------------
    #path to new simulations
    path_models = "./simulations/per_models/"

    #import simulations
    df_models = import_models(path=path_models,keys=[1])

    #transform to tidy data
    new_index_models = reindex_model_df(df_models)

    #to find number of mice at specific power, ONE MODEL
    """find_number_of_mice_at_specific_power(new_index_models,model="Listeria",plot=False)"""
    #to find number of mice at specific power, ALL MODEL
    print(new_index_models)
    find_number_of_mice_all_models(new_index_models)
    
    #plot one model
    new_index_models = new_index_models[new_index_models['infection'].isin(["Albicans","Pneumoniae","Listeria"])]
    new_index_models = new_index_models[new_index_models['threshold']=="THR_original"]
    new_index_models = new_index_models[new_index_models['N_mice']<=200]
    #figure_report_variable_N(new_index_models)
    #plot_new_index_model(new_index_models,model='All-models',N_mice=450,plot=True)
       
    #--------FIX NUMBER OF MICE---------
    #import datas
    path_fix_n = "./simulations/fix_NAnimals/per_models/"
    df_models_fix = import_simulation_from_folder(path_fix_n)

    #for better visualisation
    df_models_fix['threshold'] = df_models_fix['threshold'].apply(lambda x : x.split("_")[1])
    df_models_fix = df_models_fix[df_models_fix['infection'].isin(["Albicans","Pneumoniae","Listeria"])]
    #graph_fix_animals(df_models_fix,save=False,plot=True)

    new_folder_name = "./simulations/new_simulations/"
    df = import_simulations(new_folder_name)
    #df = df[df['Infection'].isin(['E. coli',"Listeria","S. pneumoniae","Staphylococcus aureus"])]
    print(df)
    #find_number_of_mice_all_models(df,columns_infection="Infection",V2=True)
    """sns.lineplot(df,x="N_mice",y="Power",hue="Infection",errorbar='sd')
    plt.show()"""
    """for infection in ['E. coli',"Listeria","S. pneumoniae","Staphylococcus aureus"]:
        print(infection)
        sns.regplot(df[df['Infection']==infection],x="N_mice",y="Power",logistic=True)
    plt.show()"""


    fix_n_path = "./simulations/new_simulations/fixN/"
    df_fix_N = import_simulations(fix_n_path)
    print(df_fix_N)
    sns.lineplot(df_fix_N,x='Threshold',y='Power',hue='Infection')
    plt.show()
    
    



    #----------to generate graphique for all infection VARIABLE MCIE---------------
    """infections = new_index_models['infection'].unique()
    folder_path = "./graphiques/NEW_ANALYSIS/Simulation/Variable_N_Mice/"
    for infection in infections:
        print(infection)
        path_to_save = folder_path+infection+".png"
        if infection in ["All-models",'Albicans']:
            N_mice = 450
        else:
            N_mice = 200
        plot_new_index_model(new_index_models,infection,N_mice=N_mice,save=path_to_save)"""


    #----------OLD CODE-----------
    """plot_model(df_models,name='Listeria',show=False,function_name='sigmoid',x_func_max=300,save=True,title=True)
    plot_model(df_models,name='Pneumoniae',show=False,function_name='sigmoid',x_func_max=300,save=True,title=True)
    plot_model(df_models,name='Ecoli',show=False,function_name='sigmoid',x_func_max=450,save=True,title=True)
    plot_model(df_models,name='Albicans',show=False,function_name='sigmoid',x_func_max=300,save=True,title=True)
    plot_model(df_models,name='H1N1',show=False,function_name='sigmoid',x_func_max=300,save=True,title=True)
    plot_model(df_models,name='Staphylococcus',show=False,function_name='sigmoid',x_func_max=300,save=True,title=True)
    plot_model(df_models,name='Staphylococcus-Long',show=False,function_name='sigmoid',x_func_max=200,save=True,title=True)
    plot_model(df_models,name='Staphylococcus-Short',show=False,function_name='sigmoid',x_func_max=200,save=True,title=True)
    plot_model(df_models,name='Staphylococcus-Long10',show=False,function_name='sigmoid',x_func_max=200,save=True,title=True)
    plot_model(df_models,name='Staphylococcus-Short10',show=False,function_name='sigmoid',x_func_max=200,save=True,title=True)"""
import pandas as pd
import os
import matplotlib.pyplot as plt
import seaborn as sns

def prepare_DF_p_values(path_to_folder):
    list_dir = os.listdir(path_to_folder)
    export_df = pd.DataFrame()
    list_dir = [x for x in list_dir if os.path.isfile(path_to_folder+x)]
    for name in list_dir:
        name_variable = name[:-5]
        print(name_variable)
        serie = pd.read_excel(path_to_folder+name,index_col=0,names=[name_variable])
        serie = pd.concat([pd.Series(serie.iloc[-1,0],['survival_Brut'],name=name_variable),serie[name_variable][0:-1]])
        export_df = pd.concat([export_df,serie],axis=1)
    # checking if the directory demo_folder 
    # exist or not.
    if not os.path.exists(path_to_folder+'/dataframe_concat'):
        
        # if the demo_folder directory is not present 
        # then create it.
        os.makedirs(path_to_folder+'/dataframe_concat')
    
    export_df.to_excel(path_to_folder+'/dataframe_concat/Dataframe_p_values.xlsx')


def graphique_p_values(df,col1='Censored_All',col2='NotCensored_All',name='All Experiment',to_save='./graphiques/p_values_graphiques/0.05_bonferroni',alpha=0.05,title=False):
    df[[col1,col2]].plot(figsize=(16, 9))
    if title:
        plt.title('% P_values < {} in function of sacrificed threshold\n{}'.format(alpha,name))
    plt.ylabel('Pecentage of P-values < {}'.format(alpha))
    plt.xlabel('Sacrificed threshold')
    plt.savefig(to_save)

def graphique_from_infection(df):
    sns.set_theme(style="whitegrid",palette="deep",font_scale=2)
    ax = sns.lineplot(data=df, linewidth=1.5,dashes=False)
    ax.set(ylabel='Nombre de p-values < 0.05')
    ax.set(xlabel='Seuil de perte de poids maximal')
    
    #set x ticks
    ax.set_xticks([0,6,11,16,21,26]) # <--- set the ticks first
    ax.set_xticklabels(['Originale','25%','20%','15%','10%','5%'])
    plt.show()


if __name__ == '__main__':
    print("compute graphique")

    #without H0
    NO_H0 = False
    if NO_H0:
        #models = ['E. coli', 'S. pneumoniae', 'Listeria','Staphylococcus aureus','H1N1','C. albicans']
        #models = ['Listeria','S. pneumoniae','H1N1','E. coli']
        models = ['S. pneumoniae', 'Listeria','C. albicans']
        df_list = []
        for model in models:
            print(model)
            path = './result_p_values/models/{}_2023_new_time.xlsx'.format(model)
            print(path)
            df = pd.read_excel(path,index_col=0)
            df = df.rename(columns={'':'thr',0:model})
            df['THR'] = df.index.str.split('_',expand=True).droplevel(level=0)
            df_graph = df.iloc[0:-1,:]
            df_graph = df_graph.set_index('THR')
            df_list +=[df_graph]
        result = pd.concat(df_list,axis=1)
        result = result.rename(columns={'Listeria': 'L. monocytogenes'})
        graphique_from_infection(result)

        print(result)
    
    #WITH H0
    H0 = True
    if H0:
        models = ['Candida','Listeria','Pneumo',"H1N1"]
        list_df = []
        for model in models:
            path = f"./result_p_values/WITH_H0/{model}.xlsx"
            print(path)
            df = pd.read_excel(path,index_col=0)
            list_df += [df]
        df = pd.concat(list_df,keys=models)
        df = df.reset_index().rename(columns={"level_0":"Infection","level_1":"Threshold"})
        print(df)
        infection_name = {"Candida":"C. albicans","Listeria":"L. monocytogenes","Pneumo":"S. pneumoniae"}
        df['Infection'] = df['Infection'].replace(infection_name)
        df['Threshold'] = df['Threshold'].replace({"original":"raw data"})
        infection_of_interest = ["C. albicans","L. monocytogenes","S. pneumoniae","H1N1"]
        palette = sns.color_palette("Set2")
        plt.figure(figsize=(12,8))
        sns.lineplot(df,x="Threshold",y="count_p_values",hue="Infection",palette=palette,linewidth = 3,hue_order=infection_of_interest)
        plt.xticks(rotation=45, ha='right',fontsize=16)
        plt.yticks(fontsize=16)
        plt.xlabel('Threshold of sacrifice', fontsize=18)
        plt.ylabel('# p-values < 0.05', fontsize=18)
        plt.title("", fontsize=22)
        plt.show()
    

    #----------------OLD CODE-------------------
    """path_folder_to_concat_DF ='./result_p_values/14_days/'

    prepare_DF_p_values(path_to_folder=path_folder_to_concat_DF)

    df_censored_dead = pd.read_excel('./result_p_values/14_days/dataframe_concat/Dataframe_p_values.xlsx',index_col=0)
    #dict_index = {'All':228,'Acute':136,'Chronique':57,'Training':127,'TrainingAcute':108,'TrainingChronique':13,'WithoutTraining':101}
    dict_index = {'Chronique':57}
    for index , col_name in enumerate(df_censored_dead):
        if not col_name in dict_index:
            continue
        else:
            print(index,col_name)
            index_dict = col_name.split('_',1)[1]
            total_number_p_val = dict_index[index_dict]
            df_censored_dead[col_name] = df_censored_dead[col_name]/total_number_p_val
    
    alpha = '0.05' 
    folder = '14_days'
    # checking if the directory demo_folder 
    # exist or not.
    if not os.path.exists('./graphiques/p_values_graphiques/{}'.format(folder)):
        
        # if the demo_folder directory is not present 
        # then create it.
        os.makedirs('./graphiques/p_values_graphiques/{}'.format(folder))
    for index, key in enumerate(dict_index):
        print(key)
        col1 = 'Censored_'+key
        col2 = 'NotCensored_'+key
        graphique_p_values(df_censored_dead,col1=col1,col2=col2,name=key,to_save='./graphiques/p_values_graphiques/{}/Censored_vs_Uncensored_{}'.format(folder,key),alpha=alpha)"""
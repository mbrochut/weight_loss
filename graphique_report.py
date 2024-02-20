import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os
import numpy as np
from lifelines import KaplanMeierFitter
from lifelines.statistics import multivariate_logrank_test
from PIL import Image
import imageio
from statannotations.Annotator import Annotator



def create_gif(input_folder, output_file, duration=0.2):
    images = []
    fileNames = os.listdir(input_folder)
    fileNames = fileNames[::-1]
    # Open all images in the folder
    for filename in fileNames:
        if filename.endswith(".png") or filename.endswith(".jpg"):
            filepath = os.path.join(input_folder, filename)
            img = Image.open(filepath)
            images.append(img)

    # Save as GIF using imageio
    imageio.mimsave(output_file, images, duration=duration)

def import_simulations(folder_name):
    All_names = os.listdir(folder_name)
    list_df = []
    for name in All_names:
        if name[-4:] != "xlsx":
            continue
        else:
            split_name = name.split('_')[0]
            df = pd.read_excel(folder_name+name,index_col=0)
            df['Infection'] = split_name
            list_df += [df]
    df = pd.concat(list_df,axis=0).reset_index(drop=True,inplace=False)
    return df

def import_weight_data(path_name,infection_of_interest):
    df = pd.read_excel(path_name,index_col=0)
    df['Time_infection'] = df['Time_infection'].dt.strftime("%Y-%m-%d")
    df['Time_point'] = df['Time_point'].apply(lambda x: pd.to_datetime(x.split(','),dayfirst=True,exact=False))
    df['control'] = df['exp'].apply(lambda x: 0 if x==0 else 1)
    df = df[df['Infection'].isin(infection_of_interest)]

    return df

def prepare_df_time(df,infection_to_chose=['E. coli',"Listeria","S. pneumoniae","Staphylococcus aureus"]):
    df['Control'] = df['exp'].apply(lambda x: 0 if x == 0 else 1)
    
    df['Death'] = df.apply(lambda x: 1  if x.survival_original ==1 or x.Time != x.t_origin else 0,axis=1)
    df = df[df['Infection'].isin(infection_to_chose)]
    unique = df['Threshold'].unique()
    
    thr = [float(n.split("_")[1]) if n!= "time_original" else 1 for n in unique]
    dict_change_time = dict(zip(unique,thr))
    df['Threshold'] = df['Threshold'].replace(dict_change_time)
    return df

def weight_distribution(df):
    column_survival = [n for n in df.columns.tolist() if "survival" in n]
    column_survival =  [column_survival[-1]] + column_survival[0:-1]
    data_thr_death  = df.loc[:,column_survival]
    data_thr_death[['Infection','control']] = df[['Infection','control']]
    group_sum = data_thr_death.groupby(["Infection",'control']).sum()
    count = data_thr_death.groupby(["Infection",'control']).count()
    group_sum = group_sum/count
    group_sum.reset_index(inplace=True)
    print(group_sum)
    longer = pd.wide_to_long(group_sum,"survival",['Infection',"control"],'Threshold',sep="_",suffix=".*").reset_index()
    longer = longer[longer['Infection']=="Listeria"]
    sns.barplot(longer,x='Threshold',y='survival',hue="control")
    plt.show()


def violin_plot(df,infection,max_time=10,save = False,path_to_save='',show=False):
    df_infection = df[df['Infection'] == infection]
    df_infection = df_infection[(df_infection['Time']>0) & (df_infection['Time']<=max_time)]
    plt.figure(figsize=(12,6.75))
    palette = sns.color_palette("Set2")
    palette[0] = palette[2]
    
    #ANNOTATION
    pairs = [((n,'Alive'),(n,'Dead')) for n in df_infection.Time.unique()]
    print(pairs)
    fig, ax = plt.subplots(figsize=(12,6.75))
    sns.violinplot(data=df_infection,ax=ax,x='Time',y='weight',hue='Mice',cut=0.5,palette=palette,saturation=1,inner="point",hue_order=("Alive","Dead"),legend=False)
    annotator = Annotator(ax, pairs, data=df_infection, x='Time', y='weight',hue='Mice')
    annotator.configure(test='Kruskal', text_format='star', loc='inside')
    annotator.apply_and_annotate()
    
    
    plt.xticks(rotation=0, fontsize=20)
    plt.yticks(fontsize=20)
    plt.legend([],[], frameon=False)
    plt.xlabel('Time [days]', fontsize=25)
    plt.ylabel('Weight post infection [%]', fontsize=25)
    plt.title("",fontsize=10,style='italic')
    plt.ylim(50, 135)
    if save:
        plt.savefig(path_to_save+infection+".png")

    if show:
        plt.show()

def plot_weight_by_infection(df,infection="L. monocytogenes",group="Mice",plot = False,save = False,path_to_save = "weight_over_time.png",errorbar = ("sd"),max_time=14):
    df_infection = df[df['Infection'] == infection]
    df_infection = df_infection[(df_infection['Time']<=max_time)]
    if group == "group_death":
        hue_order = ['Control_alive','Control_dead','Challenge_alive','Challenge_dead']
    elif group == "Mice":
        hue_order = ['Alive','Dead']
    else:
        hue_order = None
    plt.figure(figsize=(12,6.75))
    palette = sns.color_palette("Set2")
    palette[0] = palette[2]
    sns.lineplot(data=df_infection,x = "Time",y = "weight",hue=group,errorbar=errorbar,err_style="band",
     hue_order=hue_order,legend=False, estimator=np.median,palette=palette) #, estimator=np.median
    x_tick_label = ["Infection" if n ==0 else n for n in range(max_time+1) ]
    plt.xticks(ticks=range(0,max_time+1,1),labels=x_tick_label,rotation=0, fontsize=20)
    plt.yticks(fontsize=20)
    plt.xlabel('Time [days]', fontsize=25)
    plt.ylabel('Weight post infection [%]', fontsize=25)
    plt.title("", fontsize=10,style='italic')
    plt.subplots_adjust(top=0.9,bottom=0.2)
    plt.ylim(50, 135)
    if plot:
        plt.show()
    if save:
        plt.savefig(path_to_save+infection+"dead_vs_alive"+".png")

def remove_unwanted_mice(df):
    dead = df[df['survival_original']==1]
    dead_listeria = dead[(dead['Infection']=='Listeria')& (dead['Time'] == 'T8')]
    dead_listeria_index = dead_listeria.dropna(subset=['weight'],inplace=False).index
    infection = df[df['Infection'].isin(["Listeria","S. pneumoniae"])]
    to_high_index = infection[(infection['weight']>112) & (infection['Time'].isin(["T1","T2","T3"]))].index
    to_low_index = infection[(infection['weight']<60)].index
    #to_remove = to_high_index + to_low_index
    weight_impossible_value = to_high_index.union(to_low_index)
    index_to_remove = dead_listeria_index.union(weight_impossible_value) 
    df = df.drop(index_to_remove)

    return df
if __name__ == "__main__":
    print()
    infection_of_interest = ["C. albicans","Listeria","S. pneumoniae","H1N1","E. coli","Pseudomonas aeruginosas","Staphylococcus aureus"]
    weight_path = "./DF_for_simulation/df_2023_max14days_new_time_calculation_H0.xlsx"
    time_path = "./DF_for_simulation/time_clean.csv"
    simulation_folder = "./simulations/new_simulations/"  
    
    #graphiques to plot
    weight_histogramme = False
    violin = False
    weight_over_time = False
    time_hist = False
    log10_p_value = False
    kaplanMeier = False
    kaplanMeierTotal = False
    simulation_variable_N = False
    simulation_fix_N = True
    animated_image = False

    

    #GRAPHIQUES WEIGHT DISTRIBUTION
    if weight_histogramme:
        df_weight = import_weight_data(weight_path,infection_of_interest)
        weight_distribution(df_weight)
    
    
    #GRAPHIQUES VIOLIN PLOT BY INFECTION
    if violin:
        df = pd.read_excel("./DF_for_simulation/weigth_clean.xlsx",index_col=0)
        df = remove_unwanted_mice(df)

        status_mapping = {0: 'Alive', 1: 'Dead'}
        df["survival_original"] = df['survival_original'].replace(status_mapping)

        infection_mapping = {'Listeria':"L. monocytogenes"}
        df['Infection'] =df['Infection'].replace(infection_mapping)

        df.rename(columns={'survival_original':'Mice'},inplace=True)

        time_rename = list(np.arange(0,15,1))
        time_unique = df['Time'].unique()
        time_mapping = dict(zip(time_unique,time_rename))
        df['Time'] = df['Time'].replace(time_mapping)

        #violin_plot(df,'C. albicans',show=True)
        path_to_save = "./graphiques/NEW_ANALYSIS/violin/"
        for infection in df['Infection'].unique():
            print(infection)
            violin_plot(df,infection=infection,max_time=10,save=True,show=False,path_to_save=path_to_save)

    if weight_over_time:
        df = pd.read_excel("./DF_for_simulation/weigth_clean.xlsx",index_col=0)
        df = remove_unwanted_mice(df)
        df['Mice'] = df['survival_original']
        status_mapping = {0: 'Alive', 1: 'Dead'}
        infection_mapping = {'Listeria':"L. monocytogenes"}
        df["Mice"] = df['Mice'].replace(status_mapping)
        df['Infection'] =df['Infection'].replace(infection_mapping)
        infections = df['Infection'].unique()

        time_rename = list(np.arange(0,15,1))
        time_unique = df['Time'].unique()
        time_mapping = dict(zip(time_unique,time_rename))
        df['Time'] = df['Time'].replace(time_mapping)

        #infection_of_interest = zip(["T10","T10","T10","T10"],['C. albicans','L. monocytogenes','S. pneumoniae',"H1N1"])

        #plot_weight_by_infection(df,plot=True,max_time=10)
        for infection in df['Infection'].unique():
            print(infection)
            plot_weight_by_infection(df=df,infection=infection,save=True,path_to_save="./graphiques/NEW_ANALYSIS/weight_over_time/",max_time=10)

    if time_hist:
        df_time = pd.read_csv(time_path,index_col=0)
        df_time = prepare_df_time(df_time)
        df_time = df_time[df_time['Infection']=="Staphylococcus aureus"]
        df_time = df_time[(df_time['Threshold']>=0.23) & (df_time['Threshold']<=0.28)]
        #displot
        print(df_time)
        sns.displot(data=df_time,x="Time",col="Threshold",row='Control')
        plt.show()
    

    if log10_p_value:
        df_time = pd.read_csv(time_path,index_col=0)
        df_time = prepare_df_time(df_time)
        log_rank = df_time.groupby(['Infection','Threshold']).apply(lambda x: multivariate_logrank_test(event_durations=x.Time,groups=x.Control,event_observed=x.Death).p_value)
        log_rank = log_rank.reset_index()
        log_rank = log_rank.rename(columns={0:'p_value'})
        log_rank['p_value'] = -np.log10( log_rank['p_value'])
        log_rank['Threshold'] = log_rank['Threshold'].astype(str)
        log_rank['Threshold'] = log_rank['Threshold'].replace({"1.0":"original"})
        log_rank = log_rank.sort_values('Threshold',ascending=False)
        plt.figure(figsize=(16, 8))
        sns.lineplot(log_rank,x='Threshold',y="p_value",hue='Infection',sort=False,hue_order=['E. coli',"Listeria","S. pneumoniae","Staphylococcus aureus"])
        plt.xlabel('Weight loss threshold', fontsize=16)
        plt.ylabel('-Log10(p-value)', fontsize=16)
        plt.title("",fontsize=20)
        plt.savefig(f"./graphiques/NEW_ANALYSIS/p-values/log10_p-value.png")
        plt.close()
    
    if kaplanMeier:
        df_time = pd.read_csv(time_path,index_col=0)
        infection_of_interest = 'Pseudomonas aeruginosas'
        df_time = prepare_df_time(df_time,infection_to_chose=[infection_of_interest])
        for thr in df_time['Threshold'].unique():
            plt.figure(figsize=(16, 8))
            ax = plt.subplot(111)
            print(thr)
            kmf = KaplanMeierFitter()
            df_of_interest = df_time[df_time['Threshold']==thr]
            
            
            for name, grouped_df in df_of_interest.groupby('Control'):
                kmf.fit(grouped_df["Time"], grouped_df["Death"], label=name)
                kmf.plot_survival_function(ax=ax)

            plt.savefig(f"./graphiques/NEW_ANALYSIS/KaplanMeier/{infection_of_interest}/thr_{thr}.png")
            plt.close()


    if simulation_variable_N:
        simulation = import_simulations(simulation_folder)
        simulation['Infection'] = simulation['Infection'].replace({"Listeria":"L. monocytogenes"})
        infection_of_interest = ["C. albicans","L. monocytogenes","S. pneumoniae","H1N1"]
        simulation = simulation[simulation['Infection'].isin(infection_of_interest)]
        palette = sns.color_palette("Set2")
        plt.figure(figsize=(12,6.75))
        ax = sns.lineplot(simulation,x="N_mice",y="Power",hue="Infection",errorbar='sd',palette=palette,hue_order=infection_of_interest,linewidth = 3,legend=False)
        plt.xticks(rotation=0, fontsize=20)
        plt.yticks(fontsize=20)
        plt.xlabel('Number of mice', fontsize=25)
        plt.ylabel('Power [%]', fontsize=25)
        plt.title("",fontsize=2)
        plt.tight_layout()
        plt.savefig("./graphiques/NEW_ANALYSIS/Simulation/Variable_N_Mice/report.png")
        #plt.show()

    if simulation_fix_N:
        fix_n_path = "./simulations/new_simulations/fixN/"
        df_fix_N = import_simulations(fix_n_path)
        df_fix_N['Infection'] = df_fix_N['Infection'].replace({"Listeria":"L. monocytogenes"})
        df_fix_N = df_fix_N[df_fix_N['Infection'].isin(["L. monocytogenes","S. pneumoniae"])]
        palette = sns.color_palette("Set2")
        palette[0] = palette[1]
        palette[1] = palette[2]
        print(df_fix_N.columns)
        print(df_fix_N)
        df_fix_N['Power'] = df_fix_N['Power'].astype(float)
        plt.figure(figsize=(12,6.75))
        
        sns.lineplot(data=df_fix_N,x='Threshold',y='Power',hue='Infection',errorbar='sd',palette=palette,linewidth = 3,legend=True)
       
        thr_tick = df_fix_N['Threshold'].unique().tolist()
        thr_tick = thr_tick[::2]
        print("thr tick",thr_tick)
        plt.xticks(ticks=range(0,2*len(thr_tick),2),labels=thr_tick,rotation=45, fontsize=20,ha='right')
        plt.yticks(fontsize=20)
        plt.xlabel('Threshold', fontsize=25)
        plt.ylabel('Power [%]', fontsize=25)
        plt.title("",fontsize=2)
        plt.tight_layout()
        plt.savefig("./graphiques/NEW_ANALYSIS/Simulation/Fix_N_Mice/report.png")
        #plt.show()


    if animated_image:
        input_folder = "./graphiques/NEW_ANALYSIS/KaplanMeier/Staphylococcus aureus/thr_special"  # Replace with the path to your image folder
        output_file = "./graphiques/NEW_ANALYSIS/KaplanMeier/animation/staph_special_thr.gif"
        create_gif(input_folder, output_file)
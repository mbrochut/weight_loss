import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os
import numpy as np
from statannotations.Annotator import Annotator


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
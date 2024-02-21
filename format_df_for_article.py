import pandas as pd
import numpy as np


def keep_only_first_string_element(x):
    """
    Find the first string element in a serie containing mouse weight (float) and delete all elements after the found string
    If there is no string, do nothing
    *Arguments*
    -x: serie object (pandas object)

    *Return*
    x: serie object (pandas object)
    """
    to_keep = ['dead','sacrified','survive']
    index_first_element = x[x.isin(to_keep)].index.values
    if len(index_first_element)>1:
        x.loc[index_first_element[1]:]=np.nan
        return x
    else:
        return x

def text_formating(df):
    """
    remove/replace Unnecessary string in df
    Code cleaning to keep only specific string: "dead, sacrified, survive"
    *Arguments*
    -df: dataframe object

    *Return*
    nothing, but clean the passed df
    """
    df.loc[:,"weight_T_infection":] = df.loc[:,"weight_T_infection":].applymap(lambda s:s.lower() if isinstance(s, str) else s) #remove all uppercase in the dataframe
    df.loc[:,"weight_T_infection":] = df.loc[:,"weight_T_infection":].replace(0,np.nan) # replace value 0 by np.nan 
    df.loc[:,"weight_T_infection":] = df.loc[:,"weight_T_infection":].replace('no data',np.nan) # replace string "no data" by np.nan 
    df.loc[:,"weight_T_infection":] = df.loc[:,"weight_T_infection":].replace(r'\n',' ', regex=True) #replace all "\n" in df by spaces
    replace_dict = {r'(^.*dead.*$)|(^.*dcd.*$)':'dead'}
    df.loc[:,"weight_T_infection":] = df.loc[:,"weight_T_infection":].replace(replace_dict,regex=True) #replcae all instaces containing the world dead or dcd to "dead"
    df.loc[:,"weight_T_infection":] = df.loc[:,"weight_T_infection":].replace({'x':'dead','dec':'dead','d':'dead'}) #replcae all instaces containing the world x, dec, or d to "dead"
    df.loc[:,"weight_T_infection":] = df.loc[:,"weight_T_infection":].replace(regex={r'(^.*no weight.*$)':np.nan}) #replace world "no weight" by nan
    df.loc[:,"weight_T_infection":] = df.loc[:,"weight_T_infection":].replace(regex={r'(^.*sac.*$)': 'sacrified', 'scf': 'sacrified','moribond':'sacrified'}) #replace word containing sac/scf/moribond to sacrifice
    list_key_world =['dead','sacrified','survive',"ane"]
    df.loc[:,"weight_T_infection":] = df.loc[:,"weight_T_infection":].applymap(lambda x : x if x in list_key_world or isinstance(x,float) or isinstance(x,int) else np.nan) #only keep string contained in list_key_world, otherwise nan
    df.loc[:,"weight_T_infection":] = df.loc[:,"weight_T_infection":].apply(lambda x : keep_only_first_string_element(x),axis=1) #if a mouse have multiple string element (dead, sacrifief,survive), only keep the first one


if __name__ == "__main__":

    df = pd.read_excel("data/weight_loss_raw_data_clean.xlsx",index_col=0,parse_dates=True)
    #remove data from anesthesia
    text_formating(df)
    df = df[~df.loc[:,"weight_T_infection":"weight_T14"].apply(lambda x: x.astype(str).str.contains('ane')).any(axis=1)]
    df = df[df['Infection'].isin(['C. albicans','S. pneumoniae','Listeria',"H1N1"])]
    print(len(df))
    df['ID_Experiment'] = df['ID_Experiment'].str.split('_').str[0]+"_"+df['ID_Experiment'].str.split('_').str[1] + '_' + df['ID_Experiment'].str.split('_').str[2].str.zfill(3)

    df['rank'] = df['ID_Experiment'].rank(method='dense')

    # Create a new column with the new IDs
    serie = df['ID_Experiment'].str.split('_').str[0] + '_' + df['rank'].astype(int).astype(str)
    df.insert(loc=2,column='ID_exp',value=serie)
    df['ID_exp'] = df['ID_exp'].str.split('_').str[0]+"_"+ df['ID_exp'].str.split('_').str[1].str.zfill(3)
    df = df.drop(columns='ID_Experiment')
    # Drop the rank column if you don't need it anymore
    df = df.rename(columns={'ID_exp':"ID_Experiment"})
    df = df.drop(columns='rank')
    df.reset_index(inplace=True,drop=True)
    df.to_excel("./data/weight_loss_raw_data.xlsx")
import pandas as pd
import numpy as np
from lifelines.statistics import multivariate_logrank_test
from control_df import hyperlink_list, add_hyperlink

if __name__ == "__main__":
    print("Stats")
    T = [12,12,10,10,14,14,13,9,8,7,8,7]
    G = ["1","1","1","1","1","1","0","0","0","0","0","0"]
    E = [0,0,1,1,0,0,1,1,1,1,1,1]

    log_rank = multivariate_logrank_test(T,G,E,t_0=14)
    print("length",len(T),len(G),len(E))
    print(log_rank.summary)

    df = pd.read_excel("./cleaned_df/good_date_format_clean_2023_with_H0.xlsx",index_col=0)
    mices = pd.read_excel("./autre/Listeria_under_20.xlsx")
    links = hyperlink_list("./cleaned_df/good_date_format_clean_2023_with_H0.xlsx",column=11,sheet_name="Sheet1") #keep in memory the hyperlink of the cells
    print(df)
    print(mices)

    df_mice_under_20 = df[df['Mouse_ID'].isin(mices['Mouse_ID'].to_list())]
    print(df_mice_under_20)
    print(len(links))
    print(df_mice_under_20.index.max())
    links_with_index = [links[n] for n in df_mice_under_20.index]
    print(len(links_with_index))

    name_to_export = './autre/df_for_Didier_Listeria.xlsx'
    df_mice_under_20.to_excel(name_to_export)
    add_hyperlink(name_to_export,name_to_export,links_with_index)

import pandas as pd
print("hello")

df_original = pd.read_excel("./data/df_2023_max14days_new_time_calculation_H0.xlsx",index_col=0)
df_original = df_original[df_original['Infection'].isin(['C. albicans','S. pneumoniae','Listeria',"H1N1"])]
df_original = df_original[~df_original.loc[:,"weight_T_infection":"weight_T14"].apply(lambda x: x.astype(str).str.contains('ane')).any(axis=1)]

df_new = pd.read_excel("./data/df_for_analysis.xlsx",index_col=0)


data_original = df_original.loc[:,"weight_T_infection":].reset_index(drop=True)
data_new = df_new.loc[:,"weight_T_infection":]

are_equal = data_original.equals(data_new)

print(are_equal)  # Output: True
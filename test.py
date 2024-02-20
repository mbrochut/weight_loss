import unittest
import pandas as pd
import numpy
from pandas.testing import assert_frame_equal,assert_series_equal



# Create DataFrames
df1 = pd.read_excel("./testing/DF_for_simulation_TESTING.xlsx")
df2 = pd.read_excel("./testing/DF_for_simulation_TESTING_changeTime.xlsx")
df3 = pd.read_excel("./testing/DF_for_simulation_TESTING_changeTime2.xlsx")
df_timePLusOne = pd.read_excel("./testing/result_for_testing.xlsx",sheet_name="time_plus_one")
experiment_selected = ["ID_exp_0","ID_exp_1","ID_exp_121","ID_exp_143"]
manual_df_for_test = pd.read_excel("./testing/result_for_testing.xlsx",converters={'time_original':float,"time_0.3":float})
manual_df_for_test_plusOne = pd.read_excel("./testing/result_for_testing.xlsx",sheet_name="time_plus_one",converters={'time_original':float,"time_0.3":float})
df_home_made = pd.read_excel("./testing/test_change_time_calculation.xlsx",converters={'time_original':float,"time_0.3":float})
result_home_made = pd.read_excel("./testing/result_for_testing.xlsx",sheet_name="home_made",converters={'time_original':float,"time_0.3":float})
df_original = pd.read_excel("./testing/ORIGINAL_DF.xlsx")
serie2 = result_home_made['time_original']
# Define a test case class
class TestDataframe(unittest.TestCase):

    def test_dataframe_Tinfection_equal_TtimePoint(self):
        # Specify the columns to check

        # Slice the DataFrames to select only the specified columns
        
        df1_subset = df1[df1['ID_Experiment'].isin(experiment_selected)]
        df2_subset = df2[df2['ID_Experiment'].isin(experiment_selected)]
        df1_subset = df1_subset.loc[:,"time_0.3":"survival_original"]
        df2_subset = df2_subset.loc[:,"time_0.3":"survival_original"]

        # Use assert_frame_equal to check if the DataFrames are equal
        assert_frame_equal(df1_subset, df2_subset)

    def test_change_in_time_calculation_original(self):
        df2_exp_40 = df2[df2['ID_Experiment'] == "ID_exp_40"]
        seri1 = df2_exp_40['time_original']
        seri1.reset_index(inplace=True,drop=True)
        seri2 = manual_df_for_test['time_original']
        assert_series_equal(seri1,seri2)

    def test_change_in_time_calculation_03(self):
        df2_exp_40 = df2[df2['ID_Experiment'] == "ID_exp_40"]
        seri1 = df2_exp_40['time_0.3']
        seri1.reset_index(inplace=True,drop=True)
        seri2 = manual_df_for_test['time_0.3']
        assert_series_equal(seri1,seri2)

    def test_change_time_2(self):
        exp40 = df3[df3['ID_Experiment'] == "ID_exp_40"]
        original = exp40['time_original'].reset_index(drop=True)
        thr_03 = exp40['time_0.3'].reset_index(drop=True)
        test_original = manual_df_for_test_plusOne['time_original']
        test_03 = manual_df_for_test_plusOne['time_0.3']
        assert_series_equal(original,test_original)
    
    def test_home_made_original_time(self):
        serie1 = df_home_made['time_original']
        serie2 = result_home_made['time_original']
        assert_series_equal(serie1,serie2)

    def test_home_made_03(self):
        serie1 = df_home_made['time_0.3']
        serie2 = result_home_made['time_0.3']
        assert_series_equal(serie1,serie2)

    def test_home_made_015(self):
        serie1 = df_home_made['time_0.15']
        serie2 = result_home_made['time_0.15']
        assert_series_equal(serie1,serie2)
    
    def test_home_made_005(self):
        serie1 = df_home_made['time_0.05']
        serie2 = result_home_made['time_0.05']
        assert_series_equal(serie1,serie2)

    def test_death_original(self):
        serie1 = df_home_made['time_original']
        serie2 = result_home_made['time_original']
        assert_series_equal(serie1,serie2)
    
    def test_death_03(self):
        serie1 = df_home_made['survival_0.3']
        serie2 = result_home_made['survival_0.3']
        assert_series_equal(serie1,serie2)
    
    def test_death_015(self):
        serie1 = df_home_made['survival_0.15']
        serie2 = result_home_made['survival_0.15']
        assert_series_equal(serie1,serie2)
    
    def test_death_005(self):
        serie1 = df_home_made['survival_0.05']
        serie2 = result_home_made['survival_0.05']
        assert_series_equal(serie1,serie2)

    
    def test_original(self):
        df_to_test = pd.read_excel("./DF_for_simulation/df_2023_max14days_new_time_calculation.xlsx")
        assert_frame_equal(df_original,df_to_test)
    

    def test_time_and_death(self):
        df_to_test = pd.read_excel("./DF_for_simulation/df_2023_max14days_new_time_calculation_H0.xlsx")
        df_to_test = df_to_test.loc[:,"weight_T_infection":]
        df_original_data = df_original.loc[:,"weight_T_infection":]
        print(df_to_test.compare(df_original_data))

    def test_p_values_computing():
        print()
if __name__ == '__main__':
    unittest.main()
    
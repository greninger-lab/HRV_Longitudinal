import pandas as pd

RV_Dataframe = pd.read_csv("/Users/administrator/Desktop/HRV_SupFiles/1_AYM_All_Dataframes/Annotated_StructureWTDpulicates.csv")
RV_Dataframe2 = pd.read_csv("/Users/administrator/Desktop/HRV_SupFiles/1_AYM_All_Dataframes/Epitope_A_B_C.csv")
RV_Dataframe3 = pd.read_csv("/Users/administrator/Desktop/HRV_SupFiles/1_AYM_All_Dataframes/RV_Interpentamer_Sites.csv")
RV_Dataframe4 = pd.read_csv("/Users/administrator/Desktop/HRV_SupFiles/1_AYM_All_Dataframes/Surface_Exposure_Annotated_WTDuplicates.CSV")

RV_Dataframe = RV_Dataframe[RV_Dataframe['iSNV'] == "Yes"]
RV_Dataframe2 = RV_Dataframe2.dropna(subset=['Sample'])
RV_Dataframe4 = RV_Dataframe4[RV_Dataframe4['iSNV'] == "Yes"]

merged_df = pd.merge(RV_Dataframe, RV_Dataframe2, on=['Sample','AF','VP', 'Change','Patient','Sample ID','Day','SurfaceArea','Sequencing Method','Serotype','ChangeNum','1AYM_Position','1AYM_Position_Relative'], how='outer')
merged_df = pd.merge(merged_df, RV_Dataframe3, on=['Sample','AF','VP', 'Change','Patient','Sample ID','Day','Sequencing Method','Serotype','ChangeNum','1AYM_Position','1AYM_Position_Relative'], how='outer')
merged_df = pd.merge(merged_df, RV_Dataframe4, on=['Sample','AF','VP', 'Change','Patient','Sample ID','Day','Sequencing Method','Serotype','ChangeNum','Group','1AYM_Position','1AYM_Position_Relative'], how='outer')

merged_df.drop(columns=['iSNV_y'], inplace=True)
merged_df.drop(columns=['iSNV_x'], inplace=True)
merged_df.drop(columns=['SurfaceArea_x'], inplace=True)
merged_df.drop(columns=['SurfaceArea_y'], inplace=True)
merged_df.drop(columns=['Group_x'], inplace=True)
merged_df.drop(columns=['Group_y'], inplace=True)

#merged_df = merged_df.rename(columns={'SurfaceArea_y': 'SurfaceArea'})
#merged_df['SurfaceArea'] = merged_df['SurfaceArea_y'].copy(deep=True)
merged_df = merged_df.drop(columns=merged_df.columns[22])

merged_df.to_csv('1AYM_MergedDataframes.csv',index=False)
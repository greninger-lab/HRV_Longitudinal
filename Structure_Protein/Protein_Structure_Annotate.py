import pandas as pd

gff = pd.read_csv("/Users/administrator/Desktop/HRV_SupFiles/Structure_Protein/1AYM_full.gff", sep='\t',header = None, skiprows=[i for i in range(0,3)])

RV_Dataframe = pd.read_csv("/Users/administrator/Desktop/HRV_SupFiles/Structure_Protein/1AYM_Allignment_All.csv")

RV_Dataframe['iSNV'] = "Yes"

max_position = RV_Dataframe['1AYM_Position'].max()

positions = list(range(1, 853))

df = pd.DataFrame({'1AYM_Position': positions})

RV_Dataframe = pd.merge(RV_Dataframe, df[~df['1AYM_Position'].isin(RV_Dataframe['1AYM_Position'])], on='1AYM_Position', how='outer')

RV_Dataframe['iSNV'] = RV_Dataframe['iSNV'].fillna(value='No')

RV_Dataframe['Structure'] = ""

for i in range(len(RV_Dataframe)):
    for j in range(len(gff)):
        if(RV_Dataframe.iloc[i].loc['1AYM_Position'] >= gff.iloc[j,3] and RV_Dataframe.iloc[i].loc['1AYM_Position'] <= gff.iloc[j,4]):
            RV_Dataframe.at[i, 'Structure'] = gff.iloc[j,2]

RV_Dataframe = RV_Dataframe.drop_duplicates(subset=['1AYM_Position'])

RV_Dataframe.to_csv('Annotated_Structure.csv', index=False)
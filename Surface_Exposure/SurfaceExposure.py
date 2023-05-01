import pandas as pd
import os
import re

directory = '/Users/administrator/Desktop/HRV_SupFiles/viperdb_info/'

OneAymAll = pd.read_csv("/Users/administrator/Desktop/HRV_SupFiles/Surface_Exposure/1AYM_Allignment_All.csv")

all_dataframes = []

OneAymAll['iSNV'] = "Yes"

OneAymAll['Serotype_AB'] = OneAymAll['Serotype']

OneAymAll['Serotype_AB'] = OneAymAll['Serotype_AB'].str[:1]

for filename in os.listdir(directory):
    if filename.endswith('.csv'):
        filepath = os.path.join(directory, filename)
        df = pd.read_csv(filepath)
        all_dataframes.append(df)

merged_df = pd.concat(all_dataframes, ignore_index=True)

#D = VP4
#68
#B = VP2
#261
#C = VP3
#238
#A = VP1
#285

merged_df = merged_df.rename(columns={'sub': 'VP'})
merged_df = merged_df.rename(columns={'rid': '1AYM_Position_Relative'})

merged_df['VP'] = merged_df['VP'].replace('A', 1)
merged_df['VP'] = merged_df['VP'].replace('B', 3)
merged_df['VP'] = merged_df['VP'].replace('C', 2)
merged_df['VP'] = merged_df['VP'].replace('D', 4)

data = {'VP': ['4'] * 68 + ['2'] * 261 + ['3'] * 238 + ['1'] * 285,
        '1AYM_Position_Relative': list(range(1, 69)) + list(range(1, 262)) + list(range(1, 239)) + list(range(1, 286))}

RestOfAym = pd.DataFrame(data)

RestOfAym['VP'] = RestOfAym['VP'].astype(int)
RestOfAym['1AYM_Position_Relative'] = RestOfAym['1AYM_Position_Relative'].astype(int)

#OneAymAll = OneAymAll.reset_index(inplace=True)
OneAymAll.drop_duplicates('1AYM_Position', inplace = True)

OneAymAll = pd.merge(OneAymAll, RestOfAym, on=['VP', '1AYM_Position_Relative'], how='outer')

OneAymAll = pd.merge(OneAymAll, merged_df, on=['VP', '1AYM_Position_Relative'])

OneAymAll['iSNV'] = OneAymAll['iSNV'].fillna('No')

RestOfAym.to_csv('RestOfAym.CSV')
merged_df.to_csv('MERGED.CSV')
OneAymAll.to_csv('Surface_Exposure_Annotated.CSV',index=False)
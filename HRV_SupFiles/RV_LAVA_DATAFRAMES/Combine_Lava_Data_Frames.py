import pandas as pd
import os
import numpy as np

pd.options.mode.chained_assignment = None

metadata = pd.read_csv('/Users/administrator/Desktop/HRV_SupFiles/Metadata.csv')

directory = "/Users/administrator/Desktop/HRV_SupFiles/RV_LAVA_DATAFRAMES"

ReformatedDataframe2 = None

for filename in os.listdir(directory):

    if not filename.startswith('.') and os.path.isfile(os.path.join(directory, filename)):

        f = os.path.join(directory, filename)

        # checking if it is a file
        if os.path.isfile(f):
            print(f)

    ReformatedDataframe = pd.read_csv(f)

    #ReformatedDataframe = pd.read_csv("/Users/administrator/Desktop/RV_LAVA_DATAFRAMES/visualizationS01.csv")

    ReformatedDataframe = ReformatedDataframe[ReformatedDataframe.iloc[:,3] >= 10]

    ReformatedDataframe = ReformatedDataframe[ReformatedDataframe.iloc[:,9] >= 30]

    ReformatedDataframe = ReformatedDataframe[ReformatedDataframe.iloc[:,8] == "nonsynonymous SNV"]

    VP_list = ['VP1','VP2','VP3','VP4']

    ReformatedDataframe = ReformatedDataframe[ReformatedDataframe['Protein'].isin(VP_list)]

    ReformatedDataframe = ReformatedDataframe.loc[:, ['Sample','AF','Change','Protein']]

    ReformatedDataframe2 = pd.concat([ReformatedDataframe,ReformatedDataframe2])

ReformatedDataframe2['Sample'] = ReformatedDataframe2['Sample'].map(lambda x: x.rstrip('.fastq'))
ReformatedDataframe2['Protein'] = ReformatedDataframe2['Protein'].map(lambda x: x.lstrip('VP'))

ReformatedDataframe2 = ReformatedDataframe2.rename(columns={"Protein":"VP"})

ReformatedDataframe2 = pd.merge(ReformatedDataframe2,metadata,on=["Sample"],how='outer')

ReformatedDataframe2 = ReformatedDataframe2.drop(["Collection date","New Ct","iSNVs","Ratio of iSNVs in capsid/nonstructural genes","Ratio of NS iSNVs in capsid/nonstructural genes","Coverage (%)","NCBI SRA Accession Number","Biosample"],axis=1)

ReformatedDataframe2 = ReformatedDataframe2.dropna(subset=['AF'])

ReformatedDataframe2 = ReformatedDataframe2.drop_duplicates()

ReformatedDataframe2.to_csv('ALL_LAVA_Output.csv', index = False)

SASA = ReformatedDataframe2.loc[ReformatedDataframe2['Serotype'].isin(["A105", "A102", "A82", "A78","A39","A58","A57"])]

SASA['ChangeNum'] = SASA.loc[:,'Change']

SASA['ChangeNum'] = SASA['ChangeNum'].str.extract('(\d+)', expand=False)

AF_maxes = SASA.groupby(['Patient', 'ChangeNum']).AF.transform(max)
SASA = SASA.loc[SASA.AF == AF_maxes]

SASA['Group'] = SASA.loc[:,'Day']

condition = (SASA['Group'] >= 12) & (SASA['Group'] <= 30)

condition2 = (SASA['Group'] >= 68)

# change the values in the selected rows to "12-30"
SASA.loc[condition, 'Group'] = '12-30'

# change the values in the selected rows to "12-30"
SASA.loc[condition2, 'Group'] = '68+'

SASA.to_csv('SASA_Dataframe.csv', index = False)
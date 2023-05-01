import pandas as pd
from Bio import SeqIO

RV_Dataframe = pd.read_csv("/Users/administrator/Desktop/HRV_SupFiles/Interpentamer_Sites/1AYM_Allignment_All.csv")

interpentamer_DF = pd.read_csv("/Users/administrator/Desktop/HRV_SupFiles/Interpentamer_Sites/interpentamer-sites.csv")

for c in range(len(interpentamer_DF)):

    e = 0

    for seq_record in SeqIO.parse("/Users/administrator/Desktop/HRV_SupFiles/Interpentamer_Sites/HRV-B14_1AYM.fasta", "fasta"):
        if (e == 0):
            Epitope_Alignment = seq_record

        else:
           OneAYM_Epitope = seq_record
        e = 1

    #HRVB14
    VP4 = 69 - 1
    VP2 = 262
    VP3 = 236

    position = interpentamer_DF.iloc[c].loc['residue']

    # total position
    if (interpentamer_DF.iloc[c].loc['VP'] == 2):
        position = position + VP4
    elif (interpentamer_DF.iloc[c].loc['VP'] == 3):
        position = position + VP4 + VP2
    elif (interpentamer_DF.iloc[c].loc['VP'] == 1):
        position = position + VP4 + VP2 + VP3

    # character number with gaps
    position_GAP = position

    # Adjusts for position with gaps
    for i in range(position):
        if (Epitope_Alignment[i] == "-"):
            position_GAP = position_GAP + 1

    FourGB3_Position = position_GAP

    for i in range(position_GAP):
        if (OneAYM_Epitope[i] == "-"):
            FourGB3_Position = FourGB3_Position - 1

    if (Epitope_Alignment[position_GAP - 1] == "-" or OneAYM_Epitope[position_GAP - 1] == "-"):
        interpentamer_DF.at[c, '1AYM_Position'] = "-"
    else:
        interpentamer_DF.at[c, '1AYM_Position_Annotated'] = str(Epitope_Alignment[position_GAP - 1]) + str(
            FourGB3_Position) + str(OneAYM_Epitope[position_GAP - 1])
        interpentamer_DF.at[c, '1AYM_Position'] = FourGB3_Position

    NewFourGB3_Position = FourGB3_Position - 1

interpentamer_DF = interpentamer_DF[interpentamer_DF['1AYM_Position'] != "-"]

interpentamer_DF.to_csv("interpentamer_Position_Updated.csv", index=False)

max_position = RV_Dataframe['1AYM_Position'].max()

positions = list(range(1, 853))

df = pd.DataFrame({'1AYM_Position': positions})

RV_Dataframe = pd.merge(RV_Dataframe, df[~df['1AYM_Position'].isin(RV_Dataframe['1AYM_Position'])], on='1AYM_Position', how='outer')

RV_Dataframe['Interpentameter_Site'] = "No"

for i in range(len(RV_Dataframe)):
    for j in range(len(interpentamer_DF)):
        if(RV_Dataframe.iloc[i].loc['1AYM_Position'] == interpentamer_DF.iloc[j].loc['1AYM_Position']):
            RV_Dataframe.at[i, 'Interpentameter_Site'] = "Yes"

RV_Dataframe.to_csv("RV_Interpentamer_Sites.csv", index=False)


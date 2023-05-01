import pandas as pd
from Bio import SeqIO
import re

RV_Dataframe = pd.read_csv("/Users/administrator/Desktop/HRV_SupFiles/Epitopes/ALL_LAVA_Output.csv")
RV_Dataframe_A = pd.read_csv("/Users/administrator/Desktop/HRV_SupFiles/Epitopes/1AYM_A.csv")
Epitope_Dataframe = pd.read_csv("/Users/administrator/Desktop/HRV_SupFiles/Epitopes/Epitope_Position_Updated.csv")

RV_Dataframe = RV_Dataframe.loc[RV_Dataframe['Serotype'].isin(["B97", "B06", "C36", "C28"])]

RV_Dataframe['ChangeNum'] = RV_Dataframe.loc[:,'Change']

RV_Dataframe['ChangeNum'] = RV_Dataframe['ChangeNum'].str.extract('(\d+)', expand=False)

AF_maxes = RV_Dataframe.groupby(['Patient', 'ChangeNum']).AF.transform(max)
RV_Dataframe = RV_Dataframe.loc[RV_Dataframe.AF == AF_maxes]

RV_Dataframe['Group'] = RV_Dataframe.loc[:,'Day']

condition = (RV_Dataframe['Group'] >= 12) & (RV_Dataframe['Group'] <= 30)

condition2 = (RV_Dataframe['Group'] >= 68)

# change the values in the selected rows to "12-30"
RV_Dataframe.loc[condition, 'Group'] = '12-30'

# change the values in the selected rows to "12-30"
RV_Dataframe.loc[condition2, 'Group'] = '68+'

RV_Dataframe.to_csv('RV_Dataframe_Dataframe.csv', index = False)

FileRead = ""

VP4_1AYM = 68
VP2_1AYM = 261
VP3_1AYM = 238

for j in range(len(RV_Dataframe)):

    RV_Dataframe.index = range(len(RV_Dataframe))

    print(RV_Dataframe.iloc[j].loc['Serotype'])

    i = 0

    if(FileRead != RV_Dataframe.iloc[j].loc['Serotype']):

        fastaFilePath = "/Users/administrator/Desktop/HRV_SupFiles/Epitopes/Epitope_Alignment/B_C/" + RV_Dataframe.iloc[j].loc['Serotype'] + "_1AYM.fasta"

        for seq_record in SeqIO.parse(fastaFilePath, "fasta"):
            #print(seq_record.id)
            #print(repr(seq_record.seq))
            #print(len(seq_record))

            if(i == 0):

                SerotypeFasta = seq_record

            else:

                AYM = seq_record

            i = 1

    FileRead = RV_Dataframe.iloc[j].loc['Serotype']

    #B06
    if(SerotypeFasta.id == "B06"):
       VP4 = 69 - 1
       VP2 = 262
       VP3 = 236

    #B97
    if(SerotypeFasta.id == "B97"):
       VP4 = 69 - 1
       VP2 = 261
       VP3 = 236

    #C28
    if(SerotypeFasta.id == "C28"):
       VP4 = 67 - 1
       VP2 = 262
       VP3 = 237

    #C36
    if(SerotypeFasta.id == "C36"):
       VP4 = 67 - 1
       VP2 = 262
       VP3 = 236

    #minus one for A position because of M at begining only if VP4
    #minus one for array indexing

    #position of serotype
    position = int(re.search("\d+", RV_Dataframe.iloc[j].loc['Change'])[0])

    print(position)

    #total position
    if(RV_Dataframe.iloc[j].loc['VP'] == 2):
        position = position + VP4
    elif (RV_Dataframe.iloc[j].loc['VP'] == 3):
        position = position + VP4 + VP2
    elif (RV_Dataframe.iloc[j].loc['VP'] == 1):
        position = position + VP4 + VP2 + VP3

    print(SerotypeFasta.id)
    print(position)

    #character number with gaps
    position_GAP = position

    limitCheck = 0

    #Adjusts for position with gaps
    for i in range(position):
        if(SerotypeFasta.id == 'B06' and position > 847):
            limitCheck = 1
        elif(SerotypeFasta[i] == "-"):
            position_GAP = position_GAP + 1

    Aym_Position = position_GAP

    if(limitCheck == 0):

        for i in range(position_GAP):
            if(AYM[i] == "-"):
                Aym_Position = Aym_Position - 1

        if(AYM[position_GAP - 1] == "-"):
            RV_Dataframe.at[j, '1AYM_Position'] = "-"
        elif(limitCheck == 1):
            RV_Dataframe.at[j, '1AYM_Position'] = "-"
        else:
            RV_Dataframe.at[j,'1AYM_Position'] = str(SerotypeFasta[position_GAP - 1]) + str(Aym_Position) + str(AYM[position_GAP - 1])

        Aym_Position_Relative = Aym_Position

        #total position
        if(RV_Dataframe.iloc[j].loc['VP'] == 2):
            Aym_Position_Relative = Aym_Position_Relative - VP4_1AYM
        elif (RV_Dataframe.iloc[j].loc['VP'] == 3):
            Aym_Position_Relative = Aym_Position_Relative - VP4_1AYM - VP2_1AYM
        elif (RV_Dataframe.iloc[j].loc['VP'] == 1):
            Aym_Position_Relative = Aym_Position_Relative - VP4_1AYM - VP2_1AYM - VP3_1AYM

        RV_Dataframe.at[j,'1AYM_Position_Relative'] = Aym_Position_Relative

        #for l in range(len(RV_Dataframe_Pymol)):
        #    if(RV_Dataframe.iloc[j].loc['VP'] == RV_Dataframe_Pymol.at[l,"V Region"] and RV_Dataframe.at[j,'1AYM_Position_Relative'] == RV_Dataframe_Pymol.at[l,"Residue"]):
        #        RV_Dataframe.at[j,'SurfaceArea'] = RV_Dataframe_Pymol.at[l,"RV_Dataframe"]

RV_Dataframe = RV_Dataframe.loc[RV_Dataframe['1AYM_Position'] != '-']
RV_Dataframe = RV_Dataframe.dropna(subset=['1AYM_Position'])

RV_Dataframe.to_csv('Epitope_B_C.csv', index=False)

RV_Dataframe = pd.concat([RV_Dataframe, RV_Dataframe_A], axis=0, ignore_index=True)

RV_Dataframe['1AYM_Position'] = RV_Dataframe['1AYM_Position'].str.replace(r'\D', '')

#RV_Dataframe = pd.concat([RV_Dataframe, Epitope_Dataframe], axis=0, ignore_index=True)

RV_Dataframe['1AYM_Position'] = RV_Dataframe['1AYM_Position'].astype(int)

RV_Dataframe.to_csv('1AYM_Allignment_All.csv', index=False)

RV_Dataframe = pd.merge(RV_Dataframe, Epitope_Dataframe, on='1AYM_Position', how='outer')

RV_Dataframe.to_csv('Epitope_A_B_C.csv', index=False)


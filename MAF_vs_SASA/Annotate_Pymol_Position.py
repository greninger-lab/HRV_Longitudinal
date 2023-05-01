from Bio import SeqIO
import pandas as pd
import re
import os

SASA = pd.read_csv("/Users/administrator/Desktop/HRV_SupFiles/SASA_Dataframe.csv")
SASA_Pymol = pd.read_csv("/Users/administrator/Desktop/HRV_SupFiles/SASA_Extracted_Pymol.csv")

#total amino position
SASA['1AYM_Position'] = ""
#position relative to vp region
SASA['1AYM_Position_Relative'] = ""
SASA['SurfaceArea'] = ""
SASA['Group'] = SASA.loc[:,'Day']

condition = (SASA['Group'] >= 12) & (SASA['Group'] <= 30)
condition2 = (SASA['Group'] >= 68)

# change the values in the selected rows to "12-30"
SASA.loc[condition, 'Group'] = '12-30'

# change the values in the selected rows to "12-30"
SASA.loc[condition2, 'Group'] = '68+'

FileRead = ""

VP4_1AYM = 68
VP2_1AYM = 261
VP3_1AYM = 238

for j in range(len(SASA)):

    print(SASA.iloc[j].loc['Serotype'])

    i = 0

    if(FileRead != SASA.iloc[j].loc['Serotype']):

        fastaFilePath = "/Users/administrator/Desktop/HRV_SupFiles/Alignments/RV_A/" + SASA.iloc[j].loc['Serotype'] + "_1AYM.fasta"

        for seq_record in SeqIO.parse(fastaFilePath, "fasta"):
            #print(seq_record.id)
            #print(repr(seq_record.seq))
            #print(len(seq_record))

            print("here")

            if(i == 0):
                SerotypeFasta = seq_record

            else:

                AYM = seq_record

            i = 1

    FileRead = SASA.iloc[j].loc['Serotype']

    #A102
    if(SerotypeFasta.id == "MZ458533"):
       VP4 = 69 - 1
       VP2 = 262
       VP3 = 238

    #A105
    if(SerotypeFasta.id == "MZ542285"):
       VP4 = 69 - 1
       VP2 = 261
       VP3 = 238

    #A39
    if(SerotypeFasta.id == "MZ667419"):
       VP4 = 69 - 1
       VP2 = 265
       VP3 = 238

    #A57
    if(SerotypeFasta.id == "MZ667415"):
       VP4 = 69 - 1
       VP2 = 262
       VP3 = 238

    #A58
    if(SerotypeFasta.id == "MZ667420"):
       VP4 = 69 - 1
       VP2 = 265
       VP3 = 238

    #A78
    if(SerotypeFasta.id == "MZ667417"):
       VP4 = 69 - 1
       VP2 = 264
       VP3 = 238

    #A82
    if(SerotypeFasta.id == "MZ667414"):
       VP4 = 69 - 1
       VP2 = 260
       VP3 = 237

    #minus one for A position because of M at begining only if VP4
    #minus one for array indexing

    #position of serotype
    position = int(re.search("\d+", SASA.iloc[j].loc['Change'])[0])

    #total position
    if(SASA.iloc[j].loc['VP'] == 2):
        position = position + VP4
    elif (SASA.iloc[j].loc['VP'] == 3):
        position = position + VP4 + VP2
    elif (SASA.iloc[j].loc['VP'] == 1):
        position = position + VP4 + VP2 + VP3

    #character number with gaps
    position_GAP = position

    #Adjusts for position with gaps
    for i in range(position):
        if(SerotypeFasta[i] == "-"):
            position_GAP = position_GAP + 1

    Aym_Position = position_GAP

    for i in range(position_GAP):
        if(AYM[i] == "-"):
            Aym_Position = Aym_Position - 1

    if(AYM[position_GAP - 1] == "-"):
        SASA.at[j, '1AYM_Position'] = "-"
    else:
        SASA.at[j,'1AYM_Position'] = str(SerotypeFasta[position_GAP - 1]) + str(Aym_Position) + str(AYM[position_GAP - 1])

    Aym_Position_Relative = Aym_Position

    #total position
    if(SASA.iloc[j].loc['VP'] == 2):
        Aym_Position_Relative = Aym_Position_Relative - VP4_1AYM
    elif (SASA.iloc[j].loc['VP'] == 3):
        Aym_Position_Relative = Aym_Position_Relative - VP4_1AYM - VP2_1AYM
    elif (SASA.iloc[j].loc['VP'] == 1):
        Aym_Position_Relative = Aym_Position_Relative - VP4_1AYM - VP2_1AYM - VP3_1AYM

    SASA.at[j,'1AYM_Position_Relative'] = Aym_Position_Relative

    for l in range(len(SASA_Pymol)):
        if(SASA.iloc[j].loc['VP'] == SASA_Pymol.at[l,"V Region"] and SASA.at[j,'1AYM_Position_Relative'] == SASA_Pymol.at[l,"Residue"]):
            SASA.at[j,'SurfaceArea'] = SASA_Pymol.at[l,"SASA"]

SASA = SASA[SASA['1AYM_Position'] != "-"]

OneAYM_A = SASA

SASA = SASA[SASA['SurfaceArea'] != ""]

SASA['SurfaceArea'] = SASA['SurfaceArea'].map(lambda x: x.rstrip('%'))

OneAYM_A.to_csv('1AYM_A.csv', index = False)
SASA.to_csv('SASA_Annotated.csv', index = False)

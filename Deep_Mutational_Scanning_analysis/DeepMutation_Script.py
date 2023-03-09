import pandas as pd
from Bio import SeqIO
import re
import blosum as bl

matrix = bl.BLOSUM(62, default = 4)

pd.options.mode.chained_assignment = None

RV_A = pd.read_csv("/Users/administrator/Desktop/HRV_SupFiles/Deep_Mutation/RV_A/1AYM_A.csv")

AVG_MAX_RP = pd.read_csv("/Users/administrator/Desktop/HRV_SupFiles/Deep_Mutation/AVG_MAX_RP.csv")

AVG_MAX_RP['BlosumScore'] = ""

All_LAVA_OUTPUT = pd.read_csv("/Users/administrator/Desktop/HRV_SupFiles/Deep_Mutation/ALL_LAVA_Output.csv")

Epitopes_DF = pd.read_csv("/Users/administrator/Desktop/HRV_SupFiles/Deep_Mutation/RV_A/epitopes-list.csv")

All_LAVA_OUTPUT['ChangeNum'] = All_LAVA_OUTPUT.loc[:,'Change']

All_LAVA_OUTPUT['ChangeNum'] = All_LAVA_OUTPUT['ChangeNum'].str.extract('(\d+)', expand=False)

AF_maxes = All_LAVA_OUTPUT.groupby(['Patient', 'ChangeNum']).AF.transform(max)
All_LAVA_OUTPUT = All_LAVA_OUTPUT.loc[All_LAVA_OUTPUT.AF == AF_maxes]

All_LAVA_OUTPUT['4GB3_Position'] = ""

All_LAVA_OUTPUT['mut_virus_AVR_MFE_log2'] = ""

All_LAVA_OUTPUT['MAX_mut_virus_AVR_MFE_log2'] = ""

All_LAVA_OUTPUT['BlosumScore'] = ""

All_LAVA_OUTPUT['iSNP'] = "Yes"

RV_A = RV_A[['Sample','AF','VP','Day','Sample ID', 'Serotype', '1AYM_Position','1AYM_Position_Relative']]

RV_A['4GB3_Position'] = ""

RV_A['mut_virus_AVR_MFE_log2'] = ""

RV_A['MAX_mut_virus_AVR_MFE_log2'] = ""

RV_A['iSNP'] = "Yes"

RV_A['BlosumScore'] = ""

AVG_MAX_RP['iSNP'] = "No"

AVG_MAX_RP['Sample'] = ""

AVG_MAX_RP['VP'] = AVG_MAX_RP.loc[:,'4GB3_Position']

condition = (AVG_MAX_RP['VP'] >= 1) & (AVG_MAX_RP['VP'] <= 68)
condition2 = (AVG_MAX_RP['VP'] >= 69) & (AVG_MAX_RP['VP'] <= 331)
condition3 = (AVG_MAX_RP['VP'] >= 332) & (AVG_MAX_RP['VP'] <= 569)
condition4 = (AVG_MAX_RP['VP'] >= 570) & (AVG_MAX_RP['VP'] <= 850)

AVG_MAX_RP.loc[condition, 'VP'] = '4'
AVG_MAX_RP.loc[condition2, 'VP'] = '2'
AVG_MAX_RP.loc[condition3, 'VP'] = '3'
AVG_MAX_RP.loc[condition4, 'VP'] = '1'

AVG_MAX_RP = AVG_MAX_RP[AVG_MAX_RP.loc[:,'4GB3_Position'] != 851]

AVG_MAX_RP_A = AVG_MAX_RP
AVG_MAX_RP_B = AVG_MAX_RP
AVG_MAX_RP_C = AVG_MAX_RP

FileRead = ''

for c in range(len(Epitopes_DF)):

    e = 0

    if (FileRead != Epitopes_DF.iloc[c].loc['REF']):

        fastaFilePath = "/Users/administrator/Desktop/HRV_SupFiles/Deep_Mutation/RV_A/Epitope_Alignment/" + Epitopes_DF.iloc[c].loc['REF'] + "_1AYM.fasta"

        for seq_record in SeqIO.parse(fastaFilePath, "fasta"):
            if (e == 0):
                Epitope_Alignment = seq_record

            else:
               OneAYM_Epitope = seq_record
            e = 1

    #HRVA14
    if (Epitope_Alignment.id == "A14"):
    #if (Epitopes_DF.iloc[c].loc['REF'] == "HRV-14"):
        VP4 = 69 - 1
        VP2 = 254
        VP3 = 242
    # HRVA2
    if (Epitope_Alignment.id == "A2"):
    #if (Epitopes_DF.iloc[c].loc['REF'] == "HRV-14"):
        VP4 = 69 - 1
        VP2 = 255
        VP3 = 240

    position = Epitopes_DF.iloc[c].loc['residue']

    # total position
    if (Epitopes_DF.iloc[c].loc['VP'] == 2):
        position = position + VP4
    elif (Epitopes_DF.iloc[c].loc['VP'] == 3):
        position = position + VP4 + VP2
    elif (Epitopes_DF.iloc[c].loc['VP'] == 1):
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

    if (Epitope_Alignment[position_GAP - 1] == "-"):
        Epitopes_DF.at[c, '1AYM_Position'] = "-"
    else:
        Epitopes_DF.at[c, '1AYM_Position_Annotated'] = str(Epitope_Alignment[position_GAP - 1]) + str(
            FourGB3_Position) + str(OneAYM_Epitope[position_GAP - 1])
        Epitopes_DF.at[c, '1AYM_Position'] = FourGB3_Position

    NewFourGB3_Position = FourGB3_Position - 1

Epitopes_DF = Epitopes_DF[Epitopes_DF['1AYM_Position'] != "-"]

Epitopes_DF.to_csv("Epitope.csv")

i = 0

for seq_record in SeqIO.parse("/Users/administrator/Desktop/HRV_SupFiles/Deep_Mutation/RV_A/1AYM_4GB3.fasta", "fasta"):
    if(i == 0):
        OneAym = seq_record

    else:
        FourGB3_A = seq_record
    i = 1

for j in range(len(RV_A)):
    #character number with gaps
    position_A = int(re.search("\d+", RV_A.iloc[j].loc['1AYM_Position'])[0])

    position_GAP_A = position_A

    #Adjusts for position with gaps
    for a in range(position_A):
        if(OneAym[a] == "-"):
            position_GAP_A = position_GAP_A + 1

    while (OneAym[position_GAP_A-1] == "-"):
        position_GAP_A = position_GAP_A + 1

    FourGB3_A_position = position_GAP_A

    for b in range(position_GAP_A):
        if(FourGB3_A[b] == "-"):
            FourGB3_A_position = FourGB3_A_position - 1

    if(FourGB3_A[position_GAP_A - 1] == "-"):
        RV_A.at[j,'4GB3_Position'] = "-"
    else:
        RV_A.at[j,'4GB3_Position_Annotated'] = str(OneAym[position_GAP_A - 1]) + str(FourGB3_A_position) + str(FourGB3_A[position_GAP_A - 1])
        RV_A.at[j, '4GB3_Position'] = FourGB3_A_position

    RV_A.at[j,'mut_virus_AVR_MFE_log2'] = AVG_MAX_RP.at[FourGB3_A_position - 1,'mut_virus_AVR_MFE_log2']
    RV_A.at[j,'MAX_mut_virus_AVR_MFE_log2'] = AVG_MAX_RP.at[FourGB3_A_position - 1,'MAX_mut_virus_AVR_MFE_log2']

    AVG_MAX_RP_A = AVG_MAX_RP_A[AVG_MAX_RP_A['4GB3_Position'] != FourGB3_A_position]

    # BlosumScore
    count = 0

    for q in range(11):
        position = -6 + q + position_GAP_A
        amin = FourGB3_A[position] + OneAym[position]
        count = count + matrix[amin]

    RV_A.at[j,'BlosumScore'] = count/11

RV_A = RV_A[RV_A['4GB3_Position'] != "-"]

print(AVG_MAX_RP_A)

for k in range(len(OneAym)):

    if(OneAym[k] == '-'):
        GapCount = 0

        # Adjusts for position with gaps
        for a in range(k):
            if (FourGB3_A[a] == "-"):
                GapCount = GapCount + 1

        AVG_MAX_RP_A = AVG_MAX_RP_A[AVG_MAX_RP_A['4GB3_Position'] != k - GapCount + 1]

#print(AVG_MAX_RP_A)
AVG_MAX_RP_A.index = range(len(AVG_MAX_RP_A))

for a in range(len(AVG_MAX_RP_A)):

    count = 0
    numSlide = 11

    for q in range(11):

        # Adjusts for position with gaps
        GapCount = 0

        if(AVG_MAX_RP_A.loc[a,"4GB3_Position"] < 850):
            #print(AVG_MAX_RP_A.loc[a,"4GB3_Position"])
            for b in range(AVG_MAX_RP_A.loc[a,"4GB3_Position"]):
                if (FourGB3_A[b] == "-"):
                    GapCount = GapCount + 1

        position = -6 + q + AVG_MAX_RP_A.loc[a,"4GB3_Position"] + GapCount

        if(position < 858 and position >= 0):
            amin = FourGB3_A[position] + OneAym[position]
            count = count + matrix[amin]
        else:
            numSlide = numSlide - 1

    sumCount = (count / numSlide)

    AVG_MAX_RP_A.at[a,'BlosumScore'] = sumCount

RV_B = All_LAVA_OUTPUT.loc[All_LAVA_OUTPUT['Serotype'].isin(["B97", "B06"])]

RV_C = All_LAVA_OUTPUT.loc[All_LAVA_OUTPUT['Serotype'].isin(["C36", "C28"])]

j = 0

#No longer RV A
for l in range(2):
    if(l == 0):
        RV_Dataframe = RV_B
        RV_Dataframe.index = range(len(RV_Dataframe))
    elif (l == 1):
        RV_Dataframe = RV_C
        RV_Dataframe.index = range(len(RV_Dataframe))

    j = 0

    FileRead = ""

    for j in range(len(RV_Dataframe)):

        i = 0

        if(FileRead != RV_Dataframe.iloc[j].loc['Serotype']):

            fastaFilePath = "/Users/administrator/Desktop/HRV_SupFiles/Alignments/" + RV_Dataframe.iloc[j].loc['Serotype'] + "_4GB3.fasta"

            for seq_record in SeqIO.parse(fastaFilePath, "fasta"):

                if(i == 0):
                    SerotypeFasta = seq_record
                    if (l == 0):
                        SerotypeFastaB = SerotypeFasta
                    if (l == 1):
                        SerotypeFastaC = SerotypeFasta

                else:

                    FourGB3 = seq_record

                    if (l == 0):
                        FourGB3_B = FourGB3
                    if (l == 1):
                        FourGB3_C = FourGB3

                i = 1

        FileRead = RV_Dataframe.iloc[j].loc['Serotype']

        #B06
        if(SerotypeFasta.id == "MZ667416"):
           VP4 = 69 - 1
           VP2 = 262
           VP3 = 236

        #B97
        if(SerotypeFasta.id == "MZ667418"):
           VP4 = 69 - 1
           VP2 = 261
           VP3 = 236

        #C28
        if(SerotypeFasta.id == "MZ447875"):
           VP4 = 67 - 1
           VP2 = 262
           VP3 = 237

        #C36
        if(SerotypeFasta.id == "MZ438010"):
           VP4 = 67 - 1
           VP2 = 262
           VP3 = 236

        #minus one for A position because of M at begining only if VP4
        #minus one for array indexing

        #position of serotype
        position = int(re.search("\d+", RV_Dataframe.iloc[j].loc['Change'])[0])

        #total position
        if(RV_Dataframe.iloc[j].loc['VP'] == 2):
            position = position + VP4
        elif (RV_Dataframe.iloc[j].loc['VP'] == 3):
            position = position + VP4 + VP2
        elif (RV_Dataframe.iloc[j].loc['VP'] == 1):
            position = position + VP4 + VP2 + VP3

        #character number with gaps
        position_GAP = position

        #Adjusts for position with gaps
        for i in range(position):
            if(SerotypeFasta[i] == "-"):
                position_GAP = position_GAP + 1

        FourGB3_Position = position_GAP

        for i in range(position_GAP):
            if(FourGB3[i] == "-"):
                FourGB3_Position = FourGB3_Position - 1

        if(FourGB3[position_GAP - 1] == "-"):
            RV_Dataframe.at[j,'4GB3_Position'] = "-"
        else:
            RV_Dataframe.at[j,'4GB3_Position_Annotated'] = str(SerotypeFasta[position_GAP - 1]) + str(FourGB3_Position) + str(FourGB3[position_GAP - 1])
            RV_Dataframe.at[j,'4GB3_Position'] = FourGB3_Position

            # BlosumScore
            count = 0

            for q in range(11):
                position = -6 + q + position_GAP

                #if((FourGB3[position]) == (SerotypeFasta[position])):
                #    count = count + 1
                amin = FourGB3[position] + SerotypeFasta[position]
                count = count + matrix[amin]

            RV_Dataframe.at[j, 'BlosumScore'] = count/11

        NewFourGB3_Position = FourGB3_Position - 1

        RV_Dataframe.at[j, 'mut_virus_AVR_MFE_log2'] = AVG_MAX_RP.at[NewFourGB3_Position, 'mut_virus_AVR_MFE_log2']
        RV_Dataframe.at[j, 'MAX_mut_virus_AVR_MFE_log2'] = AVG_MAX_RP.at[NewFourGB3_Position, 'MAX_mut_virus_AVR_MFE_log2']

        if(l == 0):
            AVG_MAX_RP_B = AVG_MAX_RP_B[AVG_MAX_RP_B['4GB3_Position'] != FourGB3_Position]

        elif (l == 1):
            AVG_MAX_RP_C = AVG_MAX_RP_C[AVG_MAX_RP_C['4GB3_Position'] != FourGB3_Position]

        if (l == 0):
            RV_Dataframe_2 = AVG_MAX_RP_B
            RV_Dataframe_2.index = range(len(AVG_MAX_RP_B))
        elif (l == 1):
            RV_Dataframe_2 = AVG_MAX_RP_C
            RV_Dataframe_2.index = range(len(AVG_MAX_RP_C))

        for k in range(len(SerotypeFasta)):

            if (SerotypeFasta[k] == '-'):
                GapCount = 0

                # Adjusts for position with gaps
                for a in range(k):
                    if (FourGB3[a] == "-"):
                        GapCount = GapCount + 1

                RV_Dataframe_2 = RV_Dataframe_2[RV_Dataframe_2['4GB3_Position'] != k - GapCount + 1]

        if (l == 0):
            AVG_MAX_RP_B = RV_Dataframe_2
            AVG_MAX_RP_B.index = range(len(AVG_MAX_RP_B))
        elif (l == 1):
            AVG_MAX_RP_C = RV_Dataframe_2
            AVG_MAX_RP_C.index = range(len(AVG_MAX_RP_C))

    if (l == 0):
        RV_B = RV_Dataframe
        RV_B.index = range(len(RV_B))
    elif (l == 1):
        RV_C = RV_Dataframe
        RV_C.index = range(len(RV_C))

for d in range(2):
    if(d == 0):
        AVG_MAX_RP_Dataframe = AVG_MAX_RP_B
        FourGB3_Dataframe = FourGB3_B
        SerotypeFasta_Dataframe = SerotypeFastaB
        position_limit = 868
    if(d == 1):
        AVG_MAX_RP_Dataframe = AVG_MAX_RP_C
        FourGB3_Dataframe = FourGB3_C
        SerotypeFasta_Dataframe = SerotypeFastaC
        position_limit = 862


    for a in range(len(AVG_MAX_RP_Dataframe)):

        count = 0
        numSlide = 11

        for q in range(11):

            # Adjusts for position with gaps
            GapCount = 0

            if(AVG_MAX_RP_Dataframe.loc[a,"4GB3_Position"] < 850):
                #print(AVG_MAX_RP_A.loc[a,"4GB3_Position"])
                for b in range(AVG_MAX_RP_Dataframe.loc[a,"4GB3_Position"]):
                    if (FourGB3_Dataframe[b] == "-"):
                        GapCount = GapCount + 1

            position = -6 + q + AVG_MAX_RP_Dataframe.loc[a,"4GB3_Position"] + GapCount

            if(position < position_limit and position >= 0):
                amin = FourGB3_Dataframe[position] + SerotypeFasta_Dataframe[position]
                count = count + matrix[amin]
            else:
                numSlide = numSlide - 1

        sumCount = (count / numSlide)

        #AVG_MAX_RP_Dataframe.at[a,'BlosumScore'] = sumCount
        if (d == 0):
             AVG_MAX_RP_B.at[a,'BlosumScore'] = sumCount
        if (d == 1):
            AVG_MAX_RP_C.at[a,'BlosumScore'] = sumCount

RV_B = RV_B[RV_B['4GB3_Position'] != "-"]
RV_C = RV_C[RV_C['4GB3_Position'] != "-"]

AVG_MAX_RP_A.to_csv('AVG_MAX_RP_A.csv', index=False)
AVG_MAX_RP_B.to_csv('AVG_MAX_RP_B.csv', index=False)
AVG_MAX_RP_C.to_csv('AVG_MAX_RP_C.csv', index=False)

RV_A = pd.concat([RV_A, AVG_MAX_RP_A], axis=0, ignore_index=True)
RV_B = pd.concat([RV_B, AVG_MAX_RP_B], axis=0, ignore_index=True)
RV_C = pd.concat([RV_C, AVG_MAX_RP_C], axis=0, ignore_index=True)

RV_A.to_csv('RV_A.csv', index=False)
RV_B.to_csv('RV_B.csv', index=False)
RV_C.to_csv('RV_C.csv', index=False)

RV_Combo = pd.concat([RV_A,RV_B,RV_C], axis=0, ignore_index=True)

RV_Combo_No = RV_Combo[RV_Combo['iSNP'] != 'Yes']
RV_Combo = RV_Combo[RV_Combo['iSNP'] != 'No']

RV_Combo_No['dup_number'] = RV_Combo_No.groupby(['4GB3_Position']).cumcount()+1

RV_Combo_No = RV_Combo_No[RV_Combo_No['dup_number'] == 3]

RV_Combo_No = RV_Combo_No.drop('dup_number', axis=1)

RV_Combo = pd.concat([RV_Combo,RV_Combo_No], axis=0, ignore_index=True)

RV_Combo.to_csv('RV_Combo.csv', index=False)
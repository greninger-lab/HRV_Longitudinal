import pandas as pd

#Df = pd.read_csv("/Users/administrator/Desktop/HRV_SupFiles/mNGS_PCR/Raw/S1_Comparision_mNGS.csv")
#Df = pd.read_csv("/Users/administrator/Desktop/HRV_SupFiles/mNGS_PCR/Raw/S1_Comparision.csv")
#Df = pd.read_csv("/Users/administrator/Desktop/HRV_SupFiles/mNGS_PCR/Raw/S8_Comparison.csv")
Df = pd.read_csv("/Users/administrator/Desktop/HRV_SupFiles/mNGS_PCR/Raw/S10_Compairison.csv")

Df = Df[Df.iloc[:, 3] >= 10]

Df = Df[Df.iloc[:, 9] >= 30]

Df2 = pd.DataFrame(columns=['Sample1', 'AF_1', 'NUC_1', 'Sample2', 'AF_2', 'NUC_2'])

for i in range(len(Df)):
    for j in range(len(Df)):

        if ((Df.iloc[i,0] != Df.iloc[j,0]) and (Df.iloc[i,6] == Df.iloc[j,6])):
            print(str(Df.iloc[i,0]) + ": " + str(Df.iloc[i,3]) + " " + str(Df.iloc[i,6]) + str(Df.iloc[j,0]) + ": " + str(Df.iloc[j,3]) + " " + str(Df.iloc[j,6]))
            Df2.loc[i] = [Df.iloc[i,0], Df.iloc[i,3],Df.iloc[i,6], Df.iloc[j,0], Df.iloc[j,3], Df.iloc[j,6]]

#Df2.to_csv("BALB_BALB_2.csv", index=False)
#Df2.to_csv("NW1PCR_NW1mNGS.csv", index=False)
#Df2.to_csv("NW8PCR_NW8mNGS.csv", index=False)
Df2.to_csv("BAL10PCR_BAL10mNGS.csv", index=False)

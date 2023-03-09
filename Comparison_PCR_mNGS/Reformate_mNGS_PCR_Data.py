import pandas as pd

#Df = pd.read_csv("/Users/administrator/Desktop/HRV_SupFiles/mNGS_PCR/Raw/S11_Comparision_mNGS.csv")
#Df = pd.read_csv("/Users/administrator/Desktop/HRV_SupFiles/mNGS_PCR/Raw/S11_Comparision.csv")
#Df = pd.read_csv("/Users/administrator/Desktop/HRV_SupFiles/mNGS_PCR/Raw/S04_Comparison.csv")
Df = pd.read_csv("/Users/administrator/Desktop/HRV_SupFiles/mNGS_PCR/Raw/S03_Compairison.csv")

Df = Df[Df.iloc[:, 3] >= 10]

Df = Df[Df.iloc[:, 9] >= 30]

Df2 = pd.DataFrame(columns=['Sample1', 'AF_1', 'NUC_1', 'Sample2', 'AF_2', 'NUC_2'])

for i in range(len(Df)):
    for j in range(len(Df)):

        if ((Df.iloc[i,0] != Df.iloc[j,0]) and (Df.iloc[i,6] == Df.iloc[j,6])):
            print(str(Df.iloc[i,0]) + ": " + str(Df.iloc[i,3]) + " " + str(Df.iloc[i,6]) + str(Df.iloc[j,0]) + ": " + str(Df.iloc[j,3]) + " " + str(Df.iloc[j,6]))
            Df2.loc[i] = [Df.iloc[i,0], Df.iloc[i,3],Df.iloc[i,6], Df.iloc[j,0], Df.iloc[j,3], Df.iloc[j,6]]

#Df2.to_csv("BAL11B_BAL11B_2.csv", index=False)
#Df2.to_csv("NW11PCR_NW11mNGS.csv", index=False)
#Df2.to_csv("NW3PCR_NW3mNGS.csv", index=False)
Df2.to_csv("BAL4PCR_BAL4mNGS.csv", index=False)

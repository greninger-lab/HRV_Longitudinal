import pandas as pd

def mix_colors(number):
    # Convert the number to a value between 0 and 255
    value = int(number * 255)

    # Calculate the RGB values for purple and yellow
    purple = (128, 0, 128)
    red = (255, 0, 0)
    yellow = (255, 255, 0)

    # Mix the two colors using the value as the weight
    # For example, if value is 128, the mixed color will be halfway between purple and yellow

    mixed = tuple(int(yellow[i] * (1 - value / 255.0) + red[i] * (value / 255.0)) for i in range(3))

    hex_string = '0x%02x%02x%02x' % mixed

    return hex_string

#Peptide3d = pd.read_csv("/Users/administrator/Desktop/Old/Peptide3D_HRV/PeptideUpload.CSV")
#Peptide3d = pd.read_csv("/Users/administrator/Desktop/Old/Peptide3D_HRV/PeptideUpload.CSV")
Peptide3d = pd.read_csv("/Users/administrator/Desktop/HRV_SupFiles/MAF_vs_SASA/SASA_Annotated.csv")
#print(Peptide3d)

#Day
#Peptide3d = Peptide3d[Peptide3d.iloc[:,6] == 0]
#print(Peptide3d)
#Peptide3d = Peptide3d[Peptide3d.iloc[:,6] >= 12]
#Peptide3d = Peptide3d[Peptide3d.iloc[:,6] <= 30]

Peptide3d = Peptide3d[Peptide3d.iloc[:,6] >= 68]

Peptide3d['AF'] = Peptide3d['AF'].div(100)

Peptide3d['AF'] = Peptide3d['AF'].apply(mix_colors)

Peptide3d.index = range(len(Peptide3d))

for i in range(len(Peptide3d)):
    #Peptide3d2.iloc[i, 1] = "  ///" + str(Peptide3d.loc[i, 'VP']) + "/" + str(Peptide3d.loc[i, '1AYM_Position_Relative']) + "/"
    #Peptide3d2.iloc[i, 0] = "color " + str(Peptide3d.loc[i, 'AF'])
    Peptide3d2 = "color " + str(Peptide3d.loc[i, 'AF']) + "," + "  ///" + str(Peptide3d.loc[i, 'VP']) + "/" + str(Peptide3d.loc[i, '1AYM_Position_Relative']) + "/"

    print(Peptide3d2)

#del Peptide3d['AF']
#del Peptide3d['Day']

#Peptide3d.to_csv('PeptideHex', index = False, encoding='utf-8',header=False)

#print(Peptide3d)

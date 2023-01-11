import pandas as pd

def mix_colors(number):
    # Convert the number to a value between 0 and 255
    value = int(number * 255)

    # Calculate the RGB values for purple and yellow
    purple = (128, 0, 128)
    red = (255, 0, 0)
    yellow = (255, 255, 0)

    #black = (0, 0, 0)
    #white = (255, 255, 255)

    # Mix the two colors using the value as the weight
    # For example, if value is 128, the mixed color will be halfway between purple and yellow
    #mixed = tuple(int(purple[i] * (1 - value / 255.0) + yellow[i] * (value / 255.0)) for i in range(3))

    mixed = tuple(int(yellow[i] * (1 - value / 255.0) + red[i] * (value / 255.0)) for i in range(3))
    #mixed = tuple(int(yellow[i] * (1 - value / 255.0) + purple[i] * (value / 255.0)) for i in range(3))
    #mixed = tuple(int(black[i] * (1 - value / 255.0) + white[i] * (value / 255.0)) for i in range(3))

    #hex_string = '#%02x%02x%02x' % mixed

    hex_string = '0x%02x%02x%02x' % mixed

    return hex_string

    #return mixed

# Test the function
#print(mix_colors(0)) # Should print (128, 0, 128) (pure purple)
#print(mix_colors(0.25))
#print(mix_colors(0.5)) # Should print (191, 128, 64) (purple and yellow mix)
#print(mix_colors(0.75))
#print(mix_colors(1)) # Should print (255, 255, 0) (pure yellow)

Peptide3d = pd.read_csv("/Users/administrator/Desktop/Peptide3D_HRV/PeptideUpload.CSV")

#Day
#Peptide3d = Peptide3d[Peptide3d.iloc[:,3] == 0]

#Peptide3d = Peptide3d[Peptide3d.iloc[:,3] >= 12]
#Peptide3d = Peptide3d[Peptide3d.iloc[:,3] <= 30]

#Peptide3d = Peptide3d[Peptide3d.iloc[:,3] >= 31]
#Peptide3d = Peptide3d[Peptide3d.iloc[:,3] <= 60]

Peptide3d = Peptide3d[Peptide3d.iloc[:,3] >= 68]

Peptide3d['AF'] = Peptide3d['AF'].div(100)

Peptide3d['AF'] = Peptide3d['AF'].apply(mix_colors)

for i in range(len(Peptide3d)):
    print(i)
    Peptide3d.iloc[i,1] = "  ///" + str(Peptide3d.iloc[i, 1]) + "/" + str(Peptide3d.iloc[i, 0]) + "/"
    Peptide3d.iloc[i, 0] = "color " + str(Peptide3d.iloc[i, 2])

del Peptide3d['AF']
del Peptide3d['Day']

Peptide3d.to_csv('PeptideHex', index = False, encoding='utf-8',header=False)

print(Peptide3d)

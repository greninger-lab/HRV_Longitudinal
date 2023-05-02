Epitopes 

Epitope_Translation.py is a python script which takes, ALL_LAVA_Output.csv, 1AYM_A.csv, and Epitope_Position_Updated.csv as input. It maps HRV-B and HRV-C species to PDB:1AYM and removes iSNVs that do not have corresponding residues in the protein model alignment. The program outputs Epitope_A_B_C.csv, which is an epitope list merged to the iSNV list to show which iSNVs occurred within an epitope region. It also produces the file 1AYM_Allignment_All.csv which is a data frame of HRV iSNVs serotypes mapped to PDB:1AYM. 


#Chi Square Epitope is an r script which preforms chi square 2x2 test on epitope output  

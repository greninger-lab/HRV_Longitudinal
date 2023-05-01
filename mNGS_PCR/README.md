mNGS_PCR

Files for the mNGS vs PCR qualitative analysis are keep in the mNGS_PCR folder. Reformate_mNGS_PCR_Data.py is a python script which takes the LAVA output files for mNGS + PCR ‘pairs’ and remformats the data frame by aligning them them side by side.

#PCR vs mNGS Graphs is a collection of short R scripts used to produce scatter plots comparing pairs to one another.

Epitopes

Epitope_Translation.py is a python script which takes, ALL_LAVA_Output.csv, 1AYM_A.csv

And Epitope_Position_Updated.csv as input. It maps RVB and RVC species to 1AYM and removes iSNVs that did not map properly.

The program produces Epitope_A_B_C.csv which is an epitope list merged to the iSNV list to show which iSNVs occurred within an epitope region. It also produces the file 1AYM_Allignment_All.csv which is a data frame of RV iSNVs serotypes mapped to 1AYM.

#Chi Square Epitope is an r script which preforms chi square 2x2 test on epitope output
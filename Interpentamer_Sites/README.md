Interpentamer_Sites

Files used for the interpentamer analysis are stored in the Interpentamer_Sites folder.

RV_Interpentamer_Sites.py is a python script which takes 1AYM_Allignment_All.csv,

interpentamer-sites.csv (interpentameter site list)and HRV-B14_1AYM.fasta (reference capsid with interpentameter sites mapped to 1AYM) as input. The program checks and removes residues for interpentamer sites which do not map to 1AYM

The script produces interpentamer_Position_Updated.csv which is a list interpentameter sites mapped to 1AYM andRV_Interpentamer_Sites.csv which contains all RV iSNVs and this residue is on an interppentamer site.

#Chi Square Interpentamer was an r script to preform chi square 2x2 test on interpentameter site output.
Structure_Protein

All files used for input are keep in the Structure_Protein folder. Emboss 6.5.7 was used to predict protein structure on 1AYM, was exported as gff file.

The python program Protein_Structure_Annotate.py takes 1AYM_full.gff (Emboss prediction gff) and 1AYM_Allignment_All.csv as input.

The script produces Annotated_Structure.csv which is a dataframe of iSNVs annotated with their predicted structure as well as rest of 1AYM residues with structure annotated merged to the iSNV dataframe.

#Structure VS Residue was an R script used to compare the structures on 1AYM residues occurring iSNV vs non iSNV sites. It also preforms 2x4 chi square test

viperdb_info

Extracted ataframes from the viperDP webpage of predictured surface exposure for 1AYM are keep in viperdb_info folder.

#viperdp extract surface information is an R script from the deep mutation study ( PMC7861617) which extracts json dataframes of the ViperDB database.
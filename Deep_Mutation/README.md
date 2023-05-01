Deep mutation

Deep mutational analysis input files are stored in the Deep_mutation folder. DeepMutation_Script.py is a python script which takes ALL_LAVA_Output.csv (all lava data frames combined), 1AYM_A.csv (Aligned positions for all RVA species to 1AYM, same file as SASA_Annotated.csv from Annotate_Pymol_Position.py just renamed for organizational purposes) and AVG_MAX_RP.csv (deep mutational scanning data) as input. This program checks and updates alignment positions for 1AYM, RVB, RVC to 4GB3 and epitope alignment RVA2 and RVB4 to 1AYM. It also maps the deep mutational scanning data (PMC7861617) to the iSNVs for this RV study.

The input for this program is epitopes-list_updated.csv (list of epitope positions)

Alignments to 4GB3 in the Alignments Folder (1AYM > 4GB3, RVB and RVC to 4GB3) and epitope alignments to 1AYM in theRV_A folder Epitope_Alignment, which contains contains HRV-A2_1AYM.fasta and HRV-B14_1AYM.fasta

The output for this program is deep mutational scanning data mapped to RV_A serotypes (RV_A.csv) RVB (RV_B.csv) and RVC (RV_C.csv). It also produces RV_Combo.csv A + B + C which combines the the dataframes for all RV iSNV serotypes. This program considers any gaps in alignment for all species and removes any residues which do not map to 4GB3 regardless of whether they contain an iSNV. Deep mutational scanning data is in relation to 4GB3, a picornavirus capsid used in the previous study. Epitope list mapped to 1AYM is exported as epitopes-list_updated.csv with epitopes that did not map properly being removed.

To produce some of the graphics for the deep mutational analysis a few r scripts were used.

#Graph AF with Deep mutation scanning data was used to check and display deep mutational scanning data mapped to iSNVs as a graph comparing MFElog2 vs AF.

#Deep mutation violin plots for iSNVs produces a violin plot of MFElog2 averages for residues that did and did not contain iSNVs for each VP region. This program can take RV_A, B, C or combo dataframes.

#Graph MFE VS Blossom Score

A blossom score sliding was added to analysis to check quality of alignment but was ultimately abandoned.
Alignments 

RV A serotypes (A102, A105, A39, A57, A58, A78, A82) were mapped to 1AYM to produce 3D capsid model and for surface area accessibility analysis. RVB06, RVB97, RVC28, RVC36 and 1AYM were mapped to 4GB3 for deep mutational scanning analysis. RVB06, RVB97, RVC28 and RVC36 were also mapped to 1AYM for easier comparison between serotypes for epitope, interpentamer, protein structure and surface exposure analysis. 

 

Reference fastas for each patient were annotated using ICTV (https://ictv.global/report/chapter/picornaviridae/picornaviridae/enterovirus). VP4, VP2, VP3 and VP1 regions were extracted and translated to a protein sequences. Each species was aligned using HHpred to 1AYM and 4GB3. iSNV that feel in place where there were gaps in either of the aligned species were removed for the analysis. 

 

All alignments are stored as Fastas in the allignments folder. 

 

RV_LAVA_DATAFRAMES 

All Lava output dataframes were combined for easier analysis with all files stored in RV_LAVA_DATAFRAMES folder. Lava output files were run through Combine_Lava_Data_Frames.py to combin all lava data frames while removing duplicate mutations within same patient .Output combined file is stored as ALL_LAVA_Output.csv. 

 

For figure generation #Lava Output Graph Generation  features a short series of r scripts to use to produce iSNV graphs and coverage. The input for these scripts is the ALL_LAVA_Output.csv + path to genomecov, which is output from the LAVA pipeline. Each sample is run individually with the path to each geome coverage plot updated between graphs. 

 

To adjust the minimum depth and AF this portion of the script must be altered: 

 

viz_graph <- working_df %>% filter((working_df$syn %in% c('nonsynonymous SNV', 'synonymous SNV')) & (working_df$af >= 10) & (working_df$depth > 30)) 

 

If there are no points that fit these parameters, ‘ghost points’ need to be added for graph to work, with their position outside max xlim of graph in order to not show in final visualization. 

 

Deep mutation 

 

Deep mutational analysis input files are stored in the Deep_mutation folder. DeepMutation_Script.py  is a python script which takes ALL_LAVA_Output.csv (all lava data frames combined), 1AYM_A.csv (Aligned positions for all RVA species to 1AYM, same file as SASA_Annotated.csv from Annotate_Pymol_Position.py just renamed for organizational purposes) and AVG_MAX_RP.csv (deep mutational scanning data) as input. This program checks and updates alignment positions for 1AYM, RVB, RVC to 4GB3 and epitope alignment RVA2 and RVB4 to 1AYM. It also maps the deep mutational scanning data (PMC7861617) to the iSNVs for this RV study.  

 

The input for this program is epitopes-list_updated.csv (list of epitope positions) 

Alignments to 4GB3 in the Alignments Folder (1AYM > 4GB3, RVB and RVC to 4GB3) and epitope alignments to 1AYM in theRV_A folder Epitope_Alignment, which contains contains HRV-A2_1AYM.fasta and HRV-B14_1AYM.fasta 

 

 

The output for this program is deep mutational scanning data mapped to  RV_A serotypes (RV_A.csv) RVB (RV_B.csv) and RVC (RV_C.csv). It also produces RV_Combo.csv A + B + C which combines the the dataframes for all RV iSNV serotypes. This program considers any gaps in alignment for all species and removes any residues which do not map to 4GB3 regardless of whether they contain an iSNV. Deep mutational scanning data is in relation to 4GB3, a picornavirus capsid used in the previous study. Epitope list mapped to 1AYM is exported as epitopes-list_updated.csv with epitopes that did not map properly being removed. 

 

To produce some of the graphics for the deep mutational analysis a few r scripts were used.  

 

#Graph AF with Deep mutation scanning data was used to check and display deep mutational scanning data mapped to iSNVs as a graph comparing MFElog2 vs AF. 

 

 #Deep mutation violin plots for iSNVs produces a violin plot of MFElog2 averages for residues that did and did not contain iSNVs for each VP region. This program can take RV_A, B, C or combo dataframes. 

 

#Graph MFE VS Blossom Score  

A blossom score sliding was added to analysis to check quality of alignment but was ultimately abandoned. 

 

MAF_vs_SASA 

Files for 3d capsid model generation and surface area accessibility analysis are stored in the MAF_vs_SASA folder.  

 

Add_Colors_Pymol.py takes SASA_Annotated.csv as input and outputs iSNV residues for 1AYM with color based on a scale range between red and yellow extracted from the AF percentage. Output can also be adjusted to filter based on day range. 

 

Pymol_Inputs.txt has a brief description to generate the 1AYM capsid in pymol and color each VP region. Paste the output of Add_Colors_Pymol.py to color all the residues. 

 

SASA_Pymol.txt was produced using the surface area accessibility function in Pymol. Capsid was dehydrated prior to surface area analysis. SASA_Pymol.txt  was converted to SASA_Extracted_Pymol.csv in excel using text to column options to remove unnecessary columns. 

 

#Pymol Capsid Graph is an r script used used to add color scale to exported pymol images, and adds color guide for the VP regions. Further repositioning of elements was done in Inscape. 

 

Annotate_Pymol_Position.py takes SASA_Dataframe.csv, SASA_Extracted_Pymol.csv 

and alignments for  RV_A in the Alignments folder as input. The python program maps all iSNVs to 1AYM for RV A, removing any iSNVs that don’t map properly (gaps) and also adds surface area accessibly information to each iSNV stored as the output file SASA_Annotated.csv. 

 

#MAF vs SASA was an R script used to graph AF against surface area accessibility with iSNVs colored by group and the graph faceted V region. No VP4 region is fetured in the graph as it is mostly internal. 

 

mNGS_PCR 

Files for the mNGS vs PCR qualitative analysis are keep in the mNGS_PCR folder. Reformate_mNGS_PCR_Data.py is a python script which takes the LAVA output files for mNGS + PCR ‘pairs’ and remformats the data frame by aligning them them side by side. 

 

#PCR vs mNGS Graphs is a collection of short R scripts used to produce scatter plots comparing pairs to one another. 

 

Epitopes 

Epitope_Translation.py is a python script which takes, ALL_LAVA_Output.csv, 1AYM_A.csv 

And Epitope_Position_Updated.csv as input. It maps RVB and RVC species to 1AYM and removes iSNVs that did not map properly.  

The program produces Epitope_A_B_C.csv which is an epitope list merged to the iSNV list to show which iSNVs occurred within an epitope region. It also produces the file 1AYM_Allignment_All.csv which is a data frame of RV iSNVs serotypes mapped to 1AYM. 

 

#Chi Square Epitope is an r script which preforms chi square 2x2 test on epitope output  

 

Interpentamer_Sites 

Files used for the interpentamer analysis are stored in the Interpentamer_Sites folder. 

RV_Interpentamer_Sites.py is a python script which takes 1AYM_Allignment_All.csv, 

interpentamer-sites.csv (interpentameter site list)and HRV-B14_1AYM.fasta (reference capsid with interpentameter sites mapped to 1AYM) as input. The program checks and removes residues for interpentamer sites which do not map to 1AYM 

The script produces interpentamer_Position_Updated.csv which is a list interpentameter sites mapped to 1AYM andRV_Interpentamer_Sites.csv which contains all RV iSNVs and this residue is on an interppentamer site.  

 

#Chi Square Interpentamer was an r script to preform chi square 2x2 test on interpentameter site output. 

 

Structure_Protein 

All files used for input are keep in the Structure_Protein folder. Emboss 6.5.7 was used to predict protein structure on 1AYM, was exported as gff file. 

The python program Protein_Structure_Annotate.py takes 1AYM_full.gff (Emboss prediction gff) and 1AYM_Allignment_All.csv as input. 

The script produces Annotated_Structure.csv which is a dataframe of iSNVs annotated with their predicted structure as well as rest of 1AYM residues with structure annotated merged to the iSNV dataframe. 

 

#Structure VS Residue was an R script used to compare the structures on 1AYM residues occurring iSNV vs non iSNV sites. It also preforms 2x4 chi square test 

 

viperdb_info 

Extracted ataframes from the viperDP webpage of predictured surface exposure for 1AYM are keep in viperdb_info folder. 

 

#viperdp extract surface information is an R script from the deep mutation study (  PMC7861617) which extracts json dataframes of the ViperDB database. 

 

Surface_Exposure 

Input and output files for the surface exposure analysis are stored in the Surface_Exposure 

 Folder. The python program SurfaceExposure.py take the folder viperdb_info and the file 

1AYM_Allignment_All.csv as input, in order to map viperDP surface exposure to iSNVs. As a note the output for the residues on viperDB are annotated as such: 

#D = VP4 

#B = VP2 

#C = VP3 

#A = VP1 

 

Annotated iSNVs are stored as Surface_Exposure_Annotated.CSV. 

 

#Structure position is an r script which graphs proportional frequency of position of residue and whether an iSNV vs non iSNV occurred on that residue. It also preforms 2x4 chi square test of all positions for iSNV vs non iSNV and a 2x2 chi square test of surface vs non surface structures, iSNV non iSNV residues. 

 

1_AYM_All_Dataframes 

Files for combining the epitope analysis, interpentamer analysis, structure analysis, MAF analysis and surface exposure analysis are stored in the 1_AYM_All_Dataframes folder. 

1AYM_MergeDF.py  takesRV_Interpentamer_Sites.csv, Epitope_A_B_C.csv, Annotated_StructureWTDpulicates.csv, Surface_Exposure_Annotated_WTDuplicates.CSV as input and produces 1AYM_MergedDataframes.csv. 

 

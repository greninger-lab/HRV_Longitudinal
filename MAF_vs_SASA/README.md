MAF_vs_SASA

Files for 3d capsid model generation and surface area accessibility analysis are stored in the MAF_vs_SASA folder.

Add_Colors_Pymol.py takes SASA_Annotated.csv as input and outputs iSNV residues for 1AYM with color based on a scale range between red and yellow extracted from the AF percentage. Output can also be adjusted to filter based on day range.

Pymol_Inputs.txt has a brief description to generate the 1AYM capsid in pymol and color each VP region. Paste the output of Add_Colors_Pymol.py to color all the residues.

SASA_Pymol.txt was produced using the surface area accessibility function in Pymol. Capsid was dehydrated prior to surface area analysis. SASA_Pymol.txt was converted to SASA_Extracted_Pymol.csv in excel using text to column options to remove unnecessary columns.

#Pymol Capsid Graph is an r script used used to add color scale to exported pymol images, and adds color guide for the VP regions. Further repositioning of elements was done in Inscape.

Annotate_Pymol_Position.py takes SASA_Dataframe.csv, SASA_Extracted_Pymol.csv

and alignments for RV_A in the Alignments folder as input. The python program maps all iSNVs to 1AYM for RV A, removing any iSNVs that don’t map properly (gaps) and also adds surface area accessibly information to each iSNV stored as the output file SASA_Annotated.csv.

#MAF vs SASA was an R script used to graph AF against surface area accessibility with iSNVs colored by group and the graph faceted V region. No VP4 region is fetured in the graph as it is mostly internal.
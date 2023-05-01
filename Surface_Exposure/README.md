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
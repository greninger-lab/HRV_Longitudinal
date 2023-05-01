RV_LAVA_DATAFRAMES

All Lava output dataframes were combined for easier analysis with all files stored in RV_LAVA_DATAFRAMES folder. Lava output files were run through Combine_Lava_Data_Frames.py to combin all lava data frames while removing duplicate mutations within same patient .Output combined file is stored as ALL_LAVA_Output.csv.

For figure generation #Lava Output Graph Generation features a short series of r scripts to use to produce iSNV graphs and coverage. The input for these scripts is the ALL_LAVA_Output.csv + path to genomecov, which is output from the LAVA pipeline. Each sample is run individually with the path to each geome coverage plot updated between graphs.

To adjust the minimum depth and AF this portion of the script must be altered:

viz_graph <- working_df %>% filter((working_df$syn %in% c('nonsynonymous SNV', 'synonymous SNV')) & (working_df$af >= 10) & (working_df$depth > 30))

If there are no points that fit these parameters, ‘ghost points’ need to be added for graph to work, with their position outside max xlim of graph in order to not show in final visualization.
# Microbiome_analysis

In this folder you will find the Kraken2 reports for the samples used to the microbiome analysis in the article "Within-host rhinovirus evolution in upper and lower respiratory tract highlights capsid variability and mutation-independent compartmentalization".

*.report files = Kraken2 reports for each sample inferred using the standard prebuilt database v20221209.

Microbiome.Rmd = rMarkdown including the code for making the analysis and plots of taxonomy abundance and beta diversity. It works in RStudio with the files in the HRV-BIOM-files folder. It is recommended to review the path where the folder is downloaded and correct the path informed in the pipeline properly.

HRV-BIOM-files/BIOM.biom = BIOM table created with the Kraken2 reports.

HRV-BIOM-files/annotation_variables.csv = Table with the variables to annotate the BIOM table.

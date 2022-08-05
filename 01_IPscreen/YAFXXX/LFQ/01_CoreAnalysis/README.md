# RBP_label_free_analysis.R

Script to analyze the label free immmunoprecipitation screen experiments. 

## Input files

It requires the following input files:

- proteinGroups.txt: **Necesssary** for the analysis. It is on of the regular output files from a standard [MaxQuant](https://www.maxquant.org/) analysis. It contains the protein intensities for each experiment (raw, and LFQ and iBAQ normalized) and its associated descriptor columns (i. e., number of peptides, fasta header...). For this project, we had 12 intensity related columns for each RBP-bait; three conditions, untreated pulldown (IP), treated pulldown (RNase) and wild type (WT) in quadruplicates. 

- highlight.txt: **Not necessary** for the analysis. It contains the protein IDs you want to highlight along the graphical outputs. For this project, if provided, the IDs **must** be systematic *S. cerevisiae* IDs (i. e., YER165W for PAB1). If not provided, no ID will be highlighted. 

- volcano_contrasts.txt: **Not necessary** for the analysis. It contains two colums (left and right) in which the user can specify the direction of the contrasts in the statistical analysis as well as the volcano plot output. For instance, if you can to compare *treated* against the *control* this would go on the right and left column respectively. For this project, we were always interested in the same comparisons; for the PPI group we compared a treated pulldown (RNase, right) to a wild type (WT, left). For the RDI we compared a non treated pulldown (IP, right) to a treated pulldown (RNase, left). If not provided, a combination of all possible contrasts will be generated automatically. 

A minimal example of each file can be found in this folder.
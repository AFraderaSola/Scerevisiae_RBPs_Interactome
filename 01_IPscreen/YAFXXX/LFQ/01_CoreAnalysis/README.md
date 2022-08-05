# RBP-bait label free quantification core analysis

Script to analyze the label free immmunoprecipitation screen RBP-baits experiments. 

## Input files

It requires the following input files:

<<<<<<< HEAD
- **_proteinGroups.txt_**:  It is on of the **regular output** files from a standard [MaxQuant](https://www.maxquant.org/) analysis. It contains the **protein intensities** for each experiment (raw, and LFQ and iBAQ normalized) and its associated descriptor columns (i. e., number of peptides, fasta header...). For this project, we had 12 intensity related columns for each RBP-bait; three conditions, untreated pulldown (IP), treated pulldown (RNase) and wild type (WT) in quadruplicates. **Necesssary** for the analysis.

- **_highlight.txt_**:  It contains the **protein IDs** you want to highlight along the graphical outputs. For this project, if provided, the IDs **must** be systematic *S. cerevisiae* IDs (i. e., YER165W for PAB1). If not provided, no ID will be highlighted. **Not necessary** for the analysis.

- **_volcano_contrasts.txt_**:  It contains two colums (right and left) in which the user can specify the **direction of the contrasts** in the statistical analysis as well as the volcano plot output. For instance, if you can to compare *treated* against the *control* this would go on the right and left column respectively. For this project, we were always interested in the same comparisons; for the PPI group we compared a treated pulldown (RNase, right) to a wild type (WT, left). For the RDI we compared a non treated pulldown (IP, right) to a treated pulldown (RNase, left). If not provided, a combination of all possible contrasts will be generated automatically. **Not necessary** for the analysis.
=======
- **_proteinGroups.txt_**: **Necesssary** for the analysis. It is on of the regular output files from a standard [MaxQuant](https://www.maxquant.org/) analysis. It contains the protein intensities for each experiment (raw, and LFQ and iBAQ normalized) and its associated descriptor columns (i. e., number of peptides, fasta header...). For this project, we had 12 intensity related columns for each RBP-bait; three conditions, untreated pulldown (IP), treated pulldown (RNase) and wild type (WT) in quadruplicates. 

- **_highlight.txt_**: **Not necessary** for the analysis. It contains the protein IDs you want to highlight along the graphical outputs. For this project, if provided, the IDs **must** be systematic *S. cerevisiae* IDs (i. e., YER165W for PAB1). If not provided, no ID will be highlighted. 

- **_volcano_contrasts.txt_**: **Not necessary** for the analysis. It contains two colums (right and left) in which the user can specify the direction of the contrasts in the statistical analysis as well as the volcano plot output. For instance, if you can to compare *treated* against the *control* this would go on the right and left column respectively. For this project, we were always interested in the same comparisons; for the PPI group we compared a treated pulldown (RNase, right) to a wild type (WT, left). For the RDI we compared a non treated pulldown (IP, right) to a treated pulldown (RNase, left). If not provided, a combination of all possible contrasts will be generated automatically. 
>>>>>>> main

A minimal example of each file can be found in this folder.

## Script

The script performs the core label free quantification (LFQ) analysis. Its goal is to quantify proteins across conditions and asses signicant differencesa amongts them. This way we can determine which proteins were enriched or depleted per condition. Through the script, we cover the following basic steps:

- **_Data wrangling_**: We start by tidying up and filtering the **proteinGroups.txt** so it is ready for the analysis. In brief, we do the following steps:
  - **Add a gene name column**, so the sytematic IDs (i. e., YER165W) are readable (i. e., PAB1).
  - Filter out **the contaminants**, which are usually present in any MS experiment (i. e., keratin). A list of usual contaminants is provided by [MaxQuant](http://www.coxdocs.org/doku.php?id=maxquant:start_downloads.htm).
  - Filter out **reverse peptides**.
  - Filter peptides by **razor and unique peptide** number, which are user defined. For instance, we can decide to keep only proteinGroups with a minimum of razor+unique peptide count higher than 2 and a minimum of unique peptide count higher than 1, along all proteinGroups.txt file.
  - Filter out by **quantification event**, which is also user defined. For instance, we can decide to keep only the proteinGroups with a minimum of 2 quantification events along the replicates. This means that any protein not having an intensity value associated (not NA) to, at least, two out of four replicates will be filtered out. 
  
- **_Missing values imputation_**: A known issue of MS experiments is the **abundance of NA values** along the detected protein groups due to technical limitations. For instance, a protein group could be identified with an average intensity of 25 in three replicates and have an NA value on the fourth. Our approach to NA values is assume that the protein was **not detected due to technical limitations** not because **it was not present in protein mixture**. This way, we decide to impute an intensity value for all NA values still present **after the filtering steps**. We do so by shifting a **beta distribution** obtained from the LFQ intensity values to the **limit of quantitation**. This way, the resulting imputed values will always be random values constricted to the lowest end of our LFQ intensities.  

- **_Quality control_**: 

- **_Exploratory analysis_**:

- **_Statistical analysis_**:
  
# RBP-bait label free quantification core analysis

Script to analyze the label free immmunoprecipitation screen RBP-baits experiments. 

## Input files

It requires the following input files:

- **_proteinGroups.txt_**:  It is on of the **regular output** files from a standard [MaxQuant](https://www.maxquant.org/) analysis. It contains the **protein intensities** for each experiment (raw, and LFQ and iBAQ normalized) and its associated descriptor columns (i. e., number of peptides, fasta header...). For this project, we had 12 intensity related columns for each RBP-bait; three conditions, untreated pulldown (IP), treated pulldown (RNase) and wild type (WT) in quadruplicates. **Necesssary** for the analysis.

- **_highlight.txt_**:  It contains the **protein IDs** you want to highlight along the graphical outputs. For this project, if provided, the IDs **must** be systematic *S. cerevisiae* IDs (i. e., YER165W for PAB1). If not provided, no ID will be highlighted. **Not necessary** for the analysis.

- **_volcano_contrasts.txt_**:  It contains two colums (right and left) in which the user can specify the **direction of the contrasts** in the statistical analysis as well as the volcano plot output. For instance, if you can to compare *treated* against the *control* this would go on the right and left column respectively. For this project, we were always interested in the same comparisons; for the PPI group we compared a treated pulldown (RNase, right) to a wild type (WT, left). For the RDI we compared a non treated pulldown (IP, right) to a treated pulldown (RNase, left). If not provided, a combination of all possible contrasts will be generated automatically. **Not necessary** for the analysis.

A minimal example of each file can be found in this folder. Additionally, it requires the 00_IPScreeningMasterFile.xlsx (01_IPscreen folder) file to match the stock ID, YAFXXX, to its systematic and gene name. 

## RBP_label_free_analysis.R Script

The script performs the core label free quantification (LFQ) analysis. Its goal is to quantify proteins across conditions and asses significant differences among them. This way we can determine which proteins were enriched or depleted per condition. Through the script, we cover the following basic steps:

- **_Data wrangling_**: We start by tidying up and filtering the **proteinGroups.txt** so it is ready for the analysis. In brief, we do the following steps:
  - **Add a gene name column**, so the systematic IDs (i. e., YER165W) are readable (i. e., PAB1).
  - Filter out **the contaminants**, which are usually present in any MS experiment (i. e., keratin). A list of usual contaminants is provided by [MaxQuant](http://www.coxdocs.org/doku.php?id=maxquant:start_downloads.htm).
  - Filter out **reverse peptides**.
  - Filter peptides by **razor and unique peptide** number, which are user defined. For instance, we can decide to keep only proteinGroups with a minimum of razor+unique peptide count higher than 2 and a minimum of unique peptide count higher than 1, along all proteinGroups.txt file.
  - Filter out by **quantification event**, which is also user defined. For instance, we can decide to keep only the proteinGroups with a minimum of 2 quantification events along the replicates. This means that any protein not having an intensity value associated (not NA) to, at least, two out of four replicates will be filtered out. 
  
- **_Missing values imputation_**: A known issue of MS experiments is the **abundance of NA values** along the detected protein groups due to technical limitations. For instance, a protein group could be identified with an average intensity of 25 in three replicates and have an NA value on the fourth. Our approach to NA values is assume that the protein was **not detected due to technical limitations** not because **it was not present in protein mixture**. This way, we decide to impute an intensity value for all NA values still present **after the filtering steps**. We do so by shifting a **beta distribution** obtained from the LFQ intensity values to the **limit of quantitation**. This way, the resulting imputed values will always be random values constricted to the lowest end of our LFQ intensities.  

- **_Quality control_**: Once the values are imputed we do a series of **quality control** check:
  - Visualize the **imputed** and **original** values distribution per condition and replica, which should be **skewed to the lower** end of intensities for the imputed ones.
  - Visualize the counts of **imputed** and **original** intensity values per condition/replica.
  - Visualize **raw** and **LFQ normalized** intensities.
  - Visualize the **intensity distribution** per  conditions/value count.
 
- **_Exploratory analysis_**: With our data ready and quality control performed, we start the proper analysis by doing an exploratory analysis which will tell as the **data behaviour** and the **sample similarity**. In brief, we do the following steps:
  - Calculate the **pearson correlation coefficient** between conditions and replicas which are visualized, with and without clustering, with a **hetamap**.
  - Perform a sample **dimensionality reduction** with a principal cluster analysis (PCA). The first three components are visualized with a **scatter plot**. 

- **_Statistical analysis_**: Finally we assess differences in protein LFQ intensities between conditions by doing a [t. test](https://en.wikipedia.org/wiki/Student%27s_t-test) which are visualized with a volcano plot. 

## Output files

If the script runs properly, it will generate a series of output files that can roughly categorized into the following categories:

- **_Tables_**: Either in .txt or .csv format, a number of tables will be generated containing the analysis results. Noteworthy are the following:
  - YAFXXX_QuantifiedProteins.csv contains **all** the proteins detected in **all** conditions after the filtering steps and its associated intensity, imputed intensity, and statistical values.
  - YAFXXX_Enriched_RNase_vs_WT.csv contains the **significantly enriched** proteins from comparing a treated pulldown (RNase) to a wild type (WT). We obtain the **PPI** set for each bait from this file.
  - YAFXXX_Enriched_IP_vs_RNase.csv contains the **significantly enriched** proteins from comparing an untreated pulldown (IP) to a treated pulldown (RNase). We obtain the **RDI** set for each bait from this file.

Supplementary tables 1 and 2 from the **"RNA-dependent interactome allows network-based assignment of RBP function"** publication were generated from the output tables of this script. A table summary can be found on the 00_IPScreeningMSMasterFile.xlsx file, while the combined Quantified and Enriched tables for each RBP bait can be found at  01_CoreAnalysis_LFQ_02_QuantifiedProteins.xlsx, 01_CoreAnalysis_LFQ_03_PPI_EnrichedProteins.xlsx and 01_CoreAnalysis_LFQ_04_RNAm_EnrichedProteins.xlsx respectively (01_IPscreen folder).

- **_Plots_**: The previously mention plots (RBP_label_free_analysis.R Script section) will be generated in .pdf.

Figure 2A and 2B from the **"RNA-dependent interactome allows network-based assignment of RBP function"** publication were generated from the volcano plots obtained within this script. The script output's **plots for all RBP-baits** are included in the 01_CoreAnalysis_LFQ_01_Summary.pdf file, in the 01_IPscreen folder.
  
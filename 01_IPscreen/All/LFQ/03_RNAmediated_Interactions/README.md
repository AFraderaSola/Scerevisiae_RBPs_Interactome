# 03_RNAmediated_Interactions

This **sub-folder** contains all the scripts from the **RDI** group. 

Within this folder you can find the data and scripts to generate the following figures:

- Figure 2D

- Figure 5B

- Figure 6B

It is organized in **7 sub-folders**, usually an **Input** and **Output** folder for each script (the prefix number of the folders matches the one of the scripts)

## 00_IPsResultsFiles

This **sub-folder** contains the **quantified proteins** and the **enriched proteins** for the **RDI** groups for each YAFXXX strain.

## 01_InputFiles

This **sub-folder** contains the input files from the **01_Interactors_SGDandCensus.R** script.

## 01_OutputFiles

This **sub-folder** contains the output files from the **01_Interactors_SGDandCensus.R** script, including **Figure 2D**.

## 01_Interactors_SGDandCensus.R

This script summarizes the enriched proteins for the **RDI** group in each YAFXX strain and creates a **bar plot** showing their **inclusion** at **BioGRID** and the **RBP census**. With this script, you can generate the following figures:

- Figure 2D

## 03_InputFiles

This **sub-folder** should contain the input files from the **03_InterPro_FunctionalAnalysis.R** script, but they are too heavy for github. These are .fasta files containing the aminoacid sequences of all enriched and quantified yeast proteins. If you need them, contact me. 

## 03_OutputFiles

This **sub-folder** contains the output files from the **03_InterPro_FunctionalAnalysis.R** script.

## 03_InterPro_FunctionalAnalysis.R

This script takes the **raw results** obtained with interpro and summarize the relevant results (RNA related domains and protein families).

## 04_InputFiles

This **sub-folder** contains the input files from the **04_Full_Network_v02.R** and the **04_Individual_Networks_v02.R** scripts.

## 04_OutputFiles

This **sub-folder** contains the output files from the **04_Fulll_Network_v02.R**  and the **04_Individual_Networks_v02.R** script, including **Figure 5B** and **Figure 6B**.

## 04_Full_Network_v02.R

This script generates a **global network**, highlighting all functionalities, for the  **PPI** group in each YAFXX strain. With this script, you can generate the following figures:

- Figure 5B

## 04_Individual_Networks_v02.R

This script generates an **individual functional network** for the  **PPI** group in each YAFXX strain. With this script, you can generate the following figures:

- Figure 6B

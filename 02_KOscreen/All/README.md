# All

This **sub-folder** contains all the scripts that combines the results from all **KO strains**. 

Within this folder you can find the data and scripts to generate the following figures:

- Figure 4C, 4D and 4E

It is organized in **7 sub-folders**, usually an **Input** and **Output** folder for each script (the prefix number of the folders matches the one of the scripts)

## 00_IPsResultsFiles

This **sub-folder** contains the **quantified proteins** and the **enriched proteins** for the **PPI** and **RDI** groups for each YAFXXX strain.

## 00_KOsResultsFiles

This **sub-folder** contains the **quantified proteins** and the **enriched proteins** for the **KO-RBP** YKAFXXX strains.

## 01_OutputFiles

This **sub-folder** contains the output files from the **01_Enriched_IDs_overlap_IPs.R** script, including **Figure 4C**.

## 01_Enriched_IDs_overlap_IPs.R

This script shows the enriched proteins for the **RBP-KO** YKAFXX strains and creates a **bar plot** showing the **overlap** to the **PPI** and **RDI** immunoprecipitation screen groups. With this script, you can generate the following figures:

- Figure 4C

## 03_InputFiles

This **sub-folder** contains the input files from the **03_KEGG_Enrichment.R** script, including **Figure 4D**.

## 03_OutputFiles

This **sub-folder** contains the output files from the **03_KEGG_Enrichment.R** script, including **Figure 4D**.

## 03_KEGG_Enrichment.R

This script compares the KEGG enrichment results for the **KO** enriched proteins in each YKAFXX strain and creates a **heat map** comparing them. With this script, you can generate the following figures:

- Figure 4D

## 07_InputFiles

This **sub-folder** contains the input files from the **07_ProteinComplexAnalysisPerGroup.R** script, including **Figure 4E**.

## 07_OutputFiles

This **sub-folder** contains the output files from the **07_ProteinComplexAnalysisPerGroup.R** script, including **Figure 4E**.

## 07_ProteinComplexAnalysisPerGroup.R

This script generates the **CYC2008 complex networks** for  **KO** enriched proteins in each YKAFXX strain. With this script, you can generate the following figures:

- Figure 4E

# 01_Both

This **sub-folder** contains all the scripts that combines data from the **PPI** and **RDI** groups. 

Within this folder you can find the data and scripts to generate the following figures:

- Figure 2E

- Figure 3A, 3B and 3C

- Figure 5A and 5B

- Figure 6C, 6D and 6E

- Figure 7

It is organized in **XXX sub-folders**, usually an **Input** and **Output** folder for each script (the prefix number of the folders matches the one of the scripts)

## 00_IPsResultsFiles

This **sub-folder** contains the **quantified proteins** and the **enriched proteins** for the **PPI** and **RDI** groups for each YAFXXX strain.

## 01_OutputFiles

This **sub-folder** contains the output files from the **01_Overlap_PPIandRNAm.R** script, including **Figure 2E**.

## 01_Overlap_PPIandRNAm.R

This script compares the enriched proteins for the **PPI** and **RDI** groups in each YAFXX strain and creates a **bar plot** showing the **overlap** among them. With the script, you can generate the following figures:

- Figure 2E

## 02_InputFiles

This **sub-folder** contains the output files from the **02_KEGG_Enrichment.R** script, including **Figure 3C**.

## 02_OutputFiles

This **sub-folder** contains the output files from the **02_KEGG_Enrichment.R** script, including **Figure 3C**.

## 02_KEGG_Enrichment.R

This script compares the KEGG enrichment results for the  **PPI** and **RDI** groups in each YAFXX strain and creates a **heat map** comparing them. With the script, you can generate the following figures:

- Figure 3C
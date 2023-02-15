# 01_Both

This **sub-folder** contains all the scripts that combines data from the **PPI** and **RDI** groups. 

Within this folder you can find the data and scripts to generate the following figures:

- Figure 2E

- Figure 3A, 3B and 3C

- Figure 5A and 5B

- Figure 6C, 6D and 6E

- Figure 7

It is organized in **14 sub-folders**, usually an **Input** and **Output** folder for each script (the prefix number of the folders matches the one of the scripts)

## 00_IPsResultsFiles

This **sub-folder** contains the **quantified proteins** and the **enriched proteins** for the **PPI** and **RDI** groups for each YAFXXX strain.

## 01_OutputFiles

This **sub-folder** contains the output files from the **01_Overlap_PPIandRNAm.R** script, including **Figure 2E**.

## 01_Overlap_PPIandRNAm.R

This script compares the enriched proteins for the **PPI** and **RDI** groups in each YAFXX strain and creates a **bar plot** showing the **overlap** among them. With this script, you can generate the following figures:

- Figure 2E

## 02_InputFiles

This **sub-folder** contains the input files from the **02_KEGG_Enrichment.R** script.

## 02_OutputFiles

This **sub-folder** contains the output files from the **02_KEGG_Enrichment.R** script, including **Figure 3C**.

## 02_KEGG_Enrichment.R

This script compares the KEGG enrichment results for the  **PPI** and **RDI** groups in each YAFXX strain and creates a **heat map** comparing them. With this script, you can generate the following figures:

- Figure 3C

## 03_InputFiles

This **sub-folder** contains the input files from the **03_Individual_Networks_v02.R** script.

## 03_OutputFiles

This **sub-folder** contains the output files from the **03_Individual_Networks_v02.R** script, including **Figure 6C**.

## 03_Individual_Networks_v02.R

This script generates an **individual functional network** for the  **PPI** and **RDI** groups in each YAFXX strain. With this script, you can generate the following figures:

- Figure 6C

## 05_InputFiles

This **sub-folder** contains the input files from the **05_ClevelandDotPlot.R** script.

## 05_OutputFiles

This **sub-folder** contains the output files from the **05_ClevelandDotPlot.R** script, including **Supplementary Figure 4**.

## 05_ClevelandDotPlot.R

This script generates a **cleveland dot plot** and a **bar plot** for the  **PPI** and **RDI** groups in each YAFXX strain, showing the inclusion at **RBP census** and the **BioGRID**, respectively. With this script, you can generate the following figures:

- Supplementary Figure 4

## 06_InputFiles

This **sub-folder** contains the input files from the **06_GOMF_FigureRDIandPPI.R** and **06_InterproFigureRDIandPPI.R** scripts.

## 06_OutputFiles

This **sub-folder** contains the output files from the **06_GOMF_FigureRDIandPPI.R** and **06_InterproFigureRDIandPPI.R** scripts, including **Figure 3A** and **Figure 3B**.

## 06_GOMF_FigureRDIandPPI.R

This script generates a **bar plot** for the  **PPI** and **RDI** groups in each YAFXX strain, showing the **top5 most enriched GO molecular function terms** . With this script, you can generate the following figures:

- Figure 3B

## 06_InterproFigureRDIandPPI.R

This script generates a **bar plot** for the  **PPI** and **RDI** groups in each YAFXX strain, showing the interpro **overrepresented PFAM and SUPERFAMILY terms** . With this script, you can generate the following figures:

- Figure 3A

## 07_InputFiles

This **sub-folder** contains the input files from the **07_ProteinComplexAnalysis.R** and **07_ProteinComplexAnalysisPerGroup.R** scripts.

## 07_OutputFiles

This **sub-folder** contains the output files from the **07_ProteinComplexAnalysis.R** and **07_ProteinComplexAnalysisPerGroup.R** scripts, including **Figure 6D** and **Figure 6E**.

## 07_ProteinComplexAnalysis.R

This script generates a **dot plot** and a **heat map** for the  **PPI** and **RDI** groups in each YAFXX strain, showing the **CYC2008 complex analysis**. With this script, you can generate the following figures:

- Figure 6D
- Figure 6E

## 07_ProteinComplexAnalysisPerGroup.R

This script generates a **dot plot** and a **heat map** for the  **PPI** and **RDI** groups in each YAFXX strain, showing the **CYC2008 complex analysis**. With this script, you can generate the following figures:

- Figure 6D
- Figure 6E

## 08_InputFiles

This **sub-folder** contains the input files from the **08_ComplexNetworks.R** script.

## 08_OutputFiles

This **sub-folder** contains the output files from the **08_ComplexNetworks.R** script, including **Figure 7**.

## 08_ComplexNetworks.R

This script generates the **CYC2008 complex networks** for  **PPI** and **RDI** groups in each YAFXX strain. With this script, you can generate the following figures:

- Figure 7

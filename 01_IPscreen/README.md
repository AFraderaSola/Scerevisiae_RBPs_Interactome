# 01_IPscreen folder

This folder contains all the scripts used to analyze the MS data for the **immunoprecipitation screen** described in Figure 1B. Within this folder you can find the data and scripts to generate the following figures:

- Figure 2A, 2B, 2C, 2D and 2E

- Figure 3A, 3B, 3C, 3D and 3E

- Figure 5A, 5B, 5C and 5D

- Figure 6A, 6B, 6C and 6D

- Figure 7A, 7B, 7C and 7D

It also contains the scripts to generate Supplementary tables 1 to 8. 

It is organized into a couple sub-folders containing the analysis scripts and several files summarizing the IP screen findings. 

## All

This **sub-folder** contains a downstream analysis **including all strains**. The **input files** needed for those analysis are generated on each RBP-bait sub-folder (**YAFXXX**). The scripts to generate the previously mentioned figures are inside this folder. More details on the scripts and how to run them inside this folder.

## YAFXXX

This folder contains an example on how the main screen analysis was structured. Each **RBP-bait** had its own **sub-folder**, named after its stock strain (**YAFXXX**), containing its **Core LFQ analysis** and its **Functional Analysis**, this last one diveded already for PPI and RDI groups.

To run the Core LFQ analysis, you **need the proteinGroups.txt** (obtained from analyzing .RAW MS files with [MaxQuant](https://www.maxquant.org/)). The file should be placed at the **same folder** the analysis script is (example inside the folder). The proteinGroups file for each RBP-bait can be found on the 00_IPscreeningMasterproteinGroups.zip file. 

Additionally you can provide a **highlight.txt** and **volcano_contrasts.txt**, but this won't be necessary for the script to run. Examples of such files are provided inside the folder. 

The scripts for the **Core LFQ Analysis** and the **Functional Analysis** work in a sequential fashion; you need the results of the first to run the second.

More details on the scripts and how to run them  can be found inside each folder.

## Master Files 

First we have the MasterFiles, containing general information on the strains:

- 00_IPScreeningMasterFile.xlsx: Here you can find all the necessary information for each YAFXXX strains (i.e., gene name or selection criteria) and a summary of it data analysis (i. e., quantified proteins or enriched proteins)

- 00_IPScreeningMSMasterFile.xlsx: Here you can find the related MS information (i.e, measured date or raw file) for each YAFXXX strain.

- 00_IPscreeningMasterproteinGroups.zip: Here you can find all the protein groups for each RBP-bait.

## Core Analysis Files

Then we have the CoreAnalysis files; this files summarize the main analysis conducted on each YAFXXX strain. More details on the analysis itself can be found on each YAFXXX folder.

- 01_CoreAnalysis_LFQ_01_Summary.pdf: Here you can find a summary the CoreAnalysis figures for each YAFXXX strain. 

- 01_CoreAnalysis_LFQ_02_QuantifiedProteins.xlsx: Here you can find all quantified proteins for each YAFXXX strain and their associated statistical values.

- 01_CoreAnalysis_LFQ_03_PPI_EnrichedProteins.xlsx: Here you can find all the enriched proteins resulting from the comparison of a WT strain (YAF000) to an RNase treated YAFXXX strain. These are the proteins inside the PPI groups.

- 01_CoreAnalysis_LFQ_04_RNAm_EnrichedProteins.xlsx: Here you can find all the enriched proteins resulting from the comparison of a non-treated YAFXXX strain to an RNase treated YAFXXX strain. These are the proteins inside the RDI groups.

## Functional Analysis Files

Finally we have the functional analysis for the PPI and RDI (named RNAm in this repository) groups. Thus, for each group, we have the following files:

- 02_FunctionalAnalysis_LFQ_YYY_01_Summary.pdf: Here you can find a summary of all the figures generated during the GO (Biological process, cellular component and molecular function), Reactome and KEGG analysis for each YAFXXX strain.

- 02_FunctionalAnalysis_LFQ_YYY_02_Interpro.xlsx: Here you can find all the Interpro protein signatures associated to each YAFXXX strain and their associated statistical values.

- 02_FunctionalAnalysis_LFQ_YYY_03_UpRegulated_GO_MF.xlsx: Here you can find all the overrepresented molecular function GO terms associated to each YAFXXX strain and their associated statistical values.

- 02_FunctionalAnalysis_LFQ_YYY_04_UpRegulated_KEGG.xlsx: Here you can find all the overrepresented KEGG terms associated to each YAFXXX strain and their associated statistical values.

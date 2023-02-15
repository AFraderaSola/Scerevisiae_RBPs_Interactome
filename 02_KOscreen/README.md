# 02_KOscreen folder

This folder contains all the scripts used to analyze the MS data for the **KO screen** described in Figure 4A. Within this folder you can find the data and scripts to generate the following figures:

- Figure 4B, 4C, 4D and 4E

It also contains the scripts to generate Supplementary tables 9, 11 and 12. 

## All

This **sub-folder** contains a downstream analysis **including all strains**. The **input files** needed for those analysis are generated on each KO-RBP sub-folder (**YKAFXXX**). The scripts to generate the previously mentioned figures are inside this folder. More details on the scripts and how to run them inside this folder.

## YKAFXXX

This folder contains an example on how the main screen analysis was structured. Each **KO-RBP** had its own **sub-folder**, named after its stock strain (**YKAFXXX**), containing its **Core LFQ analysis** and its **Functional Analysis**.

To run the Core LFQ analysis, you **need the proteinGroups.txt** (obtained from analyzing .RAW MS files with [MaxQuant](https://www.maxquant.org/)). The file should be placed at the **same folder** the analysis script is (example inside the folder). The proteinGroups file for each KO-RBP can be found on the 00_KOscreeningMasterproteinGroups.zip file. 

Additionally you can provide a **highlight.txt** and **volcano_contrasts.txt**, but this won't be necessary for the script to run. Examples of such files are provided inside the folder. 

The scripts for the **Core LFQ Analysis** and the **Functional Analysis** work in a sequential fashion; you need the results of the first to run the second.

More details on the scripts and how to run them  can be found inside each folder.

## Master Files 

First we have the MasterFiles, containing general information on the strains:

- 00_KOScreeningMasterFile.xlsx: Here you can find all the necessary information for each YKAFXXX strains (i.e., gene name or selection criteria) and a summary of it data analysis (i. e., quantified proteins or enriched proteins)

- 00_KOScreeningMSMasterFile.xlsx: Here you can find the related MS information (i.e, measured date or raw file) for each YKAFXXX strain.

- 00_KOScreeningMasterproteinGroups.zip: Here you can find all the protein groups for each KO-RBP.

## Core Analysis Files

Then we have the CoreAnalysis files; this files summarize the main analysis conducted on each YKAFXXX strain. More details on the analysis itself can be found on each YKAFXXX folder.

- 01_CoreAnalysis_LFQ_01_Summary.pdf: Here you can find a summary the CoreAnalysis figures for each YKAFXXX strain. 

- 01_CoreAnalysis_LFQ_QuantifiedProteins.xlsx: Here you can find all quantified proteins for each YKAFXXX strain and their associated statistical values.

- 01_CoreAnalysis_LFQ_EnrichedProteins.xlsx: Here you can find all the enriched proteins resulting from the comparison of a WT strain (YAF000) to an RNase treated YKAFXXX strain. These are the proteins inside the PPI groups.

## Functional Analysis Files

Finally we have the functional analysis for the PPI and RDI (named RNAm in this repository) groups. Thus, for each group, we have the following files:

- 02_FunctionalAnalysis_LFQ_KO_UpRegulated_KEGG.xlsx: Here you can find all the overrepresented KEGG terms associated to each YKAFXXX strain and their associated statistical values.

- 02_FunctionalAnalysis_LFQ_KO_DownRegulated_KEGG.xlsx: Here you can find all the overrepresented KEGG terms associated to each YKAFXXX strain and their associated statistical values.

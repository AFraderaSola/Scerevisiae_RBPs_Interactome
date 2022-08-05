# PPI functional enrichment analysis

Script to perform a functional analysis, including GO, Reactome and KEGG databases, for each RBP bait.

## Input files

It requires the following input files for each RBP-bait:

- **_YAFXXX_QuantifiedProteins.csv_**:  Obtained during the **CoreAnalysis**, it contains all the quantified proteins and its associated statistical values. We need the **Protein.IDs column** which we will use as our universe (background) during the enrichment analysis. It can be found, for each RBP-bait, on the different sheets of the 01_CoreAnalysis_LFQ_02_QuantifiedProteins.xlsx file. 

- **_YAFXXX_Enriched_RNase_vs_WT.csv_**:  Obtained during the **CoreAnalysis**, it contains all the enriched proteins and its associated statistical values. We need the **Protein.IDs column** which we will compare to the universe (background) during the enrichment analysis. It can be found, for each RBP-bait, on the different sheets of the 01_CoreAnalysis_LFQ_03_PPI_EnrichedProteins.xlsx file.

Additionally, it requires the 00_IPScreeningMasterFile.xlsx (01_IPscreen folder) file to match the stock ID, YAFXXX, to its systematic and gene name. 

## LFQ_PPI_Enrichment.R Script

The script performs a functional enrichment analysis based on the [ClusterProfiler R package](https://yulab-smu.top/biomedical-knowledge-mining-book/index.html). Its goal is to find overrepresented [GO](http://geneontology.org/), [Reactome](https://reactome.org/) and [KEGG](https://www.genome.jp/kegg/) terms among the enriched PPI interactors when compared to the terms found on the quantified proteins. Through the script, we cover the following basic steps:

- **_Data wrangling_**: We retrieve the universe IDs (quantified proteins) and the target IDs (enriched proteins) from the input files and convert them to suitable identifiers. For GO terms we use the gene name (i. e., PAB1) while for Reactome and KEGG terms we use the entrezid (i. e., 856912), which we obtain from the systematic ensembl ID (i. e., YER165W).

- **_Over-representation analysis_**: The [ClusterProfiler](https://yulab-smu.top/biomedical-knowledge-mining-book/enrichment-overview.html) over-representation analysis is implemented which, in brief, calculates p-values by hypergeometric distribution from a experimentally derived list of enriched IDs (PPI set) and its background (Quantified proteins).

## Output files

If the script runs properly, it will generate a series of output files that can roughly categorized into the following categories:

- **_Tables_**: A table containing all the over-represented terms for each database (GO, Reactome and KEGG) will be generated, if any. 

Supplementary tables 4 and 7 from the **"RNA-dependent interactome allows network-based assignment of RBP function"** publication were generated from the output tables of this script. A combined table for each RBP bait GO molecular function and KEGG terms can be found at 02_FunctionalAnalysis_LFQ_PPI_03_UpRegulated_GO_MF.xlsx and  02_FunctionalAnalysis_LFQ_PPI_04_UpRegulated_KEGG.xlsx files, respectively (01_IPscreen folder).

- **_Plots_**: For each database (GO, Reactome and KEGG) the over-represented terms, if any, will be visualized as:
  - Barplot, depicting the top30 (by adjusted p-value) over-represented terms and the number of genes associated to it. Coloured by adjusted p-value
  - Dotplot, depicting the top30 (by adjusted p-value) over-represented terms and their gene ratio. Coloured by adjusted p-value. Size by number of genes per term.
  - Cnetplot, depicting the linkages of IDs and their associated terms as a network. Coloured by IDs's fold-change.
  - Heatmap, depicting a tile for each ID found on each over-represented in a particular term. Coloured by IDs's fold-change.
  - Enrichmentmap, depicting enriched terms into a network with edges connecting overlapping IDs. Coloured by IDs's fold-change.
  - GO diagram, depicting each term and its parental associated terms (only for GO). Coloured by adjusted p-value.
  - Upset plot, depicting the overlap between IDs and their associted over-represented terms. 

The script output's **plots for all RBP-baits** are included in the 02_FunctionalAnalysis_LFQ_PPI_01_Summary.pdf file, in the 01_IPscreen folder.
  
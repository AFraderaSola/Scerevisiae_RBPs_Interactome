################################
######### Libraries ############
################################

library(readr)
library(dplyr)
library(tidyr)
library(plyr)
library(grid)
library(gridExtra)
library(UpSetR)
library(RColorBrewer)
library(ggplot2)
library(pheatmap)
library(viridis)
library(readxl)
library(seqinr)
library(stringr)
library(dplyr)
library(tidyr)

################################
########## Variables ###########
################################

# Quantified proteins on YAF000

files_screen <- list.files(path = "./00_IPsResultsFiles/", pattern = "WT*.csv|RNase*.csv")

screen <- read_excel("../../../00_IPScreeningMasterFile.xlsx")

files_quantified <- list.files(path = "./00_IPsResultsFiles/", pattern = "Proteins*.csv")

################################
########## Script ##############
################################

protein_list <- c()

genes <- c()

for (i in 1:length(files_screen)) {
  
  proteins_loop <- read.csv(paste0("./00_IPsResultsFiles/",files_screen[i]))
  
  proteins_loop <- 
    proteins_loop %>%
    as_tibble() %>%
    # filter_at(vars(starts_with('enriched_')), any_vars(.)) %>%
    # glimpse() %>%
    select(Protein.IDs,
           starts_with("difference"))
  
  stock_ID <-  gsub(pattern = "_Enriched.*",replacement = "",files_screen[i])
  
  geneID <- screen[screen$Stock_Name == stock_ID,]$Gene_Name
  
  genes <- c(genes, geneID)
  
  colnames(proteins_loop)[2] <- paste0(geneID, "_", colnames(proteins_loop)[2])
  
  proteins_loop <- proteins_loop[proteins_loop[,2] > 0,]
  
  colnames(proteins_loop)[1] <- "ProteinID"
  
  colnames(proteins_loop)[2] <- gsub(pattern = "IP_RNase", replacement = "RNAmediated",x = colnames(proteins_loop[2]))
  
  colnames(proteins_loop)[2] <- gsub(pattern = "RNase_WT", replacement = "PPI",x = colnames(proteins_loop[2]))
  
  proteins_loop_list <- list(proteins_loop$ProteinID)
  
  names(proteins_loop_list) <- colnames(proteins_loop)[2]
  
  protein_list <- c(protein_list, proteins_loop_list)
  
}

## PPI and RNAmediated combined

interactors <- unique(unlist(protein_list,use.names = F))

interactors <- gsub(pattern = "-TAP",replacement = "",x = gsub(pattern = ";.*",replacement = "",x = interactors))

yeast_proteins <- read.fasta("./03_InputFiles/000000_S288C_Genome_Release_64-2-1_orf_trans_all.fasta",seqtype = "AA")

yeast_proteins_filtered <- yeast_proteins[names(yeast_proteins)%in%interactors]

yeast_proteins_cleaned <- yeast_proteins_filtered

annot <- c()

for (i in 1:length(yeast_proteins_cleaned)) {
  
  yeast_proteins_cleaned[[i]] <- yeast_proteins_filtered[[i]][!str_detect(yeast_proteins_filtered[[i]],pattern="\\*")]
  
  attr(yeast_proteins_cleaned[[i]], "Annot") <- attr(yeast_proteins_filtered[[i]],"Annot")
  
  attr(yeast_proteins_cleaned[[i]],"name") <- attr(yeast_proteins_filtered[[i]],"name")
  
  attr(yeast_proteins_cleaned[[i]],"class") <- attr(yeast_proteins_filtered[[i]],"class")
  
  annot <- c(annot, substring(attr(yeast_proteins_cleaned[[i]], "Annot"),2))
  
}


file_name <- paste0("./03_OutputFiles/", gsub(Sys.Date(),pattern = "-", replacement = ""),"_AllInteractors.fasta")

write.fasta(sequences = yeast_proteins_filtered, names = annot, file.out = file_name)

file_name <- paste0("./03_OutputFiles/", gsub(Sys.Date(),pattern = "-", replacement = ""),"_AllInteractors_Cleaned.fasta")

write.fasta(sequences = yeast_proteins_cleaned, names = annot, file.out = file_name)

## PPI

# interactors <- unique(unlist(protein_list[grep(pattern = "PPI",x = names(protein_list))],use.names = F))
# 
# interactors <- gsub(pattern = "-TAP",replacement = "",x = gsub(pattern = ";.*",replacement = "",x = interactors))
# 
# yeast_proteins <- read.fasta("./03_InputFiles/000000_S288C_Genome_Release_64-2-1_orf_trans_all.fasta",seqtype = "AA")
# 
# yeast_proteins_filtered <- yeast_proteins[names(yeast_proteins)%in%interactors]
# 
# yeast_proteins_cleaned <- yeast_proteins_filtered
# 
# annot <- c()
# 
# for (i in 1:length(yeast_proteins_cleaned)) {
#   
#   yeast_proteins_cleaned[[i]] <- yeast_proteins_filtered[[i]][!str_detect(yeast_proteins_filtered[[i]],pattern="\\*")]
#   
#   attr(yeast_proteins_cleaned[[i]], "Annot") <- attr(yeast_proteins_filtered[[i]],"Annot")
#   
#   attr(yeast_proteins_cleaned[[i]],"name") <- attr(yeast_proteins_filtered[[i]],"name")
#   
#   attr(yeast_proteins_cleaned[[i]],"class") <- attr(yeast_proteins_filtered[[i]],"class")
#   
#   annot <- c(annot, substring(attr(yeast_proteins_cleaned[[i]], "Annot"),2))
#   
# }
# 
# 
# file_name <- paste0("./03_OutputFiles/", gsub(Sys.Date(),pattern = "-", replacement = ""),"_PPIInteractors.fasta")
# 
# write.fasta(sequences = yeast_proteins_filtered, names = annot, file.out = file_name)
# 
# file_name <- paste0("./03_OutputFiles/", gsub(Sys.Date(),pattern = "-", replacement = ""),"_PPIInteractors_Cleaned.fasta")
# 
# write.fasta(sequences = yeast_proteins_cleaned, names = annot, file.out = file_name)

## RNAmediated

interactors <- unique(unlist(protein_list[grep(pattern = "RNAmediated",x = names(protein_list))],use.names = F))

interactors <- gsub(pattern = "-TAP",replacement = "",x = gsub(pattern = ";.*",replacement = "",x = interactors))

yeast_proteins <- read.fasta("./03_InputFiles/000000_S288C_Genome_Release_64-2-1_orf_trans_all.fasta",seqtype = "AA")

yeast_proteins_filtered <- yeast_proteins[names(yeast_proteins)%in%interactors]

yeast_proteins_cleaned <- yeast_proteins_filtered

annot <- c()

for (i in 1:length(yeast_proteins_cleaned)) {

  yeast_proteins_cleaned[[i]] <- yeast_proteins_filtered[[i]][!str_detect(yeast_proteins_filtered[[i]],pattern="\\*")]

  attr(yeast_proteins_cleaned[[i]], "Annot") <- attr(yeast_proteins_filtered[[i]],"Annot")

  attr(yeast_proteins_cleaned[[i]],"name") <- attr(yeast_proteins_filtered[[i]],"name")

  attr(yeast_proteins_cleaned[[i]],"class") <- attr(yeast_proteins_filtered[[i]],"class")

  annot <- c(annot, substring(attr(yeast_proteins_cleaned[[i]], "Annot"),2))

}


file_name <- paste0("./03_OutputFiles/", gsub(Sys.Date(),pattern = "-", replacement = ""),"_RNAmediatedInteractors.fasta")

write.fasta(sequences = yeast_proteins_filtered, names = annot, file.out = file_name)

file_name <- paste0("./03_OutputFiles/", gsub(Sys.Date(),pattern = "-", replacement = ""),"_RNAmediatedInteractors_Cleaned.fasta")

write.fasta(sequences = yeast_proteins_cleaned, names = annot, file.out = file_name)

## Quantified Proteins

protein_list <- c()

genes <- c()

for (i in 1:length(files_quantified)) {
  
  proteins_loop <- read.csv(paste0("./00_IPsResultsFiles/",files_quantified[i]))
  
  proteins_loop <- 
    proteins_loop %>%
    as_tibble() %>%
    # filter_at(vars(starts_with('enriched_')), any_vars(.)) %>%
    # glimpse() %>%
    select(Protein.IDs)
  
  stock_ID <-  gsub(pattern = "_Quan.*",replacement = "",files_quantified[i])
  
  geneID <- screen[screen$Stock_Name == stock_ID,]$Gene_Name
  
  genes <- c(genes, geneID)
  
  colnames(proteins_loop)[1] <- "ProteinID"
  
  proteins_loop_list <- list(proteins_loop$ProteinID)
  
  names(proteins_loop_list) <- paste0(stock_ID, "_QuantifiedProteins")
  
  protein_list <- c(protein_list, proteins_loop_list)
  
}

interactors <- unique(unlist(protein_list,use.names = F))

interactors <- gsub(pattern = "-TAP",replacement = "",x = gsub(pattern = ";.*",replacement = "",x = interactors))

yeast_proteins <- read.fasta("./03_InputFiles/000000_S288C_Genome_Release_64-2-1_orf_trans_all.fasta",seqtype = "AA")

yeast_proteins_filtered <- yeast_proteins[names(yeast_proteins)%in%interactors]

yeast_proteins_cleaned <- yeast_proteins_filtered

annot <- c()

for (i in 1:length(yeast_proteins_cleaned)) {
  
  yeast_proteins_cleaned[[i]] <- yeast_proteins_filtered[[i]][!str_detect(yeast_proteins_filtered[[i]],pattern="\\*")]
  
  attr(yeast_proteins_cleaned[[i]], "Annot") <- attr(yeast_proteins_filtered[[i]],"Annot")
  
  attr(yeast_proteins_cleaned[[i]],"name") <- attr(yeast_proteins_filtered[[i]],"name")
  
  attr(yeast_proteins_cleaned[[i]],"class") <- attr(yeast_proteins_filtered[[i]],"class")
  
  annot <- c(annot, substring(attr(yeast_proteins_cleaned[[i]], "Annot"),2))
  
}


file_name <- paste0("./03_OutputFiles/", gsub(Sys.Date(),pattern = "-", replacement = ""),"_QuantifiedProteins.fasta")

write.fasta(sequences = yeast_proteins_filtered, names = annot, file.out = file_name)

file_name <- paste0("./03_OutputFiles/", gsub(Sys.Date(),pattern = "-", replacement = ""),"_QuantifiedProteins_Cleaned.fasta")

write.fasta(sequences = yeast_proteins_cleaned, names = annot, file.out = file_name)

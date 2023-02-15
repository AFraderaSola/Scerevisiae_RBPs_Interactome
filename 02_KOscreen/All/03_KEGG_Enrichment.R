###################################
######### Libraries ###############
###################################

library(pheatmap)
library(RColorBrewer)
library(reshape2)
library(ggplot2)
library(viridis)
library(ggpubr)
library(readr)
library(ggplot2)
library(dplyr)
library(readxl)
library(scales)
library(ggthemes)
library(tidyverse)


###################################
######### Functions ###############
###################################

splitColumnBySep <- function(df, column, sep=";"){
  nameCol <- column
  names(df)[grep(column, names(df))] <- "column"
  df2 <- df %>%
    mutate(column = strsplit(column, sep)) %>%
    unnest(column)
  names(df2)[grep("column", names(df2))] <- nameCol
  return(df2)
}

###################################
######### Variables ###############
###################################

# Define files to compare

enrichment_pattern <- "_KEGG_"

enrichment_analysis <- "KEGG"

# Define heatmap style (Count, p.adjust, GeneRatio...)

style <- "GeneRatio"

################################
########## Script ##############
################################

## Up

# Define your experiment Name:

experiment <- "Up_KO"

names <- read_excel("../00_KOScreeningMasterFile.xlsx")

files_pattern <- paste0(".*",enrichment_pattern, ".*", experiment)

files_enrichment <- list.files(path = "./03_InputFiles/", pattern = files_pattern)

KEGG_Term <- c()

p.adjust <- c()

counts <- c()

df <- c()

selection <- c()

KEGG_table <- c()

for (i in 1:length(files_enrichment)) {
  
  KEGG_Term_loop <- read.csv(paste0("./03_InputFiles/", files_enrichment[i]))
  
  KEGG_Term_loop <- KEGG_Term_loop$Description
  
  
  KEGG_Term_loop <- list(KEGG_Term_loop)
  
  names(KEGG_Term_loop) <- gsub(pattern = paste0(enrichment_pattern, ".*.csv"),
                                replacement = "",files_enrichment[i])
  
  KEGG_Term <- c(KEGG_Term, KEGG_Term_loop)
  
  selection_loop <- read.csv(paste0("./03_InputFiles/", files_enrichment[i]))
  
  selection_loop <- selection_loop[,grep(pattern = style, colnames(selection_loop))]
  
  ratio <- selection_loop
  
  if (style == "GeneRatio") {
    
    ratio <- c()
    
    for (j in 1:length(selection_loop)) {
      
      loop_ratio <- selection_loop[j]
      
      loop_ratio <- as.numeric(unlist(strsplit(loop_ratio,split = "/")))
      
      loop_ratio <- loop_ratio[1] / loop_ratio [2]
      
      ratio <- c(ratio, loop_ratio)
      
    }
    
  }
  
  selection_loop <- list(ratio)
  
  names(selection_loop) <- gsub(pattern = paste0(enrichment_pattern, ".*.csv"),
                                replacement = "",files_enrichment[i])
  
  selection <- c(selection, selection_loop)
  
  df_loop <- data.frame(KEGG_Term_loop, selection_loop)
  
  colnames(df_loop) <- c(enrichment_analysis, paste0(gsub(pattern = paste0(enrichment_pattern, ".*.csv"),
                                                          replacement = "",files_enrichment[i]),
                                                     "_", style))
  
  KEGG_table_loop <- df_loop
  
  KEGG_table_loop$Stock <- gsub(pattern = ".*_",
                                replacement = "",
                                x = gsub(pattern = "_G.*",
                                         replacement = "",
                                         x = colnames(KEGG_table_loop)[2]))
  
  colnames(KEGG_table_loop)[2] <- "GeneRatio"
  
  KEGG_table_loop <- KEGG_table_loop[,c(3,1,2)]
  
  KEGG_table <- rbind(KEGG_table, KEGG_table_loop)
  
  df_loop <- list(df_loop)
  
  names(df_loop) <- gsub(pattern = paste0(enrichment_pattern, ".*.csv"),
                         replacement = "",files_enrichment[i])
  
  df <- c(df, df_loop)
  
}

## For the network plot

# KEGG_table$Gene <- KEGG_table$Stock
# 
# KEGG_table <- KEGG_table[,c(1,4,2,3)]
# 
# for (i in 1:nrow(KEGG_table)) {
# 
#   stock <- KEGG_table$Stock[i]
# 
#   gene <- names[names$Stock_Name == stock,]$Gene_Name
# 
#   KEGG_table$Gene[i] <- gene
# 
# }
# 
# write.csv(x = KEGG_table, file = "./02_OutputFiles/00000000_PPI_KEGG_table.csv",row.names = F)

## 

df_HM <- df[[1]]

for (i in 2:length(df)) {
  
  df_HM <- full_join(x = df_HM,y = df[[i]], by = enrichment_analysis)
  
}

alt_colnames <- c("KEGG")

for (i in 2:length(colnames(df_HM))) {
  
  stock <- as.character(colnames(df_HM)[i])
  
  stock <- gsub(pattern = ".*_",
                replacement = "",
                x = gsub(pattern = paste0("_",style),
                         replacement = "",
                         x = stock))
  
  gene <- names[names$Stock_Name == stock,]$Gene_Name
  
  gene <- paste0(gene,"_PPI")
  
  alt_colnames <- c(alt_colnames, gene)
  
}

colnames(df_HM) <- alt_colnames

df_HM_up <- df_HM

## RNAmediated

# Define your experiment Name:

experiment <- "Down_KO"

files_pattern <- paste0(".*",enrichment_pattern, ".*", experiment)

files_enrichment <- list.files(path = "./03_InputFiles/", pattern = files_pattern)

KEGG_Term <- c()

p.adjust <- c()

counts <- c()

df <- c()

selection <- c()

KEGG_table <- c()

for (i in 1:length(files_enrichment)) {
  
  KEGG_Term_loop <- read.csv(paste0("./03_InputFiles/", files_enrichment[i]))
  
  KEGG_Term_loop <- KEGG_Term_loop$Description
  
  
  KEGG_Term_loop <- list(KEGG_Term_loop)
  
  names(KEGG_Term_loop) <- gsub(pattern = paste0(enrichment_pattern, ".*.csv"),
                                replacement = "",files_enrichment[i])
  
  KEGG_Term <- c(KEGG_Term, KEGG_Term_loop)
  
  selection_loop <- read.csv(paste0("./03_InputFiles/", files_enrichment[i]))
  
  selection_loop <- selection_loop[,grep(pattern = style, colnames(selection_loop))]
  
  ratio <- selection_loop
  
  if (style == "GeneRatio") {
    
    ratio <- c()
    
    for (j in 1:length(selection_loop)) {
      
      loop_ratio <- selection_loop[j]
      
      loop_ratio <- as.numeric(unlist(strsplit(loop_ratio,split = "/")))
      
      loop_ratio <- loop_ratio[1] / loop_ratio [2]
      
      ratio <- c(ratio, loop_ratio)
      
    }
    
  }
  
  selection_loop <- list(ratio)
  
  names(selection_loop) <- gsub(pattern = paste0(enrichment_pattern, ".*.csv"),
                                replacement = "",files_enrichment[i])
  
  selection <- c(selection, selection_loop)
  
  df_loop <- data.frame(KEGG_Term_loop, selection_loop)
  
  colnames(df_loop) <- c(enrichment_analysis, paste0(gsub(pattern = paste0(enrichment_pattern, ".*.csv"),
                                                          replacement = "",files_enrichment[i]),
                                                     "_", style))
  
  KEGG_table_loop <- df_loop
  
  KEGG_table_loop$Stock <- gsub(pattern = ".*_",
                                replacement = "",
                                x = gsub(pattern = "_G.*",
                                         replacement = "",
                                         x = colnames(KEGG_table_loop)[2]))
  
  colnames(KEGG_table_loop)[2] <- "GeneRatio"
  
  KEGG_table_loop <- KEGG_table_loop[,c(3,1,2)]
  
  KEGG_table <- rbind(KEGG_table, KEGG_table_loop)
  
  df_loop <- list(df_loop)
  
  names(df_loop) <- gsub(pattern = paste0(enrichment_pattern, ".*.csv"),
                         replacement = "",files_enrichment[i])
  
  df <- c(df, df_loop)
  
}

## For the network plot

# KEGG_table$Gene <- KEGG_table$Stock
# 
# KEGG_table <- KEGG_table[,c(1,4,2,3)]
# 
# for (i in 1:nrow(KEGG_table)) {
#   
#   stock <- KEGG_table$Stock[i]
#   
#   gene <- names[names$Stock_Name == stock,]$Gene_Name
#   
#   KEGG_table$Gene[i] <- gene
#   
# }
# 
# write.csv(x = KEGG_table, file = "./02_OutputFiles/00000000_RNAm_KEGG_table.csv",row.names = F)

## 

df_HM <- df[[1]]

for (i in 2:length(df)) {
  
  df_HM <- full_join(x = df_HM,y = df[[i]], by = enrichment_analysis)
  
}

alt_colnames <- c("KEGG")

for (i in 2:length(colnames(df_HM))) {
  
  stock <- as.character(colnames(df_HM)[i])
  
  stock <- gsub(pattern = ".*_",
                replacement = "",
                x = gsub(pattern = paste0("_",style),
                         replacement = "",
                         x = stock))
  
  gene <- names[names$Stock_Name == stock,]$Gene_Name
  
  gene <- paste0(gene,"_RNAmediated")
  
  alt_colnames <- c(alt_colnames, gene)
  
}

colnames(df_HM) <- alt_colnames

df_HM_down <- df_HM


## Join both dataframes

df_HM <- full_join(x = df_HM_up,y = df_HM_down, by = enrichment_analysis)

df_HM <- as.data.frame(df_HM)

rownames(df_HM) <- df_HM[,1]

df_HM <- df_HM[,c(2:ncol(df_HM))]

df_HM <- as.matrix(df_HM)

colours <- c(viridis(n = 100, option = "G",direction = -1))

Enrichment <- c(rep("Up-regulated", ncol(df_HM_up[,-1])),rep("Down-regulated", ncol(df_HM_down[,-1])))

# bfc <- read.csv(file = "./03_InputFiles/00000000_BaitFunctionCriteria.csv")

# Bait <- c(bfc[bfc$Target %in% gsub(pattern = "_.*",
#                                    replacement = "",
#                                    x = colnames(df_HM)[grep(pattern = "PPI",x = colnames(df_HM))]),]$Bait_Function[order(
#                                      match(bfc[bfc$Target %in% gsub(pattern = "_.*",
#                                                                     replacement = "",
#                                                                     x = colnames(df_HM)[grep(pattern = "PPI",x = colnames(df_HM))]),]$Target,
#                                            gsub(pattern = "_.*",
#                                                 replacement = "",
#                                                 x = colnames(df_HM)[grep(pattern = "PPI",x = colnames(df_HM))])))],
#           bfc[bfc$Target %in% gsub(pattern = "_.*",
#                                    replacement = "",
#                                    x = colnames(df_HM)[grep(pattern = "RNAm",x = colnames(df_HM))]),]$Bait_Function[order(
#                                      match(bfc[bfc$Target %in% gsub(pattern = "_.*",
#                                                                     replacement = "",
#                                                                     x = colnames(df_HM)[grep(pattern = "RNAm",x = colnames(df_HM))]),]$Target,
#                                            gsub(pattern = "_.*",
#                                                 replacement = "",
#                                                 x = colnames(df_HM)[grep(pattern = "RNAm",x = colnames(df_HM))])))])

df_col <- data.frame(Enrichment)

rownames(df_col) <- colnames(df_HM)

kgf <- read.csv(file = "./03_InputFiles/00000000_KEGGGeneralFunction.csv")

Pathway <- kgf$General_Function

df_row <- data.frame(Pathway)

rownames(df_row) <- rownames(df_HM)

###
#For the network plot

# color <- brewer.pal(n = 8, name = "Set2")
# 
# KEGG <- df_row
# 
# KEGG$KEGG <- row.names(KEGG)
# 
# colnames(KEGG)[1] <- "Functionality"
# 
# KEGG <- KEGG[order(KEGG$Functionality),]
# 
# count <- KEGG %>% dplyr::count(Functionality)
# 
# KEGG$color <- c(rep(color[5],count$n[1]), #DegradationColour
#                 rep("#666666",count$n[2]), #ExportColour
#                 rep(color[7],count$n[3]), #MetabolismColour
#                 rep(color[1],count$n[4]), #RibosomeColour
#                 rep(color[4],count$n[5]), #SplicingColour
#                 rep(color[2],count$n[6])) #SynthesisColour
#                 # rep(color[6],count$n[6])) #TransportColour
# 
# write.csv(x = KEGG, file = "./02_OutputFiles/00000000_NetworkKEGG.csv",row.names = F)

####

Enrichment <- c("#e32026", "#1a79b6")
names(Enrichment) <- c("Up-regulated", "Down-regulated")

# color <- brewer.pal(n = 8, name = "Set2")
# 
# Bait <- c(color[3],color[4],color[8],
#           "#666666",
#           color[5])
# names(Bait) <- c("Capping","Splicing", "Cleavage", 
#                  "Export",
#                  "Degradation")
# 

Pathway <- c("#a6cee3",
             "#fc8d62",
             "#e5c494",
             "#b4b4b4"
             )

names(Pathway) <- c( "RNA pathway",
                     "Synthesis pathway",
                     "Metabolic pathway",
                     "Other"
                     )

anno_colors <- list(Enrichment = Enrichment,
                    Pathway = Pathway)

# df_HM <- df_HM[,c(10,1:9,25:26,20:22,11:19,23:24,37,29:36,49:51,44:45,28,38:43,46:48)]

df_HM <- df_HM[c(1,3,24,9,12,22,23,4,6:8,10,11,13,15,16,20,21,25,2,5,14,17:19,26),]

# Exclude not validated strains

df_HM <- df_HM[c(-23,-24,-25),c(-8,-9)]

filename <- paste0("./03_OutputFiles/",gsub(pattern = "-",x = Sys.Date(),replacement = ""),"_KEGG_Enrichment_Both.pdf")

gaps_col <- ncol(df_HM_up[,-1])

pheatmap(df_HM,
         col=colours,
         cluster_cols = F,
         cluster_rows = F,
         show_rownames = T,
         show_colnames = T,
         labels_col = paste0(str_to_lower(gsub(pattern = "_.*",replacement = "",x = colnames(df_HM))),"-d"),
         annotation_col = df_col,
         annotation_row = df_row,
         annotation_colors = anno_colors,
         fontsize = 16,border_color = "#b4b4b4",gaps_col = gaps_col,
         height = 11.4,width = 12,filename = filename
         )


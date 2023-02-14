#########################################
############ Libraries ##################
#########################################

library(tidyverse)
library(DOSE)
library(enrichplot)
library(ReactomePA)
library(igraph)
library(ggplot2)
library(ggparty)
library(ggraph)
library(ggnewscale)
library(plotly)
library(clusterProfiler)
library(RColorBrewer)
library(readxl)
library(viridis)
library(GGally)

options(ggrepel.max.overlaps = Inf)

set.seed(666)

#########################################
############ Variables ##################
#########################################

# Files to plot

file_network <- list.files(path = "./03_InputFiles/", pattern = "*.Table.csv")

# Plot variables

circular <- F

gene_names <- F

# Define the organism database annotation package for GO:

org <- "org.Sc.sgd.db"

layout <- "fr"

# KEGG Highlight

KEGG <- read.csv("./03_InputFiles/00000000_NetworkKEGG.csv")

IP_KEGG <- read.csv("./03_InputFiles/00000000_Both_KEGG_table.csv")

Functionality <- "Splicing"

# RBP Highliht

RBP <- read.csv("./03_InputFiles/00000000_SC_RBPs_NamedCensus.csv")

#########################################
############## Script ###################
#########################################

systematic_pattern <- c("YAL|YAR|YBL|YBR|YCL|YCR|YDL|YDR|YEL|YER|YFL|YFR|YGL|YGR|YHL|YHR|YIL|YIR|YJL|YJR|YKL|YKR|YLL|YLR|YML|YMR|YNL|YNR|YOL|YOR|YPL|YPR")

screen <- read.csv(paste0("./03_InputFiles/",file_network))

targets <- unique(screen$Target)

category <- c("Overlap","PPI","RNAmediated")

colors_interaction <- c("#7570B3","#66A61E", "#D95F02")

screen_list <- c()

category_df <- c()

for (i in 1:length(targets)) {
  
  loop_screen <- screen[screen$Target == targets[i],]
  
  for (j in 1:length(category)) {
    
    genes <- unique(loop_screen[,j])
    
    genes <- genes[!genes == ""]
    
    genes <- c(str_to_title(genes[!grepl(genes,pattern = systematic_pattern)]),genes[grepl(genes,pattern = systematic_pattern)])
    
    category_loop <- rep(category[j], length(genes))
    
    target_RBP <- rep(targets[i], length(genes))
    
    target_RBP <- str_to_title(target_RBP)
    
    colors_interaction_loop <- rep(colors_interaction[j], length(genes))
    
    category_df_loop <- data.frame(target_RBP, genes, category_loop, colors_interaction_loop)
    
    category_df <- rbind(category_df, category_df_loop)
  
  }
  
  IDs <- unique(c(loop_screen$overlap, loop_screen$PPI, loop_screen$RNAmediated))
  
  IDs <- IDs[!IDs == ""]
  
  IDs <- c(str_to_title(IDs[!grepl(IDs,pattern = systematic_pattern)]),IDs[grepl(IDs,pattern = systematic_pattern)])
  
  IDs <- list(IDs)
  
  names(IDs) <- str_to_title(targets[i])
  
  screen_list <- c(screen_list, IDs)
  
  
}

edges <- category_df[,c(1:3)]

colnames(edges) <- c("source", "target", "set")

nodes <- unique(c(edges$source, edges$target))

size <- c()

for (i in 1:length(nodes)) {
  
  if (nodes[i] %in% names(screen_list)) {
    
    loop_size <- length(unlist(screen_list[names(screen_list) == nodes[i]]))
    
  }else{
    
    loop_size <- 1
  }
  
  size <- c(size, loop_size)
}

vertices <- data.frame(nodes = nodes,
                       size = size)

vertices$Functionality <- "None"

vertices$ID <- "Prey"

vertices[vertices$nodes %in% str_to_title(targets),]$ID <- "Bait"

filter <- KEGG[KEGG$Functionality == Functionality,]$KEGG

filter <- str_to_title(IP_KEGG[IP_KEGG$KEGG %in% filter,]$Gene)

edges <- edges[edges$source %in% filter,]

filter <- unique(c(edges$source, edges$target))

vertices <- vertices[vertices$nodes %in% filter,]

filter <- KEGG[KEGG$Functionality == Functionality,]$KEGG

filter <- str_to_title(IP_KEGG[IP_KEGG$KEGG %in% filter,]$Gene)

filter2 <- as.data.frame(sort(table(c(edges$source, edges$target)), decreasing = TRUE))

filter2 <- filter2[filter2$Freq > 5,]

filter2 <- filter2$Var1

vertices[!vertices$nodes %in% filter, ]$ID <- "Prey"

vertices[vertices$ID == "Prey",]$size <- 1

vertices$nodes <- gsub(pattern = ",",replacement = "",x = vertices$nodes)

edges$target <- gsub(pattern = ",",replacement = "",x = edges$target)

# edges$set <- gsub(pattern = "Overlap",replacement = "Both",x = edges$set)
# 
# edges$set <- gsub(pattern = "RNAmediated",replacement = "RDI",x = edges$set)

write.csv(x = edges, file = paste0("./03_OutputFiles/",gsub(pattern = "-",replacement = "",x = Sys.Date()),"_SplicingEdges.csv"),quote = F,row.names = F)

write.csv(x = vertices, file = paste0("./03_OutputFiles/",gsub(pattern = "-",replacement = "",x = Sys.Date()),"_SplicingNodes.csv"),quote = F,row.names = F)

net <- graph_from_data_frame(d = edges, vertices = vertices)

E(net)$edge.color <- ifelse(E(net)$set == "Overlap", "#666666", E(net)$set)
E(net)$edge.color <- ifelse(E(net)$edge.color == "PPI","#66A61E", E(net)$edge.color)
E(net)$edge.color <- ifelse(E(net)$edge.color == "RNAmediated", "#D95F02", E(net)$edge.color)

sort(table(c(edges$source, edges$target)), decreasing = TRUE)

plot <- ggnet2(net, color = "ID", label = filter2, size = "size", edge.color = "edge.color", mode = "fruchtermanreingold",
               label.size = 8,edge.size = 0.5, 
               palette = c(
                 "Prey" = "#666666",
                 "Bait" = "#E78AC3"))+
  scale_size_discrete("Size", range = c(5, 10), breaks = seq(100, 1, -1))+
  theme_void()+
  # guides(color=guide_legend(ncol=2, byrow=TRUE,override.aes = list(size=5),order = 2),
  #        size=guide_legend(order = 1))+
  theme(legend.title = element_text(size = 16, face="bold"))+
  theme(legend.text = element_text(size = 14))+
  theme(legend.key.size = unit(1, 'cm'))+
  # guides(color = guide_legend(nrow = 2,default.unit = "cm",keywidth = 1,keyheight = 1))+
  theme(legend.position = "top", legend.box="vertical")

filename <- paste0("./03_OutputFiles/",gsub("-","",Sys.Date()),"_",Functionality,"_KEGG_Network_v02.pdf")

ggsave(plot = plot, filename, useDingbats=FALSE, height = 9.5, width = 6.6)

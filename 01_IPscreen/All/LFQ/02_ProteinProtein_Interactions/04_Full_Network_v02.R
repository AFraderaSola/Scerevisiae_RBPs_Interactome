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

files_screen <- list.files(path = "./00_IPsResultsFiles/", pattern = "*.WT.csv")

# Fold change directions

direction <- "Up"

# Plot variables

circular <- F

gene_names <- F

# Define the organism database annotation package for GO:

org <- "org.Sc.sgd.db"

# KEGG Highlight

KEGG <- read.csv("./04_InputFiles/00000000_NetworkKEGG.csv")

IP_KEGG <- read.csv("./04_InputFiles/00000000_PPI_KEGG_table.csv")

# RBP Highliht

RBP <- read.csv("./04_InputFiles/00000000_SC_RBPs_NamedCensus.csv")

#########################################
############## Script ###################
#########################################

names <- read_excel("../../../00_IPScreeningMasterFile.xlsx")

screen <- c()

IP_KEGG$Functionality <- "Random"

for (i in 1:nrow(IP_KEGG)) {
  
  IP_KEGG$Functionality[i] <- KEGG[KEGG$KEGG == IP_KEGG$KEGG[i],]$Functionality
  
}

IP_KEGG <- IP_KEGG[order(IP_KEGG$Functionality),]

IP_KEGG[IP_KEGG$Gene %in% unique(IP_KEGG[duplicated(IP_KEGG$Gene),]$Gene),]$Functionality <- "Multiple"

to_add_KEGG <- c("Multiple", "Multiple", "#666666")

KEGG <- rbind(KEGG, to_add_KEGG)

IP_KEGG$Gene <- str_to_title(IP_KEGG$Gene)

for (i in 1:length(files_screen)) {
  
  screen_loop <- read_csv(paste0("./00_IPsResultsFiles/", files_screen[i]))
  
  if (direction == "Up") {
    
    screen_loop <- screen_loop[screen_loop$difference_RNase_WT > 0,]
    
  }else{
    
    screen_loop <- screen_loop[screen_loop$difference_RNase_WT < 0,]
    
  }
  
  screen_loop_df <- screen_loop
  screen_loop <- unique(screen_loop$Gene.names)[!grepl(unique(screen_loop$Gene.names),pattern = "-TAP")]
  screen_loop <- c(stringr::str_to_title(screen_loop[!screen_loop %in% screen_loop_df$Protein.IDs]),
                   screen_loop[screen_loop %in% screen_loop_df$Protein.IDs])
  screen_loop <- list(screen_loop)
  names(screen_loop) <- gsub(pattern = "_E.*.csv",replacement = "",files_screen[i])
  screen <- c(screen, screen_loop)
  
  
}

names_screen <- names(screen)

names_gene <- c()

for (i in 1:length(names_screen)) {
  
  stock <- as.character(names_screen[i])
  
  gene <- names[names$Stock_Name == stock,]$Gene_Name
  
  names_gene <- c(names_gene, gene)
  
}


names_gene <- stringr::str_to_title(names_gene)

names(screen) <- names_gene

edges <- c()

for (i in 1:length(screen)) {
  
  loop_source <- rep(names(screen)[i], length(screen[[i]]))
  
  loop_target <- screen[[i]]
  
  loop_edges <- data_frame(source = loop_source, 
                           target = loop_target)
  
  edges <- rbind(edges, loop_edges)
  
}

nodes <- unique(c(edges$source, edges$target))

size <- c()

for (i in 1:length(nodes)) {
  
  if (nodes[i] %in% names(screen)) {
    
    loop_size <- length(unlist(screen[names(screen) == nodes[i]]))
    
  }else{
    
    loop_size <- 1
  }
  
  size <- c(size, loop_size)
}

vertices <- data.frame(nodes = nodes,
                       size = size)

vertices$Functionality <- "None"

for (i in 1:nrow(vertices)) {
  
  if (vertices$nodes[i] %in% unique(IP_KEGG$Gene)) {
    
    loop_functionality <- unique(IP_KEGG[IP_KEGG$Gene == vertices$nodes[i],]$Functionality)
    
    vertices$Functionality[i] <- loop_functionality
    
  }
}

edges$Functionality <- "Random"

for (i in 1:nrow(edges)) {
  
  edges$Functionality[i] <- vertices[vertices$nodes == edges$source[i],]$Functionality
  
}
  
  
net <- graph_from_data_frame(d = edges, vertices = vertices)

E(net)$edge.color <- ifelse(E(net)$Functionality == "None", "#b4b4b4", E(net)$Functionality)
E(net)$edge.color <- ifelse(E(net)$edge.color == "Multiple","#a157fa", E(net)$edge.color)
E(net)$edge.color <- ifelse(E(net)$edge.color == "Ribosome", "#66C2A5", E(net)$edge.color)
E(net)$edge.color <- ifelse(E(net)$edge.color == "Splicing", "#E78AC3", E(net)$edge.color)
E(net)$edge.color <- ifelse(E(net)$edge.color == "Synthesis", "#FC8D62", E(net)$edge.color)
E(net)$edge.color <- ifelse(E(net)$edge.color == "Degradation", "#A6D854", E(net)$edge.color)
E(net)$edge.color <- ifelse(E(net)$edge.color == "Export", "#666666", E(net)$edge.color)

plot <- ggnet2(net, color = "Functionality", label = unique(edges$source), size = "size", edge.color = "edge.color", mode = "fruchtermanreingold",
               label.size = 8, 
               palette = c(
                 "Multiple" = "#a157fa",
                 "Ribosome" = "#66C2A5",
                 "Splicing" = "#E78AC3",
                 "Synthesis" = "#FC8D62",
                 "Degradation" = "#A6D854",
                 "Export" = "#666666",
                 "None" = "#b4b4b4"))+
  scale_size_discrete("Size", range = c(5, 10), breaks = seq(90, 1, -5))+
  theme_void()+
  theme(legend.key.size = unit(1, 'cm'))+
  theme(legend.title = element_text(size = 14, face="bold"))+
  theme(legend.text = element_text(size = 12))+
  theme(legend.position = "top", legend.box="vertical")


name <- paste0("./04_OutputFiles/",gsub(pattern = "-",replacement = "",x = Sys.Date()),"_Full_KEGG_IncludingNone_Network_v02.pdf") 

ggsave(plot = plot, name, height = 10, width = 20)

## Excluding nodes functions

edges <- edges[!edges$Functionality == "None",]

filter <- unique(c(edges$source, edges$target))

vertices <- vertices[vertices$nodes %in% filter,]

vertices[vertices$Functionality == "None",]$size <- 1

write.csv(x = edges, file = paste0("./04_OutputFiles/",gsub(pattern = "-",replacement = "",x = Sys.Date()),"_Edges.csv"),quote = F,row.names = F)

write.csv(x = vertices, file = paste0("./04_OutputFiles/",gsub(pattern = "-",replacement = "",x = Sys.Date()),"_Nodes.csv"),quote = F,row.names = F)

net <- graph_from_data_frame(d = edges, vertices = vertices)

E(net)$edge.color <- ifelse(E(net)$Functionality == "None", "#b4b4b4", E(net)$Functionality)
E(net)$edge.color <- ifelse(E(net)$edge.color == "Multiple","#a157fa", E(net)$edge.color)
E(net)$edge.color <- ifelse(E(net)$edge.color == "Ribosome", "#66C2A5", E(net)$edge.color)
E(net)$edge.color <- ifelse(E(net)$edge.color == "Splicing", "#E78AC3", E(net)$edge.color)
E(net)$edge.color <- ifelse(E(net)$edge.color == "Synthesis", "#FC8D62", E(net)$edge.color)
E(net)$edge.color <- ifelse(E(net)$edge.color == "Degradation", "#A6D854", E(net)$edge.color)
E(net)$edge.color <- ifelse(E(net)$edge.color == "Export", "#666666", E(net)$edge.color)

plot <- ggnet2(net, color = "Functionality", label = unique(edges$source), size = "size", edge.color = "edge.color", mode = "fruchtermanreingold",
               label.size = 8, 
               palette = c(
                 "Degradation" = "#A6D854",
                 "Export" = "#666666",
                 "Multiple" = "#a157fa",
                 "Ribosome" = "#66C2A5",
                 "Splicing" = "#E78AC3",
                 "Synthesis" = "#FC8D62",
                 "None" = "#b4b4b4"))+
  scale_size_discrete("Size", range = c(5, 10), breaks = seq(90, 1, -5))+
  theme_void()+
  theme(legend.key.size = unit(1, 'cm'))+
  theme(legend.title = element_text(size = 22, face="bold"))+
  theme(legend.text = element_text(size = 20))+
  theme(legend.position = "top", legend.box="vertical")


name <- paste0("./04_OutputFiles/",gsub(pattern = "-",replacement = "",x = Sys.Date()),"_Full_KEGG_ExcludingNone_Network_v02.pdf") 

ggsave(plot = plot, name, height = 9.5, width = 20)

PPI_All_vis_nodes <- vertices

colnames(PPI_All_vis_nodes) <- c("id", "n", "Functionality")

PPI_All_vis_edges <- edges

colnames(PPI_All_vis_edges) <- c("from", "to", "Functionality")

file <- paste0("./05_OutputFiles/All_PPI_InteractiveNetwork_v02.RData")

save(PPI_All_vis_nodes, PPI_All_vis_edges, file = file)

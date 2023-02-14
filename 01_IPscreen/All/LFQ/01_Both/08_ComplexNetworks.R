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
library(geomtextpath)
library(Cairo)
library(pheatmap)
library(igraph)
library(tidyverse)
library(readr)
library(network)
library(sna)
library(ggplot2)
library(ggnet)
library(GGally)

options(ggrepel.max.overlaps = Inf)

set.seed(666)

#########################################
############ Variables ##################
#########################################

complexes <- read.csv("07_InputFiles/CYC2008_complex.csv")

complexes[complexes$Complex == "Prp9p/Prp11p/Prp21p complex(SF3a complex)",]$Complex <- "SF3a complex"

filtered_complexes <- read.csv("07_OutputFiles/20221024_Splicing_ComplexData.csv")

IP_Network <- read.csv("07_OutputFiles/20221024_Splicing_IP_Network_Cytoscape.csv")

IP_KEGG <- read.csv("./03_InputFiles/00000000_Both_KEGG_table.csv")

KEGG <- read.csv("./03_InputFiles/00000000_NetworkKEGG.csv")

network <- load("./07_InputFiles/Network_df.RData")

functionality <- "Splicing"

#########################################
############# Script ####################
#########################################

cgf <- read.csv(file = paste0("./07_InputFiles/00000000_", functionality, "_ComplexGeneralFunction.csv"))

Complex <- cgf[cgf$Functionality == functionality,]$Complex

Complex <- Complex[-6]

vertices_tosave <- c()

edges_tosave <- c()

for (i in 1:length(Complex)) {
  
  vertices <- read.csv(file = "./08_InputFiles/20221024_Splicing_IP_Network_PLUS_AllComplex_Network_Nodes.csv")
  
  edges <- read.csv(file = "./08_InputFiles/20221024_Splicing_IP_Network_PLUS_AllComplex_Network_Edges.csv")
  
  complex_loop <- Complex[i]
  
  complex_toread <- gsub(pattern = "\\/|-| ",replacement = "",x = complex_loop)
  
  vertices_detailed <- read.csv(file = paste0("./08_InputFiles/20221024_Splicing_IP_IncludedIn_", complex_toread, "_Network_PLUS_", complex_toread, "_Nodes.csv"))
  
  edges_detailed <- read.csv(file = paste0("./08_InputFiles/20221024_Splicing_IP_IncludedIn_", complex_toread, "_Network_PLUS_", complex_toread, "_Edges.csv"))
  
  # All complex Functionality Network
  
  vertices <- vertices[-2,]
  
  vertices$ID <- "Random"
  
  vertices[vertices$name %in% IP_Network$genes,]$ID <- "Prey"
  
  vertices[vertices$name %in% IP_Network$target_RBP,]$ID <- "Bait"
  
  vertices[vertices$name %in% str_to_title(unique(complexes[complexes$Complex %in% filtered_complexes$Complex,]$Name)),]$ID <- "Complex member"
  
  vertices[vertices$name %in% intersect(unique(IP_Network$genes),
                                        str_to_title(unique(complexes[complexes$Complex %in% filtered_complexes$Complex,]$Name))),]$ID <- "Prey & Complex member"
  
  vertices[vertices$name %in% intersect(unique(IP_Network$target_RBP),
                                        str_to_title(unique(complexes[complexes$Complex %in% filtered_complexes$Complex,]$Name))),]$ID <- "Bait & Complex member"
  
  vertices <-  vertices[, grep(pattern = "^name|ID",x = colnames(vertices))]
  
  edges <- edges[-2,]
  
  edges <-  edges[, grep(pattern = "^name|category",x = colnames(edges))]
  
  edges <- edges[,c(2,1)]
  
  edges$name <- gsub(pattern = " \\(interacts with\\) ", replacement = "_", x = edges$name)
  
  edges <- edges %>% separate(name,into = c("Source", "Target"),sep = "_")
  
  edges[edges$category_loop == "",]$category_loop <- "Complex interaction"
  
  
  
  net <- graph_from_data_frame(d = edges, vertices = vertices)
  
  E(net)$edge.color <- ifelse(E(net)$category_loop == "PPI", "#66A61E", E(net)$category_loop)
  E(net)$edge.color <- ifelse(E(net)$edge.color == "RNAmediated", "#D95F02", E(net)$edge.color)
  E(net)$edge.color <- ifelse(E(net)$edge.color == "Overlap", "#7570B3", E(net)$edge.color)
  E(net)$edge.color <- ifelse(E(net)$edge.color == "Complex interaction", "#b4b4b4", E(net)$edge.color)
  
  plot <- ggnet2(net, color = "ID", label = vertices[!vertices$ID == "Prey",]$name, edge.color = "edge.color",
                 label.size = 8,
                 palette = c(
                   "Prey" = "#A6CEE3",
                   "Bait" = "#1F78B4",
                   "Prey & Complex member" = "#FB9A99",
                   "Complex member" = "#E31A1C",
                   "Bait & Complex member" = "#750c0d"))+
    theme_void()+
    # scale_size_discrete("Size", range = c(5, 10), breaks = seq(100, 1, -1))+
    theme(legend.key.size = unit(1, 'cm'))+
    theme(legend.title = element_text(size = 24, face="bold"))+
    theme(legend.text = element_text(size = 22),
          plot.title = element_text(size=26, face="bold", hjust = 0.5))+
    theme(legend.position = "top", legend.box="vertical")
  
  filename <- paste0("./08_OutputFiles/",gsub("-","",Sys.Date()),"_",functionality,"_AllComplex_Network_v01.pdf")
  
  ggsave(plot = plot, filename = filename, useDingbats=FALSE, width = 20, height = 19)
  
  E(net)$edge.color <- ifelse(E(net)$edge.color == "#b4b4b4", "#A6CEE3", E(net)$edge.color)
  E(net)$edge.color <- ifelse(E(net)$edge.color == "#66A61E", "#b4b4b4", E(net)$edge.color)
  E(net)$edge.color <- ifelse(E(net)$edge.color == "#D95F02", "#b4b4b4", E(net)$edge.color)
  E(net)$edge.color <- ifelse(E(net)$edge.color == "#7570B3", "#b4b4b4", E(net)$edge.color)
  
  plot <- ggnet2(net, color = "ID", label = vertices[!vertices$ID == "Prey",]$name, edge.color = "edge.color",
                 label.size = 8,
                 palette = c(
                   "Prey" = "#CAB2D6",
                   "Bait" = "#6A3D9A",
                   "Prey & Complex member" = "#FDBF6F",
                   "Complex member" = "#FF7F00",
                   "Bait & Complex member" = "#ab5703"))+
    theme_void()+
    # scale_size_discrete("Size", range = c(5, 10), breaks = seq(100, 1, -1))+
    theme(legend.key.size = unit(1, 'cm'))+
    theme(legend.title = element_text(size = 24, face="bold"))+
    theme(legend.text = element_text(size = 22),
          plot.title = element_text(size=26, face="bold", hjust = 0.5))+
    theme(legend.position = "top", legend.box="vertical")
  
  filename <- paste0("./08_OutputFiles/",gsub("-","",Sys.Date()),"_",functionality,"_AllComplex_Network_v02.pdf")
  
  ggsave(plot = plot, filename = filename, useDingbats=FALSE, width = 20, height = 19)
  
  # Detailed complex Functionality Network
  
  # vertices_detailed <- vertices_detailed[!vertices_detailed$name == "Dur1",]
  
  vertices_detailed$ID <- "Random"
  
  vertices_detailed[vertices_detailed$name %in% IP_Network$genes,]$ID <- "Prey"
  
  vertices_detailed[vertices_detailed$name %in% IP_Network$target_RBP,]$ID <- "Bait"
  
  vertices_detailed[vertices_detailed$name %in% str_to_title(unique(complexes[complexes$Complex %in% Complex[i],]$Name)),]$ID <- "Complex member"
  
  vertices_detailed[vertices_detailed$name %in% intersect(unique(IP_Network$genes),
                                                          str_to_title(unique(complexes[complexes$Complex %in% Complex[i],]$Name))),]$ID <- "Prey & Complex member"
  
  if (length(intersect(unique(IP_Network$target_RBP),
                       str_to_title(unique(complexes[complexes$Complex %in% Complex[i],]$Name)))) > 0) {
    
    vertices_detailed[vertices_detailed$name %in% intersect(unique(IP_Network$target_RBP),
                                                            str_to_title(unique(complexes[complexes$Complex %in% Complex[i],]$Name))),]$ID <- "Bait & Complex member"
    
  }

  vertices_detailed <-  vertices_detailed[, grep(pattern = "^name|ID",x = colnames(vertices_detailed))]
  
  edges_detailed <-  edges_detailed[, grep(pattern = "^name|category",x = colnames(edges_detailed))]
  
  edges_detailed <- edges_detailed[,c(2,1)]
  
  edges_detailed$name <- gsub(pattern = " \\(interacts with\\) ", replacement = "_", x = edges_detailed$name)
  
  edges_detailed <- edges_detailed %>% separate(name,into = c("Source", "Target"),sep = "_")
  
  edges_detailed[edges_detailed$category_loop == "",]$category_loop <- "Complex interaction"
  
  net <- graph_from_data_frame(d = edges_detailed, vertices = vertices_detailed)
  
  E(net)$edge.color <- ifelse(E(net)$category_loop == "PPI", "#66A61E", E(net)$category_loop)
  E(net)$edge.color <- ifelse(E(net)$edge.color == "RNAmediated", "#D95F02", E(net)$edge.color)
  E(net)$edge.color <- ifelse(E(net)$edge.color == "Overlap", "#7570B3", E(net)$edge.color)
  E(net)$edge.color <- ifelse(E(net)$edge.color == "Complex interaction", "#b4b4b4", E(net)$edge.color)
  
  plot <- ggnet2(net, color = "ID", label = T, edge.color = "edge.color", size = "degree", size.min = 2, size.cut = T,
                 label.size = 8,
                 palette = c(
                   "Prey" = "#A6CEE3",
                   "Bait" = "#1F78B4",
                   "Prey & Complex member" = "#FB9A99",
                   "Complex member" = "#E31A1C",
                   "Bait & Complex member" = "#750c0d"))+
    theme_void()+
    # scale_size_discrete("Size", range = c(5, 10), breaks = seq(100, 1, -1))+
    theme(legend.key.size = unit(1, 'cm'))+
    ggtitle(complex_loop)+
    theme(legend.title = element_text(size = 24, face="bold"))+
    theme(legend.text = element_text(size = 22),
          plot.title = element_text(size=26, face="bold", hjust = 0.5))+
    theme(legend.position = "top", legend.box="vertical")
  
  filename <- paste0("./08_OutputFiles/",gsub("-","",Sys.Date()),"_",functionality,"_",
                     gsub(pattern = "/|-",
                          replacement = "",
                          x = gsub(pattern = " ",x = complex_loop,replacement = "")),
                     "_Network_v01.pdf")
  
  ggsave(plot = plot, filename = filename, useDingbats=FALSE, width = 20, height = 19)
  
  
  E(net)$edge.color <- ifelse(E(net)$edge.color == "#b4b4b4", "#A6CEE3", E(net)$edge.color)
  E(net)$edge.color <- ifelse(E(net)$edge.color == "#66A61E", "#b4b4b4", E(net)$edge.color)
  E(net)$edge.color <- ifelse(E(net)$edge.color == "#D95F02", "#b4b4b4", E(net)$edge.color)
  E(net)$edge.color <- ifelse(E(net)$edge.color == "#7570B3", "#b4b4b4", E(net)$edge.color)
  
  plot <- ggnet2(net, color = "ID", label = T, edge.color = "edge.color", size = "degree", size.min = 2, size.cut = T,
                 label.size = 10,
                 palette = c(
                   "Prey" = "#CAB2D6",
                   "Bait" = "#6A3D9A",
                   "Prey & Complex member" = "#FDBF6F",
                   "Complex member" = "#FF7F00",
                   "Bait & Complex member" = "#ab5703"))+
    theme_void()+
    # scale_size_discrete("Size", range = c(5, 10), breaks = seq(100, 1, -1))+
    theme(legend.key.size = unit(1, 'cm'))+
    ggtitle(complex_loop)+
    theme(legend.title = element_text(size = 24, face="bold"))+
    theme(legend.text = element_text(size = 22),
          plot.title = element_text(size=26, face="bold", hjust = 0.5))+
    theme(legend.position = "top", legend.box="vertical")
  
  filename <- paste0("./08_OutputFiles/",gsub("-","",Sys.Date()),"_",functionality,"_",
                     gsub(pattern = "/|-",
                          replacement = "",
                          x = gsub(pattern = " ",x = complex_loop,replacement = "")),
                     "_Network_v02.pdf")
  
  ggsave(plot = plot, filename = filename, useDingbats=FALSE, width = 20, height = 19)
  
  
  ## Only preys and baits in complex:
  
  vertices_complex <- vertices_detailed[grep(pattern = "Complex|Bait",x = vertices_detailed$ID),]
  
  vertices_complex <- vertices_complex[!vertices_complex$ID == "Complex member",]
  
  edges_complex <- edges_detailed[!grepl(pattern = "Complex", x = edges_detailed$category_loop),]
  
  edges_complex <- edges_complex[edges_complex$Source %in% vertices_complex[grep(pattern = "Bait",x = vertices_complex$ID),]$name,]
  
  edges_complex <- edges_complex[edges_complex$Target %in% vertices_complex$name,]
  
  edges_complex <- edges_complex[!edges_complex$Target %in% vertices_complex[vertices_complex$ID == "Bait",][!vertices_complex[vertices_complex$ID == "Bait",]$name %in% edges_complex$Source,]$name,]
  
  vertices_complex <- vertices_complex[!vertices_complex$name %in% vertices_complex[vertices_complex$ID == "Bait",][!vertices_complex[vertices_complex$ID == "Bait",]$name %in% edges_complex$Source,]$name,]
  
  net <- graph_from_data_frame(d = edges_complex, vertices = vertices_complex)
  
  E(net)$edge.color <- ifelse(E(net)$category_loop == "PPI", "#66A61E", E(net)$category_loop)
  E(net)$edge.color <- ifelse(E(net)$edge.color == "RNAmediated", "#D95F02", E(net)$edge.color)
  E(net)$edge.color <- ifelse(E(net)$edge.color == "Overlap", "#7570B3", E(net)$edge.color)
  E(net)$edge.color <- ifelse(E(net)$edge.color == "Complex interaction", "#b4b4b4", E(net)$edge.color)
  
  
  plot <- ggnet2(net, color = "ID", label = T, edge.color = "edge.color", size = 40, 
                 label.size = 12, edge.size = 4, mode = "kamadakawai",
                 palette = c(
                   "Bait" = "#CAB2D6",
                   "Prey & Complex member" = "#FDBF6F",
                   "Bait & Complex member" = "#FF7F00"))+
    theme_void()+
    # scale_size_discrete("Size", range = c(35, 40), breaks = c(seq(26, 1, -2),1))+
    theme(legend.key.size = unit(1, 'cm'))+
    ggtitle(complex_loop)+
    theme(legend.title = element_text(size = 38, face="bold"))+
    theme(legend.text = element_text(size = 36),
          plot.title = element_text(size=36, face="bold", hjust = 0.5))+
    theme(legend.position = "none", legend.box="vertical")
  
  filename <- paste0("./08_OutputFiles/",gsub("-","",Sys.Date()),"_",functionality,"_",
                     gsub(pattern = "/|-",
                          replacement = "",
                          x = gsub(pattern = " ",x = complex_loop,replacement = "")),
                     "_Exclusive_Network.pdf")
  
  ggsave(plot = plot, filename = filename, useDingbats=FALSE, width = 15, height = 13.5)
  
  vertices_complex$Complex <- complex_loop
  
  edges_complex$Complex <- complex_loop
  
  vertices_tosave <- rbind(vertices_tosave, vertices_complex)
  
  edges_tosave <- rbind(edges_tosave, edges_complex)
  
}


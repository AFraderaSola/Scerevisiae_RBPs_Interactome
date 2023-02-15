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
library(network)
library(scales)
library(sna)
library(intergraph)

options(ggrepel.max.overlaps = Inf)

set.seed(666)

#########################################
############ Variables ##################
#########################################

# Files to plot

files_screen <- list.files(path = "./00_IPsResultsFiles/", pattern = "*.RNase.csv")

# Fold change directions

direction <- "Up"

# Plot variables

circular <- F

gene_names <- F

# Define the organism database annotation package for GO:

org <- "org.Sc.sgd.db"

# KEGG Highlight

KEGG <- read.csv("./04_InputFiles/00000000_NetworkKEGG.csv")

IP_KEGG <- read.csv("./04_InputFiles/00000000_RNAm_KEGG_table.csv")

Functionality <- "Splicing"

# RBP Highliht

RBP <- read.csv("./04_InputFiles/00000000_SC_RBPs_NamedCensus.csv")

#########################################
############## Script ###################
#########################################

names <- read_excel("../../../00_IPScreeningMasterFile.xlsx")

screen <- c()

for (i in 1:length(files_screen)) {
  
  screen_loop <- read_csv(paste0("./00_IPsResultsFiles/", files_screen[i]))
  
  if (direction == "Up") {
    
    screen_loop <- screen_loop[screen_loop$difference_IP_RNase > 0,]
    
  }else{
    
    screen_loop <- screen_loop[screen_loop$difference_IP_RNase < 0,]
    
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

filter <- KEGG[KEGG$Functionality == Functionality,]$KEGG

filter <- unique(IP_KEGG[IP_KEGG$KEGG %in% filter,]$Gene)

filter <- stringr::str_to_title(filter)

screen <- screen[names(screen) %in% filter]

edges <- c()

for (i in 1:length(screen)) {
  
  loop_edges_source <- rep(names(screen[i]), length(screen[[i]]))
  
  loop_edges_target <- screen[[i]]
  
  loop_edges <- data.frame(source = loop_edges_source,
                           target = loop_edges_target)
  
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


census <- stringr::str_to_title(unique(RBP$GENENAME))

biogrid <- c()

pattern <- paste(stringr::str_to_upper(filter), collapse = "|")

files <- list.files(path = "../../../../00_Candidates/00_CandidateInteractions/BioGRID/", pattern = pattern)

for (i in 1:length(files)) {
  
  loop <- read_delim(paste0("../../../../00_Candidates/00_CandidateInteractions/BioGRID/", files[i]), "\t", escape_double = FALSE, trim_ws = TRUE,col_names = T)
  
  loop <- loop[loop$`Experimental System Type` == "physical",]
  
  loop <- unique(c(loop$`Official Symbol Interactor A`, loop$`Official Symbol Interactor B`))
  
  loop <- list(loop)
  
  biogrid <- c(biogrid, loop)
  
}

biogrid <- stringr::str_to_title(unique(unlist(biogrid, use.names = F)))

both <- stringr::str_to_title(intersect(census,biogrid))

vertices$ID <- "None"

vertices[vertices$nodes %in% biogrid,]$ID <- "BioGRID"

vertices[vertices$nodes %in% census,]$ID <- "RBP census"

vertices[vertices$nodes %in% both,]$ID <- "BioGRID & RBP census"

filter <- KEGG[KEGG$Functionality == Functionality,]$KEGG

filter <- unique(IP_KEGG[IP_KEGG$KEGG %in% filter,]$Gene)

filter <- stringr::str_to_title(filter)

vertices[vertices$nodes %in% filter,]$ID <- "Bait"

vertices$nodes <- gsub(pattern = ",",replacement = "",x = vertices$nodes)

edges$target <- gsub(pattern = ",",replacement = "",x = edges$target)

write.csv(x = edges, file = paste0("./04_OutputFiles/",gsub(pattern = "-",replacement = "",x = Sys.Date()),"_SplicingEdges.csv"),quote = F,row.names = F)

write.csv(x = vertices, file = paste0("./04_OutputFiles/",gsub(pattern = "-",replacement = "",x = Sys.Date()),"_SplicingNodes.csv"),quote = F,row.names = F)

net <- graph_from_data_frame(d = edges, vertices = vertices)

colorsA <- brewer.pal(n = 4, name = "Greys")

colorsB <- brewer.pal(n = 6, "Paired")

colors <- c(colorsA[3:4],
            colorsB[1:2],
            unique(KEGG[KEGG$Functionality == Functionality,]$color))

filter2 <- as.data.frame(sort(table(c(edges$source, edges$target)), decreasing = TRUE))

filter2 <- filter2[filter2$Freq > 4,]

filter2 <- filter2$Var1

filter3 <- c(vertices[vertices$ID == "None",]$nodes, vertices[vertices$ID == "RBP census",]$nodes, vertices[vertices$ID == "Bait",]$nodes)

plot <- ggnet2(net, size = "size", color = "ID", label = filter3, mode = "fruchtermanreingold",
               label.size = 8, edge.size = 0.4, 
               palette = c(
                 "RBP census" = colors[4],
                 "None" = colors[3],
                 "BioGRID & RBP census" = colors[2],
                 "BioGRID" = colors[1],
                 "Bait" = colors[5]))+
  scale_size_discrete("Size", range = c(5, 10), breaks = seq(75, 1, -1))+
  theme_void()+
  # guides(color=guide_legend(ncol=2, byrow=TRUE,override.aes = list(size=5),order = 2),
  #        size=guide_legend(order = 1))+
  theme(legend.title = element_text(size = 16, face="bold"))+
  theme(legend.text = element_text(size = 14))+
  # theme(legend.key.size = unit(1, 'cm'))+
  guides(color = guide_legend(nrow = 2,default.unit = "cm",keywidth = 1,keyheight = 1))+
  theme(legend.position = "top", legend.box="vertical")

filename <- paste0("./04_OutputFiles/",gsub("-","",Sys.Date()),"_",Functionality,"_KEGG_Network.pdf") 

ggsave(plot = plot, filename = filename, useDingbats=FALSE, width = 6.6, height = 9.5)


degree <- igraph::degree(net)

degree_prey <- degree[!names(degree) %in% filter]

degree_prey <- data.frame(prey = names(degree_prey),
                          degree = degree_prey,
                          Group = "RDI")

RDIdegree_prey <- degree_prey

save(RDIdegree_prey, file = "./04_OutputFiles/RDIdegree_prey.RData")

p <- ggplot(data=degree_prey, aes(x=degree, fill = Group)) +
  geom_density(aes(x = degree, y = ..density..), alpha = .4) +
  scale_fill_manual(values = c("#D95F02"))+
  theme_minimal() +
  theme(axis.text=element_text(size=22),
        axis.title=element_text(size=24,face="bold"))+
  theme(legend.text = element_text(size = 16),
        legend.title=element_text(size=18,face="bold"))+
  theme(legend.key.size = unit(1, 'cm'))+
  theme(legend.position = "top")+
  labs(fill = "Data set")+
  xlab("Degree") +
  ylab("Prey")

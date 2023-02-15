#############################################
############# Libraries #####################
#############################################

library(readr)
library(ggplot2)
library(dplyr)
library(readxl)
library(scales)
library(ggthemes)
library(extrafont)
library(viridis)
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

#############################################
############### Script ######################
#############################################

files <- list.files(path ="./00_KOsResultsFiles/", pattern = "Enriched.*\\.csv")

stock <- read_excel("../00_KOScreeningMasterFile.xlsx")

complexes <- read.csv("07_InputFiles/CYC2008_complex.csv")

comparisons <- gsub(pattern = ".*Enriched_|.csv",replacement = "",x = files)

pg <- c()

plot_df <- c()

i <- 1

for (i in 1:length(files)) {
  
  loop <- read_csv(paste0("./00_KOsResultsFiles/", files[i]))
  
  loop_up <- loop[,grep(pattern = "difference|Gene.names",x = colnames(loop))]  
  loop_up <- loop_up[unlist(list(loop_up[,2]),use.names = F) > 0,]
  
  loop_up$Strain <- gsub(x = gsub(pattern = ".csv",replacement = "",x = files[i]), pattern = "Enriched_", replacement = "")
  
  loop_up$Expression <- "Up-regulated"
  
  loop_up$Gene.names <- str_to_title(loop_up$Gene.names)
  
  loop_up <- loop_up[,c(1,3,4)]
  
  loop_down <- loop[,grep(pattern = "difference|Gene.names",x = colnames(loop))]  
  loop_down <- loop_down[unlist(list(loop_down[,2]),use.names = F) < 0,]
  
  loop_down$Strain <- gsub(x = gsub(pattern = ".csv",replacement = "",x = files[i]), pattern = "Enriched_", replacement = "")
  
  loop_down$Expression <- "Down-regulated"
  
  loop_down$Gene.names <- str_to_title(loop_down$Gene.names)
  
  loop_down <- loop_down[,c(1,3,4)]
  
  loop <- rbind(loop_down,loop_up)
  
  plot_df <- rbind(plot_df, loop)
  
}

for (i in 1:nrow(plot_df)) {
  
  exp <- plot_df$Strain[i]
  
  stock_loop <- gsub(pattern = "_.*",replacement = "",x = exp)
  
  gene <- stock[stock$Stock_Name == stock_loop,]$Gene_Name
  
  gene <- str_to_title(gene)
  
  plot_df$Strain[i] <- gene
  
}

plot_df <- plot_df[!plot_df$Strain == "Cbc2",]

plot_df <- plot_df[!plot_df$Strain == "Msl1",]

complexes <- read.csv("07_InputFiles/CYC2008_complex.csv")

complex_id <- unique(complexes$Complex)

complex_n <- complexes %>%
  group_by(Complex) %>%
  summarise(n = n())

complex_df <- c()

for (i in 1:length(unique(plot_df$Strain))) {
  
  loop_df_strain <- plot_df[plot_df$Strain == unique(plot_df$Strain)[i], ]
  
  for (j in 1:length(unique(loop_df_strain$Expression))) {
    
    loop_df_exp <- loop_df_strain[loop_df_strain$Expression == unique(loop_df_strain$Expression)[j], ]
    
    for (k in 1:length(complex_id)) {
      
      complex_loop <- complexes[complexes$Complex == complex_id[k],]
      
      loop_overlap <- length(intersect(str_to_title(complex_loop$Name), unique(loop_df_exp$Gene.names)))
      
      loop_n <- complex_n[complex_n$Complex == complex_id[k],]$n
      
      loop_ratio <- loop_overlap/loop_n
      
      loop_df <- data.frame(Strain = unique(loop_df_exp$Strain),
                            Expression = unique(loop_df_exp$Expression),
                            Complex = complex_id[k],
                            n_complex = loop_n,
                            n_network = loop_overlap,
                            ratio = loop_ratio)
      
      complex_df <- rbind(complex_df, loop_df)
      
    }
  }
}

p <- ggplot(data=complex_df, aes(x=ratio, fill= Expression)) +
  facet_wrap(~Strain)+
  geom_density(aes(x = ratio), alpha = .4) +
  scale_fill_manual(values = c("#404080","#69b3a2"))+
  theme_minimal() +
  ggtitle(paste0("n complexes = ", length(unique(complex_df$Complex))))+
  scale_x_continuous(breaks = seq(0, 1, by = 0.2))+
  theme(axis.text=element_text(size=22),
        axis.title=element_text(size=24,face="bold"))+
  theme(legend.text = element_text(size = 16),
        legend.title=element_text(size=18,face="bold"))+
  theme(legend.key.size = unit(1, 'cm'))+
  theme(legend.position = "top")+
  theme(strip.text = element_text(size = 20))+
  theme(plot.title = element_text(size = 26, face = "bold"))+
  xlab("Complex ratio") +
  ylab("N-complexes")

filename <- paste0("./07_OutputFiles/",gsub(pattern = "-",x = Sys.Date(),replacement = ""),"_ComplexDensity_byStrain.pdf")

ggsave(p, filename = filename, width = 15, height = 7.5)

p <- ggplot(data=complex_df, aes(x=ratio, fill= Strain)) +
  facet_wrap(~Expression)+
  geom_density(aes(x = ratio), alpha = .4) +
  scale_fill_viridis_d()+
  theme_minimal() +
  ggtitle(paste0("n complexes = ", length(unique(complex_df$Complex))))+
  scale_x_continuous(breaks = seq(0, 1, by = 0.2))+
  theme(axis.text=element_text(size=22),
        axis.title=element_text(size=24,face="bold"))+
  theme(legend.text = element_text(size = 16),
        legend.title=element_text(size=18,face="bold"))+
  theme(legend.key.size = unit(1, 'cm'))+
  theme(legend.position = "top")+
  theme(strip.text = element_text(size = 20))+
  theme(plot.title = element_text(size = 26, face = "bold"))+
  xlab("Complex ratio") +
  ylab("N-complexes")

filename <- paste0("./07_OutputFiles/",gsub(pattern = "-",x = Sys.Date(),replacement = ""),"_ComplexDensity_byExpression.pdf")

ggsave(p, filename = filename, width = 15, height = 7.5)

complex_df <- complex_df[complex_df$ratio != 0, ]

p <- ggplot(data=complex_df, aes(x=ratio, fill= Expression)) +
  facet_wrap(~Strain)+
  geom_density(aes(x = ratio), alpha = .4) +
  scale_fill_manual(values = c("#404080","#69b3a2"))+
  theme_minimal() +
  ggtitle(paste0("n complexes = ", length(unique(complex_df$Complex))))+
  scale_x_continuous(breaks = seq(0, 1, by = 0.2))+
  theme(axis.text=element_text(size=22),
        axis.title=element_text(size=24,face="bold"))+
  theme(legend.text = element_text(size = 16),
        legend.title=element_text(size=18,face="bold"))+
  theme(legend.key.size = unit(1, 'cm'))+
  theme(legend.position = "top")+
  theme(strip.text = element_text(size = 20))+
  theme(plot.title = element_text(size = 26, face = "bold"))+
  xlab("Complex ratio") +
  ylab("N-complexes")

filename <- paste0("./07_OutputFiles/",gsub(pattern = "-",x = Sys.Date(),replacement = ""),"_ComplexDensity_byStrain_no0.pdf")

ggsave(p, filename = filename, width = 15, height = 7.5)

p <- ggplot(data=complex_df, aes(x=ratio, fill= Strain)) +
  facet_wrap(~Expression)+
  geom_density(aes(x = ratio, y = ..density..), alpha = .4) +
  scale_fill_viridis_d()+
  theme_minimal() +
  ggtitle(paste0("n complexes = ", length(unique(complex_df$Complex))))+
  scale_x_continuous(breaks = seq(0, 1, by = 0.2))+
  theme(axis.text=element_text(size=22),
        axis.title=element_text(size=24,face="bold"))+
  theme(legend.text = element_text(size = 16),
        legend.title=element_text(size=18,face="bold"))+
  theme(legend.key.size = unit(1, 'cm'))+
  theme(legend.position = "top")+
  theme(strip.text = element_text(size = 20))+
  theme(plot.title = element_text(size = 26, face = "bold"))+
  xlab("Complex ratio") +
  ylab("N-complexes")

filename <- paste0("./07_OutputFiles/",gsub(pattern = "-",x = Sys.Date(),replacement = ""),"_ComplexDensity_byExpression_no0.pdf")

ggsave(p, filename = filename, width = 15, height = 7.5)

complex_df <- complex_df[complex_df$ratio > 0.4, ]

plot <- ggplot(complex_df, aes(fill = n_complex,y = ratio , x= Complex)) +
  facet_wrap(~Strain, scales = c("free_x"))+
  geom_bar(position=position_dodge(width = 1), stat="identity", size = 1.5)+
  scale_fill_viridis(breaks = c(5,10,15,20,25,30),option = "G",direction = -1)+
  theme_minimal()+
  theme(axis.text.x=element_text(size=12,angle = 45, hjust = 1),
        axis.text.y=element_text(size=12),
        axis.title=element_text(size=15,face="bold"))+
  theme(legend.text = element_text(size = 14),
        legend.title=element_text(size=12,face="bold"))+
  theme(legend.key.size = unit(1, 'cm'))+
  theme(strip.text = element_text(size=18,face="bold"))+
  xlab("Protein complex")+
  ylab("Ratio (n_target / n_complex)")

filename <- paste0("./07_OutputFiles/",gsub(pattern = "-",x = Sys.Date(),replacement = ""),"_Barplot.pdf")

ggsave(plot, filename = filename, width = 30, height = 15)

edges <- complex_df[,c(1,3,2,4)]

vertices <- rbind(data.frame(name = unique(edges$Strain),
                             Node = "Bait"),
                  distinct(data.frame(name = edges$Complex,
                             Node = paste0("Complex  (n = ", edges$n_complex, ")"))))
                  

net <- graph_from_data_frame(d = edges, vertices = vertices)

E(net)$edge.color <- ifelse(E(net)$Expression == "Down-regulated", "#1F78B4", E(net)$Expression)
E(net)$edge.color <- ifelse(E(net)$Expression == "Up-regulated", "#E31A1C", E(net)$edge.color)

colours <- c(viridis(n = 32, option = "G",direction = -1))

plot <- ggnet2(net, color = "Node", label = T, edge.color = "edge.color", size = 40, 
               label.size = 12, edge.size = 4,
               palette = c(
                 "Bait" = "#cab3d6",
                 "Complex  (n = 2)" = "#CEEED7FF",
                 "Complex  (n = 3)" = "#AAE1BDFF",
                 "Complex  (n = 4)" = "#7CD6AFFF",
                 "Complex  (n = 32)" = "#0B0405FF"))+
  theme_void()+
  # scale_size_discrete("Size", range = c(35, 40), breaks = c(seq(26, 1, -2),1))+
  theme(legend.key.size = unit(1, 'cm'))+
  theme(legend.title = element_text(size = 38, face="bold"))+
  theme(legend.text = element_text(size = 36))+
  theme(legend.position = "top", legend.box="vertical")

filename <- paste0("./07_OutputFiles/",gsub(pattern = "-",x = Sys.Date(),replacement = ""),"_Network.pdf")

ggsave(plot = plot, filename = filename, useDingbats=FALSE, width = 45, height = 15)

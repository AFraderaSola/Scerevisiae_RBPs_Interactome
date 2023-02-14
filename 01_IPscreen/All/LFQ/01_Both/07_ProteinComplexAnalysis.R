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

options(ggrepel.max.overlaps = Inf)

set.seed(666)

#########################################
############ Variables ##################
#########################################

complexes <- read.csv("07_InputFiles/CYC2008_complex.csv")

IP_KEGG <- read.csv("./03_InputFiles/00000000_Both_KEGG_table.csv")

KEGG <- read.csv("./03_InputFiles/00000000_NetworkKEGG.csv")

network <- load("./07_InputFiles/Network_df.RData")

functionality <- "Splicing"

#########################################
############# Script ####################
#########################################

cgf <- read.csv(file = paste0("./07_InputFiles/00000000_", functionality, "_ComplexGeneralFunction.csv"))

Complex <- cgf[cgf$Functionality == functionality,]$Complex

k <- 1

for (k in 1:length(Complex)) {
  
  target_complex <- Complex[k]
  
  filter_function <- KEGG[KEGG$Functionality == functionality,]$KEGG
  
  filter_id <- c()
  
  for (i in 1:length(filter_function)) {
    
    loop <- IP_KEGG[IP_KEGG$KEGG == filter_function[i],]$Gene
    
    filter_id <- c(filter_id, loop)
    
  }
  
  filter_id <- unique(filter_id)
  
  category_df$target_RBP <- str_to_upper(category_df$target_RBP)
  
  category_df$genes <- str_to_upper(category_df$genes)
  
  df <- category_df[category_df$target_RBP %in% filter_id,]
  
  ## For cytoscape
  
  df_cytoscape <- df
  
  df_cytoscape$genes <- str_to_title(df_cytoscape$genes)
  
  df_cytoscape$target_RBP <- str_to_title(df_cytoscape$target_RBP)
  
  filter <- c()
  
  for (i in 1:nrow(df_cytoscape)) {
    
    loop_a <- df_cytoscape$target_RBP[i]
    
    loop_b <- df_cytoscape$genes[i]
    
    if (loop_a == loop_b) {
      
      filter <- c(filter, i)
      
    }
  }
  
  df_cytoscape <- df_cytoscape[-filter,]
  
  write.csv(x = df_cytoscape, file = paste0("./07_OutputFiles/", gsub(pattern = "-",x = Sys.Date(),replacement = ""),"_",functionality, "_IP_Network_Cytoscape.csv"),quote = F,row.names = F)
  
  ##
  
  complex_id <- unique(complexes$Complex)
  
  complex_n <- complexes %>%
    group_by(Complex) %>%
    summarise(n = n())
  
  complex_df <- c()
  
  for (i in 1:length(complex_id)) {
    
    complex_loop <- complexes[complexes$Complex == complex_id[i],]
    
    loop_overlap <- length(intersect(complex_loop$Name, unique(df$genes)))
    
    loop_n <- complex_n[complex_n$Complex == complex_id[i],]$n
    
    loop_ratio <- loop_overlap/loop_n
    
    loop_df <- data.frame(Complex = complex_id[i],
                          n_complex = loop_n,
                          n_network = loop_overlap,
                          ratio = loop_ratio)
    
    complex_df <- rbind(complex_df, loop_df)
    
  }
  
  plot <- ggplot(complex_df, aes(x=ratio)) +
    geom_density(fill="#69b3a2", color="#e9ecef", alpha=0.8)+
    geom_vline(xintercept = 0.5, color = "#ad2803")+
    ggtitle(paste0("n complexes = ", length(unique(complex_df$Complex))))+
    theme_minimal()+
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=14,face="bold"))+
    theme(legend.text = element_text(size = 14),
          legend.title=element_text(size=12,face="bold"))
  
  filename <- paste0("./07_OutputFiles/",gsub(pattern = "-",x = Sys.Date(),replacement = ""),"_DensityRatioFull.pdf")
  
  ggsave(plot, filename = filename, width = 15, height = 7.5)
  
  plot <- ggplot(subset(complex_df, ratio > 0.5), aes(x=ratio)) +
    geom_density(fill="#69b3a2", color="#e9ecef", alpha=0.8)+
    ggtitle(paste0("n complexes = ", length(unique(subset(complex_df, ratio > 0.5)$Complex))))+
    theme_minimal()+
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=14,face="bold"))+
    theme(legend.text = element_text(size = 14),
          legend.title=element_text(size=12,face="bold"))
  
  filename <- paste0("./07_OutputFiles/",gsub(pattern = "-",x = Sys.Date(),replacement = ""),"_DensityFiltered.pdf")
  
  ggsave(plot, filename = filename, width = 15, height = 7.5)
  
  #### Full network ratios plot
  
  plot_df <- complex_df[complex_df$ratio > 0.5,]
  
  write.csv(x = plot_df, file = paste0("./07_OutputFiles/", gsub(pattern = "-",x = Sys.Date(),replacement = ""),
                                       "_",functionality, "_ComplexData.csv"),quote = F,row.names = F)
  
  ### For cytoscape
  
  complex_cytoscape <- complexes[complexes$Complex %in% plot_df$Complex,]
  
  complex_cytoscape <- complex_cytoscape[,c(3,2)]
  
  complex_cytoscape$Name <- str_to_title(complex_cytoscape$Name)
  
  complex_ids <- unique(complex_cytoscape$Complex)
  
  complex_cytoscape_df <- c()
  
  for (i in 1:length(complex_ids)) {
    
    loop_df <- complex_cytoscape[complex_cytoscape$Complex == complex_ids[i],]
    
    loop_complex <- loop_df$Name 
    
    for (j in 1:length(loop_complex)) {
      
      bait_loop <- rep(loop_complex[j],length(loop_complex)-1)
      
      prey_loop <- loop_complex[-j]
      
      bait_prey_loop_df <- data.frame(bait = bait_loop,
                                      prey = prey_loop,
                                      complex = rep(complex_ids[i], length(bait_loop)),
                                      n_complex = rep(length(loop_complex), length(bait_loop)))
      
      complex_cytoscape_df <- rbind(complex_cytoscape_df, bait_prey_loop_df)
    }
  }
  
  
  write.csv(x = complex_cytoscape_df, file = paste0("./07_OutputFiles/",gsub(pattern = "-",x = Sys.Date(),replacement = ""),"_", functionality, "_AllComplex_Network_Cytoscape.csv"),quote = F,row.names = F)
  
  complex_cytoscape_df[1:6,]$complex <- "SF3a complex"
  
  write.csv(x = complex_cytoscape_df[complex_cytoscape_df$complex == target_complex,], 
            file = paste0("./07_OutputFiles/", gsub(pattern = "-",x = Sys.Date(),replacement = ""),"_", 
                          functionality,
                          "_",
                          gsub(pattern = "/|-",
                               replacement = "",
                               x = gsub(pattern = " ",x = target_complex,replacement = "")), "_Network_Cytoscape.csv"),quote = F,row.names = F)
  
  #### 
  
  #### Full network complex ratios
  
  plot_df$Complex[1] <- "SF3a complex" #### Complex label abbreviation
  
  plot_df <- plot_df %>% arrange(desc(n_complex),ratio)
  
  plot_df$ComplexType <- c(rep("Splicing", 8), "Other", rep("Splicing",2), rep("Other", 8))
  
  plot_df <- plot_df %>% arrange(ComplexType, desc(n_complex),ratio)
  
  plot_df$Complex <- factor(x = plot_df$Complex, levels = plot_df$Complex)
  
  plot <- ggplot(plot_df) +
    geom_textsegment(aes(x=Complex, xend=Complex, y=n_complex, yend=n_network), color="#2e2e2e", label = round(plot_df$ratio,digits = 2), size = 5) +
    geom_point(aes(x=Complex, y=n_complex), color="#FDBF6F", size=5,) +
    geom_point(aes(x=Complex, y=n_network), color="#CAB2D6", size=5) +
    geom_point(data = subset(plot_df, ratio == 1), aes(x=Complex, y=n_complex), 
               color="#CAB2D6", fill = "#FDBF6F", shape = 21, size=5, stroke =1) +
    geom_text(data = subset(plot_df, ratio == 1), aes(x=Complex, y=n_complex), label = "1", size = 5)+
    annotate(
      xmin = c(-Inf, 9.5),
      xmax = c(9.5, Inf),
      ymin = 29.5,
      ymax = 30,
      geom = "rect",
      fill = c("#b4b4b4", "#E78AC3")
    )+
    coord_flip()+
    scale_y_continuous(breaks = seq(0, 30, by = 5))+
    theme_minimal() +
    theme(axis.text=element_text(size=16),
          axis.title=element_text(size=18,face="bold"))+
    theme(legend.title = element_text(size = 16, face="bold"))+
    theme(legend.text = element_text(size = 14))+
    theme(legend.key.size = unit(1, 'cm'))+
    theme(legend.position = "top")+
    xlab("Protein complex") +
    ylab("Number of proteins")
  
  filename <- paste0("./07_OutputFiles/",gsub(pattern = "-",x = Sys.Date(),replacement = ""),"_ComplexRatio.pdf")
  
  ggsave(plot, filename = filename, width = 10, height = 9.5)
  
  #### Bait complex ratios
  
  # Barplot
  
  targets <- unique(df$target_RBP)
  
  complex_df <- c()
  
  for (i in 1:length(complex_id)) {
    
    for (j in 1:length(targets)) {
      
      complex_loop <- complexes[complexes$Complex == complex_id[i],]
      
      loop_targets <- df[df$target_RBP == targets[j],]
      
      loop_overlap <- length(intersect(complex_loop$Name, loop_targets$genes))
      
      loop_n <- complex_n[complex_n$Complex == complex_id[i],]$n
      
      loop_ratio <- loop_overlap/loop_n
      
      loop_df <- data.frame(Target = targets[j],
                            Complex = complex_id[i],
                            n_complex = loop_n,
                            n_network = loop_overlap,
                            ratio = loop_ratio)
      
      complex_df <- rbind(complex_df, loop_df)
      
    }
  }
  
  plot_df$Complex <- as.character(plot_df$Complex)
  
  plot_df$Complex[18] <- "Prp9p/Prp11p/Prp21p complex(SF3a complex)"
  
  plot_df <- complex_df[complex_df$Complex %in% plot_df$Complex,]
  
  plot_df <- plot_df[plot_df$ratio > 0, ]
  
  plot_df$Complex[1:3] <- "SF3a complex"
  
  plot_df <- plot_df %>% group_by(Target) %>% arrange(desc(n_complex),ratio)
  
  plot_df$Complex <- factor(x = plot_df$Complex, levels = unique(plot_df$Complex))
  
  plot <- ggplot(plot_df, aes(fill = n_complex,y = ratio , x= Complex)) +
    facet_wrap(~Target, scales = c("free_x"),ncol = 2)+
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
  
  filename <- paste0("./07_OutputFiles/",gsub(pattern = "-",x = Sys.Date(),replacement = ""),"_TargetRatio.pdf")
  
  ggsave(plot, filename = filename, width = 15, height = 15)
  
  
  ### IP complexed included network
  
  bait_filter <- str_to_title(plot_df[plot_df$Complex == target_complex,]$Target)
  
  df_cytoscape <- df_cytoscape[df_cytoscape$target_RBP %in% bait_filter,]
  
  write.csv(x = df_cytoscape, file = paste0("./07_OutputFiles/", 
                                            gsub(pattern = "-",x = Sys.Date(),replacement = ""),"_", 
                                            functionality,
                                            "_IP_IncludedIn_",
                                            gsub(pattern = "/|-",
                                                 replacement = "",
                                                 x = gsub(pattern = " ",x = target_complex,replacement = "")), "_Network_Cytoscape.csv"),quote = F,row.names = F)
  
  
}

# Heatmap

plot_df$Target <- str_to_title(plot_df$Target)

targets <- unique(plot_df$Target)

df_list <- c()

for (i in 1:length(targets)) {
  
  loop_df <- plot_df[plot_df$Target == targets [i], ][,c(2,5)]
  
  colnames(loop_df)[2] <- paste0(targets[i], "_ratio")
  
  loop_df <- list(loop_df)
  
  names(loop_df) <- targets[i]
  
  df_list <- c(df_list, loop_df)
  
}

df_HM <- df_list[[1]]

for (i in 2:length(df_list)) {
  
  df_HM <- full_join(x = df_HM,y = df_list[[i]], by = "Complex")
  
}


df_HM <- as.data.frame(df_HM)

rownames(df_HM) <- df_HM[,1]

df_HM <- df_HM[,c(2:ncol(df_HM))]

df_HM <- as.matrix(df_HM)

colours <- c(viridis(n = 100, option = "G",direction = -1))

# Group <- c(rep("PPI", ncol(df_HM_ppi[,-1])),rep("RDI", ncol(df_HM_RNAmediated[,-1])))
# 

bfc <- read.csv(file = "./02_InputFiles/00000000_BaitFunctionCriteria.csv")

bfc$Target <- str_to_title(bfc$Target)

Bait <- bfc[bfc$Target %in% gsub(pattern = "_.*",
                                   replacement = "",
                                   x = colnames(df_HM)),]$Bait_Function[order(
                                     match(bfc[bfc$Target %in% gsub(pattern = "_.*",
                                                                    replacement = "",
                                                                    x = colnames(df_HM)),]$Target,
                                           gsub(pattern = "_.*",
                                                replacement = "",
                                                x = colnames(df_HM))))]

df_col <- data.frame(Bait)

rownames(df_col) <- colnames(df_HM)
 
cgf <- read.csv(file = "./07_InputFiles/00000000_Splicing_ComplexGeneralFunction.csv")

Complex <- cgf$Functionality
 
df_row <- data.frame(Complex)

rownames(df_row) <- rownames(df_HM)

display.brewer.pal(8, "Set2")

color <- brewer.pal(n = 8, name = "Set2")

Bait <- c(color[3],color[4],color[5])
names(Bait) <- c("Capping","Splicing", "Degradation")

Complex <- c(color[4],
                 "#b4b4b4")
names(Complex) <- c( "Splicing",
                         "Other")

anno_colors <- list(Bait = Bait,
                    Complex = Complex)

df_HM <- df_HM[,c(6,1:3,5,7,8,4)]

df_HM <- df_HM[c(7,11,5,13,4,3,10,2,9,1,8,16,15,18,19,12,6,14,17),]

filename <- paste0("./07_OutputFiles/",gsub(pattern = "-",x = Sys.Date(),replacement = ""),"_ComplexRatioHeatMap.pdf")

# gaps_col <- ncol(df_HM_ppi[,-1])

pheatmap(df_HM,
         col=colours,
         cluster_cols = F,
         cluster_rows = F,
         show_rownames = T,
         show_colnames = T,
         labels_col = str_to_title(gsub(pattern = "_.*",replacement = "",x = colnames(df_HM))),
         annotation_col = df_col,
         annotation_row = df_row,
         annotation_colors = anno_colors,
         fontsize = 16,
         border_color = "#b4b4b4",
         height = 9.5,width = 10,filename = filename
)


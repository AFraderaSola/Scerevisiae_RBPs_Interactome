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

groups <- c("Both", "PPI", "RDI")

ppi_targets <- c("Ist3", "Lsm2", "Dhh1", "Cbc2", "Snp1", "Sto1", "Hsh49")

rdi_targets <- c("Lsm2", "Cbc2", "Sto1", "Msl1")

#########################################
############# Script ####################
#########################################

i <- 2

for (i in 1:length(groups)) {
  
  if (groups[i] == "Both") {
    
    category_df_loop <- category_df
    
  }
  
  if (groups[i] == "PPI") {
    
    category_df_loop <- category_df[!category_df$category_loop == "RNAmediated",]
    
    category_df_loop <- category_df_loop[category_df_loop$target_RBP %in% ppi_targets,]
    
  }
  
  if (groups[i] == "RDI") {
    
    category_df_loop <- category_df[!category_df$category_loop == "PPI",]
    
    category_df_loop <- category_df_loop[category_df_loop$target_RBP %in% rdi_targets,]
    
  }
  
  
  
  cgf <- read.csv(file = paste0("./07_InputFiles/00000000_", functionality, "_ComplexGeneralFunction.csv"))
  
  filter_function <- KEGG[KEGG$Functionality == functionality,]$KEGG
  
  filter_id <- c()
  
  for (j in 1:length(filter_function)) {
    
    loop <- IP_KEGG[IP_KEGG$KEGG == filter_function[j],]$Gene
    
    filter_id <- c(filter_id, loop)
    
  }
  
  filter_id <- unique(filter_id)
  
  category_df_loop$target_RBP <- str_to_upper(category_df_loop$target_RBP)
  
  category_df_loop$genes <- str_to_upper(category_df_loop$genes)
  
  df <- category_df_loop[category_df_loop$target_RBP %in% filter_id,]
  
  complex_id <- unique(complexes$Complex)
  
  complex_n <- complexes %>%
    group_by(Complex) %>%
    summarise(n = n())
  
  complex_df <- c()
  
  for (k in 1:length(complex_id)) {
    
    complex_loop <- complexes[complexes$Complex == complex_id[k],]
    
    loop_overlap <- length(intersect(complex_loop$Name, unique(df$genes)))
    
    loop_n <- complex_n[complex_n$Complex == complex_id[k],]$n
    
    loop_ratio <- loop_overlap/loop_n
    
    loop_df <- data.frame(Complex = complex_id[k],
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
  
  filename <- paste0("./07_OutputFiles/",gsub(pattern = "-",x = Sys.Date(),replacement = ""),"_", groups[i],"_DensityRatioFull.pdf")
  
  ggsave(plot, filename = filename, width = 15, height = 7.5)
  
  plot <- ggplot(subset(complex_df, ratio > 0.5), aes(x=ratio)) +
    geom_density(fill="#69b3a2", color="#e9ecef", alpha=0.8)+
    ggtitle(paste0("n complexes = ", length(unique(subset(complex_df, ratio > 0.5)$Complex))))+
    theme_minimal()+
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=14,face="bold"))+
    theme(legend.text = element_text(size = 14),
          legend.title=element_text(size=12,face="bold"))
  
  filename <- paste0("./07_OutputFiles/",gsub(pattern = "-",x = Sys.Date(),replacement = ""),"_", groups[i],"_DensityFiltered.pdf")
  
  ggsave(plot, filename = filename, width = 15, height = 7.5)
  
  plot_df <- complex_df[complex_df$ratio > 0.5,]
  
  #### Full network complex ratios
  
  sf3a_nrow <- grep(pattern = "SF3a", plot_df$Complex)
  
  if (length(sf3a_nrow) > 0) {
    
    plot_df$Complex[sf3a_nrow] <- "SF3a complex" #### Complex label abbreviation
    
  }
  
  plot_df <- plot_df %>% arrange(desc(n_complex),ratio)
  
  plot_df$ComplexType <- "Radnom"
  
  splicing_complexes <- c("U4/U6 x U5 tri-snRNP complex", "commitment complex", "U2 snRNP complex", "U1 snRNP complex",
                          "U5 snRNP complex", "Prp19-associated complex", "U6 snRNP complex", "SF3b complex",
                          "SF3a complex", "RES complex")
  
  other_complexes <- c("THO complex", "eIF4F", "RENT complex", "transcription export complex",
                          "TRAMP complex (Air1p)", "TRAMP complex (Air2p)", "Decapping Enzyme Complex", "Gcd10p/Gcd14p complex",
                          "Ubp3p/Bre5p complex")
  
  plot_df[plot_df$Complex %in% splicing_complexes, ]$ComplexType <- "Splicing"
  
  plot_df[plot_df$Complex %in% other_complexes, ]$ComplexType <- "Other"
  
  plot_df <- plot_df %>% arrange(ComplexType, desc(n_complex),ratio)
  
  plot_df$Complex <- factor(x = plot_df$Complex, levels = plot_df$Complex)
  
  x_n <- length(plot_df[plot_df$ComplexType == "Other", ]$ComplexType)+0.5
  
  plot <- ggplot(plot_df) +
    # geom_segment(aes(x=Complex, xend=Complex, y=n_complex, yend=n_network), color="#666666") +
    # geom_point(data = subset(plot_df, ratio < 1), aes(x=Complex, y=n_complex), color="#FDBF6F", size=8,) +
    # geom_point(data = subset(plot_df, ratio < 1), aes(x=Complex, y=n_network), color="#CAB2D6", size=8) +
    # geom_point(data = subset(plot_df, ratio == 1), aes(x=Complex, y=n_complex), 
    #            color="#CAB2D6", fill = "#FDBF6F", shape = 21, size=8, stroke =1) +
    # geom_text(aes(x=Complex, y=n_complex+2.5), label = round(plot_df$ratio, digits = 2), size = 10)+
    geom_textsegment(aes(x=Complex, xend=Complex, y=n_complex, yend=n_network), color="#2e2e2e", label = round(plot_df$ratio,digits = 2), size = 15) +
    geom_point(aes(x=Complex, y=n_complex), color="#FDBF6F", size=18,) +
    geom_point(aes(x=Complex, y=n_network), color="#CAB2D6", size=18) +
    geom_point(data = subset(plot_df, ratio == 1), aes(x=Complex, y=n_complex), 
               color="#CAB2D6", fill = "#FDBF6F", shape = 21, size=18, stroke =1) +
    geom_text(data = subset(plot_df, ratio == 1), aes(x=Complex, y=n_complex), label = "1", size = 15)+
    annotate(
      xmin = c(-Inf, x_n),
      xmax = c(x_n, Inf),
      ymin = 29.5,
      ymax = 30,
      geom = "rect",
      fill = c("#b4b4b4", "#E78AC3")
    )+
    coord_flip()+
    scale_y_continuous(breaks = seq(0, 30, by = 5))+
    theme_minimal() +
    theme(axis.text=element_text(size=60),
          axis.title=element_text(size=72,face="bold"))+
    theme(legend.text = element_text(size = 72),
          legend.title=element_text(size=60,face="bold"))+
    theme(legend.key.size = unit(3, 'cm'))+
    theme(legend.position = "top")+
    xlab("Protein complex") +
    ylab("Number of proteins")
  
  filename <- paste0("./07_OutputFiles/",gsub(pattern = "-",x = Sys.Date(),replacement = ""),"_", groups[i],"_ComplexRatio.pdf")
  
  ggsave(plot, filename = filename, width = 30, height = 28.5)
  
  #### Bait complex ratios
  
  # Barplot
  
  targets <- unique(df$target_RBP)
  
  complex_df <- c()
  
  for (l in 1:length(complex_id)) {
    
    for (m in 1:length(targets)) {
      
      complex_loop <- complexes[complexes$Complex == complex_id[l],]
      
      loop_targets <- df[df$target_RBP == targets[m],]
      
      loop_overlap <- length(intersect(complex_loop$Name, loop_targets$genes))
      
      loop_n <- complex_n[complex_n$Complex == complex_id[l],]$n
      
      loop_ratio <- loop_overlap/loop_n
      
      loop_df <- data.frame(Target = targets[m],
                            Complex = complex_id[l],
                            n_complex = loop_n,
                            n_network = loop_overlap,
                            ratio = loop_ratio)
      
      complex_df <- rbind(complex_df, loop_df)
      
    }
  }
  
  plot_df$Complex <- as.character(plot_df$Complex)
  
  sf3a_nrow <- grep(pattern = "SF3a", plot_df$Complex)
  
  if (length(sf3a_nrow) > 0) {
    
    plot_df$Complex[sf3a_nrow] <- "Prp9p/Prp11p/Prp21p complex(SF3a complex)" #### Complex label abbreviation
    
  }
  
  plot_df <- complex_df[complex_df$Complex %in% plot_df$Complex,]
  
  plot_df <- plot_df[plot_df$ratio > 0, ]
  
  sf3a_nrow <- grep(pattern = "SF3a", plot_df$Complex)
  
  if (length(sf3a_nrow) > 0) {
    
    plot_df$Complex[sf3a_nrow] <- "SF3a complex" #### Complex label abbreviation
    
  }
  
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
  
  filename <- paste0("./07_OutputFiles/",gsub(pattern = "-",x = Sys.Date(),replacement = ""),"_", groups[i],"_TargetRatio.pdf")
  
  ggsave(plot, filename = filename, width = 15, height = 15)
  
  # Heatmap
  
  plot_df$Target <- str_to_title(plot_df$Target)
  
  targets <- unique(plot_df$Target)
  
  df_list <- c()
  
  for (n in 1:length(targets)) {
    
    loop_df <- plot_df[plot_df$Target == targets [n], ][,c(2,5)]
    
    colnames(loop_df)[2] <- paste0(targets[n], "_ratio")
    
    loop_df <- list(loop_df)
    
    names(loop_df) <- targets[n]
    
    df_list <- c(df_list, loop_df)
    
  }
  
  df_HM <- df_list[[1]]
  
  for (o in 2:length(df_list)) {
    
    df_HM <- full_join(x = df_HM,y = df_list[[o]], by = "Complex")
    
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
  
  cgf <- cgf[cgf$Complex %in% rownames(df_HM),]
  
  cgf <- cgf[match(rownames(df_HM),cgf$Complex),]

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
  
  if (groups[i] == "Both") {
    
    df_HM <- df_HM[,c(6,1:3,5,7,8,4)]
    
    df_HM <- df_HM[c(7,11,5,13,4,3,10,2,9,1,8,16,15,18,19,12,6,14,17),]
    
  }
  
  if (groups[i] == "PPI") {

    df_HM <- df_HM[,c(5,1,2,3,6,7,4)]

    df_HM <- df_HM[c(5,4,10,3,7,2,6,1,11,9,8),]

  }

  if (groups[i] == "RDI") {

    df_HM <- df_HM[,c(4,1:3)]

    df_HM <- df_HM[c(4,3,2,1,6,5,8,7),]

  }
  
  filename <- paste0("./07_OutputFiles/",gsub(pattern = "-",x = Sys.Date(),replacement = ""),"_", groups[i],"_ComplexRatioHeatMap.pdf")
  
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
           fontsize = 54,
           border_color = "#b4b4b4",
           height = 28.5,width = 30,filename = filename
  )
  
}


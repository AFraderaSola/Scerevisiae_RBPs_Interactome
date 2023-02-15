#############################################
############# Libraries #####################
#############################################

library(readr)
library(ggplot2)
library(dplyr)
library(readxl)
library(scales)
library(ggthemes)
library(RColorBrewer)
library(tidyverse)

#############################################
############ Variables ######################
#############################################

values <- "LFQ"

#############################################
############### Script ######################
#############################################

files_KO <- list.files(path ="./00_KOsResultsFiles/", pattern = "Enriched.*\\.csv")

files_PPI <- list.files(path = paste0("./00_IPsResultsFiles/",values,"/"), pattern = "Enriched_RNase_vs_WT")

files_RNAmediated <- list.files(path = paste0("./00_IPsResultsFiles/",values,"/"), pattern = "Enriched_IP_vs_RNase")

stock <- read_excel("../00_KOScreeningMasterFile.xlsx")

df <- c()

for (i in 1:length(files_KO)) {
  
  loop_KO <- read_csv(paste0("./00_KOsResultsFiles/", files_KO[i]))
  
  loop_KO_up <- loop_KO[,grep(pattern = "difference|Protein.IDs",x = colnames(loop_KO))]  
  
  loop_KO_up <- loop_KO[unlist(list(loop_KO_up[,2]),use.names = F) > 0,]
  
  loop_KO_up_IDs <- loop_KO_up$Protein.IDs
  
  loop_PPI <- read_csv(paste0("./00_IPsResultsFiles/",values,"/", files_PPI[i]))
  
  loop_PPI_IDs <- loop_PPI$Protein.IDs
  
  loop_PPI_IDs <- gsub(pattern = ";.*",replacement = "",x = loop_PPI_IDs)
  
  loop_RNAmediated <- read_csv(paste0("./00_IPsResultsFiles/",values,"/", files_RNAmediated[i]))
  
  loop_RNAmediated_IDs <- loop_RNAmediated$Protein.IDs
  
  loop_RNAmediated_IDs <- gsub(pattern = ";.*",replacement = "",x = loop_RNAmediated_IDs)
  
  interesect_PPIandRNAmediated <- intersect(loop_PPI_IDs, loop_RNAmediated_IDs)
  
  loop_PPI_IDs <- loop_PPI_IDs[!loop_PPI_IDs %in% interesect_PPIandRNAmediated]
  
  loop_RNAmediated_IDs <- loop_RNAmediated_IDs[!loop_RNAmediated_IDs %in% interesect_PPIandRNAmediated]
  
  intersect_KOandPPI <- intersect(loop_KO_up_IDs, loop_PPI_IDs)
  
  intersect_KOandRNAmediated <- intersect(loop_KO_up_IDs, loop_RNAmediated_IDs)
  
  intersect_KOandboth <- intersect(loop_KO_up_IDs, interesect_PPIandRNAmediated)
  
  loop_KO_up_Unique <- loop_KO_up[!loop_KO_up$Protein.IDs %in% c(intersect_KOandPPI,intersect_KOandRNAmediated,intersect_KOandboth),]
  
  loop_KO_up_Unique <- loop_KO_up_Unique %>%
    summarize(count = n())
  
  loop_KO_up_Unique$Set <- "Unique"
  
  loop_KO_up_Unique$Direction <- "UP"
  
  loop_KO_up_Unique$ID <- gsub(pattern = "_.*",replacement = "",files_KO[i])
  
  loop_KO_up_both <- loop_KO_up[loop_KO_up$Protein.IDs %in% c(intersect_KOandboth),]
  
  loop_KO_up_both <- loop_KO_up_both %>%
    summarize(count = n())
  
  loop_KO_up_both$Set <- "Both"
  
  loop_KO_up_both$Direction <- "UP"
  
  loop_KO_up_both$ID <- gsub(pattern = "_.*",replacement = "",files_KO[i])
  
  loop_KO_up_RNAmediated <- loop_KO_up[loop_KO_up$Protein.IDs %in% c(intersect_KOandRNAmediated),]
  
  loop_KO_up_RNAmediated <- loop_KO_up_RNAmediated %>%
    summarize(count = n())
  
  loop_KO_up_RNAmediated$Set <- "RNAmediated"
  
  loop_KO_up_RNAmediated$Direction <- "UP"
  
  loop_KO_up_RNAmediated$ID <- gsub(pattern = "_.*",replacement = "",files_KO[i])
  
  loop_KO_up_PPI <- loop_KO_up[loop_KO_up$Protein.IDs %in% c(intersect_KOandPPI),]
  
  loop_KO_up_PPI <- loop_KO_up_PPI %>%
    summarize(count = n())
  
  loop_KO_up_PPI$Set <- "PPI"
  
  loop_KO_up_PPI$Direction <- "UP"
  
  loop_KO_up_PPI$ID <- gsub(pattern = "_.*",replacement = "",files_KO[i])
  
  final_KO_up <- rbind(loop_KO_up_Unique,loop_KO_up_both,loop_KO_up_RNAmediated,loop_KO_up_PPI)
  
  print(sum(final_KO_up$count) == nrow(loop_KO_up))
  
  loop_KO_down <- loop_KO[,grep(pattern = "difference|Protein.IDs",x = colnames(loop_KO))]  
  
  loop_KO_down <- loop_KO[unlist(list(loop_KO_down[,2]),use.names = F) < 0,]
  
  loop_KO_down_IDs <- loop_KO_down$Protein.IDs
  
  loop_PPI <- read_csv(paste0("./00_IPsResultsFiles/",values,"/", files_PPI[i]))
  
  loop_PPI_IDs <- loop_PPI$Protein.IDs
  
  loop_PPI_IDs <- gsub(pattern = ";.*",replacement = "",x = loop_PPI_IDs)
  
  loop_RNAmediated <- read_csv(paste0("./00_IPsResultsFiles/",values,"/", files_RNAmediated[i]))
  
  loop_RNAmediated_IDs <- loop_RNAmediated$Protein.IDs
  
  loop_RNAmediated_IDs <- gsub(pattern = ";.*",replacement = "",x = loop_RNAmediated_IDs)
  
  interesect_PPIandRNAmediated <- intersect(loop_PPI_IDs, loop_RNAmediated_IDs)
  
  loop_PPI_IDs <- loop_PPI_IDs[!loop_PPI_IDs %in% interesect_PPIandRNAmediated]
  
  loop_RNAmediated_IDs <- loop_RNAmediated_IDs[!loop_RNAmediated_IDs %in% interesect_PPIandRNAmediated]
  
  intersect_KOandPPI <- intersect(loop_KO_down_IDs, loop_PPI_IDs)
  
  intersect_KOandRNAmediated <- intersect(loop_KO_down_IDs, loop_RNAmediated_IDs)
  
  intersect_KOandboth <- intersect(loop_KO_down_IDs, interesect_PPIandRNAmediated)
  
  loop_KO_down_Unique <- loop_KO_down[!loop_KO_down$Protein.IDs %in% c(intersect_KOandPPI,intersect_KOandRNAmediated,intersect_KOandboth),]
  
  loop_KO_down_Unique <- loop_KO_down_Unique %>%
    summarize(count = n())
  
  loop_KO_down_Unique$Set <- "Unique"
  
  loop_KO_down_Unique$Direction <- "DOWN"
  
  loop_KO_down_Unique$ID <- gsub(pattern = "_.*",replacement = "",files_KO[i])
  
  loop_KO_down_both <- loop_KO_down[loop_KO_down$Protein.IDs %in% c(intersect_KOandboth),]
  
  loop_KO_down_both <- loop_KO_down_both %>%
    summarize(count = n())
  
  loop_KO_down_both$Set <- "Both"
  
  loop_KO_down_both$Direction <- "DOWN"
  
  loop_KO_down_both$ID <- gsub(pattern = "_.*",replacement = "",files_KO[i])
  
  loop_KO_down_RNAmediated <- loop_KO_down[loop_KO_down$Protein.IDs %in% c(intersect_KOandRNAmediated),]
  
  loop_KO_down_RNAmediated <- loop_KO_down_RNAmediated %>%
    summarize(count = n())
  
  loop_KO_down_RNAmediated$Set <- "RNAmediated"
  
  loop_KO_down_RNAmediated$Direction <- "DOWN"
  
  loop_KO_down_RNAmediated$ID <- gsub(pattern = "_.*",replacement = "",files_KO[i])
  
  loop_KO_down_PPI <- loop_KO_down[loop_KO_down$Protein.IDs %in% c(intersect_KOandPPI),]
  
  loop_KO_down_PPI <- loop_KO_down_PPI %>%
    summarize(count = n())
  
  loop_KO_down_PPI$Set <- "PPI"
  
  loop_KO_down_PPI$Direction <- "DOWN"
  
  loop_KO_down_PPI$ID <- gsub(pattern = "_.*",replacement = "",files_KO[i])
  
  final_KO_down <- rbind(loop_KO_down_Unique,loop_KO_down_both,loop_KO_down_RNAmediated,loop_KO_down_PPI)
  
  print(sum(final_KO_down$count) == nrow(loop_KO_down))
  
  final_KO_down$count <- final_KO_down$count * -1
  
  loop_df <- rbind(final_KO_up,final_KO_down)
  
  print(sum(abs(loop_df$count)) == nrow(loop_KO))
  
  df <- rbind(df, loop_df)
  
}

for (i in 1:nrow(df)) {
  
  ID <- df$ID[i]
  
  gene <- stock[stock$Stock_Name == ID,]$Gene_Name
  
  gene <- tolower(gene)
  
  nID <- paste0(gene, "\u0394")
  
  df$ID[i] <- nID
  
}

counts <- df

counts$count <- abs(counts$count)

counts <- aggregate(counts$count, by=list(counts$ID), FUN = sum)

counts <- counts[order(counts$x),]

df$ID <- factor(x = df$ID, levels = rev(counts$Group.1))

df <- df[order(df$ID),]

df <- df[!df$ID == "cbc2Δ",]

df <- df[!df$ID == "msl1Δ",]

my_labels <- c()

for (i in 1:nrow(df)) {
  
  exp <- as.character(df$ID[i])
  
  gene <- gsub(pattern = " .*",replacement = "",x = exp)
  
  gene <- substr(gene,1,nchar(gene)-1)
  
  # expression <- bquote(.(gene)*Delta~"vs."~"WT")
  
  expression <- expr(paste(italic(!!gene),"-",Delta,sep = ""))
  
  my_labels <- c(my_labels, expression)
  
}

my_labels <- my_labels[seq(1, length(my_labels), 8)]

colors <- brewer.pal(5, "Dark2")

colors <- c(colors[5], colors[2], "#b5b5b5")

# Exclude not validated strains 

# df <- df[!df$ID == "cbc2Δ vs. WT",]
# 
# df <- df[!df$ID == "msl1Δ vs. WT",]

df <- df[!df$count == 0,]

plot <- ggplot(df, aes(fill=Set, y = count , x= ID)) +
  geom_bar(position="stack", stat="identity")+
  scale_fill_manual(labels = c("PPI",
                               "RDI",
                               "None"),
                    values = colors)+
  labs(fill = "Overlap with interactome")+
  scale_y_continuous(breaks = pretty(df$count), labels = abs(pretty(df$count)))+
  scale_x_discrete(labels = my_labels)+
  geom_hline(yintercept=0, color= "#000000") +
  coord_flip()+
  theme_minimal()+
  annotate('text', label="Upregulated", x= 13, y= 60, size = 7) +
  annotate('text', label="Downregulated", x= 13, y= -75, size = 7) +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16,face="bold"))+
  theme(legend.title = element_text(size = 16, face="bold"))+
  theme(legend.text = element_text(size = 14))+
  theme(legend.key.size = unit(1, 'cm'))+
  theme(legend.position = "top")+
  xlab("Comparison to WT")+
  ylab("Enriched interactors")+
  theme(plot.title = element_text(size = 18, face = "bold"))

filename <- paste0("./01_OutputFiles/",gsub("-","",Sys.Date()),"_EnrichedIDs_overlap_",values,"_IPs.pdf") 

ggsave(plot = plot, filename = filename, width = 10, height = 4.5)

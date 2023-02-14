################################
######### Libraries ############
################################

library(readr)
library(tidyverse)
library(grid)
library(gridExtra)
library(UpSetR)
library(RColorBrewer)
library(ggplot2)
library(pheatmap)
library(viridis)
library(readxl)

################################
########## Variables ###########
################################

# Quantified proteins on YAF000

files_screen <- list.files(path = "./00_IPsResultsFiles/", pattern = "*.csv")

screen <- read_excel("../../../00_IPScreeningMasterFile.xlsx")

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
    select(my_label,
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

genes <- unique(genes)

IDlist <- c()

df <- c()

for (i in 1:length(genes)) {
  
  list_loop <- protein_list[grep(pattern = genes[i], names(protein_list))]
  
  overlap <- intersect(list_loop[[1]], list_loop[[2]])
  
  PPI <- list_loop[grep(pattern = "PPI", names(list_loop))][[1]][!list_loop[grep(pattern = "PPI", names(list_loop))][[1]] %in% overlap]
  
  RNAmediated <- list_loop[grep(pattern = "RNAmediated", names(list_loop))][[1]][!list_loop[grep(pattern = "RNAmediated", names(list_loop))][[1]] %in% overlap]
  
  IDs <- list(overlap, PPI, RNAmediated)
  
  names(IDs) <- c("overlap", "PPI", "RNAmediated")
  
  IDs <- list(IDs)
  
  names(IDs) <- genes[i]
  
  IDlist <- c(IDlist, IDs)
  
  overlap <- length(overlap)
  
  PPI <- length(PPI)
  
  RNAmediated <- length(RNAmediated)
  
  count <- c(overlap, PPI, RNAmediated)
  
  gene <- rep(genes[i],3)
  
  set <- c("Intersection","PPI","RNAmediated")
  
  df_loop <- data.frame(gene,count,set)
  
  df <- rbind(df,df_loop)
  
  print(df_loop$count[1] + df_loop$count[3] == length(list_loop[[1]]))
  
  print(df_loop$count[1] + df_loop$count[2] == length(list_loop[[2]]))
  
}

df$gene <- str_to_title(df$gene)

totalcount <- aggregate(df$count, by=list(Category=df$gene), FUN=sum)

totalcount <- totalcount[order(totalcount$x),]

df$gene <- factor(x = df$gene, levels = rev(totalcount$Category))

colors <- c("#D95F02","#666666", "#66A61E")

df$set <- factor(x = df$set, levels = c("RNAmediated", "Intersection", "PPI"))

df_pr <- df %>%
          group_by(gene) %>%
          mutate(countT= sum(count)) %>%
          group_by(set, add=TRUE) %>%
          mutate(per=100*count/countT)

df_pr_in <- df_pr[df_pr$set == "Intersection",]

plot <- ggplot(df, aes(fill=set, y = count , x= gene)) +
  geom_bar(position="stack", stat="identity")+
  scale_fill_manual(labels = c(
                               "RNA mediated interactions",
                               "Intersection",
                               "Protein-protein interactions"
                               ),
                    values = colors)+
  labs(fill = "")+
  coord_flip()+
  scale_y_continuous(expand=c(0,2))+
  theme_minimal()+
  theme(axis.text.y = element_text(size=9),
        axis.text.x = element_text(size=16),
        axis.title=element_text(size=18,face="bold"))+
  theme(legend.title = element_text(size = 16, face="bold"))+
  theme(legend.text = element_text(size = 14))+
  theme(legend.key.size = unit(1, 'cm'))+
  theme(legend.position = "top")+
  xlab("Bait")+
  ylab("Enriched interactors")+
  theme(plot.title = element_text(size = 18, face = "bold"))

plot

filename <- paste0("./01_OutputFiles/",gsub("-","",Sys.Date()),"_OverlapPPIandRNAmediated.pdf") 

ggsave(plot = plot, filename = filename, useDingbats=FALSE, width = 20, height = 6.3)

my.l <- IDlist

new.l <- rapply(my.l, function(x) paste(x, collapse = "|"), how = "replace")

dt <- data.table::rbindlist(new.l)

dt$Target <- names(new.l)

dt <- dt %>% tidyr::separate_rows(overlap, sep = "\\|")

dt <- dt %>% tidyr::separate_rows(PPI, sep = "\\|")

dt <- dt %>% tidyr::separate_rows(RNAmediated, sep = "\\|")

dt <- data.frame(lapply(dt, function(x) {gsub("-TAP", "", x)}))

# write.csv(x = dt,file = paste0("./01_OutputFiles/",gsub("-","",Sys.Date()),"_NetworkTable.csv"),row.names = F)


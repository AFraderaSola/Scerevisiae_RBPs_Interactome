###################################
######### Libraries ###############
###################################

library(readr)
library(ggplot2)
library(ggsignif)
library(ggforce)
library(viridis)
library(tidyverse)
# library(xlsx)
library(readxl)
library(RColorBrewer)
library(ggpubr)

###################################
######### Variables ###############
###################################

census <- read.csv("./01_InputFiles/SC_RBPs_Census.csv")

census <- census$ID

files <- list.files(path = "./00_IPsResultsFiles/", pattern = "_WT.csv")

direction <- "Up"

screen <- read_excel("../../../00_IPScreeningMasterFile.xlsx")

# color_intercept <- c("#D95F02","#66A61E")

grid_colour <- "#66A61E"

###################################
########### Script ################
###################################

enriched <- c()

if (direction == "Up") {
  
  for (i in 1:length(files)) {
    
    loop <- read.csv(paste0("./00_IPsResultsFiles/",files[i]))
    
    loop <- loop[loop$difference_RNase_WT > 0,]
    
    loop <- loop$Protein.IDs
    
    loop <- list(loop)
    
    enriched <- c(enriched, loop)
    
  }
  
}else{
  
  for (i in 1:length(files)) {
    
    loop <- read.csv(paste0("./00_IPsResultsFiles/", files[i]))
    
    loop <- loop[loop$difference_RNase_WT < 0,]
    
    loop <- loop$Protein.IDs
    
    loop <- list(loop)
    
    enriched <- c(enriched, loop)
  
  }
}


files <- list.files(path = "../../../../00_Candidates/00_CandidateInteractions/SGD/", pattern = "_toHighlight.txt")

# sgd <- c()

# for (i in 1:length(files)) {
#   
#   loop <- read_delim(paste0("../../../../00_Candidates/00_CandidateInteractions/SGD/", files[i]), "\t", escape_double = FALSE, trim_ws = TRUE,col_names = F)
#   
#   loop <- loop$X1
#   
#   loop <- list(loop)
#   
#   sgd <- c(sgd, loop)
#   
# }

files <- list.files(path = "../../../../00_Candidates/00_CandidateInteractions/BioGRID/", pattern = "_BiogridInteractions.txt")

biogrid <- c()

for (i in 1:length(files)) {
  
  loop <- read_delim(paste0("../../../../00_Candidates/00_CandidateInteractions/BioGRID/", files[i]), "\t", escape_double = FALSE, trim_ws = TRUE,col_names = T)
  
  loop <- loop[loop$`Experimental System Type` == "physical",]
  
  loop <- unique(c(loop$`Systematic Name Interactor A`, loop$`Systematic Name Interactor B`))
  
  loop <- list(loop)
  
  biogrid <- c(biogrid, loop)
  
}

# overlap <- c()
# 
# for (i in 1:length(enriched)) {
#   
#   loop <- intersect(unlist(enriched[i]), unlist(sgd[i]))
#   
#   loop <- list(loop)
#   
#   overlap <- c(overlap, loop)
#   
# }

overlap_census <- c()

for (i in 1:length(enriched)) {
  
  loop <- intersect(unlist(enriched[i]), census)
  
  loop <- list(loop)
  
  overlap_census <- c(overlap_census, loop)
  
}

overlap_biogrid <- c()

for (i in 1:length(enriched)) {
  
  loop <- intersect(unlist(enriched[i]), unlist(biogrid[i]))
  
  loop <- list(loop)
  
  overlap_biogrid <- c(overlap_biogrid, loop)
  
}

IP <- gsub(pattern = "_.*",replacement = "", x = files)

df <- c()

for (i in 1:length(enriched)) {
  
  # df2 <- data.frame(unlist(enriched[i]), rep(FALSE, length(unlist(enriched[i]))), rep(FALSE, length(unlist(enriched[i]))))
  
  df2 <- data.frame(unlist(enriched[i]), rep(FALSE, length(unlist(enriched[i]))))
  
  # colnames(df2) <- c("Protein ID", "Included SGD", "Included BioGRID")
  
  colnames(df2) <- c("Protein ID", "Included BioGRID")
  
  # df2$`Included SGD` <- df2$`Protein ID` %in% unlist(overlap[i])  
  
  df2$`Included BioGRID` <- df2$`Protein ID` %in% unlist(overlap_biogrid[i])  
  
  df2$`Included BioGRID`[which(df2$`Included BioGRID` == "TRUE")] <- "TRUE.TRUE"
  
  df2$`Included BioGRID`[which(df2$`Included BioGRID` == "FALSE")] <- "FALSE.FALSE"
  
  # df2$`Included in data set` <- as.character(interaction(df2$`Included SGD`, df2$`Included BioGRID`))
  
  # df2 <- df2 %>%
  #   group_by(df2$`Included in data set`) %>%
  #   summarize(count = n())%>%
  #   mutate(perc = count/sum(count))
  
  df2 <- df2 %>%
    group_by(df2$`Included BioGRID`) %>%
    summarize(count = n())%>%
    mutate(perc = count/sum(count))
  
  df2$DataSet <- "Included at interactor database"
  
  colnames(df2)[c(1,3)] <- c("Included in data set", "Percentage of Protein IDs")
  
  df2$IP <- IP[i]
  
  df3 <- data.frame(unlist(enriched[i]), rep(FALSE, length(unlist(enriched[i]))))
  
  df3$rep.FALSE..length.unlist.enriched.i.... <- df3$unlist.enriched.i.. %in% unlist(overlap_census[i])                  
  
  colnames(df3) <- c("Protein ID", "Included in data set")
  
  df3 <- df3 %>%
    group_by(df3$`Included in data set`) %>%
    summarize(count = n())%>% 
    mutate(perc = count/sum(count))
  
  df3$DataSet <- "Included at RBP census"
  
  colnames(df3)[c(1,3)] <- c("Included in data set", "Percentage of Protein IDs")
  
  df3$IP <- IP[i]
  
  loop <- rbind(df2,df3)
  
  loop <- list(loop)
  
  df <- c(df, loop)
  
}

plot_df <- c()

for (i in 1:length(df)) {
  
  loop <- df[[i]]
  
  plot_df <- rbind(plot_df, loop)
  
}

plot_df$DataSet <- factor(plot_df$DataSet, levels = c("Included at interactor database",
                                                      "Included at RBP census"))

plot_df$Gene <- plot_df$IP

for (i in 1:nrow(plot_df)) {
  
  stock <- plot_df$IP[i]
  
  gene <- screen[screen$Stock_Name == stock,]$Gene_Name
  
  plot_df$Gene[i] <- unique(gene)
  
}

plot_df$Gene <- str_to_title(plot_df$Gene)

plot_df$Gene <- factor(x = plot_df$Gene, levels = unique(plot_df$Gene))

## GeneName Figure:

### Census

df_A <- plot_df[plot_df$DataSet == "Included at RBP census",]

# order <- df_A %>% group_by(Gene) %>% summarise(Count = sum(count))
# 
# order <- order[order(-order$Count),]
# 
# write.csv(order, file = paste0("./01_OutputFiles/", gsub(pattern = "-",replacement = "",x = Sys.Date()),"_",direction,"Regulated_CountTable.csv"),row.names = F)

order <- df_A[df_A$`Included in data set` == T,]

order <- order[order(-order$`Percentage of Protein IDs`),]

df_A$Gene <- factor(df_A$Gene, levels = unique(order$Gene))

filename <- paste0("./01_OutputFiles/",gsub("-","",Sys.Date()),"_",direction,"Regulated_CensusCounts.pdf") 

# colors <- viridis(end = 0.9, n = 2,option = "A")

# colors <- brewer.pal(n = 8, name = "Dark2")[c(1,2)]

colors <- brewer.pal(n = 4, name = "Greys")

colors <- c(colors[4],"#000000")

plot_A <- ggplot(df_A, aes(fill=`Included in data set`, y = count , x= Gene)) +
                # geom_bar(position="stack", stat="identity")+
                geom_bar(position="fill", stat="identity")+
                # scale_fill_manual(labels = c("Not included", expression(paste("Hentze, M. ", italic("et al."), " 2018"))),
                #                    values = colors)+
                scale_fill_manual(labels = c("Not included", "Included"),
                    values = colors)+
                geom_hline(yintercept=0.7, linetype="dashed", size = 1, color = "#ff8aad" )+
                geom_hline(yintercept=0.5, linetype="dashed", size = 1, color = "#990000" )+
                labs(fill = "RBP census")+
                theme_minimal()+
                theme(axis.text.x = element_text(angle = 45, hjust = 1))+
                theme(axis.text=element_text(size=16),
                      axis.title=element_text(size=18,face="bold"))+
                theme(legend.title = element_text(size = 16, face="bold"))+
                theme(legend.text = element_text(size = 14))+
                theme(legend.key.size = unit(1, 'cm'))+
                theme(legend.position = "top")+
                # theme(panel.grid.minor = element_line(linetype = 'solid', colour = grid_colour))+
                theme(panel.grid.major = element_line(linetype = 'solid', colour = grid_colour))+
                xlab("Bait")+
                ylab("Enriched interactors (%)")

plot_A

ggsave(plot = plot_A, filename = filename, useDingbats=FALSE, width = 15, height = 15)

### SGD

df_B <- plot_df[plot_df$DataSet == "Included at interactor database",]

order <- df_B[df_B$`Included in data set` == "TRUE.TRUE",]

order <- order[order(-order$`Percentage of Protein IDs`),]

diff <- setdiff(as.character(unique(df_B$Gene)),as.character(order$Gene))

df_B$Gene <- factor(df_B$Gene, levels = c(as.character(order$Gene),diff))

filename <- paste0("./01_OutputFiles/", gsub("-","",Sys.Date()),"_",direction,"Regulated_InteractionDBsCounts.pdf") 

# colors <- viridis(end = 0.7 , n = 3,option = "A")

# colors <- brewer.pal(n = 8, name = "Dark2")[c(1,2)]

colors <- brewer.pal(n = 4, name = "Greys")

colors <- c(colors[4],colors[2])

plot_B <- ggplot(df_B, aes(fill=`Included in data set`, y = count , x= Gene)) +
  # geom_bar(position="stack", stat="identity")+
  geom_bar(position="fill", stat="identity")+
  scale_fill_manual(labels = c("Not included", 
                               "Included"),
                    values = colors)+
  labs(fill = "BioGRID")+
  geom_hline(yintercept=0.5, linetype="dashed", size = 1, color = "#ff8aad" )+
  geom_hline(yintercept=0.2, linetype="dashed", size = 1, color = "#990000" )+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16,face="bold"))+
  theme(legend.title = element_text(size = 16, face="bold"))+
  theme(legend.text = element_text(size = 14))+
  theme(legend.key.size = unit(1, 'cm'))+
  theme(legend.position = "top")+
  # theme(panel.grid.minor = element_line(linetype = 'solid', colour = grid_colour))+
  theme(panel.grid.major = element_line(linetype = 'solid', colour = grid_colour))+
  xlab("Bait")+
  ylab("Enriched interactors (%)")

plot_B

ggsave(plot = plot_B, filename = filename, useDingbats=FALSE, width = 15, height = 15)

### Mirrored

order <- df_B %>% group_by(Gene) %>% summarise(Count = sum(count))

order <- order[order(-order$Count),]

df_B$count <- df_B$count*-1

df_C <- rbind(df_A, df_B)

df_C$Gene <- factor(df_C$Gene, levels = unique(order$Gene))

# colorsA <- viridis(end = 0.9, n = 2,option = "A")
# colorsB  <- viridis(end = 0.7 , n = 3,option = "A")
# 
# colors <- c(colorsA, colorsB)

colors <- brewer.pal(n = 4, name = "Greys")

# colors <- c(colors[4],  colors[2:4], "#000000")

colors <- c(colors[2],  colors[c(2,4)], "#000000")

# df_C$`Included in data set` <- factor(x = df_C$`Included in data set`, levels = c("FALSE", "TRUE.TRUE", "FALSE.TRUE","FALSE.FALSE", "TRUE"))

df_C$`Included in data set` <- factor(x = df_C$`Included in data set`, levels = c("FALSE", "FALSE.FALSE","TRUE.TRUE", "TRUE"))

yannot <- max(df_C$count)/1.25

# plot_C <- ggplot(df_C, aes(fill=`Included in data set`, y = count , x= Gene)) +
#   geom_bar(position="stack", stat="identity")+
#   scale_fill_manual(labels = c("Not included",
#                                "BioGRID & SGD",
#                                "BioGRID",
#                                "Not included",
#                                expression(paste("Hentze, M. ", italic("et al."), " 2018"))),
#                     values = colors)+
#   labs(fill = "")+
#   annotate('text', label="RBP census", x= 39, y= yannot, size = 8) +
#   annotate('text', label="Interaction databases", x= 39, y= -yannot, size = 8) +
#   geom_hline(yintercept=0, color= grid_colour) +
#   scale_y_continuous(breaks = pretty(df_C$count), labels = abs(pretty(df_C$count)))+
#   coord_flip()+
#   theme_minimal()+
#   theme(axis.text=element_text(size=14),
#         axis.title=element_text(size=16,face="bold"))+
#   theme(legend.title = element_text(size = 16, face="bold"))+
#   theme(legend.text = element_text(size = 14))+
#   theme(legend.key.size = unit(1, 'cm'))+
#   theme(legend.position = "top")+
#   xlab("IP - Protein ID")+
#   ylab("Enriched interactors")+
#   theme(plot.title = element_text(size = 18, face = "bold"))+
#   theme(
#         # panel.grid.major = element_line(linetype = 'solid', colour = "#66A61E"), 
#         panel.grid.minor = element_line(linetype = 'solid',colour = grid_colour)
#         )

plot_C <- ggplot(df_C, aes(fill=`Included in data set`, y = count , x= Gene)) +
  geom_bar(position="stack", stat="identity")+
  scale_fill_manual(labels = c("Not included",
                               "Not included",
                               "BioGRID",
                               "RBP census"),
                    values = colors)+
  labs(fill = "")+
  annotate('text', label="RBP census", x= 39, y= yannot, size = 7) +
  annotate('text', label="Interaction database", x= 39, y= -yannot, size = 7) +
  geom_hline(yintercept=0, color= grid_colour) +
  scale_y_continuous(breaks = pretty(df_C$count), labels = abs(pretty(df_C$count)))+
  coord_flip()+
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
  theme(plot.title = element_text(size = 18, face = "bold"))+
  theme(
    # panel.grid.major = element_line(linetype = 'solid', colour = "#66A61E"), 
    panel.grid.minor = element_line(linetype = 'solid',colour = grid_colour)
  )


plot_C

filename <- paste0("./01_OutputFiles/", gsub("-","",Sys.Date()),"_",direction,"Regulated_DataSetInclusions.pdf") 

ggsave(plot = plot_C, filename = filename, useDingbats=FALSE, width = 10, height = 6.3)


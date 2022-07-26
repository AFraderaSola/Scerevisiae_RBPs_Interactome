###############################################################################
#
# This script generates the candidates plot for the project RBP interactome 
# project
#
# developer: Albert Fradera-Sola <A.FraderaSola@imb.de>
# date: 2022.07.26
#
#
###############################################################################

#########################################
############ Libraries ##################
#########################################

library(ggplot2)
library(ggforce)
library(readr)
library(readxl)
library(dplyr)
library(RColorBrewer)

#########################################
############## Script ###################
#########################################

## Load and tidy up the necessary input files

# RBP census

census <- read.csv("./01_InputFiles/00_SC_RBPs_Census.csv")

# TAP commercial library

TAP <- read_excel("./00_StocksandLibraries/00_Yeast_TAP_ComercialLibrary.xlsx", 
                  sheet = "Yeast TAP coordinates")
TAP <- TAP[TAP$ORF %in% census$ID,]

# TAP (IP screen) candidates

TAP_candidates <- read_excel("./01_InputFiles/00_IPScreeningMasterFile.xlsx")
TAP_candidates <- TAP_candidates[!TAP_candidates$Gene_ID == "WT",]
TAP_candidates <- TAP_candidates[!TAP_candidates$Gene_ID == "YNL317W",]

# Delition commercial library

KO <- read_excel("./00_StocksandLibraries/00_Yeast_KO_ComercialLibrary.xlsx", 
                 sheet = "mat_a_obs")
KO <- KO[KO$`ORF name` %in% census$ID,]

# Deletion (KO screen) candidates

KO_candidates <- read_excel("./01_InputFiles/00_KOScreeningMasterFile.xlsx")
KO_candidates <- KO_candidates[-1,]


# Data formating

p_census <- census
p_census$Set <- "RBP census"
p_census$Function <- "RBP census"
p_census <- p_census %>% group_by(Set, Function) %>% tally()

p_TAP <- TAP
p_TAP$Set <- "TAP library"
p_TAP$Function <- "TAP library"
p_TAP <- p_TAP %>% group_by(Set, Function) %>% tally()

p_KO <- KO
p_KO$Set <- "KO library"
p_KO$Function <- "KO library"
p_KO <- p_KO %>% group_by(Set, Function) %>% tally()

p_TAP_candidates <- TAP_candidates
p_TAP_candidates$Set <- "IP screen candidates"
p_TAP_candidates <- p_TAP_candidates %>% group_by(Set,Function) %>%
  tally()

p_KO_candidates <- KO_candidates
p_KO_candidates$Set <- "KO screen candidates"
p_KO_candidates <- p_KO_candidates %>% group_by(Set,Function) %>%
  tally()

# Plot data frame

plot <- rbind(p_census,p_TAP,p_KO,p_TAP_candidates,p_KO_candidates)


## Supplementary Figure

plot$Set <- factor(x = plot$Set, levels = c("RBP census","TAP library", 
                                            "KO library", "IP screen candidates", 
                                            "KO screen candidates"))

plot$Function <- factor(x = plot$Function, levels = c("RBP census", "TAP library",
                                                      "KO library", "Capping", 
                                                      "Cleavage","Export", "Degradation",
                                                      "Splicing", "Polyadenylation",
                                                      "Transport"))

colorA <- brewer.pal(n = 8, name = "Set2")
colorB <- brewer.pal(n = 9, name = "Purples")
final_color <- c(colorB[c(9,6,3)], colorA[3], colorA[8],"#666666",
                 colorA[5],colorA[4],colorA[1],colorA[6])

barplot <- ggplot(plot, aes(x = Set, y = n, fill = Function))+
                geom_bar(position="stack", stat="identity")+
                facet_zoom(ylim = c(0,50),show.area = F,)+
                scale_fill_manual(labels = c(c(expression(paste("Hentze, M. ", 
                                                                italic("et al."), " 2018")),
                                               expression(paste("Giaver, G. ", 
                                                                italic("et al."), " 2002")),
                                               expression(paste("Ghaemmaghami, S. ", 
                                                                italic("et al."), " 2003")),
                                               "Capping", "Cleavage","Export", 
                                               "Degradation","Splicing","Polyadenylation", 
                                               "Transport")),
                                             values = final_color)+
                theme_bw()+
                labs(fill = "")+
                theme(axis.text.x = element_text(angle = 45, hjust = 1))+
                theme(axis.text=element_text(size=14),
                      axis.title=element_text(size=16,face="bold"))+
                theme(legend.title = element_text(size = 16, face="bold"))+
                theme(legend.text = element_text(size = 14))+
                theme(legend.key.size = unit(1, 'cm'))+
                xlab("Set")+
                ylab("Proteins")

filename <- paste0("./01_OutputFiles/",gsub("-","",Sys.Date()),"_Candidates_Selection_A.pdf")

ggsave(plot = barplot, filename = filename,width = 15, height = 7.5)

## Main Figure

p_TAP_candidates <- TAP_candidates
p_TAP_candidates$Set <- "IP screen candidates"
p_TAP_candidates <- p_TAP_candidates %>% group_by(Set,Function,Pathway_DB) %>% tally()

plot <- p_TAP_candidates
plot$Function <- factor(x = plot$Function, levels = unique(plot$Function)[c(3,6,2,4,1,7,5)])
plot$Pathway_DB <- factor(x = plot$Pathway_DB, levels = rev(c("R-SCE-72086",
                                                              "KEGG-sce03015",
                                                              "KEGG-sce03018", 
                                                              "R-SCE-927802", 
                                                              "KEGG-sce03040", 
                                                              "R-HSA-72203",
                                                              "SGD")))


colorA <- brewer.pal(n = 8, name = "Set2")
color1 <- colorRampPalette(colors = c(colorA[5], "#FFFFFF"))(6)
color2 <- colorRampPalette(colors = c(colorA[4], "#FFFFFF"))(6)
final_color <- rev(c(colorA[3],color1[c(1,3,5)], color2[c(1,3)],colorA[8]))

barplot <- ggplot(plot, aes(x = Function, y = n, fill = Pathway_DB))+
  geom_bar(position="stack", stat="identity")+
  coord_flip()+
  scale_fill_manual(values = final_color)+
  theme_minimal()+
  labs(fill = "Database")+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16,face="bold"))+
  theme(legend.title = element_text(size = 16, face="bold"))+
  theme(legend.text = element_text(size = 14))+
  theme(legend.key.size = unit(1, 'cm'))+
  xlab("Function")+
  ylab("Proteins")
 

filename <- paste0("./01_OutputFiles/",gsub("-","",Sys.Date()),"_Candidates_Selection_B.pdf")

ggsave(plot = barplot, filename = filename,width = 15, height = 7.5)
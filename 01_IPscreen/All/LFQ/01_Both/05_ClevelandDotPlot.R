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

load(file = "05_InputFiles/PPI_Census.RData")
load(file = "05_InputFiles/PPI_BioGRID.RData")
load(file = "05_InputFiles/RDI_Census.RData")
load(file = "05_InputFiles/RDI_BioGRID.RData")

PPI_color <- "#66A61E"

RDI_color <- "#D95F02"

###################################
########### Script ################
###################################

# Census Plot

Census <- PPI_Census

colnames(Census)[2] <- "Percentage_PPI"

Census$Percentage_RDI <- RDI_Census$`Percentage of Protein IDs`

Census$Gene <- str_to_title(Census$Gene)

data <- Census %>% 
  rowwise() %>% 
  mutate( mymean = mean(c(Percentage_PPI,Percentage_RDI) )) %>% 
  # arrange(Percentage_RDI) %>% 
  arrange(desc(Percentage_RDI)) %>%
  mutate(Gene=factor(Gene, Gene))

census_plot <- ggplot(data) +
  geom_segment( aes(x=Gene, xend=Gene, y=Percentage_PPI, yend=Percentage_RDI), color="darkgrey") +
  geom_point( aes(x=Gene, y=Percentage_PPI), color=PPI_color, size=5 ) +
  geom_point( aes(x=Gene, y=Percentage_RDI), color=RDI_color, size=5 ) +
  geom_hline(yintercept=0.7, linetype="dashed", size = 1, color = "#666666")+
  scale_y_continuous(labels = scales::percent)+
  coord_flip()+
  theme_minimal() +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18,face="bold"))+
  theme(legend.title = element_text(size = 16, face="bold"))+
  theme(legend.text = element_text(size = 14))+
  theme(legend.key.size = unit(1, 'cm'))+
  theme(legend.position = "top")+
  xlab("Bait")+
  ylab("Inclusion at RBP census (%)")

filename <- paste0("./05_OutputFiles/",gsub("-","",Sys.Date()),"_CensusSuppl.pdf") 

ggsave(plot = census_plot, filename = filename, useDingbats=FALSE, width = 20, height = 9.5)

# BioGRID Plot

BioGRID <- PPI_BioGRID

BioGRID$Gene <- str_to_title(BioGRID$Gene)

colnames(BioGRID)[2] <- "Percentage_PPI"

BioGRID$Percentage_RDI <- RDI_BioGRID$`Percentage of Protein IDs`

data <- BioGRID %>% 
  rowwise() %>% 
  mutate( mymean = mean(c(Percentage_PPI,Percentage_RDI) )) %>% 
  # arrange(Percentage_PPI) %>% 
  arrange(desc(Percentage_PPI)) %>%
  mutate(Gene=factor(Gene, Gene))

biogrid_plot <- ggplot(data) +
  geom_segment( aes(x=Gene, xend=Gene, y=Percentage_PPI, yend=Percentage_RDI), color="darkgrey") +
  geom_point( aes(x=Gene, y=Percentage_PPI), color=PPI_color, size=5 ) +
  geom_point( aes(x=Gene, y=Percentage_RDI), color=RDI_color, size=5 ) +
  scale_y_continuous(labels = scales::percent)+
  coord_flip()+
  theme_minimal() +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18,face="bold"))+
  theme(legend.title = element_text(size = 16, face="bold"))+
  theme(legend.text = element_text(size = 14))+
  theme(legend.key.size = unit(1, 'cm'))+
  theme(legend.position = "top")+
  xlab("Bait")+
  ylab("Inclusion at BioGRID database (%)")

filename <- paste0("./05_OutputFiles/",gsub("-","",Sys.Date()),"_BioGRIDSuppl.pdf") 

ggsave(plot = biogrid_plot, filename = filename, useDingbats=FALSE, width = 20, height = 9.5)

# Standard Bar Plot

## PPI

data <- BioGRID %>% 
  rowwise() %>% 
  mutate( mymean = mean(c(Percentage_PPI,Percentage_RDI) )) %>% 
  # arrange(Percentage_PPI) %>% 
  arrange(desc(Percentage_PPI)) %>%
  mutate(Gene=factor(Gene, Gene))


ppi_plot <- ggplot(data, aes(y = Percentage_PPI , x= Gene)) +
  geom_bar(position="stack", stat="identity", fill = PPI_color)+
  geom_hline(yintercept=0.5, linetype="dashed", size = 1, color = "#666666" )+
  scale_y_continuous(labels = scales::percent,limits = c(0,0.85))+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18,face="bold"))+
  theme(legend.title = element_text(size = 16, face="bold"))+
  theme(legend.text = element_text(size = 14))+
  xlab("Bait")+
  ylab("PPI inclusion at BioGRID (%)")

filename <- paste0("./05_OutputFiles/",gsub("-","",Sys.Date()),"_Barplot_PPI_BioGRIDSuppl.pdf") 

ggsave(plot = ppi_plot, filename = filename, useDingbats=FALSE, width = 11.62, height = 16.66)

data <- BioGRID %>% 
  rowwise() %>% 
  mutate( mymean = mean(c(Percentage_PPI,Percentage_RDI) )) %>% 
  # arrange(Percentage_PPI) %>% 
  arrange(desc(Percentage_PPI)) %>%
  mutate(Gene=factor(Gene, Gene))

data <- BioGRID %>% 
  rowwise() %>% 
  mutate( mymean = mean(c(Percentage_PPI,Percentage_RDI) )) %>% 
  # arrange(Percentage_PPI) %>% 
  arrange(desc(Percentage_RDI)) %>%
  mutate(Gene=factor(Gene, Gene))

rdi_plot <- ggplot(data, aes(y = Percentage_RDI , x= Gene)) +
  geom_bar(position="stack", stat="identity", fill = RDI_color)+
  geom_hline(yintercept=0.5, linetype="dashed", size = 1, color = "#666666" )+
  scale_y_continuous(labels = scales::percent,limits = c(0,0.85))+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18,face="bold"))+
  theme(legend.title = element_text(size = 16, face="bold"))+
  theme(legend.text = element_text(size = 14))+
  xlab("Bait")+
  ylab("RDI inclusion at BioGRID (%)")

filename <- paste0("./05_OutputFiles/",gsub("-","",Sys.Date()),"_Barplot_RDI_BioGRIDSuppl.pdf") 

ggsave(plot = rdi_plot, filename = filename, useDingbats=FALSE, width = 11.62, height = 16.66)


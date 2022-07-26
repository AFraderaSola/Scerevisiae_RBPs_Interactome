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

library(readr)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(ggpubr)

#########################################
############## Script ###################
#########################################

## Load and tidy up the necessary input files

files <- list.files(path = "./00_CandidateInteractions/BioGRID/",
                    pattern = "BiogridInteractions.txt")

counts <- c()

for (i in 1:length(files)) {
  
  interaction <- read_delim(file = paste0("00_CandidateInteractions/BioGRID/",
                                          files[i]),
                            "\t", 
                            escape_double = FALSE,
                            trim_ws = TRUE)
  
  interaction <- interaction[interaction$`Experimental System Type` == "physical",]
  
  counts_loop <- interaction %>% count(`Experimental System`)%>% 
    mutate(percent = n/sum(n))
  
  counts_loop$Gene <- gsub(pattern = "Y.*_",replacement = "",
                      x = gsub(pattern = "_B.*",replacement = "",x = files[i]))
  
  counts <- rbind(counts, counts_loop)
  
}

## Plot

order <- counts[counts$`Experimental System` == "Affinity Capture-MS",]
order <- order[order(-order$percent),]
order <- as.character(order$Gene)
order <- c(order,"HBS1")

counts$Gene <- factor(x = counts$Gene,levels = order)

counts <- counts[counts$`Experimental System` == "Affinity Capture-MS",]

hbs1 <- counts[1,]
counts <- rbind(counts,hbs1)
counts$n[40] <- 0
counts$percent[40] <- 0
counts$Gene[40] <- "HBS1"

percent <- ggplot(counts, aes(y=percent, x=Gene)) +
  geom_bar(position="stack", stat="identity", fill = "#ab5757")+
  geom_hline(yintercept=0.5, linetype="dashed", size = 1, color = "#666666")+
  scale_y_continuous(labels = scales::percent,limits = c(0,0.85))+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18,face="bold"))+
  theme(legend.key.size = unit(0.5, 'cm'))+
  theme(legend.position = "top")+
  xlab("Bait")+
  ylab("BioGRID Affinity Capture-MS evidence (%)")

ggsave(percent, filename = "./02_OutputFiles/Percent.pdf",width = 11.62,height = 16.66)
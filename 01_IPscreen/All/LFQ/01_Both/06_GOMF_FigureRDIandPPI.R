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

load(file = "06_InputFiles/20221011_RDI_GO_MF_df.RData")
load(file = "06_InputFiles/20221011_PPI_GO_MF_df.RData")

PPI_color <- "#66A61E"

RDI_color <- "#D95F02"

colour_grid <- c(RDI_color,PPI_color)

colors <- c(brewer.pal(8, "Paired")[1], "#b4b4b4")

###################################
########## Script #################
###################################

# Supp Figure

df_GO_PPI$Interactor <- "PPI"

df_GO_RDI$Interactor <- "RDI"

plot_df <- rbind(df_GO_PPI, df_GO_RDI)

plot_df <- plot_df %>%
  group_by(ID, Interactor) %>%
  summarize(Score = n())

RDIdf <- plot_df[plot_df$Interactor == "RDI",]

PPIdf <- plot_df[plot_df$Interactor == "PPI",]

RDIdf_noPPI <- RDIdf[!RDIdf$ID %in% PPIdf$ID,]

PPIdf_n0 <- RDIdf_noPPI

PPIdf_n0$Interactor <- "PPI"

PPIdf_n0$Score <- 0

PPIdf_noRDI <- PPIdf[!PPIdf$ID %in% RDIdf$ID,]

RDIdf_n0 <- PPIdf_noRDI

RDIdf_n0$Interactor <- "RDI"

RDIdf_n0$Score <- 0

plot_df <- rbind(plot_df, PPIdf_n0, RDIdf_n0)

plot_df <- plot_df[order(plot_df$Score,decreasing = T),]

plot_df$ID <- factor(plot_df$ID, levels = unique(plot_df$ID))

plot_df$Interactor <- factor(plot_df$Interactor, levels = c("RDI", "PPI"))

# plot_df$RNA <- FALSE
# 
# plot_df[grep(pattern = "RNA|nucleotide|nucleic|ribonucleoside|translation|ribosome",x = plot_df$Description),]$RNA <- TRUE
# 
# plot_df$RNA <- factor(x = plot_df$RNA, levels = c("TRUE","FALSE"))

# selected_analysis <- c("Pfam", "SUPERFAMILY")
# 
# plot_df <- plot_df[plot_df$Analysis %in% selected_analysis,]

plot <- ggplot(plot_df, aes(fill = Interactor,y = Score , x= ID, colour = Interactor)) +
  geom_bar(width = 0.6, position=position_dodge(width = 0.8), stat="identity", size = 0.8)+
  labs(fill = "Signature", colour = "Group")+
  scale_fill_manual(values = colour_grid)+
  scale_colour_manual(values = colour_grid)+
  coord_flip()+
  theme_minimal()+
  theme(axis.text=element_text(size=11),
        axis.title=element_text(size=28,face="bold"))+
  theme(legend.text = element_text(size = 24),
        legend.title=element_text(size=26,face="bold"))+
  theme(legend.key.size = unit(1, 'cm'))+
  theme(legend.position = "top")+
  theme(strip.text = element_text(size=15,face="bold"))+
  xlab("Molecular Function GO Term")+
  ylab("Baits with overrepresented term")

filename <- paste0("./06_OutputFiles/",gsub(pattern = "-",x = Sys.Date(),replacement = ""),"_SuppGOMF_Both.pdf")

ggsave(plot = plot, filename = filename, useDingbats=FALSE, width = 28, height = 13.2)


# Figure 3

plot_df_PPI <- plot_df[plot_df$Interactor == "PPI",]

plot_df_PPI <- plot_df_PPI[c(1:5),]

plot_df_RDI <- plot_df[plot_df$Interactor == "RDI",]

plot_df_RDI <- plot_df_RDI[c(1:5),]

plot_df <- rbind(plot_df_PPI,plot_df_RDI)

# plot <- ggplot(plot_df, aes(fill = RNA,y = Score , x= Description, colour = Interactor)) +
#   geom_bar(width = 0.6, position=position_dodge(width = 0.8), stat="identity", size = 1.2)+
#   labs(fill = "Term", colour = "Group")+
#   scale_fill_manual(labels = c("RNA related",
#                                "Not RNA related"
#   ),
#   values = colors)+
#   scale_color_manual(values = colour_grid)+
#   coord_flip()+
#   scale_y_continuous(expand=c(0,0.5))+
#   theme_minimal()+
#   theme(axis.text=element_text(size=13),
#         axis.title=element_text(size=15,face="bold"))+
#   theme(legend.text = element_text(size = 14),
#         legend.title=element_text(size=13,face="bold"))+
#   theme(legend.key.size = unit(1, 'cm'))+
#   theme(legend.position = "top")+
#   theme(strip.text = element_text(size=15,face="bold"))+
#   guides(fill = guide_legend(order = 2),colour = guide_legend(order = 1))+
#   xlab("Molecular Function GO Term")+
#   ylab("Baits with overrepresented term")

addition <- tibble(Description = c("ion binding",
                                       "catalytic activity",
                                       "organic cyclic compound binding",
                                       "mRNA binding"),
                       Interactor = c("RDI",
                                      "RDI",
                                      "PPI",
                                      "PPI"),
                       Score = c(0,0,0,0),
                       RNA = c("FALSE", "FALSE", "FALSE", "TRUE"))

plot_df <- rbind(plot_df, addition)

plot_df <- plot_df[order(plot_df$Score,decreasing = T),]

plot_df$Description <- factor(plot_df$Description, levels = unique(plot_df$Description))

plot_df$Interactor <- factor(plot_df$Interactor, levels = c("RDI", "PPI"))

colorado <- function(src, boulder) {
  if (!is.factor(src)) src <- factor(src)                   # make sure it's a factor
  src_levels <- levels(src)                                 # retrieve the levels in their order
  brave <- boulder %in% src_levels                          # make sure everything we want to make bold is actually in the factor levels
  if (all(brave)) {                                         # if so
    b_pos <- purrr::map_int(boulder, ~which(.==src_levels)) # then find out where they are
    b_vec <- rep("plain", length(src_levels))               # make'm all plain first
    b_vec[b_pos] <- "bold"                                  # make our targets bold
    b_vec                                                   # return the new vector
  } else {
    stop("All elements of 'boulder' must be in src")
  }
}

plot <- ggplot(plot_df, aes(fill = Interactor,y = Score , x= Description, colour = Interactor)) +
  geom_bar(width = 0.6, position=position_dodge(width = 0.8), stat="identity", size = 1.2)+
  scale_fill_manual(values = colour_grid)+
  scale_colour_manual(values = colour_grid)+
  coord_flip()+
  scale_y_continuous(expand=c(0,0.5))+
  theme_minimal()+
  labs(fill = "Group", colour = "Group")+
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18,face="bold"))+
  theme(legend.text = element_text(size = 14),
        legend.title=element_text(size=16,face="bold"))+
  theme(legend.key.size = unit(1, 'cm'))+
  theme(legend.position = "top")+
  theme(strip.text = element_text(size=18,face="bold"))+
  guides(fill = guide_legend(order = 2),colour = guide_legend(order = 1))+
  xlab("Molecular Function GO Term")+
  ylab("Baits with overrepresented term")+
  theme(axis.text.y=element_text(face=colorado(plot_df$Description, c("RNA binding", "nucleic acid binding", "mRNA binding"))))

filename <- paste0("./06_OutputFiles/",gsub(pattern = "-",x = Sys.Date(),replacement = ""),"_GOMF_Both.pdf")

ggsave(plot = plot, filename = filename, useDingbats=FALSE, width = 20, height = 4.5)


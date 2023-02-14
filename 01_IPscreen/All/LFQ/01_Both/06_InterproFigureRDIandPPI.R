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

load(file = "06_InputFiles/20221011_PPI_interpro_df.RData")
load(file = "06_InputFiles/20221011_RDI_interpro_df.RData")
load(file = "06_InputFiles/20221011_PPI_Suppinterpro_df.RData")
load(file = "06_InputFiles/20221011_RDI_Suppinterpro_df.RData")

PPI_color <- "#66A61E"

RDI_color <- "#D95F02"

colour_grid <- c(RDI_color,PPI_color)

colors <- c(brewer.pal(8, "Paired")[1], "#b4b4b4")

RDI_interpro_df$n[11] <- 17

###################################
########## Script #################
###################################

# Figure 3

PPI_interpro_df$Interactor <- "PPI"

RDI_interpro_df$Interactor <- "RDI"

plot_df <- rbind(PPI_interpro_df, RDI_interpro_df)

plot_df <- plot_df[order(plot_df$n,decreasing = T),]

plot_df$Signature <- factor(plot_df$Signature, levels = unique(plot_df$Signature))

plot_df$Interactor <- factor(plot_df$Interactor, levels = c("RDI", "PPI"))

selected_analysis <- c("Pfam", "SUPERFAMILY")

plot_df <- plot_df[plot_df$Analysis %in% selected_analysis,]

# plot <- ggplot(plot_df, aes(fill = RNA,y = n , x= Signature, colour = Interactor)) +
#   facet_wrap(~Analysis, scales = "free_y",ncol = 2)+
#   geom_bar(width = 0.6, position=position_dodge(width = 0.8), stat="identity", size = 1.2)+
#   labs(fill = "Signature", colour = "Group")+
#   scale_fill_manual(labels = c("RNA related",
#                                "Not RNA related"
#   ),
#   values = colors)+
#   scale_color_manual(values = colour_grid)+
#   coord_flip()+
#   theme_minimal()+
#   theme(axis.text=element_text(size=13),
#         axis.title=element_text(size=15,face="bold"))+
#   theme(legend.text = element_text(size = 14),
#         legend.title=element_text(size=13,face="bold"))+
#   theme(legend.key.size = unit(1, 'cm'))+
#   theme(legend.position = "top")+
#   theme(strip.text = element_text(size=15,face="bold"))+
#   xlab("Signature")+
#   ylab("Baits with overrepresented signature")

RDIdf <- plot_df[plot_df$Interactor == "RDI",]

PPIdf <- plot_df[plot_df$Interactor == "PPI",]

RDIdf_noPPI <- RDIdf[!RDIdf$`Signature accession` %in% PPIdf$`Signature accession`,]

PPIdf_n0 <- RDIdf_noPPI

PPIdf_n0$Interactor <- "PPI"

PPIdf_n0$n <- 0

PPIdf_noRDI <- PPIdf[!PPIdf$`Signature accession` %in% RDIdf$`Signature accession`,]

RDIdf_n0 <- PPIdf_noRDI

RDIdf_n0$Interactor <- "RDI"

RDIdf_n0$n <- 0

plot_df <- rbind(plot_df, PPIdf_n0, RDIdf_n0)

plot_df <- plot_df[order(plot_df$n,decreasing = T),]

plot_df[plot_df$Signature == "RNA recognition motif. (a.k.a. RRM, RBD, or RNP domain)",]$Signature <- plot_df[plot_df$Signature == "RNA recognition motif. (a.k.a. RRM, RBD, or RNP domain)",]$`Interpro description`

plot_df$Signature <- factor(plot_df$Signature, levels = unique(plot_df$Signature))

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

plot_df[plot_df$Signature == "RNA recognition motif. (a.k.a. RRM, RBD, or RNP domain)",]$Signature <- plot_df[plot_df$Signature == "RNA recognition motif. (a.k.a. RRM, RBD, or RNP domain)",]$`Interpro description`
  
plot <- ggplot(plot_df, aes(fill = Interactor,y = n , x= Signature, colour = Interactor)) +
  facet_wrap(~Analysis, scales = "free_y",ncol = 2)+
  geom_bar(width = 0.6, position=position_dodge(width = 0.8), stat="identity", size = 1.2)+
  scale_fill_manual(values = colour_grid)+
  scale_color_manual(values = colour_grid)+
  coord_flip()+
  theme_minimal()+
  labs(fill = "Group", colour = "Group")+
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18,face="bold"))+
  theme(legend.text = element_text(size = 14),
        legend.title=element_text(size=16,face="bold"))+
  theme(legend.key.size = unit(1, 'cm'))+
  theme(legend.position = "top")+
  theme(strip.text = element_text(size=18,face="bold"))+
  xlab("Signature")+
  ylab("Baits with overrepresented signature")+
  theme(axis.text.y=element_text(face=colorado(plot_df$Signature, c("RNA-binding domain, RBD"))))

filename <- paste0("./06_OutputFiles/",gsub(pattern = "-",x = Sys.Date(),replacement = ""),"_Interpro_Both.pdf")

ggsave(plot = plot, filename = filename, useDingbats=FALSE, width = 20, height = 4.5)

# Supp figure interpro

PPI_Suppinterpro_df$Interactor <- "PPI"

RDI_Suppinterpro_df$Interactor <- "RDI"

plot_df <- rbind(PPI_Suppinterpro_df, RDI_Suppinterpro_df)

RDIdf <- plot_df[plot_df$Interactor == "RDI",]

PPIdf <- plot_df[plot_df$Interactor == "PPI",]

RDIdf_noPPI <- RDIdf[!RDIdf$`Signature accession` %in% PPIdf$`Signature accession`,]

PPIdf_n0 <- RDIdf_noPPI

PPIdf_n0$Interactor <- "PPI"

PPIdf_n0$n <- 0

PPIdf_noRDI <- PPIdf[!PPIdf$`Signature accession` %in% RDIdf$`Signature accession`,]

RDIdf_n0 <- PPIdf_noRDI

RDIdf_n0$Interactor <- "RDI"

RDIdf_n0$n <- 0

plot_df <- rbind(plot_df, PPIdf_n0, RDIdf_n0)

plot_df <- plot_df[order(plot_df$n,decreasing = T),]

plot_df$Signature <- factor(plot_df$Signature, levels = unique(plot_df$Signature))

plot_df$Interactor <- factor(plot_df$Interactor, levels = c("RDI", "PPI"))

plot_df$RNA <- FALSE

RNA_filter <- c("SSF54928", "PF00076", "PS50102", "SSF50182","SM00651","PF01423","SSF54928","PF00076","SSF54211",
                "PF01138", "PS50102", "PS00178", "TIGR00308", "SSF50249", "SSF55666", "PS00039", "PS00368", 
                "G3DSA:1.10.620.20", "PF00270", "PF01248", "PF03725", "PTHR11097", "PTHR11246", "PTHR17204",
                "PTHR23338", "PTHR23409", "PTHR23409:SF20", "cd01049", "SM00322", "cd00105", "SSF50182", "SSF54791",
                "PF00013", "PF00133", "PF01423", "PF08264")

plot_df[plot_df$`Signature accession` %in% RNA_filter,]$RNA <- TRUE 

plot_df$RNA <- factor(x = plot_df$RNA, levels = c("TRUE","FALSE"))

colors <- c(brewer.pal(8, "Paired")[1], "#b4b4b4")

plot <- ggplot(plot_df, aes(fill = Interactor,y = n , x= `Signature accession`, colour = Interactor)) +
  facet_wrap(~Analysis, scales = "free_y",ncol = 3)+
  geom_bar(position=position_dodge(width = 1), stat="identity", size = 1.5)+
  labs(fill = "Signature", colour = "Group")+
  scale_fill_manual(values = colour_grid)+
  scale_color_manual(values = colour_grid)+
  coord_flip()+
  theme_minimal()+
  theme(axis.text=element_text(size=11),
        axis.title=element_text(size=28,face="bold"))+
  theme(legend.text = element_text(size = 24),
        legend.title=element_text(size=26,face="bold"))+
  theme(legend.key.size = unit(1, 'cm'))+
  theme(legend.position = "top")+
  theme(strip.text = element_text(size=22,face="bold"))+
  xlab("Signature accession number")+
  ylab("Baits with overrepresented signature")

filename <- paste0("./06_OutputFiles/",gsub(pattern = "-",x = Sys.Date(),replacement = ""),"_SuppInterpro_Both.pdf")

ggsave(plot = plot, filename = filename, useDingbats=FALSE, width = 28, height = 13.2)


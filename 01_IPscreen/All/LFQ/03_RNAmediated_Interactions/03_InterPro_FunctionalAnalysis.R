##################################
######### Libraries ##############
##################################

library(readr)
library(dplyr)
library(ggplot2)
library(viridis)
library(RColorBrewer)

##################################
######## Variables ###############
##################################

interaction <- "RNAm"

##################################
########## Script ################
##################################

if (interaction == "PPI") {
  
  pattern <- "RNase_vs_WT"
  
}else{
  
  pattern <- "IP_vs_RNase"
  
}

files_enriched <- list.files(path = "./00_IPsResultsFiles/",pattern = pattern)

pattern <- "Proteins.csv"

files_quantified <- list.files(path = "./00_IPsResultsFiles/", pattern = pattern)


if (interaction == "PPI") {
  
  interpro_enriched <- read.delim(file = "./03_InputFiles/20220429_PPIInteractors_Cleaned.fasta.tsv",
                                    header = F,
                                    sep = "\t")
  
}else{
  
  interpro_enriched <- read.delim(file = "./03_InputFiles/20220429_RNAmediatedInteractors_Cleaned.fasta.tsv",
                                    header = F,
                                    sep = "\t")
  
}


colnames(interpro_enriched) <- c("Protein Accession",
                                "Sequence MD5 digest",
                                "Sequence Length",
                                "Analysis", 
                                "Signature Accession",
                                "Signature Description",
                                "Start location",
                                "Stop location",
                                "Score",
                                "Status",
                                "Date",
                                "InterPro annotations - accession",
                                "InterPro annotations - description",
                                "GO annotations",
                                "Pathway annotations")

interpro_quantified <- read.delim(file = paste0("./03_InputFiles/20220429_QuantifiedProteins_Cleaned.fasta.tsv"),
                                header = F,
                                sep = "\t")

colnames(interpro_quantified) <- c("Protein Accession",
                                   "Sequence MD5 digest",
                                   "Sequence Length",
                                   "Analysis", 
                                   "Signature Accession",
                                   "Signature Description",
                                   "Start location",
                                   "Stop location",
                                   "Score",
                                   "Status",
                                   "Date",
                                   "InterPro annotations - accession",
                                   "InterPro annotations - description",
                                   "GO annotations",
                                   "Pathway annotations")

df_significant <- c()

for (k in 1:length(files_enriched)) {
  
  stock <- gsub(pattern = "_Q.*",replacement = "", x = files_quantified[k])
  
  enriched_ids_filter <- read.csv(file = paste0("./00_IPsResultsFiles/",files_enriched[k]),header = T)
  
  enriched_ids_filter <- enriched_ids_filter$Protein.IDs
  
  enriched_ids_filter <- gsub(pattern = "-TAP",replacement = "",x = gsub(pattern = ";.*",replacement = "",x = enriched_ids_filter))
  
  interpro_enriched_filtered <- interpro_enriched[interpro_enriched$`Protein Accession` %in% enriched_ids_filter,]
  
  quantified_ids_filter <- read.csv(file = paste0("./00_IPsResultsFiles/",files_quantified[k]),header = T)
  
  quantified_ids_filter <- quantified_ids_filter$Protein.IDs
  
  quantified_ids_filter <- gsub(pattern = "-TAP",replacement = "",x = gsub(pattern = ";.*",replacement = "",x = quantified_ids_filter))
  
  interpro_quantified_filtered <- interpro_quantified[interpro_quantified$`Protein Accession` %in% quantified_ids_filter,]
  
  analysis <- unique(interpro_enriched_filtered$Analysis)
  
  df <- c()
  
  for (i in 1:length(analysis)) {
    
    interpro_enriched_loop <- interpro_enriched_filtered[interpro_enriched_filtered$Analysis == analysis[i],]
    
    interpro_enriched_loop <- interpro_enriched_loop %>% group_by(`Protein Accession`, `Signature Accession`) %>% tally()
    
    IDsEnriched <- length(unique(interpro_enriched_loop$`Protein Accession`))
    
    interpro_quantified_loop <- interpro_quantified_filtered[interpro_quantified_filtered$Analysis == analysis[i],]
    
    interpro_quantified_loop <- interpro_quantified_loop %>% group_by(`Protein Accession`, `Signature Accession`) %>% tally()
    
    IDsQuantified <- length(unique(interpro_quantified_loop$`Protein Accession`))
    
    accession <- unique(interpro_enriched_loop$`Signature Accession`)
    
    for (j in 1:length(accession)) {
      
      a_df <- interpro_enriched_loop[interpro_enriched_loop$`Signature Accession` == accession[j],]
      
      a <- length(unique(a_df$`Protein Accession`))
      
      c_df <- interpro_enriched_loop[!interpro_enriched_loop$`Signature Accession` == accession[j],]
      
      c_df <- c_df[!c_df$`Protein Accession` %in% a_df$`Protein Accession`,]
      
      c <- length(unique(c_df$`Protein Accession`))
      
      print(a+c == IDsEnriched)
      
      b_df <- interpro_quantified_loop[interpro_quantified_loop$`Signature Accession` == accession[j],]
      
      b <- length(unique(b_df$`Protein Accession`))
      
      d_df <- interpro_quantified_loop[!interpro_quantified_loop$`Signature Accession` == accession[j],]
      
      d_df <- d_df[!d_df$`Protein Accession` %in% b_df$`Protein Accession`,]
      
      d <- length(unique(d_df$`Protein Accession`))
      
      print(b+d == IDsQuantified)
      
      pvalue <- fisher.test(matrix(c(a,b,c,d),byrow=TRUE,ncol=2),alternative="greater")$p.val
      
      estimate <- unname(fisher.test(matrix(c(a,b,c,d),byrow=TRUE,ncol=2),alternative="greater")$estimate)
      
      loop_analysis <- analysis[i]
      
      loop_accession <- accession[j]
      
      loop_description_signature <- unique(interpro_enriched_filtered[interpro_enriched_filtered$`Signature Accession` == loop_accession,]$`Signature Description`)
      
      loop_description_interpro <- unique(interpro_enriched_filtered[interpro_enriched_filtered$`Signature Accession` == loop_accession,]$`InterPro annotations - description`)
      
      interactors_in_signature <- paste0(unique(a_df$`Protein Accession`), collapse = ";")
      
      loop_df <- data.frame(stock,interactors_in_signature, loop_analysis, loop_accession, loop_description_signature, loop_description_interpro,
                            a,c,b,d,pvalue, estimate)
      
      colnames(loop_df) <- c("Pulldown stock ID","Enriched Interacors ID","Analysis", "Signature accession", "Signature description", "Interpro description",
                             "Enriched IN signature","Enriched NOT in signature", "Quantified IN signature", "Quantified NOT in signature",
                             "P-value", "Estimate")
      
      df <- rbind(df, loop_df)
      
    }
    
  }
  
  adjpvalue <- p.adjust(p = df$`P-value`, method = "BH")
  
  df$`Adjusted P-value` <- adjpvalue
 
  df_significant_loop <- df[df$`P-value` < 0.01,]
  
  write.csv(x = df_significant_loop,
            file = paste0("./03_OutputFiles/",gsub(pattern = "-",replacement = "",x = Sys.Date()),"_",
                          stock,"_",interaction, "_p01_SignificantInterPro.csv"),
            row.names = F)
  
  df_analysis_loop <- df
  
  write.csv(x = df_analysis_loop,
            file = paste0("./03_OutputFiles/",gsub(pattern = "-",replacement = "",x = Sys.Date()),"_",
                          stock,"_",interaction, "_ResultsInterPro.csv"),
            row.names = F)
  
  df_significant <- rbind(df_significant, df_significant_loop)

}


analysis_df <- df_significant

analysis <- unique(analysis_df$Analysis)

plot_df <- c()

for (l in 1:length(analysis)) {
  
  loop_analysis_df <- analysis_df[analysis_df$Analysis == analysis[l],]
  
  loop_analysis_df <- loop_analysis_df %>% group_by(`Signature accession`, `Signature description`, `Interpro description`) %>% tally()
  
  loop_analysis_df$Analysis <- analysis[l]
  
  loop_analysis_df  <- loop_analysis_df[,c(5,1:4)]
  
  plot_df <- rbind(loop_analysis_df, plot_df)
  
}



if (interaction == "PPI") {

  plot_df <- plot_df[plot_df$n > 1,]

}else{

  plot_df <- plot_df[plot_df$n > 1,]

}


plot_df <- plot_df[!plot_df$`Signature description` == "-",]

plot_df$Signature <- plot_df$`Signature description`

for (i in 1:nrow(plot_df)) {
  
  if (plot_df$`Interpro description`[i] == "-") {
    
    plot_df$`Interpro description`[i] <- plot_df$`Signature description`[i]
    
  }
  
}

plot_df[plot_df$Analysis == "SMART",]$Signature <- plot_df[plot_df$Analysis == "SMART",]$`Interpro description`

## Save plot_df for PPI + RDI together plot

RDI_Suppinterpro_df <- plot_df

save(RDI_Suppinterpro_df, file = paste0("./03_OutputFiles/",gsub("-","",Sys.Date()),"_RDI_Suppinterpro_df.RData"))

##

selected_analysis <- c("Pfam", "ProSiteProfiles", "SUPERFAMILY", "SMART")

plot_df <- plot_df[plot_df$Analysis %in% selected_analysis,]

plot_df <- plot_df[order(plot_df[,5], plot_df[,4]),]

levels <- unique(plot_df$Signature)

plot_df$Signature <- factor(x = plot_df$Signature, levels = levels)

plot_df$RNA <- TRUE

if (interaction == "PPI") {
  
  RNA_filter <- c("PF02347","PF02127","SSF101821","SSF48371","SSF47113")
  
}else{
  
  RNA_filter <- c("PF00571","PF00478","PS51371","PF00025","SM00116","SM01240","SSF54631","SSF51412","SSF52374")
  
}

plot_df[plot_df$`Signature accession` %in% RNA_filter,]$RNA <- FALSE 

plot_df$RNA <- factor(x = plot_df$RNA, levels = c("TRUE","FALSE"))

colors <- c(brewer.pal(8, "Paired")[1], "#b4b4b4")

if (interaction == "PPI") {
  
  colour_grid <- "#66A61E"
  
}else{
  
  colour_grid <- "#D95F02"
  
}


plot <- ggplot(plot_df, aes(fill = RNA,y = n , x= Signature, colour = colour_grid)) +
  facet_wrap(~Analysis, scales = "free_y",ncol = 2)+
  geom_bar(position="dodge", stat="identity", size = 1.5)+
  labs(fill = "")+
  scale_fill_manual(labels = c("RNA related",
                               "Not RNA related"
  ),
  values = colors)+
  scale_color_manual(values = colour_grid)+
  # scale_y_continuous(breaks = seq(0, 12, by = 2))+
  coord_flip()+
  theme_minimal()+
  guides(color = FALSE)+
  theme(axis.text=element_text(size=17),
        axis.title=element_text(size=19,face="bold"))+
  theme(legend.text = element_text(size = 15))+
  theme(legend.key.size = unit(1, 'cm'))+
  theme(legend.position = "top")+
  theme(strip.text = element_text(size=19,face="bold"))+
  xlab("Signature")+
  ylab("Pulldowns with enriched signature")
# theme(
#       # panel.grid.major = element_line(linetype = 'solid', colour = colour_grid), 
#       panel.grid.minor = element_line(linetype = 'solid',colour = colour_grid))

filename <- paste0("./03_OutputFiles/", gsub("-","",Sys.Date()),"_",interaction,"_p01_InterProPlot.pdf") 

ggsave(plot = plot, filename = filename, useDingbats=FALSE, width = 20, height = 9)

## Save plot_df for PPI + RDI together plot

RDI_interpro_df <- plot_df

save(RDI_interpro_df, file = paste0("./03_OutputFiles/",gsub("-","",Sys.Date()),"_RDI_interpro_df.RData"))

##

filter <- plot_df$`Signature accession`

unique_pulldowns <- length(unique(analysis_df[analysis_df$`Signature accession` %in% filter,]$`Pulldown stock ID`))

pulldowns <- analysis_df[analysis_df$`Signature accession` %in% filter,] %>% group_by(`Signature accession`) %>% tally()

pulldowns <- sum(pulldowns$n)

set <- c("Total", "Unique")

n <- c(pulldowns, unique_pulldowns)

supp_plot_df <- data.frame(set,n)

supp_plot <- ggplot(supp_plot_df, aes(y = n , x= set)) +
  geom_bar(position="dodge", stat="identity")+
  labs(fill = "")+
  theme_minimal()+
  scale_fill_viridis(discrete = T,option = "F")+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16,face="bold"))+
  theme(legend.text = element_text(size = 12))+
  theme(legend.key.size = unit(1, 'cm'))+
  theme(legend.position = "top")+
  theme(strip.text = element_text(size=16,face="bold"))+
  xlab("Pulldowns with signature")+
  ylab("Count")+
  theme(
    # panel.grid.major = element_line(linetype = 'solid', colour = "#9EB57F"), 
    panel.grid.minor = element_line(linetype = 'solid',colour = "#9EB57F"))

filename <- paste0("./03_OutputFiles/", gsub("-","",Sys.Date()),"_",interaction,"_p01_InterProSuppPlot.pdf") 

ggsave(plot = supp_plot, filename = filename, useDingbats=FALSE, width = 10, height = 10)


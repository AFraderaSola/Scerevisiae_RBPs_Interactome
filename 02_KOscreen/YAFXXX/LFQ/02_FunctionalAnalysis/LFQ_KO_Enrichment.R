#############################################
############# Libraries #####################
#############################################

library(clusterProfiler)
library(DOSE)
library(ReactomePA)
library(Cairo)
library(ggplot2)
library(ggforce)
library(DESeq2)
library(readr)
library(dplyr)
library(tidyr)
library(plyr)
library(cfpscripts)
library(topGO)
library(UpSetR)
library(GOSemSim)
library(ggupset)
library(viridis)
library(xlsx)

#############################################
############# Functions #####################
#############################################

splitColumnBySep <- function(df, column, sep=";"){
  nameCol <- column
  names(df)[grep(column, names(df))] <- "column"
  df2 <- df %>%
    mutate(column = strsplit(column, sep)) %>%
    unnest(column)
  names(df2)[grep("column", names(df2))] <- nameCol
  return(df2)
}


#############################################
############ Variables ######################
#############################################

# this.dir <- dirname(parent.frame(2)$ofile)
# setwd(this.dir)

print(getwd())

# Define your quantified protein groups path and name:

pg_quant_file <- list.files(path = "../01_CoreAnalysis/20220311_results_LFQ/",pattern = ".*QuantifiedProteins.csv")

# Define your enriched protein groups path and name:

pg_enriched_file <- list.files(path = "../01_CoreAnalysis/20220311_results_LFQ/",pattern = ".*Enriched.*.csv")

# Define your experiment Name:

experiment <- "KO"

# Define the regulation direction (any combination of "Up", "Down"):

direction_user <- c("Up","Down")

# Define the ID column name:

column_name <- "Protein.IDs"

# Define if filtering

filter <- F

# Define filter file name

filter_IDs <- "filter.txt"

# Define which ontology to use for GO terms (any combination of "BP", "CC", "MF"):

ont_user <- c("BP","CC","MF")

# Define wether you want simplified resutls or not:
## !! ATENTION !! need to take out all ggsave from plot to work. Currently (04/03/2022) it does not work when True because of that!! ##

simplified <- F

# Define the organism database annotation package for GO:

org <- "org.Sc.sgd.db"

# Define the organism name for reactome (one of "human", "rat", "mouse", "celegans", "yeast", "zebrafish", "fly"):

org_reactome <- "yeast"

# Define the organism name for KEGG (one of http://www.genome.jp/kegg/catalog/org_list.html):

org_KEGG <- "sce"

#############################################
############### Script ######################
#############################################

# IP stock name

stock <- gsub(pattern = ".*een/|/02_.*",replacement = "",x = getwd())

masterfile <- read.xlsx(file = "../../00_KOScreeningMasterFile.xlsx",header = T,sheetIndex = 1)

genename <- paste0(tolower(masterfile[masterfile$Stock_Name == gsub(pattern = ".*een/|/02_.*",replacement = "",x = getwd()),]$Gene_Name), "Δ")

for (j in 1:length(direction_user)) {
  
  direction <- direction_user[j]
  
  for (i in 1:length(ont_user)) {
    
    direction <- direction_user[j]
    
    ont <- ont_user[i]
    
    # Retrieve the IDs
    
    ## Keep the protein IDs for the universe:
    
    pg_quant <- read_csv(paste0("../01_CoreAnalysis/20220311_results_LFQ/",pg_quant_file))
    
    pg_quant <- splitColumnBySep(df = pg_quant, "Protein.IDs", sep = ";")
    
    univers_ID <- unlist(list(unique(pg_quant[,grep(pattern = column_name,x = colnames(pg_quant))])),use.names = F)
    
    ## Keep the protein IDs for the enriched:
    
    pg_enriched <- read_csv(paste0("../01_CoreAnalysis/20220311_results_LFQ/",pg_enriched_file))
    
    pg_enriched <- splitColumnBySep(df = pg_enriched, "Protein.IDs", sep = ";")
    
    if (direction == "Up") {
      
      pg_enriched <- pg_enriched[,grep(pattern = paste0("difference|", column_name),x = colnames(pg_enriched))]  
      pg_enriched <- pg_enriched[unlist(list(pg_enriched[,2]),use.names = F) > 0,]
      
      
    }else{
      
      pg_enriched <- pg_enriched[,grep(pattern = paste0("difference|", column_name),x = colnames(pg_enriched))]  
      pg_enriched <- pg_enriched[unlist(list(pg_enriched[,2]),use.names = F) < 0,]
      
    }
    
    
    enriched_ID <- unlist(list(unique(pg_enriched[,grep(pattern = column_name,x = colnames(pg_enriched))])),use.names = F)
    
    
    if (filter == T) {
      
      filter_IDs <- read_delim("filter.txt", 
                               "\t", 
                               escape_double = FALSE, trim_ws = TRUE,col_names = F)
      
      filter_IDs <- filter_IDs$X1
      
      enriched_ID <- enriched_ID[enriched_ID %in% filter_IDs]
      
    }
    
    ##Write out a table in "FASTA" format of the IDs:
    
    write.table( paste0(">", enriched_ID), "enriched_ID_FASTA.txt", append = FALSE, sep = " ", dec = ".",
                 row.names = F, col.names = F, quote = F)
    
    write.table( enriched_ID, "enriched_ID.txt", append = FALSE, sep = " ", dec = ".",
                 row.names = F, col.names = F, quote = F)
    
    #Retrieve the fold change
    
    fold_change <- as.matrix(pg_enriched)
    rownames(fold_change) <- unlist(list(fold_change[,1]),use.names = F)
    fold_change <- fold_change[,-1]
    fold_change <- apply(t(fold_change), 2,as.numeric)
    
    #Enrichment
    
    ## GO terms enrichment:
    
    genenameDeId   <- bitr(geneID = enriched_ID,OrgDb = org, fromType = "ENSEMBL",toType = "GENENAME")
    genenameUnivId <- bitr(geneID = univers_ID,OrgDb = org,fromType = "ENSEMBL",toType = "GENENAME")
    
    genenameDeId$FoldChange <- genenameDeId$ENSEMBL
    
    if (nrow(genenameDeId) > 0) {
      
      for (i in 1:nrow(genenameDeId)) {
        
        genenameDeId$FoldChange[i] <- unique(unlist(list(pg_enriched[unlist(list(pg_enriched[,1]),use.names = F) %in% genenameDeId$ENSEMBL[i],2]),use.names =  F))
      }
      
      fold_change_genename <- as.matrix(genenameDeId)
      rownames(fold_change_genename) <- genenameDeId$GENENAME
      fold_change_genename <- fold_change_genename[,-c(1:2)]
      fold_change_genename <- apply(t(fold_change_genename), 2,as.numeric)
      
    }
    
    enriched_GO <- enrichGO(gene = genenameDeId$GENENAME, 
                            OrgDb= org , 
                            keyType="GENENAME", 
                            ont= ont, 
                            universe = genenameUnivId$GENENAME,
                            pAdjustMethod = "BH", 
                            pvalueCutoff = 0.05)         
    
    ### Full
    
    if (!is.null(enriched_GO)) {
      
      if (length(enriched_GO@result[["p.adjust"]][enriched_GO@result[["p.adjust"]] < 0.05]) > 0) {
        
        file <- paste0(gsub(pattern = "-",replacement = "",x = Sys.Date()),"_", stock, "_GO_",ont, "_", direction, "_", experiment, "_table.csv")
        
        write.csv(as.data.frame(enriched_GO),
                  file= file,
                  row.names = F)
        
        file <- paste0(gsub(pattern = "-",replacement = "",x = Sys.Date()),"_", stock, "_GO_",ont, "_", direction, "_", experiment, "_barplot.pdf")
        
        barplot <- barplot(enriched_GO, showCategory=30)+
          scale_fill_viridis(begin = 0.4,end = 0.9,direction = 1)+
          ylab("Gene count")+
          xlab("GO term")+
          ggtitle(paste0(genename," ",direction, "-regulated", " GO ", ont, " barplot"))+
          theme_minimal()
        
        print(barplot)
        
        ggsave(plot = barplot, file, height = 10, width = 15)
        
        file <- paste0(gsub(pattern = "-",replacement = "",x = Sys.Date()),"_", stock, "_GO_",ont, "_", direction, "_", experiment, "_dotplot.pdf")
        
        dotplot <- dotplot(enriched_GO, showCategory=30)+
          scale_color_viridis(begin = 0.4,end = 0.9,direction = 1)+
          xlab("Gene ratio")+
          ylab("GO term")+
          ggtitle(paste0(genename," ",direction, "-regulated", " GO ", ont, " dotplot"))+
          theme_minimal()
        
        print(dotplot) 
        
        ggsave(plot = dotplot, file, height = 10, width = 15)
        
        file <- paste0(gsub(pattern = "-",replacement = "",x = Sys.Date()),"_", stock, "_GO_",ont, "_", direction, "_", experiment, "_cnetplot.pdf")
        
        cnetplot <- cnetplot(enriched_GO, foldChange = fold_change_genename, circular = F, colorEdge = T, showCategory = 30)+
          scale_color_viridis(begin = 0.4,end = 0.9,direction = 1)+
          labs(color = "FoldChange")+
          ggtitle(paste0(genename," ",direction, "-regulated", " GO ", ont, " cnetplot"))+
          theme(legend.position="top")
        
        print(cnetplot)
        
        ggsave(plot = cnetplot, file, height = 10, width = 20)
        
        file <- paste0(gsub(pattern = "-",replacement = "",x = Sys.Date()),"_", stock, "_GO_",ont, "_", direction, "_", experiment, "_heatplot.pdf")
        
        heatplot <- heatplot(enriched_GO,foldChange = fold_change_genename, showCategory = 30)+
          scale_fill_viridis(begin = 0.4,end = 0.9,direction = 1)+
          xlab("Gene name")+
          ylab("GO term")+
          ggtitle(paste0(genename," ",direction, "-regulated", " GO ", ont, " heatplot"))+
          labs(color = "FoldChange")
        
        print(heatplot)
        
        ggsave(plot = heatplot, file, height = 10, width = 15)
        
        file <- paste0(gsub(pattern = "-",replacement = "",x = Sys.Date()),"_", stock, "_GO_",ont, "_", direction, "_", experiment,  "_emapplot.pdf")
        
        if (length(enriched_GO@result[["p.adjust"]][enriched_GO@result[["p.adjust"]] < 0.05]) > 1) {
          
          emap <- enrichplot::pairwise_termsim(enriched_GO)  
          
          emapplot <- emapplot(emap,layout="nicely")+
            scale_color_viridis(begin = 0.4,end = 0.9,direction = 1)+
            ggtitle(paste0(genename," ",direction, "-regulated", " GO ", ont, " emapplot"))+
            labs(color = "FoldChange")
          
          print(emapplot)
          
          ggsave(plot = emapplot, file, height = 10, width = 15)
          
        }
        
        file <- paste0(gsub(pattern = "-",replacement = "",x = Sys.Date()),"_", stock, "_GO_",ont, "_", direction, "_", experiment,  "_GOgraph.pdf")
        
        pdf(file = file, height = 10, width = 15)
        plotGOgraph(enriched_GO)
        dev.off()
        
        print(plotGOgraph(enriched_GO))
        
        file <- paste0(gsub(pattern = "-",replacement = "",x = Sys.Date()),"_", stock, "_GO_",ont, "_", direction, "_", experiment,  "_upsetplot.pdf")
        
        upsetplot <- enrichplot::upsetplot(enriched_GO)+
          ggtitle(paste0(genename," ",direction, "-regulated", " GO ", ont, " upsetplot"))
        
        print(upsetplot)
        
        ggsave(plot = upsetplot, file, height = 10, width = 15)
        
      }
      
    }
    
    
    
    ### Simplified
    
    if (simplified) {
      
      if (length(enriched_GO@result[["p.adjust"]][enriched_GO@result[["p.adjust"]] < 0.05]) > 30) {
        
        reduced_enriched_GO <- clusterProfiler::simplify(enriched_GO, cutoff=0.7, by="p.adjust", select_fun=min)
        
        file <- paste0(stock, "_GO_",ont, "_", direction, "_", experiment, "_table_simplified.csv")
        
        write.csv(as.data.frame(reduced_enriched_GO),
                  file= file,
                  row.names = F)
        
        file <- paste0(stock, "_GO_",ont, "_", direction, "_", experiment,  "_barplot_simplified.pdf")
        
        barplot(reduced_enriched_GO, showCategory=30)+
          scale_fill_viridis(begin = 0.4,end = 0.9,direction = 1)+
          ylab("Gene count")+
          xlab("GO term")+
          theme_minimal()+
          ggsave(file, height = 10, width = 15)
        
        file <- paste0(stock, "_GO_",ont, "_", direction, "_", experiment,  "_dotplot_simplified.pdf")
        
        dotplot(reduced_enriched_GO, showCategory=30)+
          scale_color_viridis(begin = 0.4,end = 0.9,direction = 1)+
          xlab("Gene ratio")+
          ylab("GO term")+
          theme_minimal()+
          ggsave(file, height = 10, width = 15)
        
        file <- paste0(stock, "_GO_",ont, "_", direction, "_", experiment,  "_cnetplot_simplified.pdf")
        
        cnetplot(reduced_enriched_GO, foldChange = fold_change_genename, circular = F, colorEdge = T,showCategory = 30)+
          scale_color_viridis(begin = 0.4,end = 0.9,direction = 1)+
          labs(color = "FoldChange")+
          theme(legend.position="top")+
          ggsave(file, height = 10, width = 20)
        
        file <- paste0(stock, "_GO_",ont, "_", direction, "_", experiment,  "_heatplot_simplified.pdf")
        
        heatplot(reduced_enriched_GO,foldChange = fold_change_genename,showCategory = 30)+
          scale_fill_viridis(begin = 0.4,end = 0.9,direction = 1)+
          xlab("Gene name")+
          ylab("GO term")+
          labs(color = "FoldChange")+
          ggsave(file, height = 10, width = 15)
        
        file <- paste0(stock, "_GO_",ont, "_", direction, "_", experiment,  "_emapplot_simplified.pdf")
        
        if (length(enriched_GO@result[["p.adjust"]][enriched_GO@result[["p.adjust"]] < 0.05]) > 1) {
          
          emap <- enrichplot::pairwise_termsim(enriched_GO)
          
          emapplot(emap,layout="nicely",showCategory = 30)+
            scale_color_viridis(begin = 0.4,end = 0.9,direction = 1)+
            labs(color = "FoldChange")+
            ggsave(file, height = 10, width = 15)
          
        }
        
        file <- paste0(stock, "_GO_",ont, "_", direction, "_", experiment,  "_GOgraph_simplified.pdf")
        
        pdf(file = file, height = 10, width = 15)
        plotGOgraph(reduced_enriched_GO)
        dev.off()
        
        file <- paste0(stock, "_GO_",ont, "_", direction, "_", experiment,  "_upsetplot_simplified.pdf")
        
        enrichplot::upsetplot(reduced_enriched_GO)+
          ggsave(file, height = 10, width = 15)
        
      }
      
    }
    
  }
  
  ## Reactome enrichment
  
  allowed <- c("human", "rat", "mouse", "celegans", "yeast", "zebrafish", "fly")
  
  if (org_reactome %in% allowed == T) {
    
    entrezDeId   <- bitr(geneID = enriched_ID,OrgDb = org, fromType = "ENSEMBL",toType = "ENTREZID")
    entrezUnivId <- bitr(geneID = univers_ID,OrgDb = org,fromType = "ENSEMBL",toType = "ENTREZID")
    
    entrezDeId$FoldChange <- entrezDeId$ENSEMBL
    
    if (nrow(entrezDeId) > 0) {
      
      for (i in 1:nrow(entrezDeId)) {
        
        entrezDeId$FoldChange[i] <- unique(unlist(list(pg_enriched[unlist(list(pg_enriched[,1]),use.names = F) %in% entrezDeId$ENSEMBL[i],2]),use.names =  F))
      }
      
      fold_change_entrez <- as.matrix(entrezDeId)
      rownames(fold_change_entrez) <- entrezDeId$ENTREZID
      fold_change_entrez <- fold_change_entrez[,-c(1:2)]
      fold_change_entrez <- apply(t(fold_change_entrez), 2,as.numeric)
      
    }
    
    enriched_Reactome <- enrichPathway( gene = entrezDeId$ENTREZID, 
                                        organism = org_reactome , 
                                        universe = entrezUnivId$ENTREZID,
                                        pvalueCutoff = 0.05)
    
    if (!is.null(enriched_Reactome)) {
      
      if (length(enriched_Reactome@result[["p.adjust"]][enriched_Reactome@result[["p.adjust"]] < 0.05]) > 0) {
        
        file <- paste0(gsub(pattern = "-",replacement = "",x = Sys.Date()),"_", stock, "_Reactome_", direction, "_", experiment, "_table.csv")
        
        write.csv(as.data.frame(enriched_Reactome),
                  file= file,
                  row.names = F)
        
        file <- paste0(gsub(pattern = "-",replacement = "",x = Sys.Date()),"_", stock, "_Reactome_", direction, "_", experiment,  "_barplot.pdf")
        
        barplot <- barplot(enriched_Reactome, showCategory=30)+
          scale_fill_viridis(begin = 0.4,end = 0.9,direction = 1)+
          theme_minimal()+
          ylab("Gene count")+
          ggtitle(paste0(genename," ",direction, "-regulated", " Reactome barplot"))+
          xlab("Reactome pathway")
        
        print(barplot)
        
        ggsave(plot = barplot, filename = file, height = 10, width = 15)
        
        file <- paste0(gsub(pattern = "-",replacement = "",x = Sys.Date()),"_", stock, "_Reactome_", direction, "_", experiment,  "_dotplot.pdf")
        
        dotplot <- dotplot(enriched_Reactome, showCategory=30)+
          scale_color_viridis(begin = 0.4,end = 0.9,direction = 1)+
          xlab("Gene ratio")+
          ylab("Reactome pathway")+
          ggtitle(paste0(genename," ",direction, "-regulated", " Reactome dotplot"))+
          theme_minimal()
        
        print(dotplot)
        
        ggsave(filename = file, plot = dotplot, height = 10, width = 15)
        
        file <- paste0(gsub(pattern = "-",replacement = "",x = Sys.Date()),"_", stock, "_Reactome_", direction, "_", experiment,  "_cnetplot.pdf")
        
        cnetplot <- cnetplot(enriched_Reactome, foldChange = fold_change_entrez, circular = F, colorEdge = T)+
          scale_color_viridis(begin = 0.4,end = 0.9,direction = 1)+
          labs(color = "FoldChange")+
          ggtitle(paste0(genename," ",direction, "-regulated", " Reactome cnetplot"))+
          theme(legend.position="top")
        
        print(cnetplot)
        
        ggsave(filename = file, plot = cnetplot, height = 10, width = 20)
        
        file <- paste0(gsub(pattern = "-",replacement = "",x = Sys.Date()),"_", stock, "_Reactome_", direction, "_", experiment,  "_heatplot.pdf")
        
        heatplot <- heatplot(enriched_Reactome, foldChange = fold_change_entrez)+
          scale_fill_viridis(begin = 0.4,end = 0.9,direction = 1)+
          xlab("entrez ID")+
          ylab("Reactome pathway")+
          ggtitle(paste0(genename," ",direction, "-regulated", " Reactome heatplot"))+
          labs(color = "FoldChange")
        
        print(heatplot)
        
        ggsave(filename = file, plot = heatplot , height = 10, width = 15)
        
        if (length(enriched_Reactome@result[["p.adjust"]][enriched_Reactome@result[["p.adjust"]] < 0.05]) > 1) {
          
          emap <- enrichplot::pairwise_termsim(enriched_Reactome)
          
          file <- paste0(gsub(pattern = "-",replacement = "",x = Sys.Date()),"_", stock, "_Reactome_", direction, "_", experiment,  "_emapplot.pdf")
          
          emapplot <- emapplot(emap,layout="nicely")+
            scale_color_viridis(begin = 0.4,end = 0.9,direction = 1)+
            ggtitle(paste0(genename," ",direction, "-regulated", " Reactome emapplot"))+
            labs(color = "FoldChange")
          
          print(emapplot)
          
          ggsave(filename = file, plot = emapplot, height = 10, width = 15)
          
        }
        
        file <- paste0(gsub(pattern = "-",replacement = "",x = Sys.Date()),"_", stock, "_Reactome_", direction, "_", experiment,  "_upsetplot.pdf")
        
        upsetplot <- enrichplot::upsetplot(enriched_Reactome)+
          ggtitle(paste0(genename," ",direction, "-regulated", " Reactome upsetplot"))
        
        print(upsetplot)
        
        ggsave(filename = file, plot = upsetplot, height = 10, width = 15)
        
      }
      
    }
  }
  
  ## KEGG enrichment
  
  enriched_KEGG <- enrichKEGG( gene = entrezDeId$ENTREZID, 
                               organism = org_KEGG, 
                               universe = entrezUnivId$ENTREZID,
                               pvalueCutoff = 0.05,
                               keyType = "ncbi-geneid")
  
  if (!is.null(enriched_KEGG)) {
    
    
    if (length(enriched_KEGG@result[["p.adjust"]][enriched_KEGG@result[["p.adjust"]] < 0.05]) > 0) {
      
      file <- paste0(gsub(pattern = "-",replacement = "",x = Sys.Date()),"_", stock, "_KEGG_", direction, "_", experiment, "_table.csv")
      
      write.csv(as.data.frame(enriched_KEGG),
                file= file,
                row.names = F)
      
      file <- paste0(gsub(pattern = "-",replacement = "",x = Sys.Date()),"_", stock, "_KEGG_", direction, "_", experiment,  "_barplot.pdf")
      
      barplot <- barplot(enriched_KEGG, showCategory=30)+
        scale_fill_viridis(begin = 0.4,end = 0.9,direction = 1)+
        ylab("Gene count")+
        xlab("KEGG term")+
        ggtitle(paste0(genename," ",direction, "-regulated", " KEGG barplot"))+
        theme_minimal()
      
      print(barplot)
      
      ggsave(filename = file, plot = barplot, height = 10, width = 15)
      
      file <- paste0(gsub(pattern = "-",replacement = "",x = Sys.Date()),"_", stock, "_KEGG_", direction, "_", experiment,  "_dotplot.pdf")
      
      dotplot <- dotplot(enriched_KEGG, showCategory=30)+
        scale_color_viridis(begin = 0.4,end = 0.9,direction = 1)+
        xlab("Gene ratio")+
        ylab("KEGG term")+
        ggtitle(paste0(genename," ",direction, "-regulated", " KEGG dotplot"))+
        theme_minimal()
      
      print(dotplot)
      
      ggsave(filename = file, plot = dotplot, height = 10, width = 15)
      
      file <- paste0(gsub(pattern = "-",replacement = "",x = Sys.Date()),"_", stock, "_KEGG_", direction, "_", experiment,  "_cnetplot.pdf")
      
      cnetplot <- cnetplot(enriched_KEGG, foldChange = fold_change_entrez, circular = F, colorEdge = T)+
        scale_color_viridis(begin = 0.4,end = 0.9,direction = 1)+
        labs(color = "FoldChange")+
        ggtitle(paste0(genename," ",direction, "-regulated", " KEGG dotplot"))+
        theme(legend.position="top")
      
      print(cnetplot)
      
      ggsave(filename = file, plot = cnetplot, height = 10, width = 20)
      
      file <- paste0(gsub(pattern = "-",replacement = "",x = Sys.Date()),"_", stock, "_KEGG_", direction, "_", experiment,  "_heatplot.pdf")
      
      heatplot <- heatplot(enriched_KEGG,foldChange = fold_change_entrez)+
        scale_fill_viridis(begin = 0.4,end = 0.9,direction = 1)+
        xlab("entrez ID")+
        ylab("KEGG term")+
        ggtitle(paste0(genename," ",direction, "-regulated", " KEGG heatplot"))+
        labs(color = "FoldChange")
      
      print(heatplot)
      
      ggsave(filename = file, plot = heatplot, height = 10, width = 15)
      
      if (length(enriched_KEGG@result[["p.adjust"]][enriched_KEGG@result[["p.adjust"]] < 0.05]) > 1) {
        
        file <- paste0(gsub(pattern = "-",replacement = "",x = Sys.Date()),"_", stock, "_KEGG_", direction, "_", experiment,  "_emapplot.pdf")
        
        emap <- enrichplot::pairwise_termsim(enriched_KEGG)
        
        eamapplot <- emapplot(emap,layout="nicely")+
          scale_color_viridis(begin = 0.4,end = 0.9,direction = 1)+
          ggtitle(paste0(genename," ",direction, "-regulated", " KEGG emapplot"))+
          labs(color = "FoldChange")
        emapplot
        
        ggsave(filename = file, plot = emapplot , height = 10, width = 15)
        
      }
      
      file <- paste0(gsub(pattern = "-",replacement = "",x = Sys.Date()),"_", stock, "_KEGG_", direction, "_", experiment,  "_upsetplot.pdf")
      
      upsetplot <- enrichplot::upsetplot(enriched_KEGG)+
        ggtitle(paste0(genename," ",direction, "-regulated", " KEGG upsetplot"))
      
      print(upsetplot)
      
      ggsave(filename = file, plot = upsetplot, height = 10, width = 15)
      
    }
    
  }

}

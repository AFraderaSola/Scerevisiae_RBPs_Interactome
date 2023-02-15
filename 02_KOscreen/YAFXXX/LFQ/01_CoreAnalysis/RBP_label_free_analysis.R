###############################################################################
#
# This script if for analysing label free data sets and do a basic first 
# analysis. 
#
# developer: Mario Dejung <m.dejung@imb.de>
# version: 0.6.7
# date: 2019.05.29
#
# package version: 0.4.10
# package date: 2018-05-25
#
# The basic script has been optimized for the project RBP interactome project.
# developer: Albert Fradera-Sola <A.FraderaSola@imb.de>
#
# date: 2022.04.27
###############################################################################


#############################################
############# Libraries #####################
#############################################

library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(gplots)
library(dplyr)
library(ggrepel)
library(IHW)
# library(plotly)
library(coin)
# options(java.parameters = "-Xmx4096m")
options(java.parameters=c("-XX:+UseConcMarkSweepGC", "-Xmx4096m"))
library(xlsx)
# devtools::install_github('talgalili/heatmaply')
library(heatmaply)
library(rMQanalysis)

this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)
print(getwd())
rMQanalysis::copySourcedScript()

#############################################
############# Functions #####################
#############################################


Enriched2 <- function (x, y = NULL, c = 1, s0 = 1, pvalue = 0.05, ...) 
{
  if (class(x) == "data.frame") {
    if (length(x) == 2) {
      y <- x[[2]]
      x <- x[[1]]
    }
    else {
      stop("x hast to be a data frame with length 2!")
    }
  }
  results <- list()
  results$curve <- createEnrichmentCurve2(c, s0, pvalue, ...)
  results$positive <- (c/(x - s0) - y + -log10(pvalue) <= 0 & 
                         x >= s0)
  results$enriched <- results$positive 
  results$background <- !results$enriched
  class(results) <- append(class(results), "Enriched")
  return(results)
}

createEnrichmentCurve2 <-  function (c, s0, pvalue = 1, step = 0.01, xlim = 50, linetype = 2, 
          size = 0.1) 
{
  x <- seq(s0, xlim, step)
  y <- c/(x - s0) + -log10(pvalue)
  return(annotate("path", x = c(NA, x), y = c(NA, y), 
                  linetype = linetype, size = size))
}

#############################################
############# Variables #####################
#############################################

# Main variable settings --------------------------------------------------

rMQanalysis::mylog('Setting main variables', dbg_level=0)
out_dir <- paste0(gsub(pattern = "-",replacement = "",x = Sys.Date()), "_results_LFQ") # output directory for any written file
image_name <- paste0(gsub(pattern = "-",replacement = "",x = Sys.Date()), "_LFQ_Analysis.RData") # output .RData image
overwrite_out_dir <- TRUE

# Proteins with extra highlight (default is RNaseA and bait/target RBP)

masterfile <- read.xlsx(file = "../../../00_IPScreeningMasterFile.xlsx",header = T,sheetIndex = 1)
extra_highlight <- c("RNaseA", paste0(masterfile[masterfile$Stock_Name == gsub(pattern = ".*een/|/LFQ.*",replacement = "",x = getwd()),]$Gene_Name, "-TAP"))
stock <- gsub(pattern = ".*een/|/LFQ.*",replacement = "",x = getwd())

# reading and filtering the protein groups
data_dir <- '.' # path to data directory, RELATIVE to THIS MARKDOWN document! 
# Set to "." if it is in the same directory! RECOMMENDED!
pg_filename <- 'proteinGroups.txt' # The protein groups file name
razor_peptides <- 2 # minimum razor+unique peptide count for the whole protein group
unique_peptides <- 1 # minimum unique peptide count for the whole protein group
min_quant_events <- 2 # at least measured X times in any experiment replicate

# we create a link to a direct database search 
database_link <- 'ensembl' # Can be uniprot, ensembl, none or own

# write output files
pg_basename <- sub('.txt$', '', pg_filename)
pg_quant_filename <- paste0(stock, '_quantified_', pg_basename, '.txt') 
pg_volcano_filename <- paste0(stock, '_volcano_', pg_basename, '.txt')


# Imputation settings -----------------------------------------------------

impute_method <- 'beta' # can be "beta", "norm" or "original"
beta_min <- 0.002 # Check ?imputeBetaRandomNumbersPerColumn
beta_max <- 0.025
norm_modify_mean <- -4 # Check ?imputeNormRandomNumbersPerColumn
norm_sd_factor <- .2
random_seed <- 1234
original_na_replace <- 0.001
original_average_function <- 'mean'


# Contrast setup and labels -----------------------------------------------

# creating plots for different experiments
compare_to_experiment <- NULL # the standard we compare to (index of experiments). 
# Set to NULL to compare all. You can check this after executing the first few chunks and
# print(experiments) to see the order.
reverse_comparison <- TRUE # Set to TRUE to get the control on the right side of the
# volcano plots.

# using the new getLabel function. Check ?getLabel for more information.
label_column <- 'Gene.names' # set to Gene.names if it is available
label_column_regex <- '' # set to '([^;]+).*' to use as is
fallback_column <- 'Protein.IDs'
fallback_column_regex <- '([^;]+).*'
# some regular expression examples
#   '([^ ]+).*' # only the first word of the header
#   '([^;]+).*' # only the first protein ID
#   '.*\\| (NCU\\d{5}) \\|.*' # using some special string for Neurospora crassa

# Creating correlation heatmap and PCA of the LFQ intensities.
heatmap_column_regex <- '^imputed.log2.LFQ.intensity.' 
heatmap_white_area_size <- 1 # 1 == same area than blue, 2 means the white area 
# is 2 times bigger then the blue area, or the white are is 66% of the spectra
pca_column_regex <- '^imputed.log2.LFQ.intensity.' 
# use '^log2.LFQ.intensity' to create heatmap on measured without imputed values

# highlight important protein IDs. Provide a textfile in the format as highlight.txt!
highlight_file <- 'highlight.txt'


# Volcano plot setup ------------------------------------------------------

do_volcano <- TRUE # switch of volcano plot drawing.
create_volcano_pdfs <- TRUE # If you want to create PDFs
volcano_height <- 10 # you might change the pdf size together with the font size
volcano_width <- 10
font_size <- 3
volcano_threshold_linetype <- 2
volcano_threshold_linesize <- 0.5
create_volcano_enriched_files <- TRUE # write the enrichment files
# To calculate the difference between two conditions, we substract the mean or median 
# values. There are several possibilities ("median/mean2_naimp_meas_", "median/mean1_meas_" 
# and "mean/median3_less1_imp_"). Check chapter "Calculating mean and median values" to 
# get a description about the column meaning.
difference_volcano_column <- 'mean2_naimp_meas_'
# the threshold we use for enrichment
enrich_c <- .05 # increase to make it more "round" or decrease to make it more "edgy" :-)
enrich_s0 <- log2(2) # two fold change
enrich_pvalue <- 0.05

pvalue_correction <- 'no' # can be "no", "permutation" or "IHW"
ihw_test <- FALSE

ggrepel_threshold <- 10

filter_separately <- TRUE # filters before each volcano plot the quantified proteins
# for the specific experiment comparison. So subseting to at least one of the conditions
# has more or equal min_quant_events.


# Scatter plot setup ------------------------------------------------------

# Scatter plot settings to write PDF files
do_scatter <- FALSE
scatter_upper_threshold <- 1
scatter_lower_threshold <- 1
scatter_density <- .2
scatter_size <- 1.5
# median1_meas_ will plot the median measured values without imputed, labeling might fail
# median3_less1_imp_ will plot the median measured values and missing values get imputed
difference_scatter_column <- 'mean1_meas_' 
create_scatter_enriched_files <- TRUE
create_scatter_pdfs <- TRUE
scatter_height <- 8 # you might change the pdf size together with the font size
scatter_width <- 8


#############################################
############### Script ######################
#############################################

# Create output directory -------------------------------------------------

if(!file.exists(out_dir)) {
  dir.create(out_dir)
} else {
  if(!overwrite_out_dir) {
    out_dir <- paste0(out_dir, '_', format(Sys.time(), '%Y%m%d%H%M%S'))
    dir.create(out_dir)
  }
}


# Reading protein groups files --------------------------------------------

rMQanalysis::mylog('Reading protein groups', dbg_level=0)
pg <- read.delim(file.path(data_dir,pg_filename))

pg$Gene.names <-
  ifelse(grepl(".*SGDID:([^ ]*)", pg$Fasta.headers),
         sub(".* ", "",sub(" SGDID:.*", '', pg$Fasta.headers)),
         sub(".*", "NO_ID", pg$Fasta.headers))


# Renaming individual experiments -----------------------------------------

# names(pg) <- sub('RK_0308', 'Rexo4_1', names(pg))
# names(pg) <- sub('RK_0309', 'Cont_1', names(pg))
# names(pg) <- sub('RK_0310', 'Rexo4_2', names(pg))
# names(pg) <- sub('RK_0311', 'Cont_2', names(pg))
# names(pg) <- sub('RK_0312', 'Rexo4_3', names(pg))
# names(pg) <- sub('RK_0313', 'Cont_3', names(pg))
# names(pg) <- sub('RK_0314', 'Rexo4_4', names(pg))
# names(pg) <- sub('RK_0315', 'Cont_4', names(pg))
# # Check if names make sense!
# grep('^LFQ', names(pg), value=TRUE)


# Hacky way to extract some information from Fasta headers ----------------

# pg$description <- 
#   sapply(
#     strsplit(as.character(pg$Fasta.headers), ';'), function(x) {
#       paste(unname(
#         sub('^([^=|\\[]+).*', 
#             '\\1', 
#             sapply(x, function(y) {
#               paste(strsplit(y, ' ')[[1]][-1], collapse=' ')
#             }))
#       ), collapse=';')
#     })



# Merge some additional data to protein IDs -------------------------------

# rMQanalysis::mylog('merging KOG data to protein groups', dbg_level=0)
# pg_split <- splitColumnBySep(pg[c('Protein.IDs', 'id')], 'Protein.IDs')
# 
# kog_data <- xlsx::read.xlsx('annotations/JGI_zu_KOG_Daphnia_pulex.xlsx', 1)
# kog_annot <- xlsx::read.xlsx('annotations/JGI_zu_KOG_Daphnia_pulex.xlsx', 2)
# 
# dappu_pgs <- pg_split[grep('\\|', pg_split$Protein.IDs),]
# dappu_pgs$GENE.JGI_V11 <- sub('^(\\d+)\\|.*', 
#                               'GENE:JGI_V11_\\1', 
#                               dappu_pgs$Protein.IDs)
# dappu_pgs <- merge(dappu_pgs, kog_data,
#                    all.x=TRUE,
#                    sort=FALSE)
# dappu_pgs <- merge(dappu_pgs, kog_annot,
#                    by.x='KOG.ID', by.y='KOG_JGI',
#                    all.x=TRUE,
#                    sort=FALSE)
# 
# sum(!is.na(dappu_pgs$KOG.ID))
# 
# dappu_pgs <- 
#   dappu_pgs %>%
#   group_by(id) %>%
#   summarise(KOG.ID=paste(unique(KOG.ID), collapse=';'),
#             KOG.descr=paste(unique(descr), collapse=';'),
#             KOG.cat=paste(Category, collapse=';'))
# 
# pg_names <- names(pg)
# dappu_names <- names(dappu_pgs)[-1]
# 
# pg <- merge(pg, dappu_pgs,
#             all.x=TRUE,
#             sort=FALSE)
# pg <- pg[c(pg_names, dappu_names)]


rMQanalysis::mylog('Filtering protein groups', dbg_level=0)
pg_flt <- filterWholeDataset(pg)
pg_ident <- filterIdentifiedProteins(pg_flt, razor_peptides, unique_peptides)

experiments <- unique(gsub('^LFQ\\.intensity\\.(.*)_\\d+',
                           '\\1',
                           grep('^LFQ',names(pg),value=TRUE)))

if(any(grepl('\\.', experiments))) {
  experiments_save <- gsub('\\.', '', experiments)
  for(i in seq_along(experiments)) {
    names(pg) <- gsub(experiments[[i]], experiments_save[[i]], names(pg))
    names(pg_flt) <- gsub(experiments[[i]], experiments_save[[i]], names(pg_flt))
    names(pg_ident) <- gsub(experiments[[i]], experiments_save[[i]], names(pg_ident))
  }
  experiments <- unique(gsub('^LFQ\\.intensity\\.(.*)_\\d+',
                             '\\1',
                             grep('^LFQ',names(pg),value=TRUE)))
  warning('we had to remove "." from the experiment names.', call.=FALSE)
}

replicates <- unique(gsub('^LFQ\\.intensity\\..*_(\\d+)',
                          '\\1',
                          grep('^LFQ',names(pg),value=TRUE)))
rMQanalysis::mylog(
  sprintf('detected experiments: %s', paste(experiments, collapse=', ')), 
  dbg_level=0)


lfq_columns <- grep('^LFQ', names(pg_ident))
pg_ident[lfq_columns][ pg_ident[lfq_columns] == 0 ]  <- NA
pg_ident <- 
  cbind(pg_ident,
        log2=log2(pg_ident[lfq_columns]))



# Impute missing values ---------------------------------------------------

rMQanalysis::mylog('Imputation of missing values', dbg_level=0)
log2_lfq_columns <- grep('^log2.LFQ', names(pg_ident))
imputed_values <- 
  switch(impute_method,
         beta=imputeBetaRandomNumbersPerColumn(
           pg_ident[log2_lfq_columns],
           beta_min,beta_max,
           seed=random_seed),
         norm=imputeNormRandomNumbersPerColumn(
           pg_ident[log2_lfq_columns],
           modify_mean=norm_modify_mean,
           sd_factor=norm_sd_factor,
           seed=random_seed),
         original=imputeOriginalDistributionPerExperiment(
           pg_ident[log2_lfq_columns],
           na_replace=original_na_replace,
           average_function=original_average_function,
           seed=random_seed),
         stop('please set a correct impute_method!'))

pg_ident <- cbind(pg_ident,imputed_values$imputed_values)


# Missing values plot -----------------------------------------------------

rMQanalysis::mylog('Plot missing values diagnose graphs', dbg_level=0)
imputed_values$imputed_histogram$column <- 
  sub('log2.LFQ.intensity.', '', imputed_values$imputed_histogram$column)
imputed_values$imputed_histogram$column <- 
  sub('_', ' repl.:', imputed_values$imputed_histogram$column)
graph_imputed_values <- plotImputedValues(imputed_values$imputed_histogram, ncol=length(replicates)) +
  xlab('log2( LFQ intensity )') + theme_imb()
ggsave(file.path(out_dir, paste0(stock, '_imputed_values.pdf')), graph_imputed_values,
       width=8.27, height=11.69)
print(graph_imputed_values)

imputed_bar_df <- 
  imputed_values$imputed_histogram %>% 
  dplyr::group_by(column, value_type) %>% 
  dplyr::summarise(n=length(value))
graph_imputed_values_bar <- 
  ggplot(imputed_bar_df, aes(column, n, fill=value_type)) +
  geom_bar(stat='identity', position=position_dodge(), na.rm=TRUE) +
  ggtitle('comparison imputed vs. measured') +
  ylab('count') +
  theme_imb() +
  scale_fill_brewer('', palette='Set1') +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_text(angle=60, hjust=1, vjust=1))
ggsave(file.path(out_dir, paste0(stock, '_imputed_values_bar.pdf')), graph_imputed_values_bar,
       width=8.27, height=5.85)
print(graph_imputed_values_bar)

if(impute_method == 'original') {
  print(imputed_values$imputing_data[c(1,4,5)])
  imputed_data2 <- merge(imputed_values$dist_hist, imputed_values$imputing_data)
  imputed_data2$title <- sprintf('%s: %s NAs imputed with distribution', imputed_data2$experiment, imputed_data2$NA_values)
  repl_spread_graph <- ggplot(imputed_data2, aes(values)) + 
    geom_histogram(binwidth=.15, na.rm=TRUE) + facet_wrap(~ title) +
    ggtitle('replicate distribution per timepoint and # of NAs to impute') +
    theme_imb()
  ggsave(file.path(out_dir, paste0(stock, '_replicate_spread.pdf')), repl_spread_graph,
         height=8.27, width=11.69)
  print(repl_spread_graph)
}



# Compare Intensity to LFQ values -----------------------------------------

rMQanalysis::mylog('Create Intensity to LFQ plot', dbg_level=0)
intensity_comp_graph <- rMQanalysis::plotIntensityLFQComparison(pg_ident)
ggsave(file.path(out_dir, paste0(stock, '_intensity_comparison.pdf')),
       intensity_comp_graph, height=11.69, width=8.27)
print(intensity_comp_graph)



# Detect proteins to highlight --------------------------------------------


rMQanalysis::mylog('Highlight important proteins', dbg_level=0)
pg_ident$highlight <- FALSE
ids_to_highlight <- c()
if(file.exists(highlight_file)) ids_to_highlight <- read.delim(highlight_file, as.is=T)[[1]]
if(length(ids_to_highlight) > 0) {
  ids_searchstring <- paste(ids_to_highlight, collapse='|')
  matching_rows <- grep(ids_searchstring, as.character(pg_ident$Protein.IDs))
  pg_ident$highlight[matching_rows] <- TRUE
}



# Create a label column ---------------------------------------------------

rMQanalysis::mylog('Convert gene names and protein ids for labeling', dbg_level=0)
pg_ident$my_label <- 
  getLabel(pg_ident, 
           label_column, fallback_column, 
           label_column_regex, fallback_column_regex)



# Abundance graphs --------------------------------------------------------

rMQanalysis::mylog('Creating abundace graphs', dbg_level=0)
graph_abundance_rank <- rMQanalysis::plotAbundanceRank(pg_ident, 
                                                       label_column='my_label', 
                                                       subset_column='highlight')
ggsave(file.path(out_dir, paste0(stock, '_abundance_rank.pdf')),
       graph_abundance_rank, height=6, width=7)
print(graph_abundance_rank)

abundance_rank_levels <- 
  paste0(c('', 'LFQ.'), 
         rep(sub('Intensity.', 
                 '', 
                 sort(grep(
                   paste0('Intensity.', experiments, collapse='|'), 
                   names(pg), value=T)
                 )), 
             each=2))
graph_individual_abundance_rank <- 
  pg_ident %>%
  dplyr::select(my_label, highlight, 
                starts_with('LFQ.intensity'), 
                starts_with('Intensity.')) %>%
  tidyr::gather(experiment, intensity, -my_label, -highlight) %>%
  as_tibble() %>%
  mutate(intensity=ifelse(is.na(intensity), 0, intensity),
         experiment=sub('^LFQ.intensity.','LFQ.',experiment),
         experiment=sub('^Intensity.','',experiment),
         experiment=factor(experiment, 
                           levels=abundance_rank_levels)) %>%
  group_by(experiment) %>%
  dplyr::arrange(-intensity) %>%
  mutate(proteinGroups=seq(n())) %>%
  ggplot(aes(proteinGroups, log10(intensity))) +
  # geom_smooth() +
  ggtitle('Protein abundance per replicate of LFQ values') +
  geom_point(data=. %>% filter(!highlight), size=1, color='grey') +
  geom_point(data=. %>% filter(highlight), size=2, color='red') +
  geom_text_repel(data=. %>% filter(highlight), 
                  aes(label=my_label), 
                  size=2, color='black') +
  facet_wrap(~ experiment, ncol=length(replicates)*2) +
  theme(axis.text.x=element_text(angle=45, hjust=1))

ggsave(file.path(out_dir, paste0(stock, '_individual_abundance_rank.pdf')),
       graph_individual_abundance_rank, 
       height=min(8.27 * length(experiments), 49), 
       width=11.69)
print(graph_individual_abundance_rank)

# Calculating averages ----------------------------------------------------

rMQanalysis::mylog('Calculating average values', dbg_level=0)
pg_ident <- calculateAverages(as.data.frame(pg_ident))


# Filtering for minimum quantification events -----------------------------

rMQanalysis::mylog('Filtering for minimum quantification events', dbg_level=0)
pg_quant <- 
  pg_ident[apply(
    pg_ident[grep('value_count',names(pg_ident))] >= min_quant_events,
    1,
    any),]



# Protein counts graph ----------------------------------------------------

rMQanalysis::mylog('Create protein counts graph', dbg_level=0)
pg_data <- rMQanalysis::countProteinGroups(pg, pg_flt, pg_ident)
ggplot2::ggsave(file.path(out_dir, paste0(stock, '_protein_counts.pdf')), pg_data$plot,
                width=6, height=5)
print(pg_data$plot)



# Value count graph -------------------------------------------------------

rMQanalysis::mylog('Create value count graph', dbg_level=0)
graph_value_count <- createValueCountPlot(pg_ident)
ggsave(file.path(out_dir, paste0(stock,'_value_counts.png')), graph_value_count,
       width=8.27, height=11.69, dpi=200)
print(graph_value_count)



# # Write filtered protein groups file --------------------------------------
# 
# rMQanalysis::mylog('Writing protein groups file', dbg_level=0)
# write.table_imb(pg_quant,
#                 file.path(out_dir, pg_quant_filename))



# Correlation heatmap -----------------------------------------------------

rMQanalysis::mylog('Create correlation heatmap', dbg_level=0)
rMQanalysis::plotReplicateCorrelation(pg_quant, 
                                      heatmap_white_area_size, 
                                      heatmap_column_regex)()
pdf(file=file.path(out_dir, sprintf(paste0(stock, '_correlation_plot_%s.pdf'), pg_basename)),
    width=8, height=8, pointsize=8)

rMQanalysis::plotReplicateCorrelation(pg_quant, 
                                      heatmap_white_area_size, 
                                      heatmap_column_regex)()
dev.off()

pal <- c(rep("white",each=9 * heatmap_white_area_size), 
         brewer.pal(9,"Blues"))
lfq <- pg_quant[,grep(heatmap_column_regex, names(pg_quant))]
colnames(lfq) <- gsub(heatmap_column_regex,'',colnames(lfq))

group   <- as.factor(gsub("_[1-4]$","",colnames(lfq)))
lfq.avg <- t(apply(lfq,1,function(x) tapply(x,group,mean,na.rm=T)))
v <- apply(lfq,1,sd)
corr <- apply(lfq[rev(order(v)),],2,function(x) {
  apply(lfq[rev(order(v)),],2,function(y) {
    cor(x,y,use="na.or.complete")
  })
})

gplots::heatmap.2(corr,
                  trace="none",
                  Colv=T,
                  Rowv=T,
                  dendrogram="column",
                  density.info="density",
                  srtCol=45,
                  main='with clustering',
                  margins=c(10, 10), 
                  key.xlab="Pearson's correlation coefficient",
                  key.title="",
                  keysize=1.5,
                  col=pal)
pdf(file=file.path(out_dir, sprintf(paste0(stock, '_correlation_plot_cluster_%s.pdf'), pg_basename)),
    width=11.69, height=11.69, pointsize=8)

gplots::heatmap.2(corr,
                  trace="none",
                  Colv=T,
                  Rowv=T,
                  dendrogram="column",
                  density.info="density",
                  srtCol=45,
                  main='with clustering',
                  margins=c(10, 10), 
                  key.xlab="Pearson's correlation coefficient",
                  key.title="",
                  keysize=1.5,
                  col=pal)
dev.off()



# Principal component analysis --------------------------------------------

rMQanalysis::mylog('Create prinzipal component analysis', dbg_level=0)
pca_columns <- grep(pca_column_regex, colnames(pg_quant))
prcomp_list <- list(the=prcomp(t(na.omit(pg_quant[,pca_columns])), scale.=TRUE))
multi_pca_plot <- rMQanalysis::multiPCAPlot('the', prcomp_list, TRUE)

all_pca_filename <- file.path(out_dir, paste0(stock, '_all_pca_', pg_basename, '_prcomp.pdf'))
pdf(file=all_pca_filename, width=11.69, height=8.27)
grid::grid.newpage()
grid::grid.draw(multi_pca_plot)
jnk <- dev.off()

# Venn diagram ------------------------------------------------------------

if(length(experiments) <= 5) {
  rMQanalysis::mylog('Create venn diagram of quantified proteins', dbg_level=0)
  
  futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
  pg_quant %>% 
    as_tibble() %>% 
    select(id, starts_with('mean1')) %>% 
    tidyr::gather('experiment','identified', -id) %>% 
    mutate(experiment=sub('mean1_meas_','', experiment),
           identified=!is.na(identified)) %>%
    filter(identified) %>% 
    group_by(experiment) %>% 
    do(data = (.)) %>% 
    with(setNames(data, experiment) ) ->
    venn_list
  venn_1 <- VennDiagram::venn.diagram(
    lapply(venn_list, `[[`, 'id'), 
    # category.names=c('10°C', '20°C', '30°C'),
    filename=file.path(out_dir, paste0(stock, '_venn_identified.png')),
    main='Total quantified protein groups',
    width=1000, height=1000, resolution=150)
}

# Volcano plot contrasts --------------------------------------------------

volcano_contrast_file <- 'volcano_contrasts.txt'
if(!file.exists(volcano_contrast_file)) {
  if(is.null(compare_to_experiment)) {
    compare_df <- permuteNames(experiments)
  } else {
    compare_df <- permuteNames(experiments, compare_to_experiment)
  }
  if(reverse_comparison) {
    compare_df <- rev(compare_df)
  } else {
    compare_df <- compare_df
  } 
  names(compare_df) <- c('right','left')
  write.table_imb(compare_df, volcano_contrast_file)
} else {
  compare_df <- read.delim(volcano_contrast_file, 
                           colClasses=c('character', 'character'))
}

rMQanalysis::mylog('Calculate p-value and foldchange', dbg_level=0)
replicate_count <- length(replicates)
if(do_volcano & replicate_count > 2) {
  for(line in seq(nrow(compare_df))) {
    exp_a <- compare_df[line,1]
    exp_b <- compare_df[line,2] 
    regex_exp_a <- paste0('^imputed\\..+\\.', exp_a, '_.+')
    regex_exp_b <- paste0('^imputed\\..+\\.', exp_b, '_.+')
    col_a <- grep(regex_exp_a, names(pg_quant))
    col_b <- grep(regex_exp_b, names(pg_quant))
    pvalue_column <- paste('pvalue', exp_a, exp_b, sep='_')
    
    pg_quant[pvalue_column] <-
      switch(
        pvalue_correction,
        permutation=apply(pg_quant[c(col_a,col_b)],
                          1,
                          function(x) {
                            DV <- unlist(c(x[seq(length(col_a))], x[seq(length(col_b)) + length(col_a)]))
                            IV <- factor(rep(c("A", "B"), c(length(col_a), length(col_b))))
                            pvalue(oneway_test(DV ~ IV))
                          }),
        no=apply(pg_quant[c(col_a,col_b)],
                 1,
                 function(x) {
                   t.test(x[seq(length(col_a))],
                          x[seq(length(col_b)) + length(col_a)])$p.value
                 }),
        IHW={
          jnk <- apply(pg_quant[c(col_a,col_b)],
                       1,
                       function(x) {
                         t.test(x[seq(length(col_a))],
                                x[seq(length(col_b)) + length(col_a)])$p.value
                       })
          adj_pvalues(ihw(jnk ~ log2(pg_quant$Intensity), 
                          alpha=0.05, nbins=5))
        }
      ) 
    
    # pg_quant[paste0('corrected_', pvalue_column)] <- 
    #   p.adjust(pg_quant[[pvalue_column]], 'fdr')
    if(ihw_test) {
      p1 <- ggplot(pg_quant, aes_string(pvalue_column)) + 
        geom_histogram(binwidth=.01) +
        ggtitle('uncorrected')
      myp <- pg_quant[[pvalue_column]]
      p2 <- ggplot(pg_quant, aes(p.adjust(myp))) + geom_histogram(binwidth=.01) +
        coord_cartesian(ylim=c(0,100)) +
        ggtitle('BH corrected, y cut at 100')
      ihwRes <- ihw(pg_quant[[pvalue_column]] ~ log2(pg_quant$Intensity), 
                    alpha=0.1, nbins=5)
      ihwRes1 <- ihw(pg_quant[[pvalue_column]] ~ log2(pg_quant$Intensity), 
                     alpha=0.05, nbins=5)
      p3 <- ggplot(data.frame(x=adj_pvalues(ihwRes)), aes(x)) + 
        geom_histogram(binwidth=.01) +
        ggtitle('IHW corrected, alpha .1')
      p4 <- ggplot(data.frame(x=adj_pvalues(ihwRes1)), aes(x)) + 
        geom_histogram(binwidth=.01) +
        ggtitle('IHW corrected, alpha .05')
      
      p <- gridExtra::grid.arrange(p1, p2, p3, p4)
      rm(p1, p2, p3, p4)
      ggsave(file.path(out_dir, sprintf('hist_%s.pdf', pvalue_column)), 
             p, 
             width=8, height=8)    
    }
  }
  
  
  
  
  # we have to use a little work around to get the log.pvalue names in case of only one column
  log_pvalues = -log10(pg_quant[grep('^pvalue', names(pg_quant), value=TRUE)])
  colnames(log_pvalues) = paste0('log10.',names(log_pvalues))
  pg_quant <- cbind(pg_quant, log_pvalues)
  max_y_value <- max(c(unlist(pg_quant[grep('^log10.pvalue', names(pg_quant))]),
                       10^enrich_pvalue + 0.5))
} else {
  if(replicate_count < 3) {
    warning(sprintf('There are to less replicates to create Volcano plots!'))
  }
}


if((do_volcano  & length(replicates) > 2) | do_scatter) {
  for(line in seq(nrow(compare_df))) {
    exp_a <- compare_df[line,1]
    exp_b <- compare_df[line,2] 
    col_a <- paste0(difference_volcano_column, exp_a)
    col_b <- paste0(difference_volcano_column, exp_b)
    difference_column <- paste('difference', exp_a, exp_b, sep='_')
    pg_quant[difference_column] <- 
      pg_quant[col_a] - pg_quant[col_b]
  }
  
  max_x_value <- max(abs(pg_quant[grep('^difference', names(pg_quant))]), 
                     na.rm=TRUE)
  volcano_column_type <- sub('([^[:digit:]]+).*', '\\1', difference_volcano_column)
}


rMQanalysis::mylog('Create volcano plots:', dbg_level=0)

enriched_wb <- xlsx::createWorkbook()

techname <- c("IP", "RNase", "WT")

nicename <- c("RNase (-)", "RNase (+)", "Wild type")

expdf <- data.frame(techname, nicename)

color_intercept <- c("#D95F02","#66A61E")

if(do_volcano & length(replicates) > 2) {
  plotly_plots <- list()
  for(line in seq(nrow(compare_df))) {
    exp_a <- compare_df[line,1]
    exp_b <- compare_df[line,2]
    out_filename <- sprintf('%s_vs_%s',exp_a,exp_b)
    
    rMQanalysis::mylog(sprintf('(%s/%s) %s vs. %s', 
                               line, nrow(compare_df), exp_a, exp_b), dbg_level=0)
    
    if(filter_separately) {
      # This is the part which reduces the background proteins.
      val_count_a <- paste0('value_count_',exp_a) # colname for value count
      val_count_b <- paste0('value_count_',exp_b)
      pg_quant_orig <- pg_quant # saving original pg_quant, because we subset now
      filter_min_quant_events <- 
        pg_quant[val_count_a] >= min_quant_events |
        pg_quant[val_count_b] >= min_quant_events
      pg_quant <- # subset where value count in ANY of the columns to plot is above threshold
        pg_quant[filter_min_quant_events,]
      # we have to save the original pg_quant after EACH volcano plot
    }
    difference <- paste('difference', exp_a, exp_b, sep='_')
    pvalue <- paste('log10.pvalue', exp_a, exp_b, sep='_')
    enriched_col <- paste('enriched', exp_a, exp_b, sep='_')
    my_enriched <- Enriched2(pg_quant[c(difference, pvalue)],
                            c=enrich_c,
                            s0=enrich_s0,
                            pvalue=enrich_pvalue,
                            linetype=volcano_threshold_linetype,
                            size=volcano_threshold_linesize)
    pg_quant$enriched <- as.factor(my_enriched$enriched)
    
    pg_quant_orig[[enriched_col]] <- NA
    pg_quant_orig[[enriched_col]][filter_min_quant_events] <- my_enriched$enriched
    
    subtitle <- sprintf('enriched: %d; background: %d; total: %d',
                        sum(my_enriched$enriched), sum(my_enriched$background),
                        nrow(pg_quant))
    
    my_points_nohighlight <- interaction(pg_quant[!pg_quant$highlight, c('enriched', 'highlight')])
    my_points_highlight <- interaction(pg_quant[pg_quant$highlight, c('enriched', 'highlight')])
    
    niceexp_a <- expdf[expdf$techname == exp_a,]$nicename
    niceexp_b <- expdf[expdf$techname == exp_b,]$nicename
    
    plotly_string <- paste0(sprintf('log2(foldchange): %.2f, -log10(pvalue): %.2f<br />Value counts: %s; %s<br />Protein.IDs: %s<br />my_label: %s',
                                    pg_quant[[difference]], pg_quant[[pvalue]],
                                    paste(pg_quant[[paste0('value_count_', exp_a)]], 'in', niceexp_a), paste(pg_quant[[paste0('value_count_', exp_b)]], 'in', niceexp_b),
                                    sub('(.{30}).*', '\\1...', pg_quant$Protein.IDs),
                                    pg_quant$my_label),
                            switch('Gene.names' %in% names(pg_quant), 
                                   sprintf('<br />Gene.names: %s', pg_quant$Gene.names), 
                                   NULL),
                            switch('Protein.names' %in% names(pg_quant), 
                                   sprintf('<br />Protein.names: %s', pg_quant$Protein.names), 
                                   NULL),
                            switch('description' %in% names(pg_quant), 
                                   sprintf('<br />description: %s', pg_quant$description), 
                                   NULL))
    pg_quant$my_text <- plotly_string
    
    pg_quant_subset <- 
      pg_quant %>%
      rowwise() %>%
      mutate(imputation=ifelse(max(range(.data[[val_count_a]], .data[[val_count_b]])) == replicate_count &
                                 min(range(.data[[val_count_a]], .data[[val_count_b]])) == replicate_count,
                               'measured',
                               ifelse(min(range(.data[[val_count_a]], .data[[val_count_b]])) > 0,
                                      'some imputed',
                                      'ref imputed'))) %>%
      select(Protein.IDs, my_label, matches(difference), matches(pvalue), 
             my_text, highlight, enriched, imputation)
    volcanoplot <- 
      ggplot(pg_quant_subset, aes_string(x=difference, y=pvalue, text='my_text')) +   
      geom_hline(yintercept=0, color= color_intercept[line], na.rm=TRUE) +
      geom_vline(xintercept=0, color= color_intercept[line], na.rm=TRUE) +
      plot(my_enriched) + # adding the threshold line of our Enriched object
      geom_point(data=subset(pg_quant_subset, !highlight),
                 aes(alpha=my_points_nohighlight, 
                     color=my_points_nohighlight,
                     fill=my_points_nohighlight, 
                     # shape=imputation,
                     size=my_points_nohighlight), 
                 na.rm=TRUE) + 
      geom_point(data=subset(pg_quant_subset, highlight),
                 aes(alpha=my_points_highlight, 
                     color=my_points_highlight,
                     fill=my_points_highlight, 
                     # shape=imputation,
                     size=my_points_highlight), 
                 na.rm=TRUE) + 
      # we have combined ENRICHED.HIGHLIGHT
      scale_alpha_manual(values=c(`TRUE.FALSE`=.7, `FALSE.FALSE`=.4, `TRUE.TRUE`=.7, `FALSE.TRUE`=.7),
                         guide=FALSE) + 
      scale_color_manual(values=c(`TRUE.FALSE`='#000000', `FALSE.FALSE`='#b4b4b4', `TRUE.TRUE`='#e41a1c', `FALSE.TRUE`='#ff7f00'),
                         guide=FALSE) +  
      scale_fill_manual(values=c(`TRUE.FALSE`='#000000', `FALSE.FALSE`='#b4b4b4', `TRUE.TRUE`='#e41a1c', `FALSE.TRUE`='#ff7f00'),
                        guide=FALSE)+
      # scale_shape_manual(values=c(`measured`=21, `ref imputed`=25, `some imputed`=24))+
      scale_size_manual(values=c(`TRUE.FALSE`=1, `FALSE.FALSE`=.5, `TRUE.TRUE`=1.2, `FALSE.TRUE`=1),
                        guide=FALSE)
    
    color_df <- subset(pg_quant_subset, my_label %in% extra_highlight) %>% 
      as_tibble() %>% 
      mutate(side=ifelse(UQ(as.name(difference)) > 0, 'right','left')) %>% 
      group_by(side) %>% 
      arrange(desc(abs(UQ(as.name(difference)))))
    
    color <- interaction(color_df[color_df$highlight, c('enriched', 'highlight')])
    
    volcanoplot_full <- volcanoplot +
      ggtitle(sprintf('%s IP experiment\nVolcano plot of %s against %s\n%s',extra_highlight[2], niceexp_a, niceexp_b, subtitle)) +
      coord_cartesian(xlim=c(-max_x_value,max_x_value), ylim=c(-.8,max_y_value)) +
      ylab(expression(-log[10]*'(p-value)')) + 
      xlab(expression(log[2]*'(FoldChange)')) +
      annotate('text', label=niceexp_a, x=max_x_value / 2, y=-.4,size = 8) + 
      annotate('text', label=niceexp_b, x=-max_x_value / 2, y=-.4,size = 8)+
      theme_minimal()+
      # theme(panel.grid.major = element_line(linetype = 'solid', colour = colour), 
      #       panel.grid.minor = element_line(linetype = 'solid',colour = colour))+
      theme(legend.text = element_text(size = 14))+
      theme(legend.key.size = unit(1, 'cm'))+
      theme(legend.position="top")+
      theme(axis.text=element_text(size=18),axis.title=element_text(size=20,face="bold"))
    # if(sum(as.logical(pg_quant_subset$enriched)) > ggrepel_threshold) {
    #   volcanoplot_full <- volcanoplot_full + 
    #     geom_text(aes(label=my_label), size=font_size, vjust=1.2, 
    #               data=subset(pg_quant_subset, enriched == TRUE | highlight), na.rm=TRUE)
    # } else {
    volcanoplot_full <- volcanoplot_full +
      # ggrepel::geom_label_repel(data = subset(pg_quant_subset, enriched == TRUE | highlight),
      #                 aes(label=my_label),
      #                 size=1.7,
      #                 label.size = NA,, 
      # fontface = 'bold', color = 'black',
      # box.padding = 0.80, point.padding = 0.5, 
      #                 alpha = 0.6, 
      #                 label.padding=.1)
      ggrepel::geom_text_repel(data=subset(pg_quant_subset, enriched == TRUE | highlight) %>% 
                                 as_tibble() %>% 
                                 filter(!my_label %in% extra_highlight)%>%
                                 mutate(side=ifelse(UQ(as.name(difference)) > 0, 'right','left')) %>% 
                                 group_by(side) %>% 
                                 arrange(desc(UQ(as.name(pvalue))),desc(abs(UQ(as.name(difference))))) %>% 
                                 dplyr::slice(1:ggrepel_threshold)%>%
                                 filter(enriched == TRUE),
                               aes(label=my_label), 
                               size=font_size,
                               na.rm=TRUE,
                               segment.size=.3,
                               segment.alpha=.5)+
      ggrepel::geom_text_repel(data=subset(pg_quant_subset, my_label %in% extra_highlight) %>% 
                                 as_tibble() %>% 
                                 mutate(side=ifelse(UQ(as.name(difference)) > 0, 'right','left')) %>% 
                                 group_by(side) %>% 
                                 arrange(desc(abs(UQ(as.name(difference))))),
                               aes(label=my_label, color = color), 
                               size=font_size,
                               na.rm=TRUE,
                               segment.size=.3,
                               segment.alpha=.5,
                               box.padding = 1.5)
    # }
    
    print(volcanoplot_full)
    
    volcanoplot2 <- volcanoplot +
      coord_cartesian(xlim=range(labeling::extended(-max_x_value, max_x_value, m=3)), 
                      ylim=c(-.8,max_y_value)) +
      # coord_cartesian(xlim=c(-max_x_value,max_x_value), ylim=c(-.8,max_y_value)) +
      scale_x_continuous(expression(log[2]*'(foldchange)'), 
                         breaks=range(labeling::extended(-max_x_value, max_x_value, m=3)),
                         minor_breaks=labeling::extended(-max_x_value, max_x_value, m=9)) +
      scale_y_continuous(expression(-log[10]*'(p-value)'), 
                         breaks=labeling::extended(0, max_y_value, m=2, Q=c(1,5,2,4,3,6)),
                         minor_breaks=labeling::extended(0, max_y_value, m=8)) +
      theme_imb_small(9) +
      guides(shape=FALSE)
    # print(volcanoplot2)
    
    volcanoplot3 <- volcanoplot2 +
      geom_label_repel(data = subset(pg_quant_subset, highlight),
                       aes(label=my_label),
                       size=1.7,
                       label.size = NA, 
                       alpha = 0.6, 
                       label.padding=.1, 
                       fontface = 'bold', color = 'black',
                       box.padding = 0.80, point.padding = 0.35)
    # print(volcanoplot3)
    
    plotly_figure <- 
      plotly::ggplotly(volcanoplot, tooltip='text') %>%
      plotly::layout(xaxis = list(title = sprintf('difference: %s %s minus %s %s', 
                                                  volcano_column_type, exp_a, volcano_column_type, exp_b), 
                                  range = c(-max_x_value-.5,max_x_value+.5)), 
                     yaxis = list(title = "-log10(p.value)", range = c(-.8,max_y_value+.2)))
    plotly_plots[[line]] <- plotly_figure
    
    htmlwidgets::saveWidget(
      plotly_figure, 
      file=file.path(getwd(), out_dir,
                     sprintf(paste0(stock, '_interactive_volcano_%s.html'),
                             out_filename)))
    
    
    
    if(create_volcano_enriched_files) {
      pg_quant_enriched <- 
        pg_quant[as.logical(pg_quant$enriched), -length(pg_quant)]
      
      pg_quant_enriched <- 
        pg_quant_enriched[
          order(-pg_quant_enriched[[difference]] > 0, -pg_quant_enriched[[pvalue]]),]
      
      write.table_imb(pg_quant_enriched,
                      file.path(out_dir,
                                sprintf(paste0(stock, '_volcano_enriched_%s.txt'),
                                        out_filename)))
      important_enriched_columns <- 
        unique(c('Protein.IDs','my_label',"Gene.names", difference, pvalue, 
                 grep(paste0(c('value_count_', '^log2.*', '^imputed.*'), exp_a, collapse='|'), 
                      names(pg_quant_enriched), value=TRUE),
                 grep(paste0(c('value_count_', '^log2.*', '^imputed.*'), exp_b, collapse='|'), 
                      names(pg_quant_enriched), value=TRUE)))
      
      write.table_imb(pg_quant_enriched[important_enriched_columns],
                      file.path(out_dir,
                                sprintf(paste0(stock, '_volcano_enriched_reduced_%s.txt'),
                                        out_filename)))
      
      write.csv(pg_quant_enriched[important_enriched_columns],file = file.path(out_dir,sprintf(paste0(stock, '_Enriched_%s.csv'), out_filename)),row.names = F)
      
      reduced_proteinGroups_sheet <- 
        rMQanalysis::createSheetOfDF(
          enriched_wb, 
          pg_quant_enriched[important_enriched_columns], 
          out_filename,
          link=database_link)
    }
    if(create_volcano_pdfs) {
      ggsave(file.path(out_dir,sprintf(paste0(stock, '_volcano_full_%s.pdf'),out_filename)),
             volcanoplot_full,
             height=volcano_height,
             width=volcano_width)
      ggsave(file.path(out_dir,sprintf(paste0(stock, '_volcano_print_%s.pdf'),out_filename)),
             volcanoplot2,
             height=2,
             width=2)
      ggsave(file.path(out_dir,sprintf(paste0(stock, '_volcano_print_labeled_%s.pdf'),out_filename)),
             volcanoplot3,
             height=2,
             width=2)
    }
    if(filter_separately) {
      pg_quant <- pg_quant_orig 
    }
  }
  # finally, just removing the enriched column
  if(!filter_separately) {
    pg_quant['enriched'] <- NULL
  }
}
rMQanalysis::saveWorkbookMQ(enriched_wb, 
                            file.path(out_dir, paste0(stock, '_enriched_proteinGroups.xlsx')))

# Writing protein groups with enrichment values ---------------------------

rMQanalysis::mylog('Writing complete protein groups file', dbg_level=0)

write.table_imb(pg_quant,
                file.path(out_dir, pg_volcano_filename))


important_columns <- 
  c("my_label", "Protein.IDs",
    grep('^Protein.names', names(pg_quant), value=TRUE),
    grep('^Gene.names', names(pg_quant), value=TRUE),
    grep('^description', names(pg_quant), value=TRUE),
    "Peptides", "Razor...unique.peptides", "Unique.peptides", 
    "Sequence.coverage....", "Mol..weight..kDa.", "Sequence.length", 
    grep('^log2.LFQ.intensity.', names(pg_quant), value=TRUE),
    grep('^imputed.log2.LFQ.intensity.', names(pg_quant), value=TRUE),
    grep(difference_volcano_column, names(pg_quant), value=TRUE),
    grep('^value_count', names(pg_quant), value=TRUE),
    grep('^pvalue_', names(pg_quant), value=TRUE),
    grep('^log10.pvalue_', names(pg_quant), value=TRUE),
    grep('^difference_', names(pg_quant), value=TRUE),
    grep('^enriched_', names(pg_quant), value=TRUE),
    grep('^biomart.', names(pg_quant), value=TRUE))

rMQanalysis::mylog('Writing reduced protein groups file', dbg_level=0)
write.table_imb(pg_quant[important_columns],
                file.path(out_dir, paste0(stock, '_reduced_proteinGroups.txt')))

write.csv(x = pg_quant[important_columns],file = file.path(out_dir,paste0(stock, "_QuantifiedProteins.csv")),row.names = F)

tryCatch({
  excel_wb <- xlsx::createWorkbook()
  reduced_proteinGroups_sheet <- 
    rMQanalysis::createSheetOfDF(excel_wb, 
                                 pg_quant[important_columns], 
                                 'reduced proteinGroups',
                                 link=database_link)
  # rMQanalysis::highlightFilteredRows(filtered_proteins_sheet, pg$getData())
  rMQanalysis::saveWorkbookMQ(excel_wb, 
                              file.path(out_dir, paste0(stock, '_reduced_proteinGroups.xlsx')))
},
error=function(e) {
  warning('Writing of the "reduced_proteinGroups.xlsx" Excel file failed ',
          ifelse(as.numeric(format(Sys.time(), '%s')) %% 2,
                 sprintf('because it is NOT %s! ;-)', format(Sys.time() + 86400, '%A')),
                 sprintf('because it is %s! ;-)', format(Sys.time(), '%A'))
          ),
          '\nPlease format the txt within excel.',
          call. = FALSE)
}
)



# Create scatter plots ----------------------------------------------------

rMQanalysis::mylog('Creating scatter plots', dbg_level=0)
if(do_scatter) {
  if(is.null(compare_to_experiment)) {
    compare_df2 <- permuteNames(grep(difference_scatter_column, names(pg_quant), value=TRUE))
  } else {
    compare_df2 <- permuteNames(grep(difference_scatter_column, names(pg_quant), value=TRUE), 
                                compare_to_experiment)
  }
  
  checkOutsideLimits <- function(coefficients, x, y, 
                                 upper_distance=1, 
                                 lower_distance=1) {
    upper_half <- 
      coefficients[[2]] * x + coefficients[[1]] + upper_distance - y < 0 
    lower_half <- 
      coefficients[[2]] * x + coefficients[[1]] - lower_distance - y > 0 
    return(list(upper=upper_half, 
                lower=lower_half, 
                all=(upper_half | lower_half)))
  }
  difference_column_type <- sub('([^[:digit:]]+).*', '\\1', difference_scatter_column)
  for(line in seq(nrow(compare_df2))) {
    col_a <- as.character(compare_df2[line,1])
    exp_a <- sub(difference_scatter_column, '', col_a)
    col_b <- as.character(compare_df2[line,2]) 
    exp_b <- sub(difference_scatter_column, '', col_b)
    reg_coeff <- lm(pg_quant[[col_b]] ~ pg_quant[[col_a]])$coefficients
    outside <- checkOutsideLimits(reg_coeff, pg_quant[[col_a]], pg_quant[[col_b]],
                                  scatter_upper_threshold, scatter_lower_threshold)
    highlight <- nrow(na.omit(subset(pg_quant, highlight == TRUE, c(col_a, col_b)))) > 0
    my_env <- environment()
    graph <- 
      ggplot(pg_quant, aes_string(x=col_a, y=col_b), environment=my_env) +
      geom_smooth(method='lm', fullrange=TRUE, alpha=.2, na.rm=TRUE) + 
      geom_abline(intercept=reg_coeff[[1]], slope=reg_coeff[[2]], color='blue', na.rm=TRUE) +
      geom_point(alpha=scatter_density, size=scatter_size, na.rm=TRUE) +
      geom_abline(intercept=reg_coeff[[1]] + scatter_upper_threshold, na.rm=TRUE, 
                  slope=reg_coeff[[2]], color='#9999CC', lwd=.2) + 
      geom_abline(intercept=reg_coeff[[1]] - scatter_lower_threshold, na.rm=TRUE, 
                  slope=reg_coeff[[2]], color='#9999CC', lwd=.2) + 
      geom_point(data=pg_quant[outside$all,], alpha=1, color='red', na.rm=TRUE) + 
      geom_text(data=pg_quant[outside$upper,], aes(label=my_label), size=3, hjust=1.1, na.rm=TRUE) +  
      geom_text(data=pg_quant[outside$lower,], aes(label=my_label), size=3, hjust=-.1, na.rm=TRUE) +   
      theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank()) +
      theme_imb() + xlab(exp_a) + ylab(exp_b) +
      ggtitle(sprintf('%s %s vs. %s: %s upper and %s lower protein groups',
                      difference_column_type, exp_a, exp_b, sum(outside$upper, na.rm=TRUE), 
                      sum(outside$lower, na.rm=TRUE)))
    if(highlight) {
      graph <- 
        graph + 
        geom_point(alpha=1, fill='blue', color='white', pch=21, 
                   data=subset(pg_quant, highlight == TRUE), na.rm=TRUE) +
        geom_text(aes(label=my_label), size=font_size, vjust=1.3,
                  data=subset(pg_quant, highlight == TRUE), na.rm=TRUE)
    } else {
      cat('Could not label your protein IDs to highlight because of missing values. 
          Change "difference_scatter_column" to "median3_less1_imp_", in this case the NAs get imputed.\n')
    }
    suppressWarnings(print(graph))
    
    out_filename <- sprintf('%s_vs_%s', exp_a, exp_b)
    if(create_scatter_enriched_files) {
      write.table_imb(pg_quant[which(outside$upper),],
                      file.path(out_dir,sprintf('scatter_upper_enriched_%s.txt',out_filename)))
      write.table_imb(pg_quant[which(outside$lower),],
                      file.path(out_dir,sprintf('scatter_lower_enriched_%s.txt',out_filename)))
    }
    if(create_scatter_pdfs) {
      ggsave(file.path(out_dir,sprintf(paste0(stock, '_scatter_%s.pdf'),out_filename)),
             graph,
             height=scatter_height,
             width=scatter_width)
    }
  }
}


# Interactive heatmaps ----------------------------------------------------

# rMQanalysis::mylog('Creating interactive heatmap for all enriched proteins', 
#                    dbg_level=0)
# 
# pg_heatmap <- 
#   pg_quant %>%
#   as_tibble() %>%
#   filter_at(vars(starts_with('enriched_')), any_vars(.)) %>%
#   # glimpse() %>%
#   select(my_label, Protein.IDs, id, 
#          starts_with(difference_volcano_column), 
#          starts_with('value_count_')) %>%
#   mutate(short_label=sub('([^;]+).*', '\\1', my_label),
#          rowname=paste(short_label, id))
# 
# heatmaply(pg_heatmap %>% select(starts_with(difference_volcano_column)), 
#           file=file.path(out_dir, 'interactive_heatmap.html'),
#           main=sprintf('%s %.1f-fold enriched protein groups', 
#                        nrow(pg_heatmap),
#                        2^enrich_s0),
#           Colv=FALSE, 
#           labRow=pg_heatmap$rowname)
# 
# pg_heatmap2 <- 
#   pg_quant %>%
#   as_tibble() %>%
#   select(my_label, Protein.IDs, id, starts_with('value_count_')) %>%
#   mutate(short_label=sub('([^;]+).*', '\\1', my_label),
#          rowname=paste(short_label, id))
# 
# 
# heatmaply(pg_heatmap2 %>% select(starts_with('value_count_')), 
#           file=file.path(out_dir, 'interactive_heatmap_value_count.html'),
#           main=sprintf('%s proteins cluster by value count', 
#                        nrow(pg_heatmap2),
#                        2^enrich_s0),
#           Colv=FALSE, 
#           labRow=pg_heatmap2$rowname)

save.image(file = image_name)

# Session info ------------------------------------------------------------

rMQanalysis::mylog('Write session info and environment', dbg_level=0)
my_env <- c("all_pca_filename", "beta_max", "beta_min", "col_a", "col_b", 
            "compare_df", "compare_to_experiment", "create_scatter_enriched_files", 
            "create_scatter_pdfs", "create_volcano_enriched_files", "create_volcano_pdfs", 
            "data_dir", "difference", "difference_column", "difference_scatter_column", 
            "difference_volcano_column", "do_scatter", "do_volcano", "enrich_c", 
            "enrich_pvalue", "enrich_s0", "exp_a", "exp_b", "experiments", 
            "fallback_column", "fallback_column_regex", "filter_separately", 
            "font_size", "ggrepel_threshold", "graph_abundance_rank", "graph_imputed_values", 
            "graph_imputed_values_bar", "graph_value_count", "heatmap_column_regex", 
            "heatmap_white_area_size", "highlight_file", "ids_to_highlight", 
            "ihwRes", "ihwRes1", "impute_method", "imputed_bar_df", "imputed_values", 
            "intensity_comp_graph", "jnk", "label_column", "label_column_regex", 
            "lfq_columns", "line", "log_pvalues", "log2_lfq_columns", "max_x_value", 
            "max_y_value", "min_quant_events", "multi_pca_plot", "my_enriched", 
            "my_points_highlight", "my_points_nohighlight", "myp", "norm_modify_mean", 
            "norm_sd_factor", "original_average_function", "original_na_replace", 
            "out_dir", "out_filename", "p", "plotly_figure", "pca_column_regex", 
            "pca_columns", "pg", "pg_basename", "pg_data", "pg_filename", 
            "pg_flt", "pg_ident", "pg_quant", "pg_quant_filename", "pg_quant_orig", 
            "pg_volcano_filename", "prcomp_list", "pvalue", "pvalue_column", 
            "random_seed", "razor_peptides", "regex_exp_a", "regex_exp_b", 
            "replicate_count", "replicates", "reverse_comparison", "scatter_density", 
            "scatter_height", "scatter_lower_threshold", "scatter_size", 
            "scatter_upper_threshold", "scatter_width", "subtitle", "unique_peptides", 
            "val_count_a", "val_count_b", "volcano_column_type", "volcano_height", 
            "volcano_threshold_linesize", "volcano_threshold_linetype", "volcano_width"
)
# save(list=ls(), 
#      file=file.path(out_dir,'main_analysis_data.RData'))
writeLines(capture.output(sessionInfo()), 
           file.path(out_dir, paste0(stock, '_sessionInfo.txt')))

rMQanalysis::mylog(sprintf('Results are written in %s', out_dir), dbg_level=0)

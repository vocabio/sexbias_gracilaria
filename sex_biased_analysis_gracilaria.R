if(!require(ggplot2)){
  install.packages("ggplot2", lib = "E:/r_library")
  library('ggplot2')
}

if (!requireNamespace("BiocManager", quietly = TRUE)){
  install.packages("BiocManager", lib = "E:/r_library")
}

if(!require(tximport)){
  BiocManager::install("tximport", lib = "E:/r_library")
  library('tximport')
}

if(!require(DESeq2)){
  BiocManager::install("DESeq2", lib = "E:/r_library")
  library('DESeq2')
}

if(!require(topGO)){
  BiocManager::install("topGo", lib = "E:/r_library")
  library('topGo')
}

library(data.table)
library(gridExtra)
library(grid)

###DE Analysis###
create_dds <- function(row_names, condition, kallisto_results, num_samples, species_name){
  #Create DDS object from Deseq2
  coldata <- data.frame(row.names=row_names, condition)
  dds <- DESeq2::DESeqDataSetFromTximport(kallisto_results, colData=coldata, design=~condition)
  dds <- DESeq2::DESeq(dds)
  gene_universe <- rownames(dds) # a list containing all gene names needed for functional analysis
  
  #remove genes with 0 expression:
  keep <- rowSums(counts(dds)) > 1
  dds <- dds[keep,]
  res <- DESeq2::results(dds)
  
  #contrast argument switches from male vs female to female vs male
  res_female <- DESeq2::results(dds,contrast = c("condition","Female","Male")) 
  # rlogTrans avoids the distance measure being dominated by a few highly variable genes
  rld <- DESeq2::rlogTransformation(dds) 
  
  #add counts to dataframe
  resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE) 
  names(resdata)[1] <- "Gene"
  resdata_female <- merge(as.data.frame(res_female), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
  names(resdata_female)[1] <- "Gene"
  
  # Prepare data from kallisto results for tpm analysis
  abundance_df <- as.data.frame(kallisto_results$abundance, row.names = NULL)
  abundance_df <- setDT(abundance_df, keep.rownames = "Gene")[]
  abundance_df <- as.data.frame(na.omit(abundance_df))
  
  # Depending on how many samples we used, we grap different columns from the dataframe and add TPM values
  # This was done by hand here but should definitely be automated for future use
  if (num_samples == 4){
    abundance_df$MeanTPMMales <- rowMeans(subset(abundance_df, select = c(V1, V2)))
    abundance_df$MeanTPMFemales <- rowMeans(subset(abundance_df, select = c(V3, V4)))
    abundance_df$MeanTPMAll <- rowMeans(subset(abundance_df, select = c(V1, V2, V3, V4)))
  } else if (num_samples == 6){
    abundance_df$MeanTPMMales <- rowMeans(subset(abundance_df, select = c(V1, V2, V3)))
    abundance_df$MeanTPMFemales <- rowMeans(subset(abundance_df, select = c(V4, V5, V6)))
    abundance_df$MeanTPMAll <- rowMeans(subset(abundance_df, select = c(V1, V2, V3, V4, V5, V6)))
  } else if (num_samples == 8){
    abundance_df$MeanTPMMales <- rowMeans(subset(abundance_df, select = c(V1, V2, V3, V4)))
    abundance_df$MeanTPMFemales <- rowMeans(subset(abundance_df, select = c(V5, V6, V7, V8)))
    abundance_df$MeanTPMAll <- rowMeans(subset(abundance_df, select = c(V1, V2, V3, V4, V5, V6, V7, V8)))
  }
  
  # Add FoldChange column to resdata dataframe
  resdata <- merge(resdata, abundance_df[ , c("Gene", "MeanTPMMales", "MeanTPMFemales", "MeanTPMAll")], by = "Gene")
  resdata$FoldChange <- 2 ^ resdata$log2FoldChange
  resdata_female <- merge(resdata_female, abundance_df[ , c("Gene", "MeanTPMMales", "MeanTPMFemales")], by = "Gene")
  resdata_female$FoldChange <- 2 ^ resdata_female$log2FoldChange
  
  # Extract genes that are exclusively expressed in one sex
  male_linked_resdata <- subset(resdata, resdata$MeanTPMMales > 0 & resdata$MeanTPMFemales == 0 & resdata$padj <= 0.05 & resdata$log2FoldChange >= 1)
  female_linked_resdata <- subset(resdata, resdata$MeanTPMMales == 0 & resdata$MeanTPMFemales > 0 & resdata$padj <= 0.05 & resdata$log2FoldChange <= -1)
  
  
  return(list("dds" = dds, "gene_universe" = gene_universe, "res_male" = res, "rld" = rld, "resdata" = resdata,
              "male_linked_resdata" = male_linked_resdata, "female_linked_resdata" = female_linked_resdata,
              "res_female"=res_female, "resdata_female"=resdata_female, "abundance_df" = abundance_df))
}
run_dds <- function(species_name, num_samples_per_sex, abundance_file_names, rownames_species){
  working_dir <- paste("../", species_name, "/", sep="")
  condition_species <- factor(rep(c("Male", "Female"), each=num_samples_per_sex))
  files_species <- file.path(working_dir, abundance_file_names)
  tx2gene_species <- read.delim(paste("../", species_name, "/fake_tx2gene", sep=""), header = FALSE)
  species.txi.kallisto.tsv <- tximport(files_species, type = "kallisto", tx2gene = tx2gene_species, abundanceCol = "tpm" )
  head(species.txi.kallisto.tsv$counts)
  dds_results <- create_dds(rownames_species, condition_species, species.txi.kallisto.tsv, num_samples_per_sex * 2, species_name)
  
  return(list("dds_results"=dds_results, 'conditions'=condition_species))
}
write_results <- function(dds, resdata, all_table_filename, signif_table_filename, lfcThreshold=1, padj_threshold=0.05){
  write.table(resdata, file=all_table_filename, sep="\t")
  significant_DE_genes <- subset(resdata, padj<padj_threshold & abs(log2FoldChange)>lfcThreshold)
  write.table(significant_DE_genes, file=signif_table_filename, sep="\t")
}
#################
###Functional Analysis###
run_functional_analysis <- function(res_species, gene_universe, species_name){
  # File must look like this:
  # tab delimited
  # First column is the gene name that's been used in the rest of the analysis here as well
  # Second Column is a list of go terms seperated by , without white spaces
  # E.g.: GeneA GO:1234567,GO:7654321...
  # Such a file can be created from the output tsv file created by eggnog where the first column is the gene name and 10th column is the list of go terms (tab delimited)
  geneID2GO_species <- readMappings(file.choose())
  gene_selection_male <- factor(rownames(subset.data.frame(res_species, res_species$log2FoldChange >= 1 & res_species$pvalue < 0.05)))
  topGO_results_male <- create_topGO_table(gene_universe, gene_selection_male, geneID2GO_species, paste("functional_results_", species_name, "_male.tab", sep=""))
  allRes_male <- topGO_results_male$allRes
  topGO_male <- topGO_results_male$topGO
  create_revigo_table(runTest(topGO_male, algorithm = "elim", statistic = "fisher"), paste("go_term_p_value_", species_name, "_male_elim", sep=""))
  
  gene_selection_female <- factor(rownames(subset.data.frame(res_species, res_species$log2FoldChange <= 1 & res_species$pvalue < 0.05)))
  topGO_results_female <- create_topGO_table(gene_universe, gene_selection_female, geneID2GO_species, paste("functional_results_", species_name, "_female.tab", sep=""))
  allRes_female <- topGO_results_female$allRes
  topGO_female <- topGO_results_female$topGO
  create_revigo_table(runTest(topGO_female, algorithm = "elim", statistic = "fisher"),paste("go_term_p_value_", species_name, "_female_elim", sep=""))
}
create_topGO_table <- function(gene_universe, gene_selection, geneID2GO, file_name){
  gene_list <- factor(as.integer(gene_universe%in%gene_selection))
  names(gene_list) <- gene_universe
  topGO_data<- new("topGOdata",
                   description = "Functional analysis Gracilaria Chilensis", 
                   ontology = "BP",
                   allGenes = gene_list, 
                   geneSel = gene_selection,
                   annot = annFUN.gene2GO,
                   gene2GO = geneID2GO)
  topGO_data
  resultFisher <- runTest(topGO_data, algorithm = "classic", statistic = "fisher")
  resultFisher
  
  resultKS <- runTest(topGO_data, algorithm = "classic", statistic = "ks")
  resultKS_elim <- runTest(topGO_data, algorithm = "elim", statistic = "ks")
  
  summary_fisher<- summary(attributes(resultFisher)$score <= 0.05)
  numsignif <- as.integer(summary_fisher[[3]])
  
  allRes_df <- GenTable(topGO_data,
                        classicFisher = resultFisher,
                        classicKS = resultKS,
                        elimKS = resultKS_elim,
                        orderBy = "elimKS",
                        ranksOf = "classicFisher",
                        topNodes = numsignif)
  allRes_df_signif <- subset.data.frame(allRes_df, allRes_df$Significant > 0)
  write.table(allRes_df_signif, file_name, append = FALSE, sep = "\t", dec = ".",
              row.names = FALSE, col.names = TRUE)
  geneData(resultFisher)
  return(list("topGO" = topGO_data, "allRes" = allRes_df)) #Access individual values by list$topGO
}
create_revigo_table <- function(fisher_elim, file_name){
  go_p_mapping <- tibble::enframe(score(fisher_elim))
  go_p_mapping <- subset.data.frame(go_p_mapping, go_p_mapping$value < 0.05)
  go_p_mapping <- go_p_mapping[order(go_p_mapping$value),]
  write.table(go_p_mapping, file_name, append = FALSE, sep = "\t", dec = ".", row.names = FALSE, col.names = FALSE)
}
#########################
###Plotting###
create_plots <- function(dds, dds_results, rld, conditions, ma_plot_title, volcano_plot_title, heatmap_title, violin_title, lfcThreshold=1, padj_threshold=0.05){
  create_ma_plot(dds_results, lfcThreshold, padj_threshold, ma_plot_title)
  ma_plot <- recordPlot()
  volcano_data <- create_volcano_plot(na.omit(dds_results), volcano_plot_title, lfcThreshold, padj_threshold)
  volcano_plot <- recordPlot()
  create_heatmap(conditions, rld, heatmap_title)
  heatmap_plot <- recordPlot()
  create_pca_plot(rld)
  pca_plot <- recordPlot()
  #pca_plot <- DESeq2::plotPCA(rld,intgroup="condition")
  plotDispEsts(dds)
  sample_dispersion_plot <- recordPlot()
  violin_data <- create_violin_plot(dds_results, 'hochberg', violin_title)
  violin_plot <- recordPlot()
  return(list("volcano"=volcano_plot, "ma"=ma_plot, "heatmap"=heatmap_plot, "pca"=pca_plot, "sample_dispersion"=sample_dispersion_plot, "violin"=violin_plot, "violin_data"=violin_data, "volcano_data"=volcano_data))
}
create_ma_plot <- function(dds_results, lfcThreshold=1, padj_threshold=0.05, plot_title){
  require(calibrate)
  with(dds_results, plot(baseMean, log2FoldChange, pch=20, cex=0.5, log="x", main=plot_title))
  with(subset(dds_results, padj < padj_threshold, log2FoldChange > lfcThreshold), points(baseMean, log2FoldChange, col="red", pch=20, cex = 1.5))
  with(subset(dds_results, padj < padj_threshold), textxy(baseMean, log2FoldChange, labs = Gene, cex = 0.5, col = 2))
  abline(1, 0, untf = FALSE, col="green")
  abline(-1, 0, untf = FALSE, col="green")
}
create_volcano_plot <- function(dds_results, plot_title, lfcThreshold=1, padj_threshold=0.05){
  fc <- abs(dds_results$log2FoldChange)
  padj <- dds_results$padj
  dds_results$color <- ifelse(fc<lfcThreshold, "black",
                              ifelse(fc<lfcThreshold*2 & padj < 0.05, "grey",
                                     ifelse(fc<lfcThreshold*4 & padj < 0.05, "orange",
                                            ifelse(padj < 0.05, "darkgreen", "black"))))
  p <- ggplot(data=dds_results, aes(x=log2FoldChange, y=-log10(padj))) +
    geom_point(colour=dds_results$color) + theme_minimal() +
    theme(text=element_text(size=20)) + xlab("log2(FC)") + labs(title = plot_title)
  print(p)
  return(ggplotGrob(p))
}
create_heatmap <- function(conditions, rld, plot_title){
  library(gplots)
  library(RColorBrewer)
  sampleDists <- as.matrix(dist(t(SummarizedExperiment::assay(rld))))
  cols <- brewer.pal(3, "Set1")[1:length(unique(conditions))] #3 is the minimum number of colors taken from the palett
  heatmap.2(as.matrix(sampleDists), key=T, trace="none",
            col=colorpanel(100,"black", "antiquewhite3"),
            ColSideColors=cols[conditions], RowSideColors=cols[conditions],
            margin=c(10, 10), main=plot_title, cexRow = 2, cexCol = 2,
            key.title = "Color Key and Hist", keysize = 2, dendrogram = "column")
  par(cex.main=1.5)
}
create_violin_plot <- function(resdata, wilcoxon_test_type, plot_title){
  mbg_resdata <- log2(subset(resdata, resdata$log2FoldChange>= 1 & resdata$padj <= 0.05)[, c('MeanTPMMales', 'MeanTPMFemales', 'MeanTPMAll')])
  fbg_resdata <- log2(subset(resdata, resdata$log2FoldChange<= -1 & resdata$padj <= 0.05)[, c('MeanTPMMales', 'MeanTPMFemales', 'MeanTPMAll')])
  ubg_resdata <- log2(subset(resdata, resdata$log2FoldChange< 1 &
                               resdata$log2FoldChange > -1 & resdata$padj <= 0.05)[, c('MeanTPMMales', 'MeanTPMFemales', 'MeanTPMAll')])
  input_df <- data.frame(log2tpm = mbg_resdata$MeanTPMMales, type = 'MBG_males', sex='Male', genes='MBG')
  input_df <- rbind(input_df, data.frame(log2tpm = mbg_resdata$MeanTPMFemales, type = 'MBG_females', sex='Female', genes='MBG'))
  input_df <- rbind(input_df, data.frame(log2tpm = fbg_resdata$MeanTPMMales, type = 'FBG_males', sex='Male', genes='FBG'))
  input_df <- rbind(input_df, data.frame(log2tpm = fbg_resdata$MeanTPMFemales, type = 'FBG_females', sex='Female', genes='FBG'))
  input_df <- rbind(input_df, data.frame(log2tpm = ubg_resdata$MeanTPMMales, type = 'UBG_males', sex='Male', genes='UBG'))
  input_df <- rbind(input_df, data.frame(log2tpm = ubg_resdata$MeanTPMFemales, type = 'UBG_females', sex='Female', genes='UBG'))
  
  #input_df <- rbind(input_df, data.frame(log2tpm = mbg_resdata$MeanTPMAll, type = 'MBG_all'))
  #input_df <- rbind(input_df, data.frame(log2tpm = fbg_resdata$MeanTPMAll, type = 'FBG_all'))
  #input_df <- rbind(input_df, data.frame(log2tpm = ubg_resdata$MeanTPMAll, type = 'UBG_all'))
  p_value <- pairwise.wilcox.test(input_df$log2tpm, input_df$type, wilcoxon_test_type)
  print(head(input_df))
  print(p_value)
  while (!is.null(dev.list()))  dev.off()
  
  p <- ggviolin(input_df, x="sex", y="log2tpm", fill="sex",
                palette = c("#377eb8","#e41a1c"), add = "boxplot", add.params = list(fill="white"),
                short.panel.labs = FALSE) +
    facet_grid(cols=vars(genes)) +
    ggtitle(substr(plot_title, 1, nchar(plot_title)-25)) + ylab('Log2(TPM)') + xlab('') + labs(fill = '') +
    stat_compare_means(label = "p.signif", method = "wilcox.test",
                       ref.group = ".all.", pairwise=TRUE) +
    theme(text=element_text(size=20))
  print(p)
  return(ggplotGrob(p))
}
create_pca_plot <- function(rld){
  pca_data <- DESeq2::plotPCA(rld,intgroup="condition", returnData = TRUE)
  p <- ggplot(data = pca_data, aes(x=PC1, y=PC2, color=group)) + geom_point(size=4) + scale_color_brewer(palette="Set1") +
    theme(text=element_text(size=20))
  print(p)
}
create_multiplot <- function(resdata, species_name, sex){
  # Very clunky method to create side by side plots of tpm and mean expression values in bins of fold changes
  resdata <- subset(resdata, resdata$padj < 0.05)
  # Create datasets for the desired bins
  numGenes <- c(length(rownames(subset(resdata, resdata$FoldChange >= 2 & resdata$FoldChange < 4))),
                length(rownames(subset(resdata, resdata$FoldChange >= 4 & resdata$FoldChange < 10))),
                length(rownames(subset(resdata, resdata$FoldChange >= 10 & resdata$FoldChange < 20))),
                length(rownames(subset(resdata, resdata$FoldChange >= 20))))
  maleTPMs <- c(log2(mean(subset(resdata, resdata$FoldChange >= 2 & resdata$FoldChange < 4)$MeanTPMMales)),
                log2(mean(subset(resdata, resdata$FoldChange >= 4 & resdata$FoldChange < 10)$MeanTPMMales)),
                log2(mean(subset(resdata, resdata$FoldChange >= 10 & resdata$FoldChange < 20)$MeanTPMMales)),
                log2(mean(subset(resdata, resdata$FoldChange >= 20)$MeanTPMMales)))
  maleSD <- c(log10(sd(subset(resdata, resdata$FoldChange >= 2 & resdata$FoldChange < 4)$MeanTPMMales) / sqrt(numGenes[1])),
              log10(sd(subset(resdata, resdata$FoldChange >= 4 & resdata$FoldChange < 10)$MeanTPMMales) / sqrt(numGenes[2])),
              log10(sd(subset(resdata, resdata$FoldChange >= 10 & resdata$FoldChange < 20)$MeanTPMMales) / sqrt(numGenes[3])),
              log10(sd(subset(resdata, resdata$FoldChange >= 20)$MeanTPMMales) / sqrt(numGenes[4])))
  femaleTPMs <- c(log2(mean(subset(resdata, resdata$FoldChange >= 2 & resdata$FoldChange < 4)$MeanTPMFemales)),
                  log2(mean(subset(resdata, resdata$FoldChange >= 4 & resdata$FoldChange < 10)$MeanTPMFemales)),
                  log2(mean(subset(resdata, resdata$FoldChange >= 10 & resdata$FoldChange < 20)$MeanTPMFemales)),
                  log2(mean(subset(resdata, resdata$FoldChange >= 20)$MeanTPMFemales)))
  femaleSD <- c(log10(sd(subset(resdata, resdata$FoldChange >= 2 & resdata$FoldChange < 4)$MeanTPMFemales) / sqrt(numGenes[1])),
              log10(sd(subset(resdata, resdata$FoldChange >= 4 & resdata$FoldChange < 10)$MeanTPMFemales) / sqrt(numGenes[2])),
              log10(sd(subset(resdata, resdata$FoldChange >= 10 & resdata$FoldChange < 20)$MeanTPMFemales) / sqrt(numGenes[3])),
              log10(sd(subset(resdata, resdata$FoldChange >= 20)$MeanTPMFemales) / sqrt(numGenes[4])))
  # Merge datasets
  df <- data.frame(FoldChange=c("2-4", "4-10", "10-20", ">20"),
                   NumGenes=numGenes,
                   maleTPMs = maleTPMs,
                   femaleTPMs = femaleTPMs)
  df$FoldChange <- factor(df$FoldChange,levels = c("2-4", "4-10", "10-20", ">20"))
  df$maleTPMs[is.nan(df$maleTPMs)] <- 0
  df$femaleTPMs[is.nan(df$femaleTPMs)] <- 0
  
  line_chart <- ggplot(data=df, aes(x=FoldChange)) +
    geom_line(aes(y=maleTPMs, group = 1), color="blue") +
    geom_line(aes(y=femaleTPMs, group = 2), color="red") + 
    labs(y = "Log2TPM") +
    theme_minimal() + 
    theme(axis.title.x = element_blank(), axis.text.x = element_blank(), text=element_text(size=17)) + 
    ggtitle(paste(species_name, " ", sex, sep=""))
  
  bar_chart <- ggplot(data=df, aes(x=FoldChange, y=NumGenes)) +
    geom_bar(stat="identity") +
    theme_minimal() +
    theme(text=element_text(size=17)) +
    labs(y = "N. Genes")

  grid.newpage()
  res_plot <- rbind(ggplotGrob(line_chart), ggplotGrob(bar_chart), size = "last")
    
  return(res_plot)
}
##############################



###Vermiculophylla Code###
abundance_file_names_vermiculophylla <- c("33mal_abundance.tsv", "36mal_abundance.tsv", "50mal_abundance.tsv",
                          "00fem_abundance.tsv","34fem_abundance.tsv", "39fem_abundance.tsv")
rownames_vermiculophylla <- c("Male33", "Male36", "Male50", "Female00", "Female34", "Female39")
abundance_file_names_vermiculophylla <- c("03mal_abundance.tsv", abundance_file_names_vermiculophylla, "40fem_abundance.tsv")
rownames_vermiculophylla <- c("Male03", rownames_vermiculophylla, "Female40")
dds_results <- run_dds("Vermiculophylla", 4, abundance_file_names_vermiculophylla, rownames_vermiculophylla)
dds_vermiculophylla <- dds_results$dds
gene_universe_vermiculophylla <- dds_vermiculophylla$gene_universe
res_vermiculophylla <- dds_vermiculophylla$res_male
rld_vermiculophylla <- dds_vermiculophylla$rld
resdata_vermiculophylla <- dds_vermiculophylla$resdata
resdata_female_vermiculophylla <- dds_vermiculophylla$resdata_female
female_linked_resdata_vermiculophylla <- dds_vermiculophylla$female_linked_resdata
male_linked_resdata_vermiculophylla <- dds_vermiculophylla$male_linked_resdata
verm_abundance <- dds_results$dds_results$abundance_df
verm_male_specific <- subset(verm_abundance, verm_abundance$MeanTPMFemales == 0 & verm_abundance$MeanTPMMales > 0)
verm_female_specific <- subset(verm_abundance, verm_abundance$MeanTPMFemales > 0 & verm_abundance$MeanTPMMales == 0)
verm_conditions <- dds_results$conditions
plots_vermiculophylla <- create_plots(dds_vermiculophylla$dds,resdata_vermiculophylla, rld_vermiculophylla, verm_conditions,
                               "G. vermiculophylla", "G. vermiculophylla", "G. Vermiculophylla (6)", 'G. vermiculophylla gene expression variance')
#run_functional_analysis(res_vermiculophylla, gene_universe_vermiculophylla, "Vermiculophylla")
res_plot_male_verm <- create_multiplot(subset(resdata_vermiculophylla, resdata_vermiculophylla$log2FoldChange > 0), "Vermiculophylla", "Male")
res_plot_female_verm <- create_multiplot(subset(resdata_female_vermiculophylla, resdata_female_vermiculophylla$log2FoldChange > 0), "Vermiculophylla", "Female")
write_results(dds_vermiculophylla$dds, resdata_vermiculophylla, "verm_F_vs_M.tab", "verm_DEsig0.05FC2_F_vs_M.tab")
##########################
###Chilensis Code###
abundance_file_names_chilensis <- c("134mal_abundance.tsv","133mal_abundance.tsv" ,"132fem_abundance.tsv" ,"131fem_abundance.tsv")
rownames_species_chilensis <- c("Male134", "Male133", "Female132", "Female131")
dds_results <- run_dds("Chilensis", 2, abundance_file_names_chilensis, rownames_species_chilensis)
dds_chilensis <- dds_results$dds
chilensis_conditions <- dds_results$conditions
gene_universe_chilensis <- dds_chilensis$gene_universe
res_chilensis <- dds_chilensis$res_male
rld_chilensis <- dds_chilensis$rld
resdata_chilensis <- dds_chilensis$resdata
resdata_female_chilensis<- dds_chilensis$resdata_female
female_linked_resdata_chilensis <- dds_chilensis$female_linked_resdata
male_linked_resdata_chilensis <- dds_chilensis$male_linked_resdata
chil_abundance <- dds_results$dds_results$abundance_df
chil_male_specific <- subset(chil_abundance, chil_abundance$MeanTPMFemales == 0 & chil_abundance$MeanTPMMales > 0)
chil_female_specific <- subset(chil_abundance, chil_abundance$MeanTPMFemales > 0 & chil_abundance$MeanTPMMales == 0)
plots_chilensis <- create_plots(dds_chilensis$dds,resdata_chilensis, rld_chilensis, chilensis_conditions,
                               " Chilensis", "G. chilensis", "G. chilensis", 'G. chilensis gene expression variance')
#run_functional_analysis(res_chilensis, gene_universe_chilensis, "Chilensis")
res_plot_male_chil <- create_multiplot(subset(resdata_chilensis, resdata_chilensis$log2FoldChange > 0), "Chilensis", "Male")
res_plot_female_chil <- create_multiplot(subset(resdata_female_chilensis, resdata_female_chilensis$log2FoldChange > 0), "Chilensis", "Female")
write_results(dds_chilensis$dds, resdata_chilensis, "chil_F_vs_M.tab", "chil_DEsig0.05FC2_F_vs_M.tab")
####################
###Caudata Code###
abundance_file_names_caudata <- c("179mal_abundance.tsv","178mal_abundance.tsv" ,"177mal_abundance.tsv" ,"176mal_abundance.tsv", 
                          "175fem_abundance.tsv","174fem_abundance.tsv","172fem_abundance.tsv","171fem_abundance.tsv")
rownames_species_caudata <- c("Male179", "Male178", "Male177", "Male176", "Female175", "Female174", "Female173", "Female172")
dds_results <- run_dds("Caudata", 4, abundance_file_names_caudata, rownames_species_caudata)
dds_caudata <- dds_results$dds
caudata_conditions <- dds_results$conditions
gene_universe_caudata <- dds_caudata$gene_universe
res_caudata <- dds_caudata$res_male
rld_caudata <- dds_caudata$rld
resdata_caudata <- dds_caudata$resdata
resdata_female_caudata <- dds_caudata$resdata_female
female_linked_resdata_caudata <- dds_caudata$female_linked_resdata
male_linked_resdata_caudata <- dds_caudata$male_linked_resdata
caud_abundance <- dds_results$dds_results$abundance_df
plots_caudata <- create_plots(dds_caudata$dds,resdata_caudata, rld_caudata, caudata_conditions,
                               "MA Plot Caudata", "G. caudata", "Heatmap Caudata",'G. caudata gene expression variance')
#run_functional_analysis(res_caudata, gene_universe_caudata, "Caudata")
res_plot_male_caud <- create_multiplot(subset(resdata_caudata, resdata_caudata$log2FoldChange > 0), "Caudata", "Male")
res_plot_female_caud <- create_multiplot(subset(resdata_female_caudata, resdata_female_caudata$log2FoldChange > 0), "Caudata", "Female")
write_results(dds_caudata$dds, resdata_caudata, "caud_F_vs_M.tab", "caud_DEsig0.05FC2_F_vs_M.tab")
###################
###Gracilis Code###
abundance_file_names_gracilis <- c("m_A_abundance.tsv","m_B_abundance.tsv" ,"m_C_abundance.tsv" ,"m_D_abundance.tsv", 
                          "f_A_abundance.tsv","f_B_abundance.tsv","f_C_abundance.tsv","f_D_abundance.tsv")
rownames_species_gracilis <- c("MaleA", "MaleB", "MaleC", "MaleD", "FemaleA", "FemaleB", "FemaleC", "FemaleD")
dds_results <- run_dds("Gracilis", 4, abundance_file_names_gracilis, rownames_species_gracilis)
gracilis_conditions <- dds_results$conditions
dds_gracilis <- dds_results$dds
gene_universe_gracilis <- dds_gracilis$gene_universe
res_gracilis <- dds_gracilis$res_male
rld_gracilis <- dds_gracilis$rld
resdata_gracilis <- dds_gracilis$resdata
resdata_female_gracilis <- dds_gracilis$resdata_female
female_linked_resdata_gracilis <- dds_gracilis$female_linked_resdata
male_linked_resdata_gracilis <- dds_gracilis$male_linked_resdata
grac_abundance <- dds_results$dds_results$abundance_df
plots_gracilis <- create_plots(dds_gracilis$dds,resdata_gracilis, rld_gracilis, gracilis_conditions,
                               "MA Plot Gracilis", "G. gracilis", "Heatmap Gracilis", 'G. gracilis gene expression variance')
#run_functional_analysis(res_gracilis, gene_universe_gracilis, "Gracilis")
res_plot_male_grac <- create_multiplot(subset(resdata_gracilis, resdata_gracilis$log2FoldChange > 0), "Gracilis", "Male")
res_plot_female_grac <- create_multiplot(subset(resdata_female_gracilis, resdata_female_gracilis$log2FoldChange > 0), "Gracilis", "Female")
write_results(dds_gracilis$dds, resdata_gracilis, "grac_F_vs_M.tab", "grac_DEsig0.05FC2_F_vs_M.tab")
###################
###Chondrus Code###
abundance_file_names_chondrus <- c("70mal_abundance.tsv","71mal_abundance.tsv" ,"72mal_abundance.tsv", 
                                   "67fem_abundance.tsv","68fem_abundance.tsv","69fem_abundance.tsv")
rownames_species_chondrus <- c("Male70", "Male71", "Male72", "Female67", "Female68", "Female69")
dds_results <- run_dds("Chondrus", 3, abundance_file_names_chondrus, rownames_species_chondrus)
chondrus_conditions <- dds_results$conditions
dds_chondrus <- dds_results$dds
gene_universe_chondrus <- dds_chondrus$gene_universe
res_chondrus <- dds_chondrus$res_male
rld_chondrus <- dds_chondrus$rld
resdata_chondrus <- dds_chondrus$resdata
resdata_female_chondrus <- dds_chondrus$resdata_female
female_linked_resdata_chondrus <- dds_chondrus$female_linked_resdata
male_linked_resdata_chondrus <- dds_chondrus$male_linked_resdata
chon_abundance <- dds_results$dds_results$abundance_df
plots_chondrus <- create_plots(dds_chondrus$dds,resdata_chondrus, rld_chondrus, chondrus_conditions,
                               "MA Plot Chondrus", "C. crispus", "Heatmap Chondrus", 'C. crispus gene expression variance')
#run_functional_analysis(res_chondrus, gene_universe_chondrus, "Chondrus")
res_plot_male_chon <- create_multiplot(subset(resdata_chondrus, resdata_chondrus$log2FoldChange > 0), "Chondrus", "Male")
res_plot_female_chon <- create_multiplot(subset(resdata_female_chondrus, resdata_female_chondrus$log2FoldChange > 0), "Chondrus", "Female")
violin_data_chondrus <- create_violin_plot(resdata_chondrus, "hochberg", 'C. crispus gene expression variance')
write_results(dds_chondrus$dds, resdata_chondrus, "chon_F_vs_M.tab", "chon_DEsig0.05FC2_F_vs_M.tab")
###################

# Some grids used for side by soide plots
gridExtra::grid.arrange(res_plot_male_chil, res_plot_female_chil,
                        res_plot_male_caud, res_plot_female_caud,
                        res_plot_male_grac, res_plot_female_grac,
                        res_plot_male_chon, res_plot_female_chon,
                        res_plot_male_verm, res_plot_female_verm,
                        nrow=5)

gridExtra::grid.arrange(res_plot_male_chil, res_plot_female_chil,
                        res_plot_male_grac, res_plot_female_grac,
                        nrow=2)

gridExtra::grid.arrange(plots_caudata$violin_data, plots_chilensis$violin_data,
                        plots_gracilis$violin_data, plots_chondrus$violin_data,
                        nrow=2)

gridExtra::grid.arrange(plots_caudata$volcano_data, plots_chilensis$volcano_data,
                        plots_gracilis$volcano_data, plots_chondrus$volcano_data,
                        nrow=2)

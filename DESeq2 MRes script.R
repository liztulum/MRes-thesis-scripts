#R script - DESeq2 MRes script - Sample data sets for HepG2, HepaRG and MCF7 cell line dosed with Doxorubicin and Niacinamide
#Datasets analysed using DESeq2.

#load libraries relevant to this analysis

library(DESeq2)
library(tidyverse)

library(dplyr)
library(EnhancedVolcano)
library(org.Hs.eg.db)
library(ggplot2)
library(plyr)

#read the csv files containing the sample data for each cell line and chemical

counts <- read.csv("HepG2_Doxorubicin_HCl_counts.csv", row.name=1)
head(counts)

counts2 <- read.csv("HepG2_Niacinamide_counts.csv", row.name=1)
head(counts2)

counts3 <- read.csv("HepaRG_Niacinamide_counts.csv", row.name=1)
head(counts3)

counts4 <- read.csv("HepaRG_Doxorubicin_HCl_counts.csv", row.name=1)
head(counts4)

counts5 <- read.csv("MCF-7_Niacinamide_counts.csv", row.name=1)
head(counts5)

counts6 <- read.csv("MCF-7_Doxorubicin_HCl_counts.csv", row.name=1)
head(counts6)

colSums(counts)

#create list of all files to input into a loop

files <- list.files(path = ".", pattern = "_counts.csv", full.names = T)
files

#loop created to analyse all the cell lines and chemicals individually

for(i in (files)){
  print(i)
  counts_file <- i
  meta_files <- gsub("_counts.csv", "_metadata.csv", i) # create metadata csv files from the counts data files
  counts <- read.csv(counts_file, row.names = 1)
  SampleTable <- read.csv(meta_files, header = T)
  
 
  
  SampleTable$CONCENTRATION <- as.factor(SampleTable$CONCENTRATION)
  
  
  dds <- DESeqDataSetFromMatrix(countData = counts,
                                colData = SampleTable,
                                design = ~ VESSEL_ID + CONCENTRATION)
  
  dds$condition <- relevel(dds$CONCENTRATION, ref = "0")
  dds <- estimateSizeFactors(dds)
  sizeFactors(dds)
  
  
  dds <- DESeq(dds)
  res <- results(dds)
  comparisons <- resultsNames(dds)
  normalisedCounts <- counts(dds, normalized = TRUE)

    for (j in comparisons[6:12]) {
      res <- results(dds, name = j)
      padj.cutoff <- 0.05
      lfc.cutoff <- 0.58
      res_table <- res %>%
        data.frame() %>%
        rownames_to_column(var="gene") %>%
        as_tibble() #convert the results table into a tibble
      sigOE <- res_table %>%
        filter(padj < padj.cutoff & abs(log2FoldChange) > lfc.cutoff)
      #write.table(sigOE, file = paste0(j, ".txt", i)) #need to add cell line & chemical?
      write.table(sigOE, file = paste0(gsub("_counts.csv", "_DESeq", i), "_", j, ".txt"), col.names = T, row.names = T, quote = F, sep = "\t")
    }
  res_norm <- lfcShrink(dds=dds, coef=2, type="normal")

  

  save(res_norm, file = "res_norm.RData")
  
  
  png(gsub("_counts.csv", "_MA.png", i))
  ma <- plotMA(res_norm)
  print(ma)
  dev.off()
  
  png(gsub("_counts.csv", "_MA2.png", i))
  ma2 <- plotMA(res_norm, alpha = 0.05)
  print(ma)
  dev.off()
  
  
  plot_EnhancedVolcano <- EnhancedVolcano(res_norm,
                                          lab = rownames(res),
                                          x = "log2FoldChange",
                                          y = "pvalue")
 
  ggsave(plot_EnhancedVolcano, filename = gsub("_counts.csv", "_EV.png", i))
  
  normalised_counts <- normalisedCounts %>%
  data.frame() %>%
     rownames_to_column(var="gene") %>%
     as_tibble()
   
   
   top20_sigOE_genes <- res_table %>%
     arrange(padj) %>%
     pull(gene) %>%
     head(n=20)
   
   top20_sigOE_norm <- normalised_counts %>%
     filter(gene %in% top20_sigOE_genes)
   
   
   
   pivot_top20_sigOE <- top20_sigOE_norm %>%
     pivot_longer(!gene, names_to = "samplename", values_to = "normalised") 
   pivot_top20_sigOE_join <- pivot_top20_sigOE %>%
     left_join(dplyr::select(SampleTable, X, CONCENTRATION), by = c("samplename" = "X"))
   
   
   ggplot <- ggplot(pivot_top20_sigOE_join,aes(x = gene, y = normalised, color = CONCENTRATION)) +
     geom_point() +
     scale_y_log10() +
     xlab("genes") +
     ylab("normalised_counts") +
     ggtitle("Top 20 Significant DE Genes") +
     theme_bw() +
     theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
     theme(plot.title = element_text(hjust = 0.5))
   
   ggsave(ggplot, filename = gsub("_counts.csv", "_ggplot.png", i))
    rld <- vst(dds, blind=TRUE)
  pca <- plotPCA(rld, intgroup="CONCENTRATION") 
  ggsave(pca, filename = gsub("_counts.csv", "_pca.png", i))
  vst_write <- assay(vst(dds, blind=FALSE))
  write.csv(as.data.frame(vst_write), file=gsub("_counts.csv", "_results.csv", i))
  
  vst_write_copy <- vst_write # COPY SO WE DONT RUIN DATA
  colnames(vst_write_copy) <- NULL # Getting rid of colnames so we can append to data
  cols <- colnames(vst_write) # pulling sample ids
  vst_dose <- SampleTable[SampleTable$X %in% colnames(vst_write),"CONCENTRATION" ] # ordering the metadata in the order our sampleid in the normalised data then pulling concentration out
  vst_dose_dt <- data.frame(matrix(ncol=175, nrow=1)) # making the dose dataframe ready to rbind
  vst_dose_dt[1, ] <- vst_dose # setting the first row to doses
  bmd_input <- rbind(vst_dose_dt, data.frame(vst_write_copy)) # combining dose row with gene normalised data
  rownames(bmd_input)[1] <- "Dose" # getting dose in the rownames
  toprow <- data.frame(matrix(ncol = 175, nrow=1))  # new data  frame to get the 'SampleID' bit in
  toprow[1, ] <- cols # set first row to sample ID. This will be our actual colnames name but we dont set it here
  rownames(toprow) <- "SampleID" # corner value
  final_bmd_input <- rbind(toprow, bmd_input)
  write.table(final_bmd_input, file = gsub("_counts.csv", "_BMDInput.txt", i), col.names = F, row.names =T, quote = F, sep = "\t") # write withour colnames as these are current X1. X2


  
  print(paste(i, "finished"))
}




















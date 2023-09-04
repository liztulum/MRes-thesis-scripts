library(glue)
library(purrr)
library(dplyr)
library(UpSetR)



upset_gene_intersect <- function(cell_line,chemical, conc1, conc2, conc3, conc4, conc5, conc6, conc7){
  
  # check if there are 6 or 7 concentrations available and make list of all the concentrations
  
    
    tsv_files <- c(file.path( glue::glue("{cell_line}_{chemical}_DESeq_CONCENTRATION_{conc1}_vs_0.txt")),
                   file.path( glue::glue("{cell_line}_{chemical}_DESeq_CONCENTRATION_{conc2}_vs_0.txt")),
                   file.path( glue::glue("{cell_line}_{chemical}_DESeq_CONCENTRATION_{conc3}_vs_0.txt")),
                   file.path( glue::glue("{cell_line}_{chemical}_DESeq_CONCENTRATION_{conc4}_vs_0.txt")),
                   file.path( glue::glue("{cell_line}_{chemical}_DESeq_CONCENTRATION_{conc5}_vs_0.txt")),
                   file.path( glue::glue("{cell_line}_{chemical}_DESeq_CONCENTRATION_{conc6}_vs_0.txt")),
                   file.path( glue::glue("{cell_line}_{chemical}_DESeq_CONCENTRATION_{conc7}_vs_0.txt")))
    
    
    conc_list <- c(conc1,conc2, conc3, conc4, conc5, conc6, conc7)
    

   
  
  # make a list of dataframes containing the data from the concentration with 1 or more degs
  dfs <- list()

  for (i in 1:length(tsv_files)){

    num_rows <- nrow(read.csv(tsv_files[i], sep = "\t"))
    
    conc <- toString(conc_list[i])

    if (num_rows > 0) {
    dfs[[conc]] <- read.csv(tsv_files[i], sep = "\t")
    }
  }


  # check how many of the concentrations have degs, based on the number of data frames in the list
  if (length(dfs) == 7){


    #combine data frames, given all the concentrations have 1 or more degs
    combined <- purrr::reduce(list(data.frame(gene = dfs[[1]]$gene, conc1 = 1),
                                   data.frame(gene = dfs[[2]]$gene, conc2 = 1),
                                   data.frame(gene = dfs[[3]]$gene, conc3 = 1),
                                   data.frame(gene = dfs[[4]]$gene, conc4 = 1),
                                   data.frame(gene = dfs[[5]]$gene, conc5 = 1),
                                   data.frame(gene = dfs[[6]]$gene, conc6 = 1),
                                   data.frame(gene = dfs[[7]]$gene, conc7 = 1)), full_join)
    combined[is.na(combined)] <- 0

    # give columns more meaningful names
    names(combined) <- c("genes",
                         names(dfs)[1],
                         names(dfs)[2],
                         names(dfs)[3],
                         names(dfs)[4],
                         names(dfs)[5],
                         names(dfs)[6],
                         names(dfs)[7])

  } else if (length(dfs) == 6){

    #combine data frames, given 6 of the concentrations have 1 or more degs
    combined <- purrr::reduce(list(data.frame(gene = dfs[[1]]$gene, conc1 = 1),
                                   data.frame(gene = dfs[[2]]$gene, conc2 = 1),
                                   data.frame(gene = dfs[[3]]$gene, conc3 = 1),
                                   data.frame(gene = dfs[[4]]$gene, conc4 = 1),
                                   data.frame(gene = dfs[[5]]$gene, conc5 = 1),
                                   data.frame(gene = dfs[[6]]$gene, conc6 = 1)), full_join)
    combined[is.na(combined)] <- 0

    # give columns more meaningful names
    names(combined) <- c("genes",
                         names(dfs)[1],
                         names(dfs)[2],
                         names(dfs)[3],
                         names(dfs)[4],
                         names(dfs)[5],
                         names(dfs)[6])

  } else if (length(dfs) == 5){

    #combine data frames, given 5 of the concentrations have 1 or more degs
    combined <- purrr::reduce(list(data.frame(gene = dfs[[1]]$Probe, conc1 = 1),
                                   data.frame(gene = dfs[[2]]$Probe, conc2 = 1),
                                   data.frame(gene = dfs[[3]]$Probe, conc3 = 1),
                                   data.frame(gene = dfs[[4]]$Probe, conc4 = 1),
                                   data.frame(gene = dfs[[5]]$Probe, conc5 = 1)), full_join)
    combined[is.na(combined)] <- 0

    # give columns more meaningful names
    names(combined) <- c("genes",
                         names(dfs)[1],
                         names(dfs)[2],
                         names(dfs)[3],
                         names(dfs)[4],
                         names(dfs)[5])

  } else if (length(dfs) == 4){

    #combine data frames, given 4 of the concentrations have 1 or more degs
    combined <- purrr::reduce(list(data.frame(gene = dfs[[1]]$gene, conc1 = 1),
                                   data.frame(gene = dfs[[2]]$gene, conc2 = 1),
                                   data.frame(gene = dfs[[3]]$gene, conc3 = 1),
                                   data.frame(gene = dfs[[4]]$gene, conc4 = 1)), full_join)
    combined[is.na(combined)] <- 0

    # give columns more meaningful names
    names(combined) <- c("genes",
                         names(dfs)[1],
                         names(dfs)[2],
                         names(dfs)[3],
                         names(dfs)[4])


  } else if (length(dfs) == 3){

    #combine data frames, given 3 of the concentrations have 1 or more degs
    combined <- purrr::reduce(list(data.frame(gene = dfs[[1]]$gene, conc1 = 1),
                                   data.frame(gene = dfs[[2]]$gene, conc2 = 1),
                                   data.frame(gene = dfs[[3]]$gene, conc3 = 1)), full_join)
    combined[is.na(combined)] <- 0

    # give columns more meaningful names
    names(combined) <- c("genes",
                         names(dfs)[1],
                         names(dfs)[2],
                         names(dfs)[3])


  } else if (length(dfs) == 2){

    #combine data frames, given 2 of the concentrations have 1 or more degs
    combined <- purrr::reduce(list(data.frame(gene = dfs[[1]]$gene, conc1 = 1),
                                   data.frame(gene = dfs[[2]]$gene, conc2 = 1)), full_join)
    combined[is.na(combined)] <- 0

    # give columns more meaningful names
    names(combined) <- c("genes",
                         names(dfs)[1],
                         names(dfs)[2])


  } else {

    # the chemical either has only one or no concentration with degs
    # an upset plot is not needed
    print("There are not enough data sets containing data!")
  }

  # check the chemical has at least 2 or more concentrations with degs
  if (length(dfs) > 1) {

    # get  vector names of concentrations being compared from the combined data frame column headers
    # needed to set concentrations displayed order in the upset plot
    upset_sets <- names(combined)[-1]

    # create upset plot
    plot <- UpSetR::upset(combined, nsets = length(names(combined)),keep.order = T, sets = upset_sets)

    # export upset plot
   png(file.path(glue::glue("{cell_line}_{chemical}_shared_genes_upset.png")), width = 1300, height = 800, res = 100)
   print(plot)
   dev.off()

  }
}

upset_gene_intersect("HepaRG", "Niacinamide", 8000, 1600, 320, 64, 12.8, 2.56, 0.512)
upset_gene_intersect("HepG2", "Niacinamide", 60000, 12000, 2400, 480, 96, 19.2, 3.84)
upset_gene_intersect("MCF-7", "Niacinamide", 60000, 12000, 2400, 480, 96, 19.2, 3.84)
upset_gene_intersect("HepaRG", "Doxorubicin_HCl", 1, 0.2, 0.04, 0.008, 0.0016, 0.00032, "6.4e.05")
upset_gene_intersect("HepG2", "Doxorubicin_HCl", 1, 0.2, 0.04, 0.008, 0.0016, 0.00032, "6.4e.05")
upset_gene_intersect("MCF-7", "Doxorubicin_HCl", 1, 0.2, 0.04, 0.008, 0.0016, 0.00032, "6.4e.05")

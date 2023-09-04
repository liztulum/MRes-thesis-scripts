library(logger)
library(glue)
library(tibble)

#' Create output directory if it doesn't exist. Warn if existing directory is not empty
#' @param path directory path to create (string)
#' @importFrom logger log_debug log_warn
#'
#' @export
create_dir <- function(path) {
  if (!dir.exists(path)) {
    logger::log_debug("Creating output directory {dQuote(path)}")
    dir.create(path, recursive = T)
  } else {
    logger::log_debug("Output directory {dQuote(path)} is existing")
    if (length(dir(path = path, all.files = FALSE)) > 0) {
      logger::log_warn("Output directory {dQuote(path)} is not empty, files may be overwritten")
    }
  }
  return(path)
}

#' BMD and pathway data checks and error handling
#'
#' 
#' @param data variable to check intended to be either bmd data or pathway data from BMDExpress2
#' @param function_name name of the function running this check function. This is used in error messages
#
#' @importFrom logger log_fatal
#' @importFrom glue glue
check_data <- function(data, function_name) {
  list(
    type = function() {
      if (!class(data) == "data.frame") {
        logger::log_fatal("data passed to function '{function_name}' is not a data.frame")
        stop(glue::glue("Please supply a data.frame object to function '{function_name}' with columns from BMDExpress2 results BMD export."))
      }
    },
    
    cols = function(keep_cols) {
      is_missing <- !(keep_cols %in% colnames(data))
      if (any(is_missing)) {
        missing_cols <- keep_cols[is_missing]
        logger::log_fatal("data passed to function '{function_name}' is missing column(s) {missing_cols}")
        stop(glue::glue("Please supply a data.frame object to function '{function_name}' with columns from BMDExpress2 results BMD export."))
      }  
    },
    
    empty = function() {
      if (nrow(data) == 0) {
        e <- "data passed to function '{function_name}' has no data."
        logger::log_fatal(e)
        stop(e)
      }
    }, 
    analysis_length = function() {
      if (length(unique(data$Analysis)) != 1) {
        logger::log_fatal("data passed to function '{function_name}' has more than one analysis group.")
        stop("pathway_filtered passed to function '{function_name}' has more than one analysis group. Please pass data with one unique value of column 'Analysis' only")
      }
    }
  )
}

#' Write filtered data set
#'
#' Write a data frame of data to csv with a table of filter parameters used writtten above data.
#' 
#' @param filtered_data Dataframe of data to be written
#' @param filter_info Dataframe of filter parameters to be written above the filtered_data
#' @param output_file_name the file name (included '.csv') to be written
#' @param output_dir the directory where the file should be written.
#' @importFrom logger log_fatal log_info
#' @importFrom glue glue
write_filtered_data <- function(filtered_data, output_file_name, filter_info = NULL, output_dir = ".") {
  # filter parameters
  filter_by <- tibble::lst(filtered_data, 
                           filter_info)
  
  for (f in names(filter_by)) {
    filter_value <- filter_by[[f]]
    if(!is.null(filter_value)) {
      if(!is.data.frame(filter_value)) {
        logger::log_fatal("{f} '{filter_value}' passed to function 'write_filtered_data' is not a dataframe")
        stop(glue::glue("{filter_value} passed to function 'write_filtered_data' is not a dataframe. Please supply a dataframe"))
      }
    }
  }
  
  create_dir(output_dir)
  
  # Write filter parameters table and the filtered data to a csv file by appending tables on top another (with blank dataframe for spacing)
  file_path <- file.path(output_dir, output_file_name)
  write.table(filter_info, file_path, col.names = T, sep = ", ", append = F, row.names = F)
  suppressWarnings(write.table(data.frame(), file_path, col.names = T, sep = ",", append = T))
  suppressWarnings(write.table(filtered_data, file_path, col.names = T, sep = ",", append = T, row.names = F))
  
  logger::log_info("Filtered data written to {file_path}.")
}


#' Filter Significant Gene BMDs
#'
#' Filtering for significant can include filtering out bmds which exceed the highest tested concentration, filtering based on the fit-p value and on the BMDU and BMDL ratio.
#' This function can be run on different analysis groups (cell line - chemical groups) but it must be noted that different analysis groups likely have different highest concentrations and therefore require a different highest_conc_filter value. If this this the case the function should be run on a subset of the input data for one analysis.
#' This function will write the returned data.frame of filtered gene level BMDs to csv file at {output_dir} with file name {output_prefix}_gene_bmd_filtered.csv
#' @param bmd_data Dataframe subset of BMDExpress gene bmd data. Dataframe should only include one analysis (cell line - chemical) group.
#' Must include the following columns as a minimum
#'  - Analysis
#'  - Probe ID
#'  - Entrez Gene IDs
#'  - Genes Symbols
#'  - Best Model
#'  - Best BMD
#'  - Best BMDL
#'  - Best BMDU
#'  - Best fitPValue
#'  - Best fitLogLikelihood
#'  - Best AIC
#'  - Best BMD/BMDL
#'  - Best BMDU/BMDL
#'  - Best BMDU/BMD
#'  - Max Fold Change
#'  - Max Fold Change Absolute Value
#' @param highest_conc_filter Numeric parameter indicating the value to filter 'Best BMD'. Data are filtered 'Best BMD' <= highest_conc_filer. To not use the filter set highest_conc_filter as NULL.
#' @param bmdl_bmdu_filter Numeric parameter indicating the value to filter 'Best BMDU/BMDL'. Data are filtered 'Best BMDU/BMDL' < bmdl_bmdu_filter. To not use the filter set bmdl_bmdu_filter as NULL. Default = 40
#' @param fitP_filter Numeric parameter indicating the value to filter on 'Best fitPValue'. Data are filtered 'Best fitPValue' > fitP_filter. To not use the filter set fitP_filter as NULL. Default = 0.1
#' @param gene_csv_output_prefix Character string to use as the prefix the the outputed csv of filtered data. Recommended is the analysis group from BMDExpress. If no csv is required, use NULL (default).
#' @param gene_csv_output_dir The directory path to where the output csv of filtered data should be written to. This directory should already exist.

#' @importFrom logger log_fatal log_info
#' @importFrom glue glue
#' @return dataframe of filtered gene BMD data
#' @export
filter_gene_bmds <- function(bmd_data, highest_conc_filter = NULL, bmdl_bmdu_filter = 40, fitP_filter = 0.1, gene_csv_output_prefix = NULL, gene_csv_output_dir = ".") {
  keep_cols <- c(
    "Analysis", # chosen as the most relavent data columns to save in csvs
    "Probe ID",
    "Entrez Gene IDs",
    "Genes Symbols",
    "Best Model",
    "Best BMD",
    "Best BMDL",
    "Best BMDU",
    "Best fitPValue",
    "Best fitLogLikelihood",
    "Best AIC",
    "Best BMD/BMDL",
    "Best BMDU/BMDL",
    "Best BMDU/BMD",
    "Max Fold Change",
    "Max Fold Change Absolute Value"
  )
  ## Parameter Error Handling
  # bmd_data
  check <- check_data(bmd_data, "filter_gene_bmds")
  check$type()
  check$cols(keep_cols)
  check$empty()
  
  # filter parameters
  filter_by <- tibble::lst(fitP_filter, 
                           bmdl_bmdu_filter, 
                           highest_conc_filter)
  
  for (f in names(filter_by)) {
    filter_value <- filter_by[[f]]
    if(!is.null(filter_value)) {
      if(!is.numeric(filter_value)) {
        logger::log_fatal("{f} '{filter_value}' passed to function 'filter_gene_bmds' is not numeric")
        stop(glue::glue("{filter_value} passed to function 'filter_gene_bmds' is not numeric. Please supply numeric value or NULL."))
      }
    }
  }
  
 
  FC_cols <- colnames(bmd_data)[grep("FC Dose Level", colnames(bmd_data))]
  keep_cols <- c(keep_cols, FC_cols)
  
  bmd_filtered <- bmd_data[bmd_data$`Best BMD` != "none" & bmd_data$`Best BMD` != "NaN", keep_cols] # get rid of nones and remove unneeded columns
  
  for (i in 1:ncol(bmd_filtered)) {
    if (!any(is.na(suppressWarnings(as.numeric(bmd_filtered[, i]))))) { # check whether columns can be converted to numeric, if so convert. (NAs appear when non numeric characters are attempted to be converted so we check for NAs to determine suitability) 
      bmd_filtered[, i] <- as.numeric(bmd_filtered[, i])
    }
  }
  if (!is.null(highest_conc_filter)) {
    bmd_filtered <- bmd_filtered[bmd_filtered[, "Best BMD"] <= highest_conc_filter, ]
  } else {
    highest_conc_filter <- "None"
  }
  if (!is.null(bmdl_bmdu_filter)) {
    bmd_filtered <- bmd_filtered[bmd_filtered$`Best BMDU/BMDL` < bmdl_bmdu_filter, ]
  } else {
    bmdl_bmdu_filter <- "None"
  }
  if (!is.null(fitP_filter)) {
    bmd_filtered <- bmd_filtered[bmd_filtered$`Best fitPValue` > fitP_filter, ]
  } else {
    fitP_filter <- "None"
  }
  
  if (!is.null(gene_csv_output_prefix)) {
    # Get table of filter parameters used
    filters_used <- data.frame(c(
      paste0("Best BMD <= ", highest_conc_filter),
      paste0("Best BMDU/BMDL < ", bmdl_bmdu_filter),
      paste0("Best fitPValue > ", fitP_filter)
    ))
    colnames(filters_used) <- "These data are filtered using the following column filters:"
    write_filtered_data(filtered_data = bmd_filtered,
                        filter_info = filters_used,
                        output_file_name = paste0(gene_csv_output_prefix, "_gene_bmd_filtered.csv"),
                        output_dir = gene_csv_output_dir) 
    }
  return(bmd_filtered)
}

#' Calculate gene level PoDs from BMD data
#'
#' Using data from the filter_gene_bmds function, PoDs are calculated via 2 method:
#' 1. Average of BMD/BMDLs of 20 genes with the largest fold change.
#' 2. Average of the 25th and 75th percentile of BMD/BMDLs
#' This function assumed all genes inputted to the function have significant BMD/BMDLs and that only one analysis group is passed (eg one cell line- chemical)
#' @param bmd_filtered output of the filter_gene_bmds. Data needs columns "Probe ID", "Max Fold Change Absolute Value" and "Best BMD" or "Best BMDL"
#' @param bmd_param "BMD" or "BMDL" indicating which value to base calculations around.
#' @param bmdExpress_input_dir directory path where BMD input files for BMDExpress are located
#' @importFrom logger log_info
#' @importFrom glue glue
#' @importFrom stats quantile
#' 
#' @return dataframe of gene level PoDs for one analysis group
#' @export
calculate_gene_PoD <- function(bmd_filtered, bmd_param, bmdExpress_input_dir) {
  # Parameter Error Handling
  # bmd_param
  check_cols <- c("BMD", "BMDL")
  if (!toupper(bmd_param) %in% check_cols) {
    e <- glue::glue("bmd_param '{bmd_param}' passed to function 'filter_gene_bmds' is not {paste(check_cols, collapse = ' or ')}")
    logger::log_fatal(e)
    stop(e)
  }
  
  col_names <- c(
    "Analysis", # minimum columns needed for filtering
    "Best BMDU/BMDL",
    "Max Fold Change Absolute Value",
    paste0("Best ", bmd_param)
  )
  check <- check_data(bmd_filtered, "calculate_gene_PoD")
  check$type()
  check$cols(col_names)
  check$empty()
  check$analysis_length()
  bmd_filtered[, paste0("Best ", bmd_param)] <- as.numeric(bmd_filtered[, paste0("Best ", bmd_param)])
  
  if(nrow(bmd_filtered) >= 20) {
    top20FC <- bmd_filtered[order(bmd_filtered$`Max Fold Change Absolute Value`, decreasing = T)[1:20], paste0("Best ", bmd_param)] # get BMD/BMDL values for the genes with 20 highest FC
    PoD_top20FC <- mean(top20FC)
    analysis_name <- unique(bmd_filtered$Analysis)
  } else {
    PoD_top20FC <- mean(bmd_filtered[, paste0("Best ", bmd_param)])
    analysis_name <- paste0(unique(bmd_filtered$Analysis), " *(", nrow(bmd_filtered), " bmds)")
  }
  PoD_percentile_25 <- quantile(bmd_filtered[, paste0("Best ", bmd_param)], c(.25))
  PoD_percentile_75 <- quantile(bmd_filtered[, paste0("Best ", bmd_param)], c(.75))
  PoD_percentile <- mean(bmd_filtered[bmd_filtered[, paste0("Best ", bmd_param)] >= PoD_percentile_25 & 
                                        bmd_filtered[, paste0("Best ", bmd_param)] >= PoD_percentile_75,
                                      paste0("Best ", bmd_param)])
  lowest_probe <-  bmd_filtered[order(bmd_filtered[, paste0("Best ", bmd_param)], decreasing = F)[1],] 
  lowest_probe_BMD_param <-  lowest_probe[, paste0("Best ", bmd_param)] # get lowest BMD/BMDL Conc
  lowest_probe <-  lowest_probe[, "Probe ID"] # get probe name with lowest bmd param
  
  lowest_conc <- get_highest_lowest_conc_value(analysis = analysis_name, bmdExpress_input_dir = bmdExpress_input_dir, get_value = "Lowest" )
  bmd_filtered_above_1st_conc <- bmd_filtered[bmd_filtered[, paste0("Best ", bmd_param)] >= lowest_conc,]
  lowest_bmd_filtered_above_1st_conc <- bmd_filtered_above_1st_conc[order(bmd_filtered_above_1st_conc[, paste0("Best ", bmd_param)], decreasing = F)[1],]
  lowest_bmd_filtered_above_1st_conc_BMD_param <- lowest_bmd_filtered_above_1st_conc[, paste0("Best ", bmd_param)] # get lowest BMD/BMDL Conc above lowest tested concentration
  lowest_bmd_filtered_above_1st_conc_probe <- lowest_bmd_filtered_above_1st_conc[, "Probe ID"] # get lowest BMD/BMDL Probe above lowest tested concentration
  
  PoD_percentiles <- t(as.data.frame(quantile(bmd_filtered[, paste0("Best ", bmd_param)], c(.05, 0.10, .2, .3, .4, .5, .6 ,.7, .8, .9, .95))))
  
  gene_PoDs <- data.frame(analysis_name,
                          PoD_top20FC,
                          PoD_percentile,
                          "",
                          lowest_probe_BMD_param,
                          lowest_probe,
                          lowest_bmd_filtered_above_1st_conc_BMD_param,
                          lowest_bmd_filtered_above_1st_conc_probe,
                          stringsAsFactors = F)
  colnames(gene_PoDs) <- c("Analysis",
                           paste0("Avg of 20 ", bmd_param, "s with highest FC"),
                           paste0("Avg of 25th-75th percentile of ", bmd_param, "s"),
                           "-",
                           paste0("Lowest probe ", bmd_param),
                           paste0("Lowest probe ID"),
                           paste0("Lowest probe ", bmd_param, "after lowest dose"),
                           paste0("Lowest probe ID after lowest dose")
  )
  gene_PoDs <- cbind(gene_PoDs, PoD_percentiles)
  row.names(gene_PoDs) <- analysis_name
  return(gene_PoDs)
}

#' Filter Significant Pathway BMDs
#'
#' BMDExpress2 Pathway level BMD results can be extracted using bmdexpress2-cmd export --bm2-file {bm2_file_name.bm2} --analysis-group categorical --output-file-name {filename.txt}. These data should be filtered using this function for significance.
#' Filtering for significant pathway level BMDs include:
#' 1. Filtering by number of total genes (with a dose response) in the pathway.
#' 2. Filtering by the number of significant dose responsive genes found in the pathway. NOTE: significant genes to filter pathway BMDs are defined in the BMDExpress2 configuration file not by this R package (NOT filter_gene_bmds()).
#' 3. Filtering by Fishers 2-tail p value.
#' This function can be run on multiple analysis groups if all analysis groups should be filtered the same.
#' This function will write the returned data.frame of filtered pathway level BMDs to csv file at {path_csv_output_dir} with file name {path_csv_output_prefix}_path_bmd_filtered.csv

#'
#' @param path_data Dataframe of BMDExpress pathway bmd data. Minimum required columns:
#' - Input Genes
#' - Genes That Passed All Filters
#' - Fisher's Exact Two Tail
#' @param min_total_genes Numeric parameter indicating the value to filter 'Input Genes'. Data are filtered 'Input Genes' >= min_total_genes To not use the filter set min_total_genes as NULL. Default = 3
#' @param min_sig_genes Numeric parameter indicating the value to filter 'Genes That Passed All Filters'. Data are filtered 'Genes That Passed All Filters' <= min_sig_genes. To not use the filter set min_sig_genes as NULL. Default = 2
#' @param fishers_p_val Numeric parameter indicating the value to filter on 'Fisher's Exact Two Tail'. Data are filtered 'Fisher's Exact Two Tail' < fishers_p_val To not use the filter set fishers_p_val as NULL. Default = 0.1
#' @param path_csv_output_prefix Character string to use as the prefix the the outputed csv of filtered data. Recommended is the analysis group from BMDExpress. If no csv is required, use NULL (default).
#' @param path_csv_output_dir The directory path to where the output csv of filtered data should be written to. This directory should already exist.
#' 
#' @importFrom glue glue
#' @importFrom logger log_fatal log_info
#' @return dataframe of filtered pathway BMD data
#' @export
filter_pathway_bmds <- function(path_data, min_total_genes = 3, min_sig_genes = 2, fishers_p_val = 0.1, path_csv_output_prefix = NULL, path_csv_output_dir = ".") {
  # Parameter Error Handling
  # path_data
  col_names <- c("Analysis", "Input Genes", "Genes That Passed All Filters", "Fisher's Exact Two Tail") # minimum columns for filtering
  check <- check_data(path_data, "filter_pathway_bmds")
  check$type()
  check$cols(col_names)
  check$empty()
  

  # filter parameters
  filter_by <- tibble::lst(min_total_genes, 
                           min_sig_genes, 
                           fishers_p_val)
  
  for (f in names(filter_by)) {
    filter_value <- filter_by[[f]]
    if(!is.null(filter_value)) {
      if(!is.numeric(filter_value)) {
        logger::log_fatal("{f} '{filter_value}' passed to function 'filter_pathway_bmds' is not numeric")
        stop(glue::glue("'{filter_value}' passed to function 'filter_pathway_bmds' is not numeric. Please supply numeric value or NULL."))
      }
    }
  }
  
  # path_csv_output_dir
  create_dir(path_csv_output_dir)

  path_data <- path_data[path_data$`Genes That Passed All Filters` != 0, ] # remove blank pathway enrichment rows (as no genes present)
  for (i in 1:ncol(path_data)) {
    if (!any(is.na(suppressWarnings(as.numeric(path_data[, i]))))) {
      path_data[, i] <- as.numeric(path_data[, i])
    }
  }
  
  filtered_path_data <- path_data
  if (!is.null(min_total_genes)) {
    filtered_path_data <- filtered_path_data[filtered_path_data$`Input Genes` >= min_total_genes, ]
  } else {
    min_total_genes <- "None"
  }
  if (!is.null(min_sig_genes)) {
    filtered_path_data <- filtered_path_data[filtered_path_data$`Genes That Passed All Filters` >= min_sig_genes, ]
  } else {
    min_sig_genes <- "None"
  }
  if (!is.null(fishers_p_val)) {
    filtered_path_data <- filtered_path_data[filtered_path_data$`Fisher's Exact Two Tail` < fishers_p_val, ]
  } else {
    fishers_p_val <- "None"
  }
  
  if (!is.null(path_csv_output_prefix)) {
    # Get a table of filter values
    filters_used <- data.frame(c(
      paste0("Input Genes >= ", min_total_genes),
      paste0("Genes That Passed All Filters >= ", min_sig_genes),
      paste0("Fisher's Exact Two Tail < ", fishers_p_val)
    ))
    colnames(filters_used) <- "These data are filtered using the following column filters:"
    write_filtered_data(filtered_data = filtered_path_data,
                        filter_info = filters_used,
                        output_file_name = paste0(path_csv_output_prefix, "_pathway_bmd_filtered.csv"),
                        output_dir = path_csv_output_dir) 
  }
  return(filtered_path_data)
}

#' Calculate pathway level PoDs from BMD data
#'
#' Using data from the filter_pathway_bmds function, PoDs are calculated via 3 methods.
#' 1. Average of BMD/BMDLs of 20 pathways with the lowest Fishers 2-tail p values.
#' 2. Average of BMD/BMDLs of 20 pathways with the lowest BMD/BMDLs mean values.
#' 3. Taking the lowest pathways BMD/BMDLs mean value.
#'
#' This function assumed all pathways inputted to the function have significant BMD/BMDLs and that only one analysis group is passed (eg one cell line- chemical)
#'
#' @param pathway_filtered output of the filter_pathway_bmds. Minimum columns required:
#' - Analysis
#' - GO/Pathway/Gene Set Name
#' - GO/Pathway/Gene Set ID
#' - Fisher's Exact Two Tail
#' - Mean BMD or Mean BMDL (dependent on the bmd_param parameter)
#'
#' @param bmd_param "BMD" or "BMDL" indicating which value to base calculations around.
#' @param bmdExpress_input_dir directory path where BMD input files for BMDExpress are located
#' @importFrom logger log_info
#' @importFrom stats quantile
#' 
#' @return dataframe of pathway level PoDs for one analysis group
#' @export
calculate_pathway_PoD <- function(pathway_filtered, bmd_param, bmdExpress_input_dir) {
  # Parameter Error Handling
  # bmd_param
  if (toupper(bmd_param) != "BMD" & toupper(bmd_param) != "BMDL") {
    logger::log_fatal("bmd_param '{bmd_param}' passed to function 'calculate_gene_PoD' is not 'BMD' or 'BMDL'.")
    stop("bmd_param passed to function 'calculate_gene_PoD' is not 'BMD' or 'BMDL'.")
  } else {
    bmd_param <- toupper(bmd_param)
  }
  # pathway_filtered
  #replace version 2.3 columns anems with version 2.0 column names
  colnames(pathway_filtered) <- gsub("GO/Pathway/Gene Set/Gene Name", "GO/Pathway/Gene Set Name", colnames(pathway_filtered))
  colnames(pathway_filtered) <- gsub("GO/Pathway/Gene Set/Gene ID", "GO/Pathway/Gene Set ID", colnames(pathway_filtered))
  
  col_names <- c("Analysis", "GO/Pathway/Gene Set Name", "GO/Pathway/Gene Set ID", "Fisher's Exact Two Tail", paste0(bmd_param, " Mean"))
  check <- check_data(pathway_filtered, "calculate_pathway_PoD")
  check$type()
  check$cols(col_names)
  check$empty()
  check$analysis_length()
  
  if(nrow(pathway_filtered) >= 20) {
    min20bmd <- pathway_filtered[order(pathway_filtered[, paste0(bmd_param, " Mean")], decreasing = F)[1:20], paste0(bmd_param, " Mean")] # Get 20 lowest BMD/BMDL values
    minbmd <- min(min20bmd) # get the lowest BMD/BMDL value (complete.cases in line above retains orignal order of data so min has to be used in this line)
    min20p <- pathway_filtered[order(pathway_filtered$`Fisher's Exact Two Tail`, decreasing = F)[1:20], paste0(bmd_param, " Mean")] # Get BMD/BMDL values for the pathways with the 20 lowest fishers 2 tail p value
    analysis_name <- unique(pathway_filtered$Analysis)
    } else {
    min20bmd <- pathway_filtered[, paste0(bmd_param, " Mean")] # Get 20 lowest BMD/BMDL values
    minbmd <- min(min20bmd) # get the lowest BMD/BMDL value (complete.cases in line above retains orignal order of data so min has to be used in this line)
    min20p <- pathway_filtered[, paste0(bmd_param, " Mean")]
    analysis_name <- paste0(unique(pathway_filtered$Analysis), " *(", length(min20bmd), " bmds)")
  }
  
  PoD_min20bmd <- mean(min20bmd)
  PoD_minbmd <- minbmd
  PoD_min20p <- mean(min20p)

  lowest_conc <- get_highest_lowest_conc_value(analysis = analysis_name, bmdExpress_input_dir = bmdExpress_input_dir, get_value = "Lowest")
  path_filtered_above_1st_conc <- pathway_filtered[pathway_filtered[, paste0(bmd_param, " Mean")] >= lowest_conc,]
  lowest_path_filtered_above_1st_conc <- path_filtered_above_1st_conc[order(path_filtered_above_1st_conc[, paste0(bmd_param, " Mean")], decreasing = F)[1],]
  lowest_path_filtered_above_1st_conc_BMD_param <- lowest_path_filtered_above_1st_conc[, paste0(bmd_param, " Mean")] # get lowest BMD/BMDL Conc above lowest tested concentration
  lowest_path_filtered_above_1st_conc_path <- lowest_path_filtered_above_1st_conc[,"GO/Pathway/Gene Set Name"] # get lowest BMD/BMDL pathway above lowest tested concentration
  PoD_percentiles <- t(as.data.frame(quantile(pathway_filtered[, paste0(bmd_param, " Mean")], c(.05, 0.10, .2, .3, .4, .5, .6 ,.7, .8, .9, .95))))
  
  pathway_PoDs <- data.frame(analysis_name,
                             PoD_min20bmd,
                             PoD_minbmd,
                             PoD_min20p,
                             "",
                             lowest_path_filtered_above_1st_conc_BMD_param,
                             lowest_path_filtered_above_1st_conc_path,
                             stringsAsFactors = F)
  colnames(pathway_PoDs) <- c("Analysis",
                  paste0("Avg of 20 lowest pathway ", bmd_param, "s"),
                  paste0("The lowest pathway ", bmd_param), 
                  paste0("Avg of 20 pathway ", bmd_param, "s with lowest 2-tail fisher P values"),
                  "-",
                  paste0("Lowest pathway ", bmd_param, "after lowest dose"),
                  paste0("Lowest pathway after lowest dose")
  )
  pathway_PoDs <- cbind(pathway_PoDs, PoD_percentiles)                
  rownames(pathway_PoDs) <- analysis_name
  return(pathway_PoDs)
}


#' Determine the highest concentration filter to be used based on a combination of highest_conc_filter and bmdExpress_input_dir parameters
#'
#' This function's set the rules of high_conc_filter is determined for use in function calculate_PoDs_from_BMDExpress2.
#' - Rule 1: When a numeric value is given to highest_conc_filter, then this value will always be used.
#' - Rule 2: When both highest_conc_filter and bmdExpress_input_dir are NULL then the final highest concentration wont be filtered and therefore final value set to NULL.
#' - Rule 3: When highest_conc_filter is NULL but bmdExpress_input_dir is a valid file path to text files of the input files for BMDExpress2, the highest concentration value is set to the highest value found in the files which matches the BMDExpress2 analysis group.
#'
#' @param analysis BMDExpress2 analysis group name. This shall be used to match bmdExpress_input_dir files if present.
#' @param bmdExpress_input_dir The directory path for the BMDExpress2 input text file of which highest concentration shall be extracted. Can also be NULL in cases of rule 1 and 2.
#' @param highest_conc_filter The parameter value passed to calculate_PoDs_from_BMDExpress2. Can be numeric value for rule 1 or NULL for rules 2 and 3.
#' @param get_value "Highest" or "Lowest" for concentration to be used
#' 
#' @importFrom logger log_info
#' @importFrom utils read.delim
#' @return The final determined highest_conc_filter value
get_highest_lowest_conc_value <- function(analysis, bmdExpress_input_dir, highest_conc_filter = NULL, get_value = "Highest") {
  get_value_accepted <- c("Highest", "Lowest")
  if(!get_value %in% get_value_accepted){
    logger::log_fatal("get_value passes to get_highest_lowest_conc_value is not either of {get_value_accepted}. Value passed: {get_value}.")
    stop("Incorrect get_value parameter.")
  }
  if (is.null(highest_conc_filter) | get_value == "Lowest") {
    if (!is.null(bmdExpress_input_dir)) {
      input_list <- list.files(path = bmdExpress_input_dir)
      analysis_input <- input_list[unlist(lapply(
        X = gsub(".txt", "", input_list),
        FUN = grepl,
        x = analysis
      ))] # Get BMDExpress2 input file from bmdExpress_input_dir which has the part of the analysis string in - The original data
      
      if(length(analysis_input) == 0) {
        logger::log_fatal("BMDExpress input file for {analysis} cannot be found at {bmdExpress_input_dir}. This means a concentration to filter cannot be determined and analysis is terminated.")
        stop("Analysis group not found in provide BMDExpress input files.")
      }
      
      input_file_path  <- file.path(bmdExpress_input_dir, analysis_input)
      logger::log_debug("{get_value} concentration for {analysis} is pulled from {input_file_path}")
      input_data <- read.delim(input_file_path,
                               header = F, skip = 1, nrows = 1)[1, -1]
      if(get_value == "Highest") {
        chosen_filter <- max(input_data) # Get max concentration for the BMD input file from the second row of the text file.
      } else {
        input_data_not0 = input_data[input_data!=0]
        chosen_filter <- min(input_data_not0) # Get min concentration for the BMD input file from the second row of the text file.
      }
    } else {
      chosen_filter <- NULL
    }
  } else {
    chosen_filter <- highest_conc_filter
  }
  if(get_value == "Highest") {
    logger::log_info("For {analysis} the chosen highest_conc_filter is {chosen_filter}.")
  } else {
    logger::log_info("For {analysis} the chosen lowest concentration is {chosen_filter}.")
  }
  return(chosen_filter)
}

#' Write PoD file
#'
#' PoD's shall be written to a csv file at {pod_csv_output_prefix} with file name {pod_csv_output_dir}_PoDs.csv. File will include filters applied to the data and gene and pathway level PoDs
#
#' @param defined_conc concentration to be used (TODO Jade to add more detail here)
#' @param bmd_param "BMD" or "BMDL" indicating which value to base calculations around
#' @param collated_gene_PoDs Calculated gene level PoDs (data to be included in csv)
#' @param collated_pathway_PoDs Calculated pathway level PoDs (data to be included in csv)
#' @param pod_csv_output_prefix Character string to use as the prefix the the outputed csv of PoDs calculated at both the gene and pathway level. Recommended is the analysis group from BMDExpress. If no csv is required, use NULL (default).
#' @param pod_csv_output_dir The directory path to where the output csv of PoDs calculated at both the gene and pathway level should be written to. This directory should already exist. If NULL no file will be written.

#' @inheritParams filter_gene_bmds
#' @inheritParams filter_pathway_bmds
#' 
#' @importFrom logger log_info
write_PoDs_to_file <- function(pod_csv_output_prefix,
                               defined_conc,
                               highest_conc_filter,
                               bmd_param, 
                               bmdl_bmdu_filter,
                               fitP_filter,
                               min_total_genes,
                               min_sig_genes,
                               fishers_p_val,
                               pod_csv_output_dir,
                               collated_gene_PoDs,
                               collated_pathway_PoDs) {
    
    filters_used_gene <- data.frame(c(
      paste0("PoDs based on ", bmd_param),
      paste0("Best ", bmd_param, " <= ", defined_conc),
      paste0("Best BMDU/BMDL < ", bmdl_bmdu_filter),
      paste0("Best fitPValue > ", fitP_filter)
    ))
    colnames(filters_used_gene) <- "Gene BMDExpress2 results data are filtered using the following column filters:"
    
    
    if (is.null(min_total_genes)) min_total_genes <- "None"
    if (is.null(min_sig_genes)) min_sig_genes <- "None"
    if (is.null(fishers_p_val)) fishers_p_val <- "None"
    
    filters_used_path <- data.frame(c(
      paste0("Input Genes >= ", min_total_genes),
      paste0("Genes That Passed All Filters >= ", min_sig_genes),
      paste0("Fisher's Exact Two Tail < ", fishers_p_val)))
    colnames(filters_used_path) <- "Pathway BMDExpress2 results are filtered using the following column filters:"
    
    create_dir(pod_csv_output_dir)
    
    file_path <- file.path(pod_csv_output_dir, paste0(pod_csv_output_prefix,"_", bmd_param, "_PoDs.csv"))
    write.table(filters_used_gene, file_path, col.names = T, sep = ",", row.names = F)
    suppressWarnings(write.table(data.frame(), file_path, col.names = T, sep = ",", append = T))
    suppressWarnings(write.table(filters_used_path, file_path, col.names = T, sep = ",", append = T, row.names =F))
    suppressWarnings(write.table(data.frame(), file_path, col.names = T, sep = ",", append = T))
    suppressWarnings(write.table(collated_gene_PoDs, file_path, col.names = T, sep = ",", append = T, row.names = F))
    suppressWarnings(write.table(data.frame(), file_path, col.names = T, sep = ",", append = T))
    suppressWarnings(write.table(collated_pathway_PoDs, file_path, col.names = T, sep = ",", append = T, row.names = F))
    logger::log_info("PoD file successfully written to {file_path}.")
  
}
#' Calculate gene level PoDs from BMDExpress2 output for a group of analysis
#' 
#' Loop through each analysis group and call get_highest_lowest_conc_value, filter_gene_bmds and calculate_gene_PoD and then collate and output all gene PoDs.
#'
#' @param gene_bmd_file_path File path for the gene level bmd text file {filename.txt}, directly from bmdexpress2-cmd export --bm2-file {bm2_file_name.bm2} --analysis-group bmd --output-file-name {filename.txt}.
#' @param bmdExpress_input_dir Directory path for where the BMDExpress2 input files are. These files are read in for each analysis to get the highest concentration used to filter gene level BMDs by. This will override any value used set by highest_conc_filter.
#' @inheritParams filter_gene_bmds
#' @inheritParams calculate_gene_PoD
#' 
#' @importFrom logger log_info
#' @importFrom utils read.delim
#' @return dataframe gene level PoDs for all analysis groups in the BMD files.
run_gene_level_BMD_analysis <- function(gene_bmd_file_path,
                                        gene_csv_output_dir,
                                        bmdExpress_input_dir = NULL,
                                        highest_conc_filter = NULL,
                                        bmd_param = "BMDL",
                                        bmdl_bmdu_filter = 40,
                                        fitP_filter = 0.1
                                        ) {
  # gene_bmd_file_path
  if(!file.exists(gene_bmd_file_path)) {
    logger::log_fatal("gene_bmd_file_path '{gene_bmd_file_path}' passed to function 'calculate_PoDs_from_BMDExpress2' has no such file or directory")
    stop("gene_bmd_file_path passed to function 'calculate_PoDs_from_BMDExpress2' does not exist")
  }
  # bmdExpress_input_dir
  if (!is.null(bmdExpress_input_dir)) {
    if (!dir.exists(bmdExpress_input_dir)) {
      logger::log_fatal("bmdExpress_input_dir '{bmdExpress_input_dir}' passed to function 'calculate_PoDs_from_BMDExpress2' does not exist")
      stop("bmdExpress_input_dir passed to function 'calculate_PoDs_from_BMDExpress2' does not exist")
    }
  }
  gene_data <- as.data.frame(read.delim(gene_bmd_file_path, sep = "\t", header = F, stringsAsFactors = F, skip = 1)) # header= T throws errors due to NA columns
  # set column names
  gene_data_cols <- as.character(gene_data[1, ])
  colnames(gene_data) <- gene_data_cols
  gene_data <- gene_data[-1, ]
  gene_data <- gene_data[, !apply(is.na(gene_data), 2, all)] # remove NA columns
  
  collated_gene_PoDs <- data.frame(matrix(nrow = 0, ncol = 19))

  for (analysis in unique(gene_data$Analysis)) {
    analysis_short_name <- gsub("_williams_0.05_NOMTC_foldfilter1.5_BMD", "", analysis)
    analysis_bmd <- gene_data[gene_data$Analysis == analysis, ]
    defined_conc <- get_highest_lowest_conc_value(analysis, bmdExpress_input_dir, highest_conc_filter)
    if (nrow(analysis_bmd) > 0) {
      bmd_filtered <- filter_gene_bmds(
        bmd_data = analysis_bmd,
        highest_conc_filter = defined_conc,
        bmdl_bmdu_filter = bmdl_bmdu_filter,
        fitP_filter = fitP_filter,
        gene_csv_output_prefix = analysis_short_name,
        gene_csv_output_dir = gene_csv_output_dir
      )
      if (nrow(bmd_filtered) > 0) {
        gene_PoDs <- calculate_gene_PoD(bmd_filtered = bmd_filtered, bmd_param = bmd_param, bmdExpress_input_dir = bmdExpress_input_dir)
        collated_gene_PoDs <- rbind(collated_gene_PoDs, gene_PoDs)
      } else {
        collated_gene_PoDs[paste0(analysis, "  *(0 bmds)"),] <- c(paste0(analysis, " *(0 bmds)"),  rep(NA, 2), "",  rep(NA, 15)) # NA for when no pathway data is returned so no PoD can be calculated. 
        
      }
    }
  }
  for (i in 1:ncol(collated_gene_PoDs)) {
    if (!any(is.na(suppressWarnings(as.numeric(collated_gene_PoDs[!is.na(collated_gene_PoDs[,i]), i]))))) { # check whether columns can be converted to numeric, if so convert. (NAs appear when non numeric characters are attempted to be converted so we check for NAs to determine suitability) 
      collated_gene_PoDs[, i] <- as.numeric(collated_gene_PoDs[, i])
    }
  }
  return(collated_gene_PoDs)
  
}
#' Calculate pathway level PoDs from BMDExpress2 output for a group of analysis
#' 
#' Loop through each analysis group and call filter_pathway_bmds and calculate_pathway_PoD and then collate and output all pathway PoDs.
#'
#' @param path_bmd_file_path File path for the pathway level bmd text file {filename.txt}, directly from bmdexpress2-cmd export --bm2-file {bm2_file_name.bm2} --analysis-group categorical --output-file-name {filename.txt}.
#' @param highest_conc_filter The parameter value passed to calculate_PoDs_from_BMDExpress2. Can be numeric value for rule 1 or NULL for rules 2 and 3.
#' @inheritParams filter_pathway_bmds
#' @inheritParams calculate_pathway_PoD
#' @importFrom logger log_fatal
#' @importFrom utils read.delim
#' @return dataframe pathway level PoDs for all analysis groups in the BMD files.
#' @export
run_pathway_level_BMD_analysis <- function(path_bmd_file_path,
                                           highest_conc_filter,
                                           bmd_param = "BMDL",
                                           min_total_genes = 3,
                                           min_sig_genes = 2,
                                           fishers_p_val = 0.1,
                                           path_csv_output_dir,
                                           bmdExpress_input_dir) {
  
  # path_bmd_file_path
  if(!file.exists(path_bmd_file_path)) {
    logger::log_fatal("path_bmd_file_path '{path_bmd_file_path}' passed to function 'calculate_PoDs_from_BMDExpress2' has no such file or directory")
    stop("path_bmd_file_path passed to function 'calculate_PoDs_from_BMDExpress2' does not exist")
  }
  
  path_data <- as.data.frame(read.delim(path_bmd_file_path, sep = "\t", header = F, stringsAsFactors = F, skip = 1)) # header= T throws errors due to NA columns
  # set column names
  path_data_cols <- as.character(path_data[1, ])
  colnames(path_data) <- path_data_cols
  path_data <- path_data[-1, ]
  path_data <- path_data[, !apply(is.na(path_data), 2, all)] # remove na columns
  
  collated_pathway_PoDs <- data.frame(matrix(nrow = 0, ncol = 18), stringsAsFactors = F)
  colnames(collated_pathway_PoDs) <- c("Analysis",
                                       "Avg of 20 lowest pathway BMDLs",
                                       "The lowest pathway BMDL",
                                       "Avg of 20 pathway BMDLs with lowest 2-tail fisher P values",
                                       "-",
                                       "Lowest pathway BMDLafter lowest dose",
                                       "Lowest pathway after lowest dose",
                                       "5%",
                                       "10%",
                                       "20%",
                                       "30%",
                                       "40%",
                                       "50%",
                                       "60%",
                                       "70%",
                                       "80%",
                                       "90%",
                                       "95%")
  for (analysis in unique(path_data$Analysis)) {
    analysis_short_name <- gsub("_williams_0.05_NOMTC_foldfilter1.5_BMD_WT_Human_REACTOME_true_true_pval0.1_ratio40", "", analysis)
    analysis_bmd <- path_data[path_data$Analysis == analysis, ]
    analysis_pathway_filtered <- filter_pathway_bmds(
      path_data = analysis_bmd,
      min_total_genes = min_total_genes,
      min_sig_genes = min_sig_genes,
      fishers_p_val = fishers_p_val,
      path_csv_output_prefix = analysis_short_name, 
      path_csv_output_dir = path_csv_output_dir
    )
    if (nrow(analysis_pathway_filtered) > 0) {
      path_PoDs <- calculate_pathway_PoD(pathway_filtered = analysis_pathway_filtered, bmd_param = bmd_param, bmdExpress_input_dir = bmdExpress_input_dir)
      collated_pathway_PoDs <- rbind(collated_pathway_PoDs, path_PoDs)
    } else {
      collated_pathway_PoDs[paste0(analysis, " *(0 bmds)"),] <- c(paste0(analysis, " *(0 bmds)"), rep(NA, 3), "", rep(NA,13)) # NA for when no pathway data is returned so no PoD can be calculated. 
    }
  }
  for (i in 1:ncol(collated_pathway_PoDs)) {
    if (!any(is.na(suppressWarnings(as.numeric(collated_pathway_PoDs[!is.na(collated_pathway_PoDs[,i]), i]))))) { # check whether columns can be converted to numeric, if so convert. (NAs appear when non numeric characters are attempted to be converted so we check for NAs to determine suitability) 
      collated_pathway_PoDs[, i] <- as.numeric(collated_pathway_PoDs[, i])
    }
  }
  return(collated_pathway_PoDs)
}

#' Calculate gene and pathway level PoDs from BMDExpress2 output
#'
#' Point of departures are calculated from BMDExpress2's gene results for each analysis group (BMDExpress2 input file). Gene level BMDs can be extracted with bmdexpress2-cmd export --bm2-file {bm2_file_name.bm2} --analysis-group bmd --output-file-name {filename.txt}.
#' Pathway level BMDs can be extracted with bmdexpress2-cmd export --bm2-file {bm2_file_name.bm2} --analysis-group categorical --output-file-name {filename.txt}.
#' PoDs are calculated from a subset of the genes/pathways which are deemed as significant based on given filter parameters and the chemical PoD are defined using different methods:
#'
#' Gene Level PoDs
#' 1. Average of BMD/BMDLs of 20 genes with the largest fold change.
#' 2. Average of the 25th and 75th percentile of BMD/BMDLs
#'
#' Pathway Level PoDs
#' 3. Average of BMD/BMDLs of 20 pathways with the lowest Fishers 2-tail p values.
#' 4. Average of BMD/BMDLs of 20 pathways with the lowest BMD/BMDLs mean values.
#' 5. Taking the lowest pathways BMD/BMDLs mean value.
#'
#' This function loops through each analysis group (cell line- chemical) in the BMDExpress2 text files and run functions 'filter_gene_bmds', 'calculate_gene_PoD', 'filter_pathway_bmds' and 'calculate_pathway_PoD' and collate results.
#' Each BMDExpress analysis will have filtered bmd and pathway csv exported. Output files will have BMDExpress2 analysis name with appended with "_bmd_filtered" and "_pathway_bmd_filtered".Where BMDExpress2 bmd analysis names have "_williams_0.05_NOMTC_foldfilter1.5_BMD" and pathway analysis names have "_williams_0.05_NOMTC_foldfilter1.5_BMD_WT_Human_REACTOME_true_true_pval0.1_ratio40" in them, this text will be removed from the file output name to reduce characters, otherwise the whole BMDExpress analysis name will be seen. 
#' 
#' PoD's shall be written to a csv file at {pod_csv_output_prefix} with file name {pod_csv_output_dir}_PoDs.csv
#'
#' @param gene_bmd_file_path File path for the gene level bmd text file {filename.txt}, directly from bmdexpress2-cmd export --bm2-file {bm2_file_name.bm2} --analysis-group bmd --output-file-name {filename.txt}.
#' @param path_bmd_file_path File path for the pathway level bmd text file {filename.txt}, directly from bmdexpress2-cmd export --bm2-file {bm2_file_name.bm2} --analysis-group categorical --output-file-name {filename.txt}.
#' @param bmdExpress_input_dir Directory path for where the BMDExpress2 input files are. These files are read in for each analysis to get the highest concentration used to filter gene level BMDs by. This will override any value used set by highest_conc_filter.
#' @param bmd_param "BMD" or "BMDL" indicating which value to base calculations around.
#' @inheritParams filter_gene_bmds
#' @inheritParams filter_pathway_bmds
#' @inheritParams write_PoDs_to_file
#' 
#' @return list of 2 dataframes, gene_PoDs =  gene level PoDs and pathway_PoDs = pathway level PoDs for all analysis groups in the BMD files.
#' @export
calculate_PoDs_from_BMDExpress2 <- function(gene_bmd_file_path,
                                            path_bmd_file_path,
                                            bmdExpress_input_dir = NULL,
                                            highest_conc_filter = NULL,
                                            bmd_param = "BMDL",
                                            bmdl_bmdu_filter = 40,
                                            fitP_filter = 0.1,
                                            min_total_genes = 3,
                                            min_sig_genes = 2,
                                            fishers_p_val = 0.1,
                                            gene_csv_output_dir = ".",
                                            path_csv_output_dir = ".",
                                            pod_csv_output_prefix = "",
                                            pod_csv_output_dir = ".") {
 
  collated_gene_PoDs <- run_gene_level_BMD_analysis(gene_bmd_file_path = gene_bmd_file_path,
                                                    bmdExpress_input_dir = bmdExpress_input_dir,
                                                    highest_conc_filter = highest_conc_filter,
                                                    bmd_param = bmd_param,
                                                    bmdl_bmdu_filter = bmdl_bmdu_filter,
                                                    fitP_filter = fitP_filter,
                                                    gene_csv_output_dir = gene_csv_output_dir)
  collated_pathway_PoDs <- run_pathway_level_BMD_analysis(path_bmd_file_path = path_bmd_file_path,
                                                          bmd_param = bmd_param,
                                                          highest_conc_filter = highest_conc_filter,
                                                          min_total_genes = min_total_genes,
                                                          min_sig_genes = min_sig_genes,
                                                          fishers_p_val = fishers_p_val,
                                                          path_csv_output_dir = path_csv_output_dir,
                                                          bmdExpress_input_dir = bmdExpress_input_dir)
  if (!is.null(pod_csv_output_dir)) {
    if (is.null(bmdExpress_input_dir) & is.null(highest_conc_filter)) {
      conc <- "None"
    } else if (!is.null(bmdExpress_input_dir) & is.null(highest_conc_filter)) {
      conc <- "Highest tested concentration"
    } else {
      conc <- highest_conc_filter
    }
  
  write_PoDs_to_file(pod_csv_output_prefix,
                     conc,
                     highest_conc_filter,
                     bmd_param, 
                     bmdl_bmdu_filter,
                     fitP_filter,
                     min_total_genes,
                     min_sig_genes,
                     fishers_p_val,
                     pod_csv_output_dir,
                     collated_gene_PoDs,
                     collated_pathway_PoDs)
  
  }
  
  return(list(gene_PoDs = collated_gene_PoDs,
              pathway_PoDs = collated_pathway_PoDs))
}

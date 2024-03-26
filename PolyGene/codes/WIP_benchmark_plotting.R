library(tidyverse)
library(pROC)
library(here)
library(data.table)

#load the benchmarking function
source(str_c(here(), "/PolyGene/codes/benchmark_PolyNet_scores.R"))

#read in the Mendelian ground truths file
gt <- str_c(here(),"/data/freund_2018_monogenic_to_complex_ground_truth.csv")

#make the gene rankings for the best conditions
best_rwr_filename = str_c(here(),"/data/BestRWR_Rankings",
                          "RWR_restart",
                          0.7,
                          "_softmax",
                          TRUE,
                          "_nSeeds",
                          str_replace_na(500, "NA"),
                          "_HubGeneAdjust",
                          FALSE,
                          ".csv")
pops_filename = str_c(here(),"/data/POPs_MAGMA_0kb.txt")
old_rwr_filename = str_c(here(),"/data/original_output/PPI_STRING_MAGMA.txt")
old_pn_full_filename = str_c(here(),"/data/original_output/Combo_PPI_CoexpressDB_STRING_WtdMean.txt")


###output the rankings for the best RWR conditions
#best_rwr_rankings <- run_RWR(gene_scores = fread(str_c(here(), "/data/MAGMA_v108_GENE_0_ZSTAT.txt")),
#                             restart_prob = 0.7,
#                             thresh_score = NA,
#                             softmax = TRUE,
#                             n_seed_genes = 500,
#                             adjust_hub_genes = FALSE,
#                             graph = create_string_graph(edge_threshold = 400))
#fwrite(best_rwr_rankings, best_rwr_filename)

#make a parameter file
file_info_df <- tribble(
  ~filename, ~label, ~restart_prob, ~softmax, ~n_seeds, ~adj_hub,
  pops_filename, "PoPs", NA, NA,NA,NA,
  best_rwr_filename, "Best_PolyNetSTRING", 0.7,TRUE,500,FALSE,
  old_rwr_filename, "Old_PolyNetSTRING", 0.5, FALSE,11000,FALSE,
  old_pn_full_filename, "Old_PolyNetSTRING+CoEx", 0.5, FALSE,11000,FALSE
)

#' Performing gene and phenotype intersection filtering for downstream benchmarking
#' of input gene scoring files
#'
#' This function preprocesses benchmarking data for overlap analysis by finding the intersection of
#' overlapping phenotypes and genes and filtering the input datasets to that intersection
#'
#' @param file_info_df A data frame containing the columns filename, label.
#' Filename is the path to the gene scores to be benchmarked.
#' Label is what to label this dataset (e.g. in plotting).
#' 
#' @return A list containing a processed dataframe for each row of the "file_info_df"
#' dataframe. The names of each list element will be the provided label.
#' 
#' @examples
#'file_info <- tribble(
#'  ~filename, ~label,
#'  "/data/POPs_MAGMA_0kb.txt", "PoPs",
#'  "/data/original_output/PPI_STRING_MAGMA.txt", "Best_PolyNetSTRING",
#')
#' benchmarking_data_overlap_preprocessing(file_info)
#'
#' @import tidyverse
#' @import data.table
#' 
benchmarking_data_overlap_preprocessing <- function(file_info_df) {
  library(tidyverse)
  library(data.table)
  
  # Step 1: Read all files and find the intersection of overlapping_phenos and overlapping_genes
  overlapping_phenos <- character(0)
  overlapping_genes <- character(0)
  
  for (i in 1:nrow(file_info_df)) {
    filename <- file_info_df$filename[i]
    file_data <- fread(filename)
    
    if ("V1" %in% colnames(file_data)) {
      file_data <- file_data %>% rename("gene_symbol" = "V1")
    }
    
    if (i == 1) {
      overlapping_phenos <- colnames(file_data)
      overlapping_genes <- file_data$gene_symbol
    } else {
      overlapping_phenos <- intersect(overlapping_phenos, colnames(file_data))
      overlapping_genes <- intersect(overlapping_genes, file_data$gene_symbol)
    }
  }
  
  # Step 2: Loop through each file again and filter based on overlapping genes and overlapping phenos
  result_list <- list()
  
  for (i in 1:nrow(file_info_df)) {
    filename <- file_info_df$filename[i]
    file_label <- file_info_df$label[i]
    
    file_data <- fread(filename)
    
    if ("V1" %in% colnames(file_data)) {
      file_data <- file_data %>% rename("gene_symbol" = "V1")
    }
    
    # Filter file_data to include only overlapping genes and phenotypes
    file_data_filtered <- file_data %>%
      select(gene_symbol, all_of(overlapping_phenos)) %>%
      filter(gene_symbol %in% overlapping_genes)
    
    # Store the processed dataframe in the result list with label as the list element name
    result_list[[file_label]] <- file_data_filtered
  }
  
  return(result_list)
}


#' Perform benchmarking for lists of data
#'
#' This function performs benchmarking for lists of data,
#' computing AUC and Fisher's Exact Enrichment for genes above
#' the specified threshold
#'
#' @param data_list A list of datasets for benchmarking. Typically the output
#' of benchmarking_data_overlap_preprocessing
#' @param gt_file Ground truth data that specifies true positives and true negatives.
#' str_c(here(),"/data/freund_2018_monogenic_to_complex_ground_truth.csv") by default.
#' @param grouping_column How to group the phenotypes. "mendelian_disease_group"
#' by default. Could also be "complex_trait" 
#' @param percentile_threshold Percentile threshold for enrichment calculation. 0.95 by default.
#' @param label Label for the dataset being processed.
#' @param restart_prob Optional: Restart probability of the input datasets. NA by default.
#' @param softmax Optional: Softmax TRUE/FALSE for the input datasets. NA by default.
#' @param n_seeds Optional: Number of seeds for benchmarking. NA by default.
#' @param adj_hub Optional: Hub Adjustment TRUE/FALSE. NA by default.
#' @param thresholds Threshold values for sensitivity and specificity calculations.
#'
#' @return A list containing benchmarking results for each dataset.
#'
#' @examples
#' benchmarking_for_lists(data_list = pp_out, gt_file = gt, grouping_column = "mendelian_disease_group", 
#'                       percentile_threshold = 0.95, label = NA, restart_prob = NA, softmax = NA, 
#'                       n_seeds = NA, adj_hub = NA, thresholds = c(0, 0.01, 0.05, 0.2, 0.3, 0.4, 0.5, 
#'                       0.6, 0.7, 0.8, 0.9, 0.95, 0.99, 1))
#'
#' @importFrom dplyr mutate
#'
benchmarking_for_lists <- function(data_list = pp_out,
                         gt_file = str_c(here(),"/data/freund_2018_monogenic_to_complex_ground_truth.csv"),
                         grouping_column = "mendelian_disease_group",
                         percentile_threshold = 0.95,
                         label = NA,
                         restart_prob = NA,
                         softmax = NA,
                         n_seeds = NA,
                         adj_hub = NA,
                         thresholds = c(0, 0.01, 0.05,
                                        0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,
                                        0.95, 0.99, 1)) {
  
  result_list <- list() 
  for (i in seq_along(data_list)) {
    
    #get the label
    if (!is.na(label)) {
      current_label <- label
    } else {
      current_label <- names(data_list)[[i]]
    }
    data <- data_list[[i]]
    
    bm_results <- compute_sensitivity_specificity(data = data,
                                                  ground_truth_df = gt_file,
                                                  thresholds = thresholds)
    
    bm_summary <- summarize_benchmarks(restart_prob = restart_prob, softmax = softmax, n_seeds = n_seeds, adj_hub = adj_hub,
                                       benchmark_filename = bm_results,
                                       grouping_column = grouping_column,
                                       percentile_threshold = percentile_threshold) %>%
      mutate(label = current_label)
    
    result_list[[current_label]] <-  bm_summary
  }
  return(result_list)
}

#' Plot comparison of benchmarking summary (i.e. AUC and Enrichment OR)
#'
#' This function plots the comparison of benchmarking performance metrics 
#' (AUC or Enrichment OR) across a list of input dataframes
#'
#' @param bm_summary_list A list of benchmarking summary data frames that contain
#' the fisher_or and/or auc columns and their CIs.
#' @param y_axis_column The column name for the y-axis. "auc" by defualt.
#' @param x_axis_column The group column to plot on the x-axis.
#' "mendelian_disease_group by default"
#'
#' @return A ggplot object representing the comparison plot.
#'
#' @examples
#' plot_comparison(bm_summary_list = test_bm_summary_list, y_axis_column = "auc", 
#'                 x_axis_column = "mendelian_disease_group")
#'
#' @import ggplot2
#' @importFrom dplyr bind_rows
#' @importFrom ggplot2 geom_point geom_errorbar labs scale_color_discrete theme_bw theme
#'
plot_comparison <- function(bm_summary_list = test_bm_summary_list,
                            y_axis_column = "auc",
                            x_axis_column = "mendelian_disease_group") {
  
  combined_summary <- do.call(bind_rows, bm_summary_list)
  
  # Determine CI column names based on the selected y-axis column
  ci_lower_col <- paste0(y_axis_column, "_lower_ci")
  ci_upper_col <- paste0(y_axis_column, "_upper_ci")
  
  ggplot(combined_summary,
         aes_string(x = x_axis_column, y = y_axis_column, color = "label")) +
    geom_point(position = position_dodge(width = 0.4)) +
    geom_errorbar(aes_string(ymin = ci_lower_col, ymax = ci_upper_col),
                  position = position_dodge(width = 0.4),
                  width = 0.3) +
    labs(title = paste("Compare", y_axis_column), "Performance Across Methods",
         x = "Mendelian Phenotype Group") +
    scale_color_discrete(name = "Algorithm") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 75, hjust = 1))
}

#' Perform benchmarking, visualization, and combine results
#'
#' This function performs benchmarking of multiple gene ranking dataframes against
#' a ground truth dataset, plots a comparison of each dataset with respect to
#' AUC and Fisher's Exact Enrichment above some threshold, and combines the results
#' into a single dataframe based on an input dataframe of filenames and labels.
#'
#' @param file_info_df A dataframe containing information about the files to be processed
#' by benchmarking_data_overlap_preprocessing(). The dataframe should include filenames and corresponding labels.
#' @param gt_file Ground truth data file. str_c(here(),"/data/freund_2018_monogenic_to_complex_ground_truth.csv") by default.
#' @param threshold Threshold for percentile calculation. 0.95 by default.
#' @param thresholds Threshold values for benchmarking calculations to report. c(0.95) by default.
#' @param y_axis_columns c("fisher_or", "auc") by default.
#' @param x_axis_column "mendelian_disease_group"
#'
#' @return A list containing the following elements:
#' \item{fisher_or_plot}{A ggplot object representing the comparison plot for Fisher OR.}
#' \item{auc_plot}{A ggplot object representing the comparison plot for AUC.}
#' \item{combined_df}{A combined dataframe containing benchmarking results}
#'
#' @examples
#' benchmark_plot_and_results(file_info_df = file_info_df,
#' gt_file = str_c(here(),"/data/freund_2018_monogenic_to_complex_ground_truth.csv"),
#' threshold = 0.95,
#' thresholds = c(0.95),
#' y_axis_columns = c("fisher_or", "auc"),
#' x_axis_column = "mendelian_disease_group")
#'
#' @import dplyr
#' @import ggplot2
benchmark_plot_and_results <- function(file_info_df,
                                       gt_file = str_c(here(),"/data/freund_2018_monogenic_to_complex_ground_truth.csv"),
                                       threshold = 0.95,
                                       thresholds = c(threshold),
                                       y_axis_columns = c("fisher_or", "auc"),
                                       x_axis_column = "mendelian_disease_group") {
    
    # Step 1: Data Preprocessing
    pp_out <- benchmarking_data_overlap_preprocessing(file_info_df)
    
    # Step 2: Benchmarking
    bm_summary_list <- benchmarking_for_lists(data_list = pp_out,
                                    gt_file = gt_file,
                                    grouping_column = x_axis_column,
                                    percentile_threshold = threshold)
    
    # Step 3: Visualization
    plots <- list()
    for(metric in y_axis_columns) {
      plot_comparison <- plot_comparison(bm_summary_list = bm_summary_list,
                                         y_axis_column = metric,
                                         x_axis_column = "mendelian_disease_group")
      plots[[metric]] <- plot_comparison
    }
    
    # Step 4: Combine and return
    
    combined_df <- do.call(bind_rows, bm_summary_list) %>%
      select(-restart_prob, -softmax, -n_seeds, -adj_hub) %>%
      left_join(file_info_df, by = "label") %>%
      select(label, everything())
    
    return(list(fisher_or_plot = plots[["fisher_or"]],
                auc_plot = plots[["auc"]],
                combined_df = combined_df))
}
  

# Example usage
output_list <- benchmark_plot_and_results(file_info_df)
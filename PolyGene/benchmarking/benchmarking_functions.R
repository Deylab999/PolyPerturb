library(here)
library(tidyverse)
library(data.table)

source(str_c(here(), "/PolyGene/PolyNet/create_STRING_graph_object.R"))
source(str_c(here(), "/PolyGene/PolyNet/run_RWR.R"))
source(str_c(here(), "/PolyGene/PolyNet/run_pagerank.R"))

#' Compute sensitivity, specificity, precision, recall, Fisher's exact test p-value, and odds ratio.
#'
#' This function computes performance metrics including sensitivity, specificity, precision, recall,
#' Fisher's exact test p-value, and odds ratio for each phenotype at different percentile thresholds.
#'
#' @param data The data containing "gene_symbol" column, then all phenotypes and scores.
#' Metrics will be based on percentiles and therefore the input metric does not matter. 
#' @param ground_truth The ground truth data containing "gene_symbol", "phenotype",
#' and "ground_truth" labels. ground_truth==1 is a true positive relationship
#' between the gene and the phenotype and ground_truth==0 is a true negative relationship.
#' @param thresholds A vector of thresholds for defining positive predictions.
#' These are the percentile of PolyNet score (e.g. threshold=0.9 means that the top 10% of
#' genes by PolyNet score are evaluated)
#'
#' @return A lon-format data frame containing performance metrics for each phenotype at different thresholds.
#' This dataframe contains everything required to reconstruct the confusion matrix, 
#' sensitivity, specificity, precision, and recall at each threshold,
#' as well as AUC with DeLong CIs (auc_lower_ci, auc_upper_ci) for each phenotype.
#'
#' @import tidyverse
#' @import pROC
#'
compute_sensitivity_specificity <- function(data,
                                            ground_truth_df,
                                            thresholds = c(0, 0.01, 0.05,
                                                           0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,
                                                           0.95, 0.99, 1)) {
  library(tidyverse)
  library(data.table)
  library(pROC)
  
  # Check if 'data' is a filepath or a dataframe
  if (is.character(data)) {
    data <- fread(data)
  }
  
  if (is.character(ground_truth_df)) {
    ground_truth_df <- fread(ground_truth_df)
  }
  
  if("V1" %in% colnames(data)){
    data <- data %>% rename("gene_symbol" = "V1")
  }
  
  # Calculate percentile ranks for each phenotype
  percentile_ranks <- data %>%
    mutate(across(-gene_symbol, ~na_if(., NA) %>%
                    replace_na(mean(.)) %>%
                    percent_rank())) %>%
    pivot_longer(-gene_symbol, names_to = "phenotype", values_to = "percentile_rank")
  
  # Merge with true positives and true negatives
  merged_data <- ground_truth_df %>%
    inner_join(percentile_ranks, by = c("gene_symbol", "phenotype")) %>%
    drop_na()
  
  # Helper function to compute sensitivity, specificity, precision, and recall at many thresholds
  compute_metrics <- function(data_merged_with_gt = merged_data, threshold = 0.9) {
    contingency_tables <- data_merged_with_gt %>%
      group_by(phenotype) %>%
      mutate(predicted_label = ifelse(percentile_rank >= threshold, 1, 0)) %>%
      summarise(
        n_ground_truth_positives = sum(ground_truth == 1),
        n_ground_truth_negatives = sum(ground_truth == 0),
        n_predicted_positives = sum(predicted_label == 1),
        n_predicted_negatives = sum(predicted_label == 0),
        n_true_positives = sum(ground_truth == 1 & predicted_label == 1),
        n_true_negatives = sum(ground_truth == 0 & predicted_label == 0),
        n_false_positives = sum(ground_truth == 0 & predicted_label == 1),
        n_false_negatives = sum(ground_truth == 1 & predicted_label == 0),
        sensitivity = sum(ground_truth == 1 & predicted_label == 1) / sum(ground_truth == 1),
        specificity = sum(ground_truth == 0 & predicted_label == 0) / sum(ground_truth == 0),
        precision = sum(ground_truth == 1 & predicted_label == 1) / sum(predicted_label == 1),
        recall = sensitivity,
        n_genes_above_thresh = n_predicted_positives,
        n_genes_total = sum(predicted_label == 0 | predicted_label == 1)
      ) %>%
      ungroup() %>%
      mutate(
        # Compute Fisher's exact test p-value and odds ratio
        fisher_p_value = mapply(function(tp, fp, tn, fn) {
          fisher.test(matrix(c(tp, fp, fn, tn), nrow = 2))$p.value
        }, n_true_positives, n_false_positives, n_true_negatives, n_false_negatives),
        odds_ratio = mapply(function(tp, fp, tn, fn) {
          (tp / (tp + fn)) / (fp / (fp + tn))
        }, n_true_positives, n_false_positives, n_true_negatives, n_false_negatives)
      ) %>%
      mutate(threshold := threshold)
    
    # Compute AUC for each phenotype
    auc_list <- list()
    for (p in unique(data_merged_with_gt$phenotype)) {
      subset_df <- data_merged_with_gt %>%
        filter(phenotype == p)
      roc_obj <- invisible(roc(
        data = subset_df,
        response = "ground_truth",
        predictor = "percentile_rank",
        ci = TRUE,
        direction = "<"
      ))
      roc_auc <- as.numeric(auc(roc_obj))
      auc_lower_ci <- roc_obj$ci[1]
      auc_upper_ci <- roc_obj$ci[3]
      
      # Estimate standard error
      z_value <- qnorm(0.975)  # for 95% confidence interval
      auc_se <- (auc_upper_ci - auc_lower_ci) / (2 * z_value)
      
      auc_list[[p]] <- c(roc_auc, auc_lower_ci, auc_upper_ci, auc_se)
    }
    
    auc_df <- data.frame(
      phenotype = names(auc_list),
      auc = unlist(lapply(auc_list, `[`, 1)),
      auc_lower_ci = unlist(lapply(auc_list, `[`, 2)),
      auc_upper_ci = unlist(lapply(auc_list, `[`, 3)),
      auc_se = unlist(lapply(auc_list, `[`, 4))
    )
    
    # Merge AUC and contingency table data
    metric_summary <- left_join(contingency_tables, auc_df, by = "phenotype")
    
    return(metric_summary)
  }
  
  # Compute metrics for each threshold and phenotype
  performance_metrics <- map_dfr(thresholds,
                                 ~compute_metrics(data = merged_data,
                                                  threshold = .x))
  performance_metrics <- merged_data %>%
    select(-gene_symbol, -percentile_rank, -ground_truth) %>%
    distinct() %>%
    left_join(performance_metrics, by="phenotype", relationship = "many-to-many")
    
  return(performance_metrics)
}

#' Run benchmarking starting from RWR with specified parameters
#' 
#' @param restart_prob The restart probability to use
#' @param softmax TRUE/FALSE. If TRUE, use softmax normalization of seed genes.
#' @param outfile_name Output filepath for benchmarking data
#' @param ground_truth A data frame containing ground truth information.
#' Mendelian ground truths are used by default (i.e. from freund_2018_monogenic_to_complex_ground_truth.csv)
#' @param graph A graph object representing the network. Will create a STRING graph wtih edge_threshold of 400
#' if an igraph object is not provided
#' @return NULL. The function writes benchmark results to the specified CSV file in outfile_name.


# Define a function to run compute_sensitivity_specificity()
# across a grid of parameters
run_benchmarking_with_parameters <- function(restart_prob = 0.8,
                                             softmax = TRUE,
                                             n_seeds = 100,
                                             adj_hub = FALSE,
                                             outfile_name = str_c(here(),
                                                                  "/PolyGene/benchmarking/PolyNet/benchmark.csv"),
                                             ground_truth = fread(str_c(here(),
                                                                        "/data/freund_2018_monogenic_to_complex_ground_truth.csv")),
                                             graph = create_string_graph(edge_threshold = 400)) {
  
  #run RWR with the set parameters
  graph_gene_rankings <- run_RWR(restart_prob = restart_prob,
                                 softmax = softmax,
                                 n_seed_genes = n_seeds,
                                 adjust_hub_genes = adj_hub,
                                 graph = graph,
                                 )
  #calculate the benchmarks
  result <- compute_sensitivity_specificity(data = graph_gene_rankings,
                                            ground_truth_df = gt)
  #write out the benchmarks to a .csv file
  fwrite(result, file = outfile_name)
}


#' Summarizes benchmark data from a given file.
#'
#' @param restart_prob Numeric, restart probability to use in RWR.
#' 0=hub genes prioritized, 1=seed genes prioritized. Default is 0.5.
#' @param softmax Logical, indicating whether SOFTMAX normalization of the seed nodes should be used
#' prior to RWR. A way of weighting the seed genes.
#' @param n_seeds Numeric, filter to n top genes to use as seeds in RWR
#' @param adj_hub Logical, indicating whether to adjust for hub genes.
#' @param percentile_threshold Numeric, the threshold at which to evaluate the enrichment ORs.
#' Default is 0.95, corresponding to calculating enrichment ORs at genes in the top 5% of the score.
#' @param grouping_column Character, a column to summarize the AUCs and ORs by. This is done
#' byt inverse variance weighting. "mendelian_disease_group" by default, but another useful group is the "complex_trait" column.
#' @param benchmark_filename Character, path to the benchmarking file output by run_benchmarking_with_parameters()
#' @return A data frame containing summarized AUC benchmark results, both overall and split by group.
#' @import data.table
#' @import dplyr
#' @importFrom purrr where
#' @importFrom tidyr distinct
#' @importFrom tibble as_tibble
#' @importFrom dplyr across group_by mutate summarise ungroup select

summarize_benchmarks <- function(restart_prob = NA,
                                 softmax = NA,
                                 n_seeds = NA,
                                 adj_hub = NA,
                                 percentile_threshold = 0.95, #look at the top 5% of prioritized genes
                                 grouping_column = "mendelian_disease_group", 
                                 benchmark_filename = str_c(here(),
                                                            "/PolyGene/benchmarking/PolyNet/",
                                                            "Mendelian_Freund_2018_Benchmarking/",
                                                            "MAGMA_0kb/STRING_PolyNet/",
                                                            "RWR_restart0.7_softmaxTRUE_nSeeds500_HubGeneAdjustFALSE.csv")){
  bm <- benchmark_filename
  #read in the benchmarking file from disc
  if (is.character( benchmark_filename)) {
    bm <- fread(benchmark_filename)
  }
  
  #add Fisher's exact SE
  bm <- bm %>%
    mutate(log_odds_ratio = log(odds_ratio),
           log_odds_ratio_se = sqrt(1 / (n_true_positives + n_false_positives) + 
                                      1 / (n_false_positives + n_true_negatives) + 
                                      1 / (n_true_positives + n_false_negatives) + 
                                      1 / (n_false_negatives + n_true_negatives)))
  
  
  #drop the columns to reconstruct the confusion matrix to simplify the output
  cols_to_keep <- data.frame(colnames = colnames(bm)) %>%
    filter(!str_detect(colnames, "^n_") |
             (colnames %in% c("n_genes_above_thresh",
                              "n_true_positives",
                              "n_genes_total"))) %>%
    pull(colnames)
  
  #grouped summary
  grouped_summary <- bm %>%
    filter(threshold == percentile_threshold) %>%
    distinct() %>%
    group_by(!!sym(grouping_column), threshold) %>%
    summarize(auc = weighted.mean(auc, 1 / auc_se^2),
              auc_se = sqrt(1 / sum(1 / auc_se^2)),
              fisher_or = exp(weighted.mean(log_odds_ratio, 1 / log_odds_ratio_se^2)),
              fisher_se = sqrt(1 / sum(1 / log_odds_ratio_se^2))) %>%
    ungroup() %>%
    transmute(
      !!sym(grouping_column),
      threshold,
      restart_prob = restart_prob,
      softmax = softmax,
      n_seeds = n_seeds,
      adj_hub = adj_hub,
      fisher_or,
      fisher_se,
      fisher_or_lower_ci = exp(log(fisher_or) - 1.96 * fisher_se),
      fisher_or_upper_ci = exp(log(fisher_or) + 1.96 * fisher_se),
      fisher_p_value = 2 * pnorm(-abs(log(fisher_or) / fisher_se)),
      auc,
      auc_se,
      auc_lower_ci = auc - 1.96 * auc_se,
      auc_upper_ci = auc + 1.96 * auc_se
    ) %>%
    arrange(desc(fisher_or))
  
  #compute the overall summary for each phenotype
  # at the selected percentile threshold
  overall_summary <- grouped_summary %>%
    group_by(threshold) %>%
    summarize(auc = weighted.mean(auc, 1 / auc_se^2),
              auc_se = sqrt(1 / sum(1 / auc_se^2)),
              fisher_or = exp(weighted.mean(log(fisher_or), 1 / fisher_se^2)),
              fisher_se = sqrt(1 / sum(1 / fisher_se^2))) %>%
    ungroup() %>%
    transmute(
      !!grouping_column := "Overall",
      threshold,
      restart_prob = restart_prob,
      softmax = softmax,
      n_seeds = n_seeds,
      adj_hub = adj_hub,
      fisher_or,
      fisher_se,
      fisher_or_lower_ci = exp(log(fisher_or) - 1.96 * fisher_se),
      fisher_or_upper_ci = exp(log(fisher_or) + 1.96 * fisher_se),
      fisher_p_value = 2 * pnorm(-abs(log(fisher_or) / fisher_se)),
      auc,
      auc_se,
      auc_lower_ci = auc - 1.96 * auc_se,
      auc_upper_ci = auc + 1.96 * auc_se
    ) %>%
    arrange(desc(fisher_or))
    
  #combine
  full_summary <- bind_rows(overall_summary, grouped_summary)
  return(full_summary)
}
  
#' Generate AUC vs Parameter Plots
#'
#' This function generates plots of AUC (Area Under the Curve) against specified parameters for different phenotypes.
#'
#' @param data A dataframe containing the AUC values along with the parameters and phenotype information.
#' @param params A character vector specifying the names of the parameters to be plotted.
#' @param phenotype_column The name of the column in \code{data} containing the phenotype information.
#'
#' @return A list containing the generated ggplot2 plots, with each plot corresponding to one of the specified parameters.
#'
#' @details This function iterates over each parameter specified in \code{params} and generates a separate plot for each one.
#' The plots show the AUC values against the parameter values for different phenotypes.
#'
#' @examples
#' \dontrun{
#' # Example usage
#' plots <- generate_auc_param_plots(data = ct_out_summary,
#'                                   params = c("restart_prob", "softmax", "n_seeds", "adj_hub"),
#'                                   phenotype_column = "complex_trait")
#'                                   }
#'
#' @import ggplot2
#' @import dplyr
#' @importFrom stringr str_c
#' @importFrom rlang sym syms !!!

generate_auc_param_plots <- function(data = ct_out_summary,
                                     params = c("restart_prob", "softmax", "n_seeds", "adj_hub"),
                                     phenotype_column = "complex_trait") {
  # Convert parameters to factors for better visualization
  for(param in params) {
    data[[param]] <- as.factor(data[[param]])
  }
  
  # Initialize list to store plots
  plots <- list()
  
  # Iterate over each parameter and generate plots
  for(param in params) {
    # Construct other_params by combining all parameters except the current one
    other_params <- params[params != param]
    tmp_plot_df <- data %>% mutate(other_params = str_c(!!sym(phenotype_column),
                                                        !!!syms(other_params)))
    
    # Create plot
    plot <- ggplot(tmp_plot_df,
                   aes_string(x = param, y = "auc", color = phenotype_column)) +
      geom_point() +
      geom_line(aes(group = other_params)) +
      labs(title = paste("AUC vs", param, "for Different Phenotypes"),
           x = param, y = "AUC") +
      theme_bw()
    
    # Store plot in the list
    plots[[param]] <- plot
  }
  
  # Return the list of plots
  return(plots)
}

#' Generate AUC vs Parameter Plots
#'
#' This function generates plots of AUC (Area Under the Curve) against specified parameters for different phenotypes.
#'
#' @param data A dataframe containing the AUC values along with the parameters and phenotype information.
#' @param params A character vector specifying the names of the parameters to be plotted.
#' @param phenotype_column The name of the column in \code{data} containing the phenotype information.
#' @param  y_variable The name of the metric column name you'd like plottet on the Y-axis. "auc" is default, but "fisher_or" can
#' also be useful.
#' @return A list containing the generated ggplot2 plots, with each plot corresponding to one of the specified parameters.
#'
#' @details This function iterates over each parameter specified in \code{params} and generates a separate plot for each one.
#' The plots show the AUC values against the parameter values for different phenotypes.
#'
#' @examples
#' \dontrun{
#' # Example usage
#' plots <- generate_auc_param_plots(data = ct_out_summary,
#'                                   params = c("restart_prob", "softmax", "n_seeds", "adj_hub"),
#'                                   phenotype_column = "complex_trait")
#'                                   }
#'
#' @import ggplot2
#' @import dplyr
#' @importFrom stringr str_c
#' @importFrom rlang sym syms !!!

generate_param_plots <- function(data,
                                 y_variable = "auc",
                                 params = c("restart_prob", "softmax", "n_seeds", "adj_hub"),
                                 phenotype_column = "complex_trait") {
  # Convert parameters to factors for better visualization
  for(param in params) {
    data[[param]] <- as.factor(data[[param]])
  }
  
  # Initialize list to store plots
  plots <- list()
  
  # Iterate over each parameter and generate plots
  for(param in params) {
    # Construct other_params by combining all parameters except the current one
    other_params <- params[params != param]
    tmp_plot_df <- data %>% mutate(other_params = str_c(!!sym(phenotype_column),
                                                        !!!syms(other_params)))
    
    # Create plot
    plot <- ggplot(data = filter(tmp_plot_df, !!sym(phenotype_column) != "Overall"),
                   aes_string(x = param, y = y_variable, color = phenotype_column)) +
      geom_point() +
      geom_line(aes(group = other_params)) +
      labs(title = paste(y_variable, "vs", param, "for Different Phenotypes"),
           x = param, y = y_variable) +
      theme_bw()
    
    # Add thicker black line if complex_trait is Overall
    if ("Overall" %in% unique(tmp_plot_df[[phenotype_column]])) {
      plot <- plot +
        geom_point(data = filter(tmp_plot_df, !!sym(phenotype_column) == "Overall"),
                   color = "black", size = 4) +
        geom_line(data = filter(tmp_plot_df, !!sym(phenotype_column) == "Overall"),
                  aes(group = other_params), color = "black", size = 1.5)
    }
    
    # Store plot in the list
    plots[[param]] <- plot
  }
  
  # Return the list of plots
  return(plots)
}

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
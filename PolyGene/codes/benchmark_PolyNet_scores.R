library(here)
library(tidyverse)
library(data.table)

source(str_c(here(), "/PolyGene/codes/create_string_graph_object.R"))
source(str_c(here(), "/PolyGene/codes/ppi_string_RWR.R"))
source(str_c(here(), "/PolyGene/codes/run_pagerank.R"))

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
  library(pROC)
  
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
    
    #compute AUC for each phenotype
    auc_list <- list()
    for (p in unique(data_merged_with_gt$phenotype)) {
      
      subset_df <- data_merged_with_gt %>%
        filter(phenotype==p)
      
      roc_obj <- invisible(roc(
        data = subset_df,
        response = "ground_truth",
        predictor = "percentile_rank",
        ci=TRUE, 
        direction = "<"
      ))
      
      roc_auc <- as.numeric(auc(roc_obj))
      auc_lower_ci <- roc_obj$ci[1]
      auc_upper_ci <- roc_obj$ci[3]
      
      auc_list[[p]] <- c(roc_auc, auc_lower_ci, auc_upper_ci)
    }
    
    auc_df <- data.frame(
      phenotype = names(auc_list),
      auc = unlist(lapply(auc_list, `[`, 1)),
      auc_lower_ci = unlist(lapply(auc_list, `[`, 2)),
      auc_upper_ci = unlist(lapply(auc_list, `[`, 3))
    )
    
    #Create a dataframe to re-generate confusion matrix and plot ROC
    metric_summary <- data_merged_with_gt %>%
      group_by(phenotype) %>%
      mutate(predicted_label = ifelse(percentile_rank >= threshold, 1, 0)) %>%
      summarise(
        n_ground_truth_positives = sum(ground_truth == 1),
        TP = sum(ground_truth == 1 & predicted_label == 1),
        TN = sum(ground_truth == 0 & predicted_label == 0),
        FP = sum(ground_truth == 0 & predicted_label == 1),
        FN = sum(ground_truth == 1 & predicted_label == 0),
        sensitivity = sum(ground_truth == 1 & predicted_label == 1) / sum(ground_truth == 1),
        specificity = sum(ground_truth == 0 & predicted_label == 0) / sum(ground_truth == 0),
        precision = sum(ground_truth == 1 & predicted_label == 1) / sum(predicted_label == 1),
        recall = sensitivity
      ) %>%
      ungroup() %>%
      mutate(threshold = threshold) %>%
      left_join(auc_df, by="phenotype")
    
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

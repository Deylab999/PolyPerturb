library(here)
library(tidyverse)
library(data.table)

source(str_c(here(), "/PolyGene/codes/create_string_graph_object.R"))
source(str_c(here(), "/PolyGene/codes/ppi_string_RWR.R"))
source(str_c(here(), "/PolyGene/codes/run_pagerank.R"))

#create some data
graph <- create_string_graph(edge_threshold = 400)
test_rwr <- run_RWR(graph=graph, softmax = TRUE, restart_prob = 0.5)

#read in the ground truth data
gt <- fread(str_c(here(),"/data/freund_2018_monogenic_to_complex_ground_truth.csv")) %>%
  select(-mendelian_disease_group)

#' Compute sensitivity, specificity, precision, recall, Fisher's exact test p-value, and odds ratio.
#'
#' This function computes performance metrics including sensitivity, specificity, precision, recall,
#' Fisher's exact test p-value, and odds ratio for each phenotype at different percentile thresholds.
#'
#' @param data The data containing "gene_symbol" column, then all phenotypes and scores.
#' Metrics will be based on percentiles and therefore the input metric does not matter. 
#' @param ground_truth The ground truth data containing "gene_symbol", "phenotype",
#' and "ground_truth" labels. ground_truth== 1 is a true positive relationship
#' and ground_truth==0 is a true negative relationship
#' @param thresholds A vector of thresholds for defining positive predictions.
#'
#' @return A data frame containing performance metrics for each phenotype at different thresholds.
#'
#' @examples
#' compute_sensitivity_specificity()
#'
#' @import tidyverse
#' @import pROC
#'
compute_sensitivity_specificity <- function(data = test_rwr,
                                            ground_truth_df = gt,
                                            thresholds = c(0.9, 0.95, 0.99, 0.995)) {
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
    left_join(percentile_ranks, by = c("gene_symbol", "phenotype")) %>%
    drop_na()
  
  # Helper function to compute sensitivity, specificity, precision, and recall at many thresholds
  compute_metrics <- function(data_merged_with_gt = merged_data, threshold = 0.9) {
    metric_summary <- data_merged_with_gt %>%
      group_by(phenotype) %>%
      mutate(predicted_label = ifelse(percentile_rank >= threshold, 1, 0)) %>%
      summarise(
        n_ground_truths = n(),
        TP = sum(ground_truth == 1),
        TN = sum(ground_truth == 0),
        sensitivity = sum(ground_truth == 1 & predicted_label == 1) / sum(ground_truth == 1),
        specificity = sum(ground_truth == 0 & predicted_label == 0) / sum(ground_truth == 0),
        precision = sum(ground_truth == 1 & predicted_label == 1) / sum(predicted_label == 1),
        recall = sensitivity
      ) %>%
      ungroup() %>%
      mutate(threshold = threshold)
    return(metric_summary)
  }
  
  # Compute metrics for each threshold and phenotype
  performance_metrics <- map_dfr(thresholds,
                                 ~compute_metrics(data = merged_data,
                                                  threshold = .x))
  return(performance_metrics)
}

test_out <- compute_sensitivity_specificity()
test_summary <- test_out %>%
  group_by(threshold) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE))
test_summary 
  


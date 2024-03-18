library(tidyverse)
library(pROC)
library(here)
library(data.table)


#load the benchmarking function
source(str_c(here(), "/PolyGene/codes/benchmark_PolyNet_scores.R"))

#read in the Mendelian ground truths file
gt <- fread(str_c(here(),"/data/freund_2018_monogenic_to_complex_ground_truth.csv"))

#create the STRING graph to use for RWR
string_graph <- create_string_graph(edge_threshold = 400)

# Define the grid of parameters
parameter_grid <- expand.grid(restart_prob = c(0.2, 0.8),
                              softmax = c(TRUE, FALSE)) %>%
  mutate(outfile_filename=str_c(here(),
                                "/PolyGene/benchmarking/PolyNet/",
                                "Mendelian_Disease_Gene_Benchmarks_",
                                "RWR_restartprob",
                                restart_prob,
                                "_softmax",
                                softmax,
                                ".csv"))


# Apply the function to each combination of parameters specified in the parameter_grid df
output_files <- pmap(parameter_grid,
                     ~run_rwr_benchmarking_with_parameters(restart_prob = ..1,
                                                           softmax = ..2,
                                                           outfile_name = ..3,
                                                           ground_truth = gt,
                                                           graph = string_graph))


# Create a function that reads the benchmarking output file from disc and then
# summarizes the mean AUC overall and by phenotype.

summarize_benchmarks <- function(restart_prob = 0.8,
                                 softmax = TRUE,
                                 grouping_column = "complex_trait", 
                                 benchmark_filename = "/home/robertg1/PolyPerturb/PolyGene/benchmarking/PolyNet/Mendelian_Disease_Gene_Benchmarks_RWR_restartprob0.8_softmaxTRUE.csv"){
  
  #read in the benchmarking file from disc
  bm <- fread(benchmark_filename)
  
  #compute the overall summary
  overall_summary <- bm %>%
    group_by(threshold) %>%
    summarise(across(where(is.numeric), mean, na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(complex_trait = "Overall",
           restart_prob = restart_prob,
           softmax = softmax) %>%
    select(complex_trait, restart_prob, softmax, auc, auc_lower_ci, auc_upper_ci) %>%
    distinct()
  
  #grouped summary
  grouped_summary <- bm %>%
    group_by(complex_trait, threshold) %>%
    summarise(across(where(is.numeric), mean, na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(restart_prob = restart_prob,
           softmax = softmax) %>%
    select(complex_trait, restart_prob, softmax, auc, auc_lower_ci, auc_upper_ci) %>%
    distinct() %>%
    arrange(desc(auc))
  
  #combine
  full_summary <- bind_rows(overall_summary, grouped_summary)
  return(full_summary )
}

# Apply the benchmarking summarizing function to get a final dataframe
out_summary <- pmap_df(parameter_grid,
                     ~summarize_benchmarks(restart_prob = ..1,
                                           softmax = ..2,
                                           benchmark_filename = ..3,
                                           grouping_column = "complex_trait"))
out_summary %>%
  filter(complex_trait == "Overall" | complex_trait == "Rheumatoid Arthritis") %>%
  arrange(complex_trait)



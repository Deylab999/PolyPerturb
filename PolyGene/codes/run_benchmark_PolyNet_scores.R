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
parameter_grid <- expand.grid(restart_prob = c(0, 0.2, 0.5, 0.8, 1),
                              softmax = c(TRUE, FALSE),
                              n_seeds = c(100, 200, 500, 1000),
                              adj_hub = c(TRUE, FALSE)) %>%
  mutate(outfile_filename=str_c(here(),
                                "/PolyGene/benchmarking/PolyNet/",
                                "Mendelian_Benchmarks_",
                                "RWR_restartprob",
                                restart_prob,
                                "_softmax",
                                softmax,
                                "_nSeeds",
                                n_seeds,
                                "_HubGeneAdjust",
                                adj_hub,
                                ".csv"))


# Apply the function to each combination of parameters specified in the parameter_grid df
pmap(parameter_grid,
     ~run_benchmarking_with_parameters(restart_prob = ..1,
                                       softmax = ..2,
                                       n_seeds = ..3,
                                       adj_hub = ..4,
                                       outfile_name = ..5,
                                       ground_truth = gt,
                                       graph = string_graph))


# Apply the benchmarking summarizing function to get a final dataframe
out_summary <- pmap_df(parameter_grid,
                     ~summarize_benchmarks(restart_prob = ..1,
                                           softmax = ..2,
                                           n_seeds = ..3,
                                           adj_hub = ..4,
                                           benchmark_filename = ..5,
                                           grouping_column = "mendelian_disease_group"))

# Plot the impact of each of the parameters on AUC
plots <- generate_auc_param_plots(data = out_summary,
                                  params = c("restart_prob", "softmax", "n_seeds", "adj_hub"),
                                  phenotype_column = "mendelian_disease_group")

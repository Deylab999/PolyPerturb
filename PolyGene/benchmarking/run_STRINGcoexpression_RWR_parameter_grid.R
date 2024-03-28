library(tidyverse)
library(pROC)
library(here)
library(data.table)

#load the benchmarking function
source(str_c(here(), "/PolyGene/benchmarking/benchmarking_functions.R"))

#read in the Mendelian ground truths file
gt <- fread(str_c(here(),"/data/freund_2018_monogenic_to_complex_ground_truth.csv"))

#create the STRING graph to use for RWR
string_graph <- create_STRING_coexp_graph(edge_threshold = 400)

# Define the grid of parameters to generate results for
parameter_grid <- expand.grid(restart_prob = c(0, 0.2, 0.5, 0.7, 0.8, 0.9, 1),
                              softmax = c(TRUE, FALSE),
                              n_seeds = c(500),
                              adj_hub = c(FALSE)) %>%
  mutate(outfile_filename=str_c(here(),
                                "/PolyGene/benchmarking/PolyNet/",
                                "Mendelian_Freund_2018_Benchmarking/",
                                "MAGMA_0kb/STRINGcoexpression_PolyNet/",
                                "RWR_restart",
                                restart_prob,
                                "_softmax",
                                softmax,
                                "_nSeeds",
                                str_replace_na(n_seeds, "NA"),
                                "_HubGeneAdjust",
                                adj_hub,
                                ".csv"))


# Apply the function to each combination of parameters specified in the parameter_grid df
# NOTE: THIS TAKES ~2 HOURS WITH CURRENT PARAM GRID
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
out_summary <- out_summary %>%
  mutate(n_seeds = str_replace_na(n_seeds, "7921")) %>%
  mutate(n_seeds = as.numeric(n_seeds))


#####PLOTS######
# Plot the impact of each of the parameters on AUC
plots <- generate_param_plots(data = out_summary,
                              y_variable="auc",
                              params = c("restart_prob", "softmax", "n_seeds", "adj_hub"),
                              phenotype_column = "mendelian_disease_group")
plots[[1]]

plots_restart_0.7 <- generate_param_plots(data = filter(out_summary, restart_prob==0.7),
                                          y_variable="auc",
                                          params = c("restart_prob", "softmax", "n_seeds", "adj_hub"),
                                          phenotype_column = "mendelian_disease_group")
plots_restart_0.7[[4]]


plots_best_param <- generate_param_plots(data = filter(out_summary, softmax==FALSE & n_seeds==500 & adj_hub==FALSE),
                                         y_variable="auc",
                                         params = c("restart_prob", "softmax", "n_seeds", "adj_hub"),
                                         phenotype_column = "mendelian_disease_group")
plots_best_param[[1]]
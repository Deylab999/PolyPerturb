library(tidyverse)
library(pROC)
library(here)
library(data.table)

#load the benchmarking function
source(str_c(here(), "/PolyGene/benchmarking/benchmarking_functions.R"))

#make the gene rankings for the best conditions
best_rwr_STRINGppi_filename = str_c(here(),"/data/BestPolyNetPPI_Rankings",
                          "RWR_restart",
                          0.7,
                          "_softmax",
                          TRUE,
                          "_nSeeds",
                          str_replace_na(500, "NA"),
                          "_HubGeneAdjust",
                          FALSE,
                          ".csv")
best_rwr_STRINGcoex_filename = str_c(here(),"/data/BestPolyNetCoExpression_Rankings",
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

###output the rankings for the best PolyNet STRING RWR conditions
best_rwr_rankings_string <- run_RWR(gene_scores = fread(str_c(here(), "/data/MAGMA_v108_GENE_0_ZSTAT.txt")),
                             restart_prob = 0.7,
                             thresh_score = NA,
                             softmax = TRUE,
                             n_seed_genes = 500,
                             adjust_hub_genes = FALSE,
                             graph = create_string_graph(edge_threshold = 400))
fwrite(best_rwr_rankings_string, best_rwr_STRINGppi_filename)

###output the rankings for the best PolyNet STRING CoExpression RWR conditions
best_rwr_rankings_coex <- run_RWR(gene_scores = fread(str_c(here(), "/data/MAGMA_v108_GENE_0_ZSTAT.txt")),
                             restart_prob = 0.5,
                             thresh_score = NA,
                             softmax = FALSE,
                             n_seed_genes = 500,
                             adjust_hub_genes = FALSE,
                             graph = create_STRING_coexp_graph(edge_threshold = 400))
fwrite(best_rwr_rankings_coex, best_rwr_STRINGcoex_filename)

#make a parameter file
file_info_df <- tribble(
  ~filename, ~label, ~restart_prob, ~softmax, ~n_seeds, ~adj_hub,
  pops_filename, "PoPs", NA, NA,NA,NA,
  best_rwr_STRINGppi_filename, "Best_PolyNetSTRING_PPI", 0.7,TRUE,500,FALSE,
  #best_rwr_STRINGcoex_filename, "Best_PolyNetSTRING_CoExpression", 0.7,TRUE,500,FALSE,
  #old_rwr_filename, "Old_PolyNetSTRING", 0.5, FALSE,11000,FALSE,
  #old_pn_full_filename, "Old_PolyNetSTRING+CoEx", 0.5, FALSE,11000,FALSE
)

###run the benchmarking
mend_grouped_output <- benchmark_plot_and_results(file_info_df = file_info_df,
                                                  gt_file = str_c(here(),"/data/freund_2018_monogenic_to_complex_ground_truth.csv"),
                                                  threshold = 0.95,
                                                  thresholds = c(threshold),
                                                  y_axis_columns = c("fisher_or", "auc"),
                                                  x_axis_column = "mendelian_disease_group")

mend_grouped_output[[1]]
mend_grouped_output[[2]]

ct_grouped_output <- benchmark_plot_and_results(file_info_df = file_info_df,
                                                  gt_file = str_c(here(),"/data/freund_2018_monogenic_to_complex_ground_truth.csv"),
                                                  threshold = 0.95,
                                                  thresholds = c(threshold),
                                                  y_axis_columns = c("fisher_or", "auc"),
                                                  x_axis_column = "complex_trait")
ct_grouped_output[[1]]
ct_grouped_output[[2]]






library(tidyverse)
library(pROC)
library(here)
library(data.table)

#load the benchmarking function
source(str_c(here(), "/PolyGene/benchmarking/benchmarking_functions.R"))

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

#make a parameter file
file_info_df <- tribble(
  ~filename, ~label, ~restart_prob, ~softmax, ~n_seeds, ~adj_hub,
  pops_filename, "PoPs", NA, NA,NA,NA,
  best_rwr_filename, "Best_PolyNetSTRING", 0.7,TRUE,500,FALSE,
  old_rwr_filename, "Old_PolyNetSTRING", 0.5, FALSE,11000,FALSE,
  old_pn_full_filename, "Old_PolyNetSTRING+CoEx", 0.5, FALSE,11000,FALSE
)

###output the rankings for the best RWR conditions
best_rwr_rankings <- run_RWR(gene_scores = fread(str_c(here(), "/data/MAGMA_v108_GENE_0_ZSTAT.txt")),
                             restart_prob = 0.7,
                             thresh_score = NA,
                             softmax = TRUE,
                             n_seed_genes = 500,
                             adjust_hub_genes = FALSE,
                             graph = create_string_graph(edge_threshold = 400))
fwrite(best_rwr_rankings, best_rwr_filename)

###run the benchmarking
mend_grouped_output <- benchmark_plot_and_results(file_info_df)




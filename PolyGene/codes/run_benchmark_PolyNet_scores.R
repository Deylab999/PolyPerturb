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

# Define the grid of parameters to generate results for
parameter_grid <- expand.grid(restart_prob = c(0, 0.2, 0.5, 0.7, 0.8, 0.9, 1),
                              softmax = c(TRUE, FALSE),
                              n_seeds = c(NA, 100, 500),
                              adj_hub = c(TRUE, FALSE)) %>%
  mutate(outfile_filename=str_c(here(),
                                "/PolyGene/benchmarking/PolyNet/",
                                "Mendelian_Freund_2018_Benchmarking/",
                                "MAGMA_0kb/STRING_PolyNet/",
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
  mutate(n_seeds = str_replace_na(n_seeds, "14891")) %>%
  mutate(n_seeds = as.numeric(n_seeds))


#####PLOTS######
# Plot the impact of each of the parameters on AUC
plots <- generate_param_plots(data = out_summary,
                              y_variable="auc",
                              params = c("restart_prob", "softmax", "n_seeds", "adj_hub"),
                              phenotype_column = "mendelian_disease_group")
plots[[1]]

plots_restart_0.7 <- generate_param_plots(data = filter(out_summary, restart_prob==0.7),
                                          y_variable="fisher_or",
                                          params = c("restart_prob", "softmax", "n_seeds", "adj_hub"),
                                          phenotype_column = "mendelian_disease_group")
plots_restart_0.7[[4]]

plots_softmaxTRUE <- generate_param_plots(data = filter(out_summary, restart_prob==0.7 & softmax==TRUE),
                                          y_variable="fisher_or",
                                          params = c("restart_prob", "softmax", "n_seeds", "adj_hub"),
                                          phenotype_column = "mendelian_disease_group")
plots_softmaxTRUE[[3]]

plots_restart_nseeds500 <- generate_param_plots(data = filter(out_summary, restart_prob==0.7 & n_seeds==500),
                                                    y_variable="fisher_or",
                                                    params = c("restart_prob", "softmax", "n_seeds", "adj_hub"),
                                                    phenotype_column = "mendelian_disease_group")
plots_restart_nseeds500[[2]]

plots_best_param <- generate_param_plots(data = filter(out_summary, softmax==TRUE & n_seeds==500 & adj_hub==FALSE),
                                         y_variable="fisher_or",
                                         params = c("restart_prob", "softmax", "n_seeds", "adj_hub"),
                                         phenotype_column = "mendelian_disease_group")
plots_best_param[[1]]

##########################################################################################
######################run benchmarking on PoPs scores &  #################################
################## compare to STRING PolyNet with suggested params  ######################
##########################################################################################

# read in the PoPs data
pops_raw <- fread(str_c(here(),
                        "/data/POPs_MAGMA_0kb.txt")) %>%
  rename(gene_symbol = V1)

# read in the PolyNet scores for the best set of paramters
rwr_best_param <- fread(str_c(here(),
                                    "/PolyGene/benchmarking/PolyNet/",
                                    "Mendelian_Freund_2018_Benchmarking/",
                                    "MAGMA_0kb/STRING_PolyNet/",
                                    "RWR_restart",
                                    "0.7",
                                    "_softmax",
                                    "TRUE",
                                    "_nSeeds",
                                    "500",
                                    "_HubGeneAdjust",
                                    "FALSE",
                                    ".csv"))


#find overlap of genes and phenotypes in pops vs. rwr
overlapping_phenos <- intersect(colnames(pops_raw),
                                unique(rwr_best_param$phenotype))

#only use overlapping genes between STRING and PoPs to benchmark
string_graph <- create_string_graph(edge_threshold = 400)
overlapping_genes <- intersect(as.character(V(string_graph)$name),
                               pops_raw$gene_symbol)

#filter pops data to overlapping genes and phenos
pops_raw <- pops_raw %>%
  select(gene_symbol, all_of(overlapping_phenos)) %>%
  filter(gene_symbol %in% overlapping_genes)

#run pops benchmarking
pops_bm <- compute_sensitivity_specificity(data = rwr_best_param,
                                             ground_truth_df = gt,
                                             thresholds = c(0, 0.01, 0.05,
                                                            0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,
                                                            0.95, 0.99, 1))
fwrite(pops_bm, str_c(here(),
                            "/PolyGene/benchmarking/PolyNet/",
                            "Mendelian_Freund_2018_Benchmarking/",
                            "MAGMA_0kb/",
                            "POPs_MAGMA_0kb.txt"))

pops_bm_summary <- summarize_benchmarks(restart_prob = 0,
                                        softmax = 0,
                                        n_seeds = 0,
                                        adj_hub = 0,
                                        benchmark_filename = str_c(here(),
                                                                   "/PolyGene/benchmarking/PolyNet/",
                                                                   "Mendelian_Freund_2018_Benchmarking/",
                                                                   "MAGMA_0kb/",
                                                                   "POPs_MAGMA_0kb.txt"),
                                        grouping_column = "mendelian_disease_group") %>%
  mutate(algo = "POPs_MAGMA_0kb.txt") %>%
  select(-restart_prob, -softmax, -n_seeds, -adj_hub)

#run old PolyGene score benchmarking
old_PolyNet <- fread(str_c(here(),
                        "/data/original_output/PPI_STRING_MAGMA.txt")) %>%
  rename(gene_symbol = V1) %>%
  select(gene_symbol, any_of(overlapping_phenos)) %>%
  filter(gene_symbol %in% overlapping_genes)

old_pn_bm <- compute_sensitivity_specificity(data = old_PolyNet,
                                           ground_truth_df = gt,
                                           thresholds = c(0, 0.01, 0.05,
                                                          0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,
                                                          0.95, 0.99, 1))
fwrite(old_pn_bm, str_c(here(),
                      "/PolyGene/benchmarking/PolyNet/",
                      "Mendelian_Benchmarks_",
                      "old_PPI_STRING_MAGMA.txt"))

old_pn_bm_summary <- summarize_benchmarks(restart_prob = 0,
                                        softmax = 0,
                                        n_seeds = 0,
                                        adj_hub = 0,
                                        benchmark_filename = str_c(here(),
                                                                   "/PolyGene/benchmarking/PolyNet/",
                                                                   "Mendelian_Benchmarks_",
                                                                   "old_PPI_STRING_MAGMA.txt"),
                                        grouping_column = "mendelian_disease_group") %>%
  mutate(algo = "old_PPI_STRING_MAGMA.txt") %>%
  select(-restart_prob, -softmax, -n_seeds, -adj_hub)


####combine everything together
best_params_pops_benchmark <- rwr_best_param %>%
  filter(phenotype %in% overlapping_phenos) %>%
  group_by(!!sym("mendelian_disease_group"), threshold) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(algo = "PolyNet_STRING_PPI") %>%
  select(all_of(colnames(pops_bm_summary))) %>%
  bind_rows(., pops_bm_summary, old_pn_bm_summary) 

overall_best <- best_params_pops_benchmark %>%
  filter(algo == "PolyNet_STRING_PPI") %>%
  group_by(algo) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(mendelian_disease_group = "Overall")

best_params_pops_benchmark <- best_params_pops_benchmark %>%
  bind_rows(overall_best ) %>%
  arrange(mendelian_disease_group)

#plot the comparison
ggplot(best_params_pops_benchmark,
       aes(x = mendelian_disease_group, y = auc, color = algo)) +
  geom_point(position = position_dodge(width = 0.2)) +
  geom_errorbar(aes(ymin = auc_lower_ci, ymax = auc_upper_ci), 
                position = position_dodge(width = 0.2), 
                width = 0.2) +
  labs(title = "Compare AUC for PoPs and Best Parameters of STRING PolyNet",
       x = "Mendelian Phenotype Group",
       y = "AUC") +
  scale_color_discrete(name = "Algorithm") +
  theme_bw() +
  scale_y_continuous(breaks = seq(0.5, 1, by = 0.05), labels = seq(0.5, 1, by = 0.05)) +
  theme(axis.text.x = element_text(angle = 75, hjust = 1))


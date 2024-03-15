library(here)
library(tidyverse)

source(str_c(here(), "/PolyGene/codes/create_string_graph_object.R"))
source(str_c(here(), "/PolyGene/codes/ppi_string_RWR.R"))
source(str_c(here(), "/PolyGene/codes/run_pagerank.R"))

#create some data
graph <- create_string_graph(edge_threshold = 400)
test_rwr <- run_RWR(graph=graph, softmax = TRUE, restart_prob = 0.5)

true_positives <- bind_cols(gene_symbol = sample(test_rwr$gene_symbol, 20),
                            phenotype = sample(colnames(test_rwr)[2:247], 20)) %>%
  mutate(ground_truth = 1)
true_negatives <- bind_cols(gene_symbol = sample(test_rwr$gene_symbol, 20),
                            phenotype = sample(colnames(test_rwr)[2:247], 20)) %>%
  mutate(ground_truth = 0)

ground_truth <- bind_rows(true_positives, true_negatives) 

compute_sensitivity_specificity <- function(data = test_rwr,
                                            ground_truth = ground_truth) {
  library(tidyverse)
  
    # Calculate percentile ranks for each phenotype
  percentile_ranks <- data %>%
    mutate(across(-gene_symbol, ~na_if(., NA) %>% replace_na(mean(.)) %>% percent_rank())) %>%
    pivot_longer(-gene_symbol, names_to = "phenotype", values_to = "percentile_rank")
    
    # Merge with true positives and true negatives
    merged_data <- percentile_ranks %>%
      left_join(ground_truth, by = c("gene_symbol", "phenotype"))
    
    
    # Calculate sensitivity and specificity for each phenotype
    sensitivity_specificity <- percentile_ranks %>%
      group_by(phenotype) %>%
      summarize(sensitivity = sum(ground_truth == 1 & percentile_rank >= quantile(percentile_rank, 0.8)) / sum(ground_truth == 1),
                specificity = sum(ground_truth == 0 & percentile_rank < quantile(percentile_rank, 0.2)) / sum(ground_truth == 0))
    
    
    return(sensitivity_specificity)
  }
  
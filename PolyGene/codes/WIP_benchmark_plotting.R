library(tidyverse)
library(pROC)
library(here)
library(data.table)

#load the benchmarking function
source(str_c(here(), "/PolyGene/codes/benchmark_PolyNet_scores.R"))

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


###output the rankings for the best RWR conditions
#best_rwr_rankings <- run_RWR(gene_scores = fread(str_c(here(), "/data/MAGMA_v108_GENE_0_ZSTAT.txt")),
#                             restart_prob = 0.7,
#                             thresh_score = NA,
#                             softmax = TRUE,
#                             n_seed_genes = 500,
#                             adjust_hub_genes = FALSE,
#                             graph = create_string_graph(edge_threshold = 400))
#fwrite(best_rwr_rankings, best_rwr_filename)

#make a parameter file
file_info_df <- tribble(
  ~filename, ~label, ~restart_prob, ~softmax, ~n_seeds, ~adj_hub,
  pops_filename, "PoPs", NA, NA,NA,NA,
  best_rwr_filename, "Best_PolyNetSTRING", 0.7,TRUE,500,FALSE,
  old_rwr_filename, "Old_PolyNetSTRING", 0.5, FALSE,11000,FALSE,
  old_pn_full_filename, "Old_PolyNetSTRING+CoEx", 0.5, FALSE,11000,FALSE
)

# Step 1: Data Preprocessing
data_preprocessing <- function(file_info_df) {
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

pp_out <- data_preprocessing(file_info_df)

# Step 2: Benchmarking
benchmarking <- function(data_list = pp_out,
                         gt_file = gt,
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

# Example usage with a list containing a single dataframe
test_bm_summary_list <- benchmarking(data_list = pp_out)

# Step 3: Visualization
plot_comparison <- function(bm_summary_list = test_bm_summary_list,
                            y_axis_column = "fisher_or",
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

benchmark_plot_and_results <- function(file_info_df,
                                       gt_file = gt,
                                       threshold = 0.95,
                                       thresholds = c(threshold),
                                       y_axis_columns = c("fisher_or", "auc"),
                                       x_axis_column = "mendelian_disease_group") {
    
    # Step 1: Data Preprocessing
    pp_out <- data_preprocessing(file_info_df)
    
    # Step 2: Benchmarking
    bm_summary_list <- benchmarking(data_list = pp_out,
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
  

# Example usage
output_list <- benchmark_plot_and_results(file_info_df)





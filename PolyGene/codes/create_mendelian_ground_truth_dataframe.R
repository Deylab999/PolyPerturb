#' Create Mendelian ground truth mapping file
#'
#' This function creates a Mendelian ground truth mapping file based on provided input files.
#' It identifies "True Negatives" as genes that are in the ground truth file but are not present for any other member of the larger group.
#'
#' @param mend_gt_map_file Character string specifying the file path for the Mendelian to GWAS phenotype mapping CSV file.
#' @param mend_gt_genes_file Character string specifying the file path for the CSV file containing monogenic disease genes by
#' Mendelian disease group from Freund et al (2018).
#' @return A data frame containing the Mendelian ground truth mapping to use in sensitivity/specificity benchmarking.
#'
mendelian_ground_truth <- function(mend_gt_map_file = str_c(here(), "/data/mendelian_to_gwas_pheno_mapping.csv"),
                                   mend_gt_genes_file = str_c(here(),"/data/freund_2018_monogenic_disease_genes.csv")) {
  
  require(data.table)
  require(dplyr)
  require(here)
  
  mend_gt_map <- fread(mend_gt_map_file)
  mend_gt_genes <- fread(mend_gt_genes_file)
  
  mend_gt <- mend_gt_map %>%
    drop_na() %>%
    left_join(mend_gt_genes, by = "mendelian_disease_group", relationship = "many-to-many") %>%
    select(gene_symbol, phenotype, mendelian_disease_group) %>%
    drop_na()
  
  gt_final <- data.frame()
  for (p in unique(mend_gt$phenotype)) {
    
    positive_df <- mend_gt %>%
      filter(phenotype == p) %>%
      mutate(ground_truth = 1)
    
    group <- positive_df %>%
      pull(mendelian_disease_group) %>%
      unique()
    
    # get a list of all genes not represented as true positives
    # in the same phenotype group
    negative_genes <- mend_gt %>%
      filter(!gene_symbol %in% positive_df$gene_symbol) %>%
      pull(gene_symbol) %>%
      unique()
    
    #define true negatives as those genes that are in the negative_genes list
    negative_df <- bind_cols(
      gene_symbol = negative_genes,
      phenotype = rep(p, length(negative_genes))
    ) %>%
      mutate(mendelian_disease_group = group[1],
             ground_truth = 0)
    
    gt_final <- rbind(gt_final, positive_df, negative_df)
  }
  
  return(gt_final)
}

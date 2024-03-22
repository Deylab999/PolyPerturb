library(here)
library(tidyverse)

source(str_c(here(), "/PolyGene/codes/create_string_graph_object.R"))

#' @param gene_scores  matrix of gene scores with genes along the rows and gene sets along the columns
#' @param restart_prob The probability of restart for RWR method (default is 0.8, which appears to 
#' perform well for STRING).
#' @param thresh_score A minimum score required for a gene to be considered a seed gene
#' @param softmax When FALSE, all seed genes are weighted equally. When TRUE, seed genes are weighted by softmax normalization
#' @param n_seed_genes The number of seed genes to use (default is NA, which means all genes are used as seeds. If
#' softmax=TRUE, then softmax weighting will be applied across all genes.)
#' @param graph iGraph object to run RWR on. Default is to re-compute STRING graph using create_string_graph()
#' function, although the function is faster if a pre-created iGraph object is added.
#' @param adjust_hub_genes Strong correction for bias towards prioritizing "hub" genes.
#' When TRUE, subtract off the RWR values when the restart parameter is set to 0 and all overlapping genes
#' are used as seed genes. NOTE: This is very likely to result in negative values. FALSE by default.
#' @results:  A list with two entries : 'score' : The prioritization score between 0 and 1 based on the PPI RWR.
#' We recommend thresholding it
#' to top X% genes. 'score_sd': Standard deviation of the score.
#' @import dnet, igraph, here, data.table, tidyverse
#'

run_RWR = function(gene_scores = fread(str_c(here(), "/data/MAGMA_v108_GENE_0_ZSTAT.txt")),
                   restart_prob = 0.8,
                   thresh_score = NA,
                   softmax = TRUE,
                   n_seed_genes = NA,
                   adjust_hub_genes = FALSE,
                   graph = create_string_graph(edge_threshold = 400)){ 

  library(here)
  library(tidyverse)
  library(dnet)
  library(data.table)
  library(igraph)
  
  #####################  Subset the gene aliases to proteins that occur in the graph  #######################
  
  #load in the graph
  string_graph = graph
  
  #deal with no gene column name in input gene score file
  if("V1" %in% colnames(gene_scores)){
    gene_scores <- gene_scores %>% rename("gene_symbol" = "V1")
  }
  
  #throw some warnings if something seems off
  if(!'gene_symbol' %in% colnames(gene_scores)){
    warning("A column named `gene_symbol` was not provided in the input gene_scores data.")
  }
  if(is.null(colnames(gene_scores))){
    warning("The names of gene sets in columns of gene_scores not provided. We assign them names Pheno1, Pheno2....")
    colnames(gene_scores) = c("gene_symbol", str_c("Pheno", 1:ncol(gene_scores)))
  }

  common_genes = intersect(gene_scores$gene_symbol, V(string_graph)$name)
  print(str_c("There are ",
              length(common_genes),
              " genes in common between the input gene scores and nodes in the STRING graph"))
  
  if(length(common_genes) < 100){
    stop("The number of genes in the network matched in the gene scores matrix is too few (<100).
    Either check for match discepancy between the network gene names and gene_scores input genes or
    include more genes in the gene scores file.")
  }
  
    #############  Create seed genes for RWR  method    ###################################
  
  seed_genes_df <- gene_scores %>%
    pivot_longer(!gene_symbol, names_to = "phenotype", values_to = "score") %>%
    filter(gene_symbol %in% common_genes)
    
  if(!is.na(thresh_score)) {
    seed_genes_df <- seed_genes_df %>%
      filter(score >= thresh_score) #filter by a minimum score threshold
  }
  
  if(!is.na(n_seed_genes)) {
    seed_genes_df <- seed_genes_df %>%
      group_by(phenotype) %>%
      top_n(n_seed_genes, wt = score) %>% #filter to the specified number of seed genes
      ungroup()
  }

  seed_genes_df <- seed_genes_df %>%
    group_by(phenotype) %>%
    mutate(min_score = min(score, na.rm=TRUE)) %>% #ensure the score is always positive
    mutate(score = ifelse(min_score < 0, 
                          score + abs(min_score),
                          score)) %>%
    select(-min_score) %>%
    mutate(personalization_vector = exp(score)/sum(exp(score), na.rm=TRUE)) %>% #SOFTMAX to weight seeds
    ungroup() %>%
    arrange(phenotype, desc(score))
  
  seed_mat <- seed_genes_df %>%
    select(-score) %>%
    pivot_wider(names_from = phenotype,
                values_from = personalization_vector,
                values_fill = 0)
  tmp_rownames <- seed_mat$gene_symbol
  seed_mat <- seed_mat %>% select(-gene_symbol) %>% as.matrix(seed_mat)
  rownames(seed_mat) <- tmp_rownames
  
  # Make sure that all nodes in the graph
  # have a corresponding row in the seed matrix
  all_genes <- V(string_graph)$name
  seed_mat <- rbind(seed_mat,
                    matrix(0, nrow = length(all_genes) - nrow(seed_mat),
                           ncol = ncol(seed_mat),
                           dimnames = list(setdiff(all_genes, rownames(seed_mat)),
                                           colnames(seed_mat))))
  
  #retain the ability to weight all seed genes equally
  if(softmax == FALSE){
    seed_mat <- ifelse(seed_mat != 0, 1, 0)
  }
  
  #############  Run RWR    ################################################
  set.seed(893)
  PTmatrix <- dRWR(g=string_graph, normalise="laplacian", setSeeds=seed_mat,
                     restart=restart_prob, parallel=TRUE)
    rownames(PTmatrix) = as_ids(V(string_graph))
    colnames(PTmatrix) = colnames(seed_mat)
    
  if(adjust_hub_genes == TRUE){
    all_genes_seeds <- matrix(rep(1/length(common_genes), length(common_genes)), 
                              nrow = length(common_genes),
                              dimnames = list(common_genes, NULL))
    hub_genes_rwr <- dRWR(g=string_graph, normalise="laplacian", setSeeds=all_genes_seeds,
                     restart=0, parallel=TRUE)
    rownames(hub_genes_rwr) = as_ids(V(string_graph))
    adjusted_PTmatrix <- apply(as.matrix(PTmatrix), 2,
                               function(col) col - as.matrix(hub_genes_rwr))
    adjusted_PTmatrix <- matrix(unlist(adjusted_PTmatrix),
                                nrow = nrow(PTmatrix),
                                ncol = ncol(PTmatrix), byrow = FALSE)
    rownames(adjusted_PTmatrix) = as_ids(V(string_graph))
    colnames(adjusted_PTmatrix) = colnames(seed_mat)
    PTmatrix = adjusted_PTmatrix
  }
    
  out_rwr_df <- as.data.frame(as.matrix(PTmatrix)) %>%
    rownames_to_column(var = "gene_symbol")
  return(out_rwr_df)
}

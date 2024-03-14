library(here)
library(tidyverse)

source(str_c(here(), "/PolyGene/codes/create_string_graph_object.R"))

#' @param gene_scores  matrix of gene scores with genes along the rows and gene sets along the columns
#' @param damping_factor The damping factor, which controls diffusion in the graph. Some similarities to the restart probability in RWR.
#' @param thresh_score An upper threshold for the scores in the gene score matrix.
#' @param softmax When FALSE, all seed genes are weighted equally. When TRUE, seed genes are weighted by softmax normalization
#' @param n_seed_geens The number of seed genes to use. (default is 100).
#' @param graph iGraph object to run RWR on. Default is to re-compute STRING graph using create_string_graph()
#' function, although the function is faster if a pre-created iGraph object is added.
#' @param adjust_hub_genes Strong correction for bias towards prioritizing "hub" genes.
#' When TRUE, subtract off the PageRank values when the restart parameter is set to 0 and all overlapping genes
#' are used as seed genes. NOTE: This is very likely to result in negative values. FALSE by default.
#' @results:  A list with two entries : 'score' : The prioritization score between 0 and 1 based on the PPI RWR.
#' We recommend thresholding it
#' to top X% genes. 'score_sd': Standard deviation of the score.
#' @import igraph, here, data.table, tidyverse
#'

run_pagerank = function(gene_scores = fread(str_c(here(), "/data/MAGMA_v108_GENE_0_ZSTAT.txt")),
                   damping_factor = 0.75,
                   thresh_score = 5,
                   softmax = TRUE,
                   n_seed_geens=100,
                   adjust_hub_genes = FALSE,
                   graph = create_string_graph(edge_threshold = 400)){ 
  
  library(here)
  library(tidyverse)
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
  
  #############  Create seed genes for PageRank  method    ###################################
  
  seed_genes_df <- gene_scores %>%
    pivot_longer(!gene_symbol, names_to = "phenotype", values_to = "score") %>%
    filter(gene_symbol %in% common_genes) %>%
    group_by(phenotype) %>%
    filter(score >= thresh_score) %>% #filter by a minimum threshold
    top_n(n_seed_geens, wt = score) %>% #filter to the specified number of seed genes
    mutate(personalization_vector = exp(score)/sum(exp(score))) %>% #SOFTMAX to weight seeds
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
  
  #############  Run PageRank  #############################################
  set.seed(893)
  pr_matrix <- matrix(NA, nrow = nrow(seed_mat), ncol = ncol(seed_mat))
  for (i in 1:ncol(seed_mat)) {
    pr_scores <- page_rank(graph,
                           damping = damping_factor,
                           personalized = seed_mat[, i],
                           weights=NA, #removes edge weights
                           )
    
    #do hub gene adjustment if specified
    if (adjust_hub_genes) {
      hub_genes_pr <- page_rank(graph, damping = damping_factor)
      pr_scores$vector <- pr_scores$vector - hub_genes_pr$vector
    }
    
    pr_matrix[, i] <- pr_scores$vector
    row_names <- names(pr_scores$vector)
  }
  colnames(pr_matrix) = colnames(seed_mat)
  rownames(pr_matrix) <- row_names
  
  out_pr_df <- as.data.frame(as.matrix(pr_matrix)) %>%
    rownames_to_column(var = "gene_symbol")
  return(out_pr_df)
}

#' @param gene_scores  matrix of gene scores with genes along the rows and gene sets along the columns
#' @param restart_prob The probability of restart for RWR method (default is 0.5).
#' @param thresh_score An upper threshold for the scores in the gene score matrix.
#' @param softmax When FALSE, all seed genes are weighted equally. When TRUE, seed genes are weighted by softmax normalization
#' @param NUM_SUBGENES The number of genes to keep in each sub-sample per run for a gene score. If NULL, we take 80% of genes per gene
#'                     score, thresholded below at 100.
#' @results:  A list with two entries : 'score' : The prioritization score between 0 and 1 based on the PPI RWR. We recommend thresholding it
#' to top X% genes. 'score_sd': Standard deviation of the score.
#' @import dnet, igraph, here, data.table, tidyverse
#'
#'

ppi_string_RWR = function(gene_scores = fread(str_c(here(), "/data/MAGMA_v108_GENE_0_ZSTAT.txt")),
                   restart_prob = 0.5,
                   thresh_score = 5,
                   softmax = TRUE,
                   NUM_RUNS = 5,
                   NUM_SUBGENES=100){ 

  library(here)
  library(tidyverse)
  library(dnet)
  library(data.table)
  library(igraph)
  
  #############################  Create a String database graph ##############################################
  
  #get the STRING to gene symbol mapping
  string_alias <- str_c(here(),
                        "/STRING/9606.protein.aliases.v12.0.txt.gz")
  
  alias_df <- fread(string_alias) %>%
    filter(source == "BioMart_HUGO") %>%
    select(string_id = `#string_protein_id`, symbol = alias) %>%
    distinct()
  
  #create dataframe of all pairwise edges to be included in the graph
  string_local <- str_c(here(),
                        "/STRING/9606.protein.physical.links.v12.0.txt.gz")
  
  string_df <- fread(string_local) %>%
    filter(combined_score >= 400) %>%
    left_join(alias_df, by=c("protein1" = "string_id"), relationship = "many-to-many") %>%
    rename(symbol1 = symbol) %>%
    left_join(alias_df, by=c("protein2" = "string_id"), relationship = "many-to-many") %>%
    rename(symbol2 = symbol) %>%
    select(symbol1, symbol2, combined_score) %>%
    filter(symbol1 != "NA") %>%
    filter(symbol2 != "NA")
    
  #create an empty graph and fill it in from the dataframe
  string_graph <- graph_from_data_frame(string_df, directed = FALSE)

  #####################  Subset the gene aliases to proteins that occur in the graph  #######################

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
    filter(gene_symbol %in% common_genes) %>%
    group_by(phenotype) %>%
    filter(score >= thresh_score) %>% #filter by a minimum threshold
    top_n(NUM_SUBGENES, wt = score) %>% #filter to the top NUM_SUBGENES
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
  
  #retain the ability to weight all seed genes equally
  if(softmax==FALSE){
    seed_mat <- ifelse(seed_mat != 0, 1, 0)
  }
  
  #############  Run RWR    ################################################
  PTmatrix <- dRWR(g=string_graph, normalise="laplacian", setSeeds=seed_mat,
                     restart=restart_prob, parallel=TRUE)
    rownames(PTmatrix) = as_ids(V(string_graph))
    colnames(PTmatrix) = colnames(seed_mat)
    
  score_sd = apply(PTmatrix, 1, sd)
  outlist = list("score" =  as.matrix(PTmatrix), "score_sd" = score_sd)
  return(outlist)
}

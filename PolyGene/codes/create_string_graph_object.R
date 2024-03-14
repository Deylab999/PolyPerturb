#' @param string_network_df_path  path to a text file that contains all functional relationships between genes.
#' See https://string-db.org/cgi/download for files. STRING/9606.protein.physical.links.v12.0.txt.gz is default.
#' @param string_gene_names_path path to a text file that contains aliases for the STRING gene IDs.
#' This is used to map from STRING IDs to HUGO Gene Symbols. STRING/9606.protein.aliases.v12.0.txt.gz is default.
#' @param edge_threshold Minimum threshold to draw an edge in the graph between two genes (default is 400)
#' @results: An iGraph object that represents the STRING graph to be used by the RWR function.
#' @import igraph, here, data.table, tidyverse

create_string_graph = function(string_network_df_path = str_c(here(),
                                                              "/STRING/9606.protein.physical.links.v12.0.txt.gz"),
                               string_alias_path = str_c(here(),
                                                         "/STRING/9606.protein.aliases.v12.0.txt.gz"),
                               edge_threshold = 400){
  library(here)
  library(tidyverse)
  library(data.table)
  library(igraph)
  
  #############################  Create a String database graph ##############################################
  
  #get the STRING to gene symbol mapping
  alias_df <- fread(string_alias_path) %>%
    filter(source == "BioMart_HUGO") %>%
    select(string_id = `#string_protein_id`, symbol = alias) %>%
    distinct()
  
  #create dataframe of all pairwise edges to be included in the graph
  string_df <- fread(string_network_df_path) %>%
    filter(combined_score >= edge_threshold) %>%
    left_join(alias_df, by=c("protein1" = "string_id"), relationship = "many-to-many") %>%
    rename(symbol1 = symbol) %>%
    left_join(alias_df, by=c("protein2" = "string_id"), relationship = "many-to-many") %>%
    rename(symbol2 = symbol) %>%
    select(symbol1, symbol2, combined_score) %>%
    filter(symbol1 != "NA") %>%
    filter(symbol2 != "NA")
  
  #create an empty graph and fill it in from the dataframe
  string_graph <- graph_from_data_frame(string_df, directed = FALSE)
  return(string_graph)
}

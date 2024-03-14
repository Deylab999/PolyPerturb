source(str_c(here(), "/PolyGene/codes/create_string_graph_object.R"))
source(str_c(here(), "/PolyGene/codes/ppi_string_RWR.R"))
source(str_c(here(), "/PolyGene/codes/run_pagerank.R"))

graph <- create_string_graph(edge_threshold = 400)

test_rwr <- run_RWR(graph=graph, softmax = FALSE)
test_ppr <- run_pagerank(graph=graph, softmax = FALSE)

matched_colnames  <- colnames(test_rwr)[!colnames(test_rwr) %in% c("gene_symbol", "method")]


correlation_results <- lapply(matched_colnames, function(col_name) {
  cor_result <- cor(test_rwr[[col_name]], test_ppr[[col_name]], method = "spearman")
  return(data.frame(Column = col_name, Spearman_Correlation = cor_result))
})
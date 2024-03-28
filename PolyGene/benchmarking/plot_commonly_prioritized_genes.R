pops_filename = str_c(here(),"/data/POPs_MAGMA_0kb.txt")
best_rwr_STRINGppi_filename = str_c(here(),"/data/BestPolyNetPPI_Rankings",
                                    "RWR_restart",
                                    0.7,
                                    "_softmax",
                                    TRUE,
                                    "_nSeeds",
                                    str_replace_na(500, "NA"),
                                    "_HubGeneAdjust",
                                    FALSE,
                                    ".csv")

rank_original <- #fread(pops_filename) %>%
  fread("~/PolyPerturb/data/original_output/PPI_STRING_MAGMA.txt") %>%
  rename(gene_symbol = V1) %>%
  mutate_at(vars(-gene_symbol), funs(rank(-., ties.method = "min"))) %>%
  mutate(mean_rank = rowMeans(select(., -gene_symbol))) %>%
  mutate(median_rank = apply(select(., -gene_symbol), 1, median) ) %>%
  mutate(n_traits_rank_top10 = rowSums(select(., -gene_symbol) <= 10)) %>%
  mutate(percent_traits_top10 = round(n_traits_rank_top10 / 120, 2)) %>%
  arrange(desc(percent_traits_top10)) %>%
  mutate(gene_symbol = factor(gene_symbol, levels = gene_symbol)) %>%
  head(20)

top_genes <- rank_original %>%
  head(20) %>%
  pull(gene_symbol)

gt <- fread(str_c(here(),"/data/freund_2018_monogenic_to_complex_ground_truth.csv"))

gt_most_common <- 
  
  tmp <- gt %>%
  filter(ground_truth==1) %>%
  filter(gene_symbol=="EGFR")
  select(mendelian_disease_group, gene_symbol) %>%
  distinct() %>%
  mutate(top_hub = ifelse(gene_symbol %in% top_genes, 1, 0)) %>%
  group_by(mendelian_disease_group) %>%
  summarise(n_top_hub = sum(top_hub),
            n_total = n(),
            percent = n_top_hub / n_total) %>%
  arrange(desc(percent)) %>%
  ungroup()


top_genes_df <- rank_original %>%
  filter(gene_symbol %in% top_genes) %>%
  select(-gene_symbol) %>%
  t()

colnames(top_genes_df) <- top_genes

plot_topgenes <- ggplot(rank_original, aes(x = gene_symbol, y = percent_traits_top10)) +
  geom_bar(stat = "identity") +  # Create the barplot
  geom_text(aes(label = n_traits_rank_top10), vjust = -0.5) + 
  labs(x = "Top 20 Genes by STRING PolyNet Ranking Across Traits",
       y = "% Traits with Gene Ranked in Top 10",
       title = "PoPs Commonly Prioritized Genes") +
  theme(axis.text.x = element_text(angle = 75, vjust = 0.5)) +
  ylim(0,1)
plot_topgenes

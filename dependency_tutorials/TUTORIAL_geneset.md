
============================================================================
Construct gene sets for tissues and connect them to GWAS traits
============================================================================

- /n/groups/price/kushal/GeneRank/code/Tissue_Trait_specific
  - ABC_gene_score_by_tissue.R : ABC (intergenic) and ABC (genic+intergenic) top genes for different tissues
  - EDS_gene_score_by_tissue.R : Probabilistic EDS score top genes for different tissues
  - SEG_gene_score_by_tissue.R : SEG top genes for different tissues
  - PPI_string_tissue_genes.R : construct PPI-enhancer driven genes as combination of the ABC, EDS and SEG genes for each tissue.

The resulting gene sets are saved in /n/groups/price/kushal/GeneRank/data/Gene_Scores/tissue
    - arranged by tissue: blood, brain, adipose, gut, lung, liver, heart

Gene set -> (ABC+Roadmap) S2G strategy -> S-LDSC -> postprocess 

The postprocessed output are in /n/groups/price/kushal/GeneRank/output/Tissue_results_Jul12

    - Postprocessed results to (Trait x Tissue-specific gene sets rankings):
      - /n/groups/price/kushal/GeneRank/code/Tissue_Trait_specific/trait_specific_scores.R
      Output: /n/groups/price/kushal/GeneRank/output/Trait_Specific_Genes

single cell L2 level modules from kushal/singlecellLDSC project are processed similarly and the postprocessed files are at:

- /n/groups/price/kushal/GeneRank/output/scRResults_Jul12

  - Postprocessed results to (Trait x	singlecell-specific	gene sets rankings):
    - /n/groups/price/kushal/GeneRank/code/Tissue_Trait_specific/trait_specific_scores_sc.R
    Output: /n/groups/price/kushal/GeneRank/output/Trait_Specific_Genes_singlecell


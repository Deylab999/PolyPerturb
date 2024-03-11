###### We present a tutorial on MAGMA  and how to get gene prioritization for different traits  ###################

Installation:
Download and Install MAGMA. Please download the static Linux version (2nd option) from here:
https://ctg.cncr.nl/software/magma [RUN ALL MAGMA codes from kkd14 as that is where I have installed the MAGMA static version]

Next run ./magma to see if the function is working.

MAGMA analysis pipeline:

Step 1: Define the S2G strategy for MAGMA. We take a 50kb window around the gene to compute the averaging in MAGMA. 
Code: 1_magma_annotate.sh . We use the gene locations, and EUR PLINK files available with the MAGMA software.

Step 2: Define the gene-level enrichment analysis in MAGMA
Code: 2_gene_level_magma.sh. This script uses the additional script `gene_magma.sh` and `sumstats_to_pfile.R` that converts the 
sumstats Z scores to p-values and then uses S2G strategy from Step 1 to compute gene level enrichment.

Step 3: Define gene set level enrichment analysis in MAGMA
For any new class of gene sets, make the first column the column of genes (Entrez) and the other columns representing continuous Z scores
representing the candidate membership of the gene in each score. 
Code: 3_genecovar_level_magma.sh. Builds on the gene level enrichments from Step 2 and additional continuous covariates created using the 
continuous Z scores fo the different classes.

2 gene sets we investigated:

DEPICT pathways, Mouse Phenotype, KEGG, PPI : /n/groups/price/kushal/GeneRank/data/GPL570-GPL96-GPL1261-GPL1355TermGeneZScores-MGI_MF_CC_RT_IW_BP_KEGG_z_z_Entrez.txt

single cell CTS gene sets : /n/groups/price/kushal/GeneRank/data/Karthik_Kushal_scLDSC_z_z_Entrez.txt

####################  Postprocessing the MAGMA gene covariates results  ##########################################

code: magma_postprocess_covar.R (we postprocess the MAGMA covariate results fore enrichment analysis of the gene scores)

We save the postprocessed significant traits and gene sets pairs in

/n/groups/price/kushal/GeneRank/data/MAGMA_Karthik_Kushal_scLDSC_postprocess.txt

/n/groups/price/kushal/GeneRank/data/MAGMA_MGI_MF_CC_RT_IW_BP_KEGG_CELL_postprocess.txt



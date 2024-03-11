

#######################  Katie Gieger-Schuller, Basak et al  (Aviv Regev data)  #########################

Data location: /n/groups/price/karthik/scdata/karthik/
1032 guides ~ 20k genes (mouse)
Guides are chosen as genes that are associated with ubiquitin ligase action pathway (some of them are directly linked 
to ubiquitin ligase protein, others are functionally associated with it)

Codebase: /n/groups/price/kushal/GeneRank/Gene_Programs/Basak2022

Codes: 

basak_perturb_seq.R
 Takes input the original full Perturb-seq data 
 Input:/n/groups/price/karthik/scdata/karthik/adata-hash-features_singlets_SingleKO_08122020_PerGENE.h5ad
 Creates metadata output for guides and individuals and shows how such data can be read and explored in R

basak_sigbeta_genescores.R
 Takes in the MIMOSCA beta and p-values matrices 
	Input: /n/groups/price/karthik/scdata/karthik/ME_LMBetaCoefsALL.csv ( beta coefficients ) 
	       /n/groups/price/karthik/scdata/karthik/ME_LMPValuesALL.csv ( p values )
	Creates two types of gene lists for each guide: 
	(i) Genes with p-value < 0.05 and also among top 100 interms of magnitude of the beta coefficient. This gene score
	    is saved under
	    /n/groups/price/kushal/GeneRank/Gene_Programs/Basak2022/Gene_Scores/Basak2022_Betasigs
	(ii) Genes with p-value < 0.05
	     Output saved under
	     /n/groups/price/kushal/GeneRank/Gene_Programs/Basak2022/Gene_Scores/Basak2022_Betasigs2/

basak_PMD_genescores.R
  Takes in the MIMOSCA beta and p-values matrices 
  	Input: /n/groups/price/karthik/scdata/karthik/ME_LMBetaCoefsALL.csv ( beta coefficients )
               /n/groups/price/karthik/scdata/karthik/ME_LMPValuesALL.csv ( p values )
	First the code assigns all beta values with p > 0.05 to 0, then performs a PMD matrix factorization on the 
	coefficient data matrix and focus it only on the factor programs for genes. (6685 genes, not 1032 guides)
	For each factor program, we only consider genes that have factor loading > 0.10.
	Also considered a bigger threshold 0.25, but 0.10 was better
	Also tried K=50
	The gene scores are saved in 
	/n/groups/price/kushal/GeneRank/Gene_Programs/Basak2022/Gene_Scores/Basak2022_PMD_K25_thresh_10

basak_fgsea_immune.R
  Takes in the gene scores from 
  	Input: /n/groups/price/kushal/GeneRank/Gene_Programs/Basak2022/Gene_Scores/Basak2022_Betasigs
  	and the PolyGene scores
  	Input: /n/groups/price/kushal/GeneRank/Genes_by_X/PolyGene_plus/PolyGene_PPI_PoPS_justMAGMA_Mean.txt
  	also you can try using PoPS and MAGMA 
  	and peform PolyGene program score for the Perturb-seq programs and save as
  	Output: /n/groups/price/kushal/GeneRank/Gene_Programs/Basak2022/Basak2022_Betasigs_immune.rda

basak_fgsea_PMD_immune.R
   Takes in the gene scores from	
        Input: /n/groups/price/kushal/GeneRank/Gene_Programs/Basak2022/Gene_Scores/Basak2022_PMD_K25_thresh_10
        and the PolyGene scores
        Input: /n/groups/price/kushal/GeneRank/Genes_by_X/PolyGene_plus/PolyGene_PPI_PoPS_justMAGMA_Mean.txt
        also you can try using PoPS and MAGMA	
        and peform PolyGene program score for	the Perturb-seq	programs and save as
        Output: /n/groups/price/kushal/GeneRank/Gene_Programs/Basak2022/Basak2022_PMD_K25_thresh_10_immune.rda


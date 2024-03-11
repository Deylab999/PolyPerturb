This is a README demonstrating the codes and data relevant to the PolyGene paper and project 
(Dey, Jagadeesh et al 2021 in prep). The README is structured based on different gene score
prioritization approaches as explined below

1. MAGMA/: directory containing scripts to perform MAGMA gene-level analysis (de Leeuw et al 2015 PLoS Comp biol). 
   	   g1000_eur*: The LD ref panel from 1000 Genomes distributed with the MAGMA software
	   Download and Install MAGMA. Please download the static Linux version (2nd option) from here:
	   https://ctg.cncr.nl/software/magma [RUN ALL MAGMA codes from kkd14 as that is where I have installed 
	   the MAGMA static version]
	   cd /n/groups/price/kushal/GeneRank/MAGMA.
	   Next run ./magma to see if the function is working.
	   NCBI37.3.gene.loc: gene nomenclature file that comes with MAGMA.
	   codes/: subdirectory containing code to run MAGMA gene-level analyses. 
	   	   1_magma_annotate.sh: Create a set of loci mapped to a gene for different gene windows:
		   			gene +- 0kb, 5kb, 10kb and 50kb. The script saves loci output in LOC_CELL
					namely LOC_CELL_{0kb, 5kb, 10kb, 50kb}. The output files in these folders
					comprise of 
					step1.genes.annot files that ae of following form:
					..................................................
					gene_entrez chr:start:end  rsids-mapped-to-gene
					79501     1:69091:70008	rs140739101	 rs200505207	rs527512746 rs200676709
					..................................................
		   1_magma_annotate_abc.R(.sh): Create a list of rsids mapped to each gene based on Activity-By-Contact
		   				approach for different tissues. Uses GRanges to map rsids from 1000G bim
						files used in S-LDSC and gene nomenclature file above: NCBI37.3.gene.loc
						to assign SNPs matched to genes.
						The results are saved in folder ~/MAGMA/ABC_LOC, with subdirectories ABC_LOC_{Tissue}. 
		   1_magma_annotate_abcUroadmap.R(.sh): Same as	above function but instead of ABC uses  Roadmap.
                   				     	    	  	   The results are saved in ~/MAGMA/Roadmap_LOC.
		   1_magma_annotate_abcUroadmap.R(.sh): Same as above function but instead of ABC uses union of ABC and Roadmap.
		   					The results are saved in ~/MAGMA/ABC_Roadmap_LOC.  
							
                   2_gene_level_magma.sh: MAGMA gene-level analysis based on loc file:
		   			  Input: /n/groups/price/kushal/GeneRank/MAGMA/LOC_CELL_0kb
					  Output: /n/groups/price/kushal/GeneRank/MAGMA/GENE_CELL_0
				          This code calls the workhorse script: gene_magma.sh.
					  The gene_magma.sh converts .sumstats file to MAGMA input .sumstats.pfile file 
					  and then use that a input to run MAGMA, and save output files in Output.
					  .genes.out, .genes.raw files. .genes.raw file contains MAGMA gene-level LD info.
					  .genes.out contains gene-level prioritization scores.
					  .........................................................................
					  GENE       CHR      START       STOP  NSNPS  NPARAM       N        ZSTAT            P
					  148398       1     859993     879961      5       2  207661      -1.4819      0.93081
					  26155        1     879583     894679     11       3  207661     -0.31466      0.62349
					  ..........................................................................
					  The MAGM P nd ZSTAT scores are used for downstream analysis. 
		  2_gene_level_magma_enhancer.sh: same as 2_gene_level_magma.sh but instead of 0kb or window-based MAGMA, we use the 
		  				  ABC-U-Roaadmap. 
						  Input: /n/groups/price/kushal/GeneRank/MAGMA/ABC_Roadmap_LOC/ABC_Roadmap_LOC_*
						  Output: /n/groups/price/kushal/GeneRank/MAGMA/ABC_Roadmap_GENE_CELL/ABC_Roadmap_*
		  2_gene_level_exopro.sh: same as 2_gene_level_magma.sh	but instead of 0kb or window-based MAGMA, we use the
		  			  ExoPro from Steven's where the loc files were provided by Steven Gazal.
					  Input: /n/groups/price/kushal/GeneRank/MAGMA/LOC_EXOPRO
					  Output: /n/groups/price/kushal/GeneRank/MAGMA/EXOPRO_GENE_CELL 
		  3_magma_process_genes_window.R: Aggregate MAGMA (window based) Z and P statistics from different traits into a matrix 
		  				  Input: /n/groups/price/kushal/GeneRank/MAGMA/GENE_CELL_*
						  Output: /n/groups/price/kushal/GeneRank/MAGMA/output/MAGMA_v108_*_ZSTAT
						  	  /n/groups/price/kushal/GeneRank/MAGMA/output/MAGMA_v108_*_PSTAT
		  3_magma_process_genes_abc_roadmap.R: Same as above but using MAGMA ABC-U-Roadmap instead of window-based.
		  				      Input: /n/groups/price/kushal/GeneRank/MAGMA/ABC_Roadmap_GENE_CELL/
						      Output: /n/groups/price/kushal/GeneRank/MAGMA/output/MAGMA_v108_*{ZSTAT, PSTAT}
		  
		  All MAGMA files are copied in the folder: 
		      /n/groups/price/kushal/GeneRank/Genes_by_X/MAGMA-v108/MAGMA_v108_GENE_{0, 5, 10, 50}_{PSTAT, ZSTAT}.txt
		      /n/groups/price/kushal/GeneRank/Genes_by_X/MAGMA-v108-ABCURoadmap/MAGMA_v108_ABC_Roadmap_{Tissue}_{PSTAT, ZSTAT}.txt
		      Linear combination of window and enhancer-gene MAGMA in /n/groups/price/kushal/GeneRank/Genes_by_X/extras/MAGMA-v108_combo_ABC_try
		      created by the code ~/MAGMA/codes/combo_magma_window_abcroad.R

2. TWAS/: directory comprising of codes to generate gene prioritization for a trait based on TWAS associations 
   	  fusion-twas-master: codebase for performing TWAS analysis downloaded from FUSION/TWAS webpage (Gusev et al).
	  TWAS-Hub: a hub of TWAS association dat files downloaded from TWAS/FUSION webpage for different traits and cell types combos.
	            .dat files for each trait compises of TWAS Z scores for genes for different cell and types aggregated together.
          twas_traits.txt: names of all TWAS coded traits 
	  TWAS_data_process_for_trait.R: For each trait, a matrix of TWAS Z scores for genes and tissue, saved as .txt files in the folder
	  				 /n/groups/price/kushal/GeneRank/TWAS/TWAS-Hub/Zmats
					 Same done also for cis h2 as estimated by TWAS GCTA and saved in 
					 /n/groups/price/kushal/GeneRank/TWAS/TWAS-Hub/HSQ
					 However H2 does not depend on trait, so the matrix will be same for all traits, will be different 
					 only for different genes and tissues.
	  TWAS_postprocess.R: takes max TWAS abs z score across cell types for each gene and trait pair. This Gene-trait score is saved 
	  		      in the path  /n/groups/price/kushal/GeneRank/TWAS/TWAS-Hub.
			      This Genes_by_X file is then moved to /n/groups/price/kushal/GeneRank/Genes_by_X/TWAS.
			      Genes_by_X_TWAS_abs_Z.txt 
			      Genes_by_X_TWAS_HSQ.txt
	  /n/groups/price/kushal/GeneRank/Genes_by_X/TWAS/MAGMA_TWAS_matched_traits.txt: matched traits between MAGMA and TWAS.
	  download_twas.sh: code to download all TWAS association statistics from online.
	  TWAS_models_traits.txt: description of the TWAS data paths along with information about the trait. 

	  All TWAS gene-trait files are coped in the folder:
	  /n/groups/price/kushal/GeneRank/Genes_by_X/TWAS/Genes_by_X_*

3. ACAT/: Instead of MAGMA aggregatio method, use ACAT-V aggregtion strategy for SNPs mapped to the gene by gene windows
          and enhancer-gene links.
	  codes
	  1_acatv_genescore.R(.sh): Similar to MAGMA, computes gene score using input
	  			    /n/groups/price/kushal/GeneRank/ACAT/LOC_CELL/LOC_CELL_{0, 50, ..., ALL}
				    for gene-windows (0, 5, 10, 50kb and ABC+Roadmap in different tissues)
				    (data tranported from MAGMA loc files)
				    and output the results in 
				    /n/groups/price/kushal/GeneRank/ACAT/GENE_CELL
	  2_combine_atac_gscores.R: code to combine the results from previouss step to create a gene-trait score matrix 
	  			    for different linking strategies.
				    Input: /n/groups/price/kushal/GeneRank/ACAT/GENE_CELL
				    Output: /n/groups/price/kushal/GeneRank/Genes_by_X/ACAT/*_{TSTAT, PSTAT}
	  NCBI37.3.gene.loc : gene nomenclature file as MAGMA for comparison

	  All ACAT gene-trait files are saved in the folder:
	  /n/groups/price/kushal/GeneRank/Genes_by_X/ACAT/*_{TSTAT, PSTAT}

4. Gazal-h2: /n/groups/price/kushal/GeneRank/Gazalh2/beta2: beta2 estimates fro PolyFun+SuSIE (Omer's fine-mapping approach)
   	     codes/: cd /n/groups/price/kushal/GeneRank/Gazalh2/codes
	     	     create_sgscore_abc_roadmap.R: creates sgscore files mapping SNPs to Genes as proposed by Gazal:  
		     				   SGscores are saved in /n/groups/price/kushal/GeneRank/Gazalh2/ABC_Roadmap.
						   A sgscore file has the following structure

						   ```
						   SNP     GENE    SGSCORE INFO
						   rs7089857       ZMYND11 1       ABC-U-Roadmap
						   rs11253280      ZMYND11 1       ABC-U-Roadmap
						   rs4881493       ZMYND11 1       ABC-U-Roadmap
						   rs4881494       ZMYND11 1       ABC-U-Roadmap
						   ```
	             gazal_geneh2g_calc_abc_roadmap.R(.sh): scripts to perform sum of beta2 from 
		     					    Input: /n/groups/price/kushal/GeneRank/Gazalh2/beta2
							    to a gene using the SGscore files in 
							    SGscores are saved in /n/groups/price/kushal/GeneRank/Gazalh2/ABC_Roadmap
							    to score a gene. The gene scores are then saved in 
							    Output: /n/groups/price/kushal/GeneRank/Gazalh2/h2g_generank
							    	    (results saved forr different ABC+Roadmap enhancer-gene links to compare with cS2G)
		    gazalh2g_annotations.R(.sh): creates SNP level annotations for different SGScores. 
		    gazal_process_genes.R: Combines Gazal-h2 results across multiple traits to generate a gene-trait matrix in 
		    			   /n/groups/price/kushal/GeneRank/Gazalh2/processed/Gazalh2_*

		All Gazal-h2 files are stored in: /n/groups/price/kushal/GeneRank/Genes_by_X/Gazal-h2-ABCURoadmap 
		This folder also includes the file Gazalh2_ABC_Roadmap_cS2G.txt based on the cS2G that Steven uses.


5. PoPS/: directory containing scripts to run the PoPS (Weeks et al medRxiv 2020) approach.
   	  codes/: Codebase to run PoPS predictive scoring of genes.
	  pops.sample_features.magma.py: Feature pre-selection step in PoPS performed using MAGMA-0kb scores from 
	  				 /n/groups/price/kushal/GeneRank/MAGMA/GENE_CELL_0 created in `MAGMA` step above.
					 We use 57k PoPS features from file:
					 /n/groups/price/kushal/GeneRank/POPS/supplemental_data/PoPS.features.txt.gz
					 gene loci: /n/groups/price/kushal/extras/NCBI37.3.ensembl.gene.loc (same as MAGMA)
					 Output: /n/groups/price/kushal/GeneRank/POPS/features/MAGMA_0kb
	pops.sample_features.atacv.py: similar file as pops.sample_features.magma.py, but modified to truncate the TSTAT scores from ACAT-V as
				       they can be crazy big. (-80 to 80). Also the input files are slightly different. 
				       /n/groups/price/kushal/GeneRank/ACAT/GENE_CELL/GENE_CELL_0kb
	pops.sample_features.magma_saige.py: similar to pops.sample_features.magma.py and pops.sample_features.atacv.py, more similar to latter
					     in terms of input structure.
					     /n/groups/price/kushal/GeneRank/POPS/supplemental_data/MAGMA_SAIGE (for .gscore files)
					     /n/groups/price/kushal/GeneRank/POPS/supplemental_data/MAGMA_SAIGE_eqwtd (for .gscore files; this is what we use)
					     We use the same gene-level LD file as MAGMA from 
					     /n/groups/price/kushal/GeneRank/MAGMA/GENE_CELL_0/
					     as MAGMA-SAIGE is a linear combination of MAGMA and SAIGE, and SAIGE gene-level LD is 0, and are 
					     mostly uncorrelated with MAGMA. 
	features/: directory where all pre-selected features from different PoPS model fits are	saved in per chr basis.
	pooled_features: directory where all pre-selected features from different PoPS model fits are saved (aggregating over features/ files)
	pops.combine_chunk_features.R: combines per chr  pre-selected features in features/ to get a single file per trait in pooled_features/

	pops.predict_scores.magma.py: PoPS predictive model for MAGMA using pre-selected features from pooled_features/  
   	    			      and PoPS features and a set of control features (pe-designated)
				      /n/groups/price/kushal/GeneRank/POPS/pooled_features/MAGMA_0kb/*
				      /n/groups/price/kushal/GeneRank/POPS/supplemental_data/PoPS.features.txt.gz
				      /n/groups/price/kushal/GeneRank/POPS/supplemental_data/control.features
				      Output: /n/groups/price/kushal/GeneRank/POPS/predicted_scores/ (per chr predicted scores)
	pops.predict_scores.atacv.py: Same as pops.predict_scores.magma.py but using ACAT-V, and so main difference is different
				      gene-LD from MAGMA and ACAT gene scores; 
				      gene_result: /n/groups/price/kushal/GeneRank/ACAT/GENE_CELL/GENE_CELL_0kb
				      gene_raw: /n/groups/price/kushal/GeneRank/MAGMA/GENE_CELL_0
	pops.predict_scores.magma_saige.py: Same as pops.predict_scores.magma.py but using MAGMA-SAIGE instead of MAGMA. 
	pops.munge_predicted_scores.R: combine the predicted PoPS scores in /n/groups/price/kushal/GeneRank/POPS/predicted_scores/ to
				       generate gene-trait matrix in /n/groups/price/kushal/GeneRank/Genes_by_X/POPS-plus.
	You can tune the scripts to run PoPS on ABC+Roadmap, I found that this did not work very well compared to standard 0kb, you can check
	the results in POPs_MAGMA_ABC_Roadmap_ALL.txt in /n/groups/price/kushal/GeneRank/Genes_by_X/POPS-plus
	pops.qmatch_magma.R: Match the quantiles of the PoPs scores with MAGMA scores. 
	
	All PoPS gene-trait scores are available at: /n/groups/price/kushal/GeneRank/Genes_by_X/POPS-plus
	Particular relevance 
	/n/groups/price/kushal/GeneRank/Genes_by_X/POPS-plus/POPS_0kb_qmatched_MAGMA.txt (PoPS on MAGMA).
	/n/groups/price/kushal/GeneRank/Genes_by_X/POPS-plus/POPs_MAGMA_SAIGE_eqwtd_0kb_qmatched_MAGMA.txt (PoPS on MAGMA-S).


6. SAIGE_GENE/: UKBB-200KWES-SAIGE-GENE-03222021-PriceLab.txt
   		SAIGE-GENE results delivered by Wei Zhou and others from Lee lab.
		SAIGE: scriptss to perform SAIGE analysis on traits.
		All SAIGE-GENE scores are saved in: /n/groups/price/kushal/GeneRank/Genes_by_X/SAIGE_GENE/SAIGEGENE_SKATO_*{ZSTAT,PSTAT}
		codes/:
				process_saige.R: process SAIGE-GENE scores for {Burden, SKAT, SKAT-O} for different MAF cut-offs. 
				corr_saige.R: correlation between SAIGE-GENE scores at different cut-offs.
				process_magma_saige.R: Generate two types of MAGMA+SAIGE scores - 1) average of MAGMA and SAIGE-GENE at three MAF cutoffs, 
				       and 2) average of MAGMA (0.5) and flat average of SAIGE-GENE (0.167 each) (MAGMA_SAIGE_eqwtd),
				        all scores were quantile normalized with non-negative MAGMA. 
	        		create_magma_saige_gscore.R: create MAGMA+SAIGE (eqwtd) gscore files for running PoPS.
				cor_magma_saige.R: correlation of MAGMA and the three SAIGE-GENE scores.
				control_gene_scores.R: generate control gene scores by permuting MAGMA genes. The control gene scores are saved in
				       /n/groups/price/kushal/GeneRank/Genes_by_X/Control/
				saige_gene_qmatched_magma.R: SAIGE-GENE scores at different cut-offs quantile matched to MAGMA gene score.
	
7. PolyGene/:  codes/: codebase to perform PolyGene gene prioritization analysis.
   	       ppi_genenet_RWR.R: compute PPI-Net when input edgelist graph and a gene score file provided as seed genes. 
	       ppi_string_RWR.R: compute PPI-Net when input gene score as seed and using latest STRING network as underlying network.
	       create_ppi_string_scores.R: Creates 3 PPI-Net scores corresponding to STRING network and 3 types of seed genes
	       				   1) combination of PoPS and MAGMA
					   2) just MAGMA
					   3) just PoPS
					   The output are saved in 
					   /n/groups/price/kushal/GeneRank/Genes_by_X/Network_plus/PPI_STRING_{MAGMA_PoPS, MAGMA, PoPS}.txt
	       create_ppi_coexpressdb_scores.R: Creates 3 PPI-Net scores corresponding to STRING network and 3 types of seed genes
	       			           1) combination of PoPS and MAGMA
                                           2) just MAGMA
                                           3) just PoPS
					   and 3 types of edgelist edge weight cut-offs (0.5, 0.75, 0.25)
					   The output are saved in 
					   PPI_CoExpressDB_MAGMA_{MAGMA_PoPS, MAGMA, PoPS}{_NA, _25, _75}.txt
	       create_ppi_{depmap, humap, dorothea}_scores.R: same as above two but using gene networks from DepMap, HuMap and DoRothEA.
	       combine_ppi_networks.R: various linear combos of different PPI-Nets trained on PoPS+MAGMA
	       			       most relevant combination output:
	       			       /n/groups/price/kushal/GeneRank/Genes_by_X/Network_plus/Combo_PPI_CoexpressDB_STRING_WtdMean.txt
				       which is 0.2*coexpressdb + 0.8*string
	       combine_ppi_networks_justmagma.R: various linear combos of different PPI-Nets trained on MAGMA
	       					 which is 0.2*coexpressdb + 0.8*string
						 /n/groups/price/kushal/GeneRank/Genes_by_X/Network_plus/Combo_PPI_CoexpressDB_STRING_WtdMean_justMAGMA.txt
						 Combo_PPI_CoexpressDB_STRING_WtdMean_justMAGMA.txt copied also as Net_MAGMA_STRING_CoexprressDB.txt
	       create_ppi_coexpressdb_magmasaige_scores.R: CoExpressDB PPI-Net (edgelist 0.5 and up) with seed genes determined by 
	       						   MAGMA_v108_SAIGE_GENE_ZSTAT_equiweighted_May2021.txt (MAGMA+SAIGE)
							   (PPI-CoexpressDB-MAGMA-SAIGE-pooled), output saved as
							   /n/groups/price/kushal/GeneRank/Genes_by_X/Network_plus/PPI_CoExpressDB_MAGMA_SAIGE_pooled.txt
							   +
							   PPI-Net on a combination of {MAGMA, SAIGE-1e-02, SAIGE-1e-03, SAIGE-1e-04}
							   (PPI-CoexpressDB-MAGMA-SAIGE-indiv), output saved as
							   /n/groups/price/kushal/GeneRank/Genes_by_X/Network_plus/PPI_CoExpressDB_MAGMA_SAIGE_indiv.txt
	       create_ppi_string_magmasaige_scores.R: same as above function but using the STRING network
	       					      outputs saved as:
						      /n/groups/price/kushal/GeneRank/Genes_by_X/Network_plus/PPI_STRING_MAGMA_SAIGE_pooled.txt
						      /n/groups/price/kushal/GeneRank/Genes_by_X/Network_plus/PPI_STRING_MAGMA_SAIGE_indiv.txt
	       combine_ppi_networks_magmasaige.R: similar weighted combination of CoexpressDB and STRING PPI-Net scores for pooled and indiv combos
	      					 separately. The two corresponding output are 
						 /n/groups/price/kushal/GeneRank/Genes_by_X/Network_plus/Net_MAGMA_SAIGE_STRING_CoexprressDB_pooled.txt
						 /n/groups/price/kushal/GeneRank/Genes_by_X/Network_plus/Net_MAGMA_SAIGE_STRING_CoexprressDB_indiv.txt

	       The three PPI-Net files from above scripts used for downstream polygene calculation are 
	       /n/groups/price/kushal/GeneRank/Genes_by_X/Network_plus/
	       Net_MAGMA_STRING_CoexprressDB.txt, 
	       Net_MAGMA_SAIGE_STRING_CoexprressDB_pooled.txt, 
	       Net_MAGMA_SAIGE_STRING_CoexprressDB_indiv.txt

	       create_polygene_scores.R: Take max and mean of the PoPs score and PPI-Net score generated from above steps. 
	       				 Quantile matched with MAGMA (+ve)
					 Input:
					 /n/groups/price/kushal/GeneRank/Genes_by_X/POPS-plus/POPs_MAGMA_0kb_qmatched_MAGMA.txt
					 /n/groups/price/kushal/GeneRank/Genes_by_X/Network_plus/Net_MAGMA_STRING_CoexprressDB.txt
					 Output:
	       				 /n/groups/price/kushal/GeneRank/Genes_by_X/PolyGene_plus/PolyGene_PPI_PoPS_justMAGMA_Mean.txt
					 /n/groups/price/kushal/GeneRank/Genes_by_X/PolyGene_plus/PolyGene_PPI_PoPS_justMAGMA_Max.txt

   	       create_polygene_magma_saige.R: Take max and mean of the PoPs-S score and PPI-Net-S score generated from above steps.
	       				      Quantile matched with MAGMA.
					      Input:
					      /n/groups/price/kushal/GeneRank/Genes_by_X/POPS-plus/POPs_MAGMA_SAIGE_eqwtd_0kb_qmatched_MAGMA.txt
					      /n/groups/price/kushal/GeneRank/Genes_by_X/Network_plus/Net_MAGMA_SAIGE_STRING_CoexprressDB_{pooled, indiv}.txt
					      Output:
					      /n/groups/price/kushal/GeneRank/Genes_by_X/PolyGene_plus/PolyGene_PPI_PoPS_justMAGMA_SAIGE_pooled.txt
					      /n/groups/price/kushal/GeneRank/Genes_by_X/PolyGene_plus/PolyGene_PPI_PoPS_justMAGMA_SAIGE_indiv.txt
		Benchmarker: A set of codes to geenerate predictive PolyGene score and perform comparative Benchmarker (100kb) and Benchmarker (Enhancer-gene)
			     of the PolyGene-predict and the PoPs scores for different traits (PoPs seems to outperform PolyGene in this comparison).
			     
8. Metric_Genes/: a directory containing all relevant metric gene sets and codebase to perform enrichment in metric gene sets.
   		  data/: subdirectory containing all relevant metric gene sets
		  	 essentiality: collection of all essential genes
			 mendelian: different Mendelian gene sets
			 drug_targets: All drug target gene sets 
		  code/: create_drug_targets.R: crate drug targets for all disease, immune related diseases and brain related diseases.
		  	 test_*.R: test enrichment of different gene scores for various metric gene sets relevant to the trait.

				    
	  
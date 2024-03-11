

#######################################  PPI-CoexpressDB-MAGMA-SAAIGE-pooled   #####################################################################

library(preprocessCore)
library(dnet)
library(igraph)
source("/n/groups/price/kushal/Enhancer_MasterReg/Dey_Enhancer_MasterReg/GSSG/code/calc_PPI_scores/ppi_genenet_RWR.R")

gg_edgelist = read.table("/n/groups/price/kushal/extras/edgelists/coexpressdb_edgelist_50", header=F)
tabb1 = read.table("/n/groups/price/kushal/GeneRank/Genes_by_X/MAGMA-v108/MAGMA_v108_SAIGE_GENE_ZSTAT_equiweighted_May2021.txt")
tabb1[is.na(tabb1)] = 0
tabb2 = read.table("/n/groups/price/kushal/GeneRank/Genes_by_X/MAGMA-v108/MAGMA_v108_GENE_0_ZSTAT.txt")
tabb2[is.na(tabb2)] = 0

common_genes = Reduce(intersect, list(rownames(tabb1), rownames(tabb2)))
common_traits = Reduce(intersect, list(colnames(tabb1), colnames(tabb2)))

ccmat = c()
for(traitname in common_traits){

  dff = cbind(tabb1[match(common_genes, rownames(tabb1)),traitname],
              tabb2[match(common_genes, rownames(tabb2)),traitname])
  dff[dff>25] = 25
  dff[dff < 0] = 0
  rownames(dff) = common_genes
  qdff = normalize.quantiles.use.target(dff,dff[,2],copy=TRUE,subset=NULL)
  rownames(qdff) = common_genes
  set.seed(1)
  tt = ppi_genenet_RWR(gg_edgelist, cbind(as.matrix(qdff)[,1], as.matrix(qdff)[,1]), NUM_RUNS = 10)
  vecc = rep(0, length(common_genes))
  vecc[match(intersect(names(tt$score), common_genes),
             common_genes)] = tt$score[match(intersect(names(tt$score), common_genes), names(tt$score))]
  qvecc = as.vector(normalize.quantiles.use.target(as.matrix(vecc), qdff[,2], copy=T, subset=NULL))
  ccmat = cbind(ccmat, qvecc)
  cat("We are at trait:", traitname, "\n")
}

rownames(ccmat) = common_genes
colnames(ccmat) = common_traits
write.table(ccmat, file = paste0("/n/groups/price/kushal/GeneRank/Genes_by_X/Network_plus/PPI_CoExpressDB_MAGMA_SAIGE_pooled.txt"),
            row.names = T, col.names = T, sep = "\t", quote=F)


#######################################  PPI-CoexpressDB-MAGMA-SAIGE-indiv   #####################################################################

library(preprocessCore)
library(dnet)
library(igraph)
source("/n/groups/price/kushal/Enhancer_MasterReg/Dey_Enhancer_MasterReg/GSSG/code/calc_PPI_scores/ppi_genenet_RWR.R")

gg_edgelist = read.table("/n/groups/price/kushal/extras/edgelists/coexpressdb_edgelist_50", header=F)
tabb1 = read.table("/n/groups/price/kushal/GeneRank/Genes_by_X/SAIGE_GENE/SAIGEGENE_SKATO_1E-02_ZSTAT.txt")
tabb2 = read.table("/n/groups/price/kushal/GeneRank/Genes_by_X/SAIGE_GENE/SAIGEGENE_SKATO_1E-03_ZSTAT.txt")
tabb3 = read.table("/n/groups/price/kushal/GeneRank/Genes_by_X/SAIGE_GENE/SAIGEGENE_SKATO_1E-04_ZSTAT.txt")
tabb4 = read.table("/n/groups/price/kushal/GeneRank/Genes_by_X/MAGMA-v108/MAGMA_v108_GENE_0_ZSTAT.txt")
tabb4[is.na(tabb4)] = 0

common_genes = Reduce(intersect, list(rownames(tabb1), rownames(tabb2), rownames(tabb3),
                                      rownames(tabb4)))
colnames(tabb1)[which(colnames(tabb1) == "UKB_460K.Multiple_sclerosis")] = "PASS_Multiple_sclerosis"
colnames(tabb1)[which(colnames(tabb1) == "UKB_460K.disease_T1D")] = "PASS_Type_1_Diabetes"
colnames(tabb1)[which(colnames(tabb1) == "UKB_460K.Hypertension")] = "UKB_460K.disease_HYPERTENSION_DIAGNOSED"
colnames(tabb1)[which(colnames(tabb1) == "UKB_460K.stroke")] = "PASS_Stroke_Malik2018"
colnames(tabb1)[which(colnames(tabb1) == "UKB_460K.ischaemic_stroke")] = "PASS_IschemicStroke_Malik2018"
colnames(tabb1)[which(colnames(tabb1) == "UKB_460K.mental_SCZ")] = "PASS_Schizophrenia_Ruderfer2018"
colnames(tabb1)[which(colnames(tabb1) == "UKB_460K.Bipolar")] = "PASS_BipolarDisorder_Ruderfer2018"
colnames(tabb1)[which(colnames(tabb1) == "UKB_460K.Lupus")] = "PASS_Lupus"
colnames(tabb1)[which(colnames(tabb1) == "UKB_460K.celiac")] = "PASS_Celiac"
colnames(tabb1)[which(colnames(tabb1) == "UKB_460K.IBD")] = "PASS_IBD_deLange2017"
colnames(tabb1)[which(colnames(tabb1) == "UKB_460K.Crohns_disease")] = "PASS_CD_deLange2017"
colnames(tabb1)[which(colnames(tabb1) == "UKB_460K.Ulcerative_colitis")] = "PASS_UC_deLange2017"
colnames(tabb1)[which(colnames(tabb1) == "UKB_460K.Rheumatoid_arthritis")] = "PASS_Rheumatoid_Arthritis"
colnames(tabb1)[which(colnames(tabb1) == "UKB_460K.Anorexia")] = "PASS_Anorexia"
colnames(tabb1)[which(colnames(tabb1) == "UKB_460K.insomnia")] = "PASS_Insomnia_Jansen2019"
colnames(tabb1)[which(colnames(tabb1) == "UKB_460K.Autism")] = "PASS_Autism"
colnames(tabb1)[which(colnames(tabb1) == "UKB_460K.high_cholesterol")] = "UKB_460K.biochemistry_Cholesterol"


colnames(tabb2) = colnames(tabb1)
colnames(tabb3) = colnames(tabb1)
common_traits = intersect(colnames(tabb1), colnames(tabb4))

ccmat = c()
for(traitname in common_traits){

  dff = cbind(tabb1[match(common_genes, rownames(tabb1)),traitname],
              tabb2[match(common_genes, rownames(tabb2)),traitname],
              tabb3[match(common_genes, rownames(tabb3)),traitname],
              tabb4[match(common_genes, rownames(tabb4)),traitname])
  dff[dff>25] = 25
  dff[dff < 0] = 0
  rownames(dff) = common_genes
  qdff = normalize.quantiles.use.target(dff,dff[,4],copy=TRUE,subset=NULL)
  rownames(qdff) = common_genes
  set.seed(1)
  tt = ppi_genenet_RWR(gg_edgelist, as.matrix(qdff), NUM_RUNS = 10)
  vecc = rep(0, length(common_genes))
  vecc[match(intersect(names(tt$score), common_genes),
             common_genes)] = tt$score[match(intersect(names(tt$score), common_genes), names(tt$score))]
  qvecc = as.vector(normalize.quantiles.use.target(as.matrix(vecc), qdff[,2], copy=T, subset=NULL))
  ccmat = cbind(ccmat, qvecc)
  cat("We are at trait:", traitname, "\n")
}

rownames(ccmat) = common_genes
colnames(ccmat) = common_traits
write.table(ccmat, file = paste0("/n/groups/price/kushal/GeneRank/Genes_by_X/Network_plus/PPI_CoExpressDB_MAGMA_SAIGE_indiv.txt"),
            row.names = T, col.names = T, sep = "\t", quote=F)

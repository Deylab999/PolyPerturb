library(preprocessCore)
library(dnet)
library(igraph)
library(STRINGdb)

#######################################  PPI-STRING-MAGMA-PoPS   #####################################################################


source("/n/groups/price/kushal/Enhancer_MasterReg/Dey_Enhancer_MasterReg/GSSG/code/calc_PPI_scores/ppi_string_RWR.R")

tabb1 = read.table("/n/groups/price/kushal/GeneRank/Genes_by_X/POPS-plus/POPS_0kb_qmatched_MAGMA.txt")
tabb2 = read.table("/n/groups/price/kushal/GeneRank/Genes_by_X/MAGMA-v108/MAGMA_v108_GENE_0_ZSTAT.txt")
tabb2[is.na(tabb2)] = 0
colnames(tabb1)[which(colnames(tabb1) == "PASS_Bipolar_Disorder")] = "PASS_BipolarDisorder_Ruderfer2018"
colnames(tabb1)[which(colnames(tabb1) == "PASS_ChildOnsetAsthma_Ferreira2019")] = "UKB_460K.disease_ASTHMA_DIAGNOSED"

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
  tt = ppi_string_RWR(as.matrix(qdff), NUM_RUNS = 5)
  vecc = rep(0, length(common_genes))
  vecc[match(intersect(names(tt$score), common_genes),
             common_genes)] = tt$score[match(intersect(names(tt$score), common_genes), names(tt$score))]
  qvecc = as.vector(normalize.quantiles.use.target(as.matrix(vecc), qdff[,2], copy=T, subset=NULL))
  ccmat = cbind(ccmat, qvecc)
  cat("We are at trait:", traitname, "\n")
}

rownames(ccmat) = common_genes
colnames(ccmat) = common_traits
write.table(ccmat, file = paste0("/n/groups/price/kushal/GeneRank/Genes_by_X/Network_plus/PPI_STRING_MAGMA_PoPS.txt"),
            row.names = T, col.names = T, sep = "\t", quote=F)



#######################################  PPI-STRING-MAGMA   #####################################################################


source("/n/groups/price/kushal/Enhancer_MasterReg/Dey_Enhancer_MasterReg/GSSG/code/calc_PPI_scores/ppi_string_RWR.R")

tabb1 = read.table("/n/groups/price/kushal/GeneRank/Genes_by_X/POPS-plus/POPS_0kb_qmatched_MAGMA.txt")
tabb2 = read.table("/n/groups/price/kushal/GeneRank/Genes_by_X/MAGMA-v108/MAGMA_v108_GENE_0_ZSTAT.txt")
tabb2[is.na(tabb2)] = 0
colnames(tabb1)[which(colnames(tabb1) == "PASS_Bipolar_Disorder")] = "PASS_BipolarDisorder_Ruderfer2018"
colnames(tabb1)[which(colnames(tabb1) == "PASS_ChildOnsetAsthma_Ferreira2019")] = "UKB_460K.disease_ASTHMA_DIAGNOSED"

common_genes = Reduce(intersect, list(rownames(tabb1), rownames(tabb2)))
common_traits = Reduce(intersect, list(colnames(tabb1), colnames(tabb2)))

ccmat = c()
for(traitname in common_traits){

  dff = cbind(tabb2[match(common_genes, rownames(tabb1)),traitname],
              tabb2[match(common_genes, rownames(tabb2)),traitname])
  dff[dff>25] = 25
  dff[dff < 0] = 0
  rownames(dff) = common_genes
  qdff = normalize.quantiles.use.target(dff,dff[,2],copy=TRUE,subset=NULL)
  rownames(qdff) = common_genes
  set.seed(1)
  tt = ppi_string_RWR(as.matrix(qdff), NUM_RUNS = 5)
  vecc = rep(0, length(common_genes))
  vecc[match(intersect(names(tt$score), common_genes),
             common_genes)] = tt$score[match(intersect(names(tt$score), common_genes), names(tt$score))]
  qvecc = as.vector(normalize.quantiles.use.target(as.matrix(vecc), qdff[,2], copy=T, subset=NULL))
  ccmat = cbind(ccmat, qvecc)
  cat("We are at trait:", traitname, "\n")
}

rownames(ccmat) = common_genes
colnames(ccmat) = common_traits
write.table(ccmat, file = paste0("/n/groups/price/kushal/GeneRank/Genes_by_X/Network_plus/PPI_STRING_MAGMA.txt"),
            row.names = T, col.names = T, sep = "\t", quote=F)


#######################################  PPI-STRING-PoPS   #####################################################################


source("/n/groups/price/kushal/Enhancer_MasterReg/Dey_Enhancer_MasterReg/GSSG/code/calc_PPI_scores/ppi_string_RWR.R")

tabb1 = read.table("/n/groups/price/kushal/GeneRank/Genes_by_X/POPS-plus/POPS_0kb_qmatched_MAGMA.txt")
tabb2 = read.table("/n/groups/price/kushal/GeneRank/Genes_by_X/MAGMA-v108/MAGMA_v108_GENE_0_ZSTAT.txt")
tabb2[is.na(tabb2)] = 0
colnames(tabb1)[which(colnames(tabb1) == "PASS_Bipolar_Disorder")] = "PASS_BipolarDisorder_Ruderfer2018"
colnames(tabb1)[which(colnames(tabb1) == "PASS_ChildOnsetAsthma_Ferreira2019")] = "UKB_460K.disease_ASTHMA_DIAGNOSED"

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
  tt = ppi_string_RWR(as.matrix(cbind(qdff[,1], qdff[,1])), NUM_RUNS = 5)
  vecc = rep(0, length(common_genes))
  vecc[match(intersect(names(tt$score), common_genes),
             common_genes)] = tt$score[match(intersect(names(tt$score), common_genes), names(tt$score))]
  qvecc = as.vector(normalize.quantiles.use.target(as.matrix(vecc), qdff[,2], copy=T, subset=NULL))
  ccmat = cbind(ccmat, qvecc)
  cat("We are at trait:", traitname, "\n")
}

rownames(ccmat) = common_genes
colnames(ccmat) = common_traits
write.table(ccmat, file = paste0("/n/groups/price/kushal/GeneRank/Genes_by_X/Network_plus/PPI_STRING_PoPS.txt"),
            row.names = T, col.names = T, sep = "\t", quote=F)


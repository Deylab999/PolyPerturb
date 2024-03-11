library(dorothea)
library(bcellViper)
library(dplyr)
library(viper)
library(preprocessCore)
library(igraph)
library(data.table)

coexpressdb_cor = data.frame(fread("/n/groups/price/kushal/extras/edgelists/coexpressdb_cor", sep = ",", header=T))
matt = coexpressdb_cor[,-1]
entrez_ids = coexpressdb_cor[,1]
gene_nomenclature = read.delim("/n/groups/price/kushal/extras/edgelists/coexpressdb_gene_nomenclature", header=T)
common_ids = intersect(gene_nomenclature$gene_id, entrez_ids)

matt2 = matt[match(common_ids, entrez_ids), match(common_ids, entrez_ids)]
rownames(matt2) = gene_nomenclature$gene[match(common_ids, gene_nomenclature$gene_id)]
colnames(matt2) = rownames(matt2)
diag(matt2) = 0

idx = which(matt2 > 0.75, arr.ind = T)

outdf = cbind.data.frame(rownames(matt2)[idx[,1]], colnames(matt2)[idx[,2]], 1)

write.table(outdf, file = "/n/groups/price/kushal/extras/edgelists/coexpressdb_edgelist_75",
            row.names = F, col.names = F, sep = "\t", quote=F)



#######################################  PPI-CoexpressDB-MAGMA-PoPS   #####################################################################

library(preprocessCore)
library(dnet)
library(igraph)
source("/n/groups/price/kushal/Enhancer_MasterReg/Dey_Enhancer_MasterReg/GSSG/code/calc_PPI_scores/ppi_genenet_RWR.R")

gg_edgelist = read.table("/n/groups/price/kushal/extras/edgelists/coexpressdb_edgelist_75", header=F)
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
  tt = ppi_genenet_RWR(gg_edgelist, as.matrix(qdff), NUM_RUNS = 5)
  vecc = rep(0, length(common_genes))
  vecc[match(intersect(names(tt$score), common_genes),
             common_genes)] = tt$score[match(intersect(names(tt$score), common_genes), names(tt$score))]
  qvecc = as.vector(normalize.quantiles.use.target(as.matrix(vecc), qdff[,2], copy=T, subset=NULL))
  ccmat = cbind(ccmat, qvecc)
  cat("We are at trait:", traitname, "\n")
}

rownames(ccmat) = common_genes
colnames(ccmat) = common_traits
write.table(ccmat, file = paste0("/n/groups/price/kushal/GeneRank/Genes_by_X/Network_plus/PPI_CoExpressDB_MAGMA_PoPS_75.txt"),
            row.names = T, col.names = T, sep = "\t", quote=F)


#######################################  PPI-CoexpressDB-MAGMA   #####################################################################

library(preprocessCore)
library(dnet)
library(igraph)
source("/n/groups/price/kushal/Enhancer_MasterReg/Dey_Enhancer_MasterReg/GSSG/code/calc_PPI_scores/ppi_genenet_RWR.R")

gg_edgelist = read.table("/n/groups/price/kushal/extras/edgelists/coexpressdb_edgelist_75", header=F)
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
  tt = ppi_genenet_RWR(gg_edgelist, as.matrix(qdff), NUM_RUNS = 5)
  vecc = rep(0, length(common_genes))
  vecc[match(intersect(names(tt$score), common_genes),
             common_genes)] = tt$score[match(intersect(names(tt$score), common_genes), names(tt$score))]
  qvecc = as.vector(normalize.quantiles.use.target(as.matrix(vecc), qdff[,2], copy=T, subset=NULL))
  ccmat = cbind(ccmat, qvecc)
  cat("We are at trait:", traitname, "\n")
}

rownames(ccmat) = common_genes
colnames(ccmat) = common_traits
write.table(ccmat, file = paste0("/n/groups/price/kushal/GeneRank/Genes_by_X/Network_plus/PPI_CoExpressDB_MAGMA_75.txt"),
            row.names = T, col.names = T, sep = "\t", quote=F)


#######################################  PPI-CoexpressDB-PoPS   #####################################################################

library(preprocessCore)
library(dnet)
library(igraph)
source("/n/groups/price/kushal/Enhancer_MasterReg/Dey_Enhancer_MasterReg/GSSG/code/calc_PPI_scores/ppi_genenet_RWR.R")

gg_edgelist = read.table("/n/groups/price/kushal/extras/edgelists/coexpressdb_edgelist_75", header=F)
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
  tt = ppi_genenet_RWR(gg_edgelist, as.matrix(cbind(qdff[,1], qdff[,1])), NUM_RUNS = 5)
  vecc = rep(0, length(common_genes))
  vecc[match(intersect(names(tt$score), common_genes),
             common_genes)] = tt$score[match(intersect(names(tt$score), common_genes), names(tt$score))]
  qvecc = as.vector(normalize.quantiles.use.target(as.matrix(vecc), qdff[,2], copy=T, subset=NULL))
  ccmat = cbind(ccmat, qvecc)
  cat("We are at trait:", traitname, "\n")
}

rownames(ccmat) = common_genes
colnames(ccmat) = common_traits
write.table(ccmat, file = paste0("/n/groups/price/kushal/GeneRank/Genes_by_X/Network_plus/PPI_CoExpressDB_PoPS_75.txt"),
            row.names = T, col.names = T, sep = "\t", quote=F)


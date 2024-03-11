options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
print(args)
numchr <- as.numeric(toString(args[1]))

library(dplyr)
library(preprocessCore)
library(igraph)
library(data.table)

source("/n/groups/price/kushal/Enhancer_MasterReg/Dey_Enhancer_MasterReg/GSSG/code/calc_PPI_scores/ppi_genenet_RWR.R")

gg_edgelist = read.table("/n/groups/price/kushal/extras/edgelists/coexpressdb_edgelist_50", header=F)
tabb1 = read.table("/n/groups/price/kushal/GeneRank/Genes_by_X/POPS-plus/POPS_0kb_qmatched_MAGMA.txt")
tabb2 = read.table("/n/groups/price/kushal/GeneRank/Genes_by_X/MAGMA-v108/MAGMA_v108_GENE_0_ZSTAT.txt")
tabb2[is.na(tabb2)] = 0
colnames(tabb1)[which(colnames(tabb1) == "PASS_Bipolar_Disorder")] = "PASS_BipolarDisorder_Ruderfer2018"
colnames(tabb1)[which(colnames(tabb1) == "PASS_ChildOnsetAsthma_Ferreira2019")] = "UKB_460K.disease_ASTHMA_DIAGNOSED"

common_genes = Reduce(intersect, list(rownames(tabb1), rownames(tabb2)))
common_traits = Reduce(intersect, list(colnames(tabb1), colnames(tabb2)))

gene_nomenclature = read.table("/n/groups/price/kushal/extras/Homo_sapiens.GRCh37.gene.nomenclature.txt", header=T)

#numchr=22

common_chr_genes = intersect(gene_nomenclature$gene_name[which(gene_nomenclature$seqnames == numchr)], common_genes)

ccmat = c()
for(traitname in common_traits){

  dff = cbind(tabb2[match(common_genes, rownames(tabb1)),traitname],
              tabb2[match(common_genes, rownames(tabb2)),traitname])
  dff[dff>25] = 25
  dff[dff < 0] = 0
  rownames(dff) = common_genes
  qdff = normalize.quantiles.use.target(dff,dff[,2],copy=TRUE,subset=NULL)
  rownames(qdff) = common_genes
  qdff[match(common_chr_genes, common_genes), ] = 0
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
write.table(ccmat, file = paste0("/n/groups/price/kushal/GeneRank/Genes_by_X/Benchmaker_plus/PPI_CoExpressDB_MAGMA_50_chr", numchr, ".txt"),
            row.names = T, col.names = T, sep = "\t", quote=F)


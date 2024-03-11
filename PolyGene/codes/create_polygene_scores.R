tabb1 = read.table("/n/groups/price/kushal/GeneRank/Genes_by_X/POPS-plus/POPs_MAGMA_0kb_qmatched_MAGMA.txt")
tabb1[tabb1 < 0] = 0
colnames(tabb1)[which(colnames(tabb1) == "PASS_Bipolar_Disorder")] = "PASS_BipolarDisorder_Ruderfer2018"
colnames(tabb1)[which(colnames(tabb1) == "PASS_ChildOnsetAsthma_Ferreira2019")] = "UKB_460K.disease_ASTHMA_DIAGNOSED"

tabb2 = read.table("/n/groups/price/kushal/GeneRank/Genes_by_X/Network_plus/Net_MAGMA_STRING_CoexprressDB.txt")

common_traits = colnames(tabb2)
common_genes = rownames(tabb2)

ccmat1 = c()
ccmat2 = c()

for(traitname in common_traits){
  dff = cbind(tabb1[, traitname], tabb2[, traitname])
  dff[dff > 25] = 25
  vecc = apply(dff, 1, max)
  vecc2 = apply(dff, 1, mean)
  qvecc = as.vector(normalize.quantiles.use.target(as.matrix(vecc), dff[,1], copy=T, subset=NULL))
  qvecc2 = as.vector(normalize.quantiles.use.target(as.matrix(vecc2), dff[,1], copy=T, subset=NULL))
  ccmat1 = cbind(ccmat1, qvecc)
  ccmat2 = cbind(ccmat2, qvecc2)
  cat("We are at trait:", traitname, "\n")
}

rownames(ccmat1) = common_genes
colnames(ccmat1) = common_traits

rownames(ccmat2) = common_genes
colnames(ccmat2) = common_traits


write.table(ccmat1, file = paste0("/n/groups/price/kushal/GeneRank/Genes_by_X/PolyGene_plus/PolyGene_PPI_PoPS_justMAGMA_Max.txt"),
            row.names = T, col.names = T, sep = "\t", quote=F)
write.table(ccmat2, file = paste0("/n/groups/price/kushal/GeneRank/Genes_by_X/PolyGene_plus/PolyGene_PPI_PoPS_justMAGMA_Mean.txt"),
            row.names = T, col.names = T, sep = "\t", quote=F)


################################################### use PPI (MAGMA, PoPS)   ##############################################################


tabb1 = read.table("/n/groups/price/kushal/GeneRank/Genes_by_X/POPS-plus/POPS_0kb_qmatched_MAGMA.txt")
tabb1[tabb1 < 0] = 0
colnames(tabb1)[which(colnames(tabb1) == "PASS_Bipolar_Disorder")] = "PASS_BipolarDisorder_Ruderfer2018"
colnames(tabb1)[which(colnames(tabb1) == "PASS_ChildOnsetAsthma_Ferreira2019")] = "UKB_460K.disease_ASTHMA_DIAGNOSED"

tabb2 = read.table("/n/groups/price/kushal/GeneRank/Genes_by_X/Network_plus/Combo_PPI_CoexpressDB_STRING_WtdMean.txt")

common_traits = colnames(tabb2)
common_genes = rownames(tabb2)

ccmat1 = c()
ccmat2 = c()

for(traitname in common_traits){
  dff = cbind(tabb1[, traitname], tabb2[, traitname])
  dff[dff > 25] = 25
  vecc = apply(dff, 1, max)
  vecc2 = apply(dff, 1, mean)
  qvecc = as.vector(normalize.quantiles.use.target(as.matrix(vecc), dff[,1], copy=T, subset=NULL))
  qvecc2 = as.vector(normalize.quantiles.use.target(as.matrix(vecc2), dff[,1], copy=T, subset=NULL))
  ccmat1 = cbind(ccmat1, qvecc)
  ccmat2 = cbind(ccmat2, qvecc2)
  cat("We are at trait:", traitname, "\n")
}

rownames(ccmat1) = common_genes
colnames(ccmat1) = common_traits

rownames(ccmat2) = common_genes
colnames(ccmat2) = common_traits


write.table(ccmat1, file = paste0("/n/groups/price/kushal/GeneRank/Genes_by_X/PolyGene_plus/PolyGene_PPI_PoPS_Max.txt"),
            row.names = T, col.names = T, sep = "\t", quote=F)
write.table(ccmat2, file = paste0("/n/groups/price/kushal/GeneRank/Genes_by_X/PolyGene_plus/PolyGene_PPI_PoPS_Mean.txt"),
            row.names = T, col.names = T, sep = "\t", quote=F)



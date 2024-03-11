

tabb1 = read.table("/n/groups/price/kushal/GeneRank/Genes_by_X/Network_plus/PPI_CoExpressDB_MAGMA_PoPS.txt")
tabb2 = read.table("/n/groups/price/kushal/GeneRank/Genes_by_X/Network_plus/PPI_Dorothea_MAGMA_PoPS.txt")
tabb3 = read.table("/n/groups/price/kushal/GeneRank/Genes_by_X/Network_plus/PPI_Depmap_MAGMA_PoPS.txt")
tabb4 = read.table("/n/groups/price/kushal/GeneRank/Genes_by_X/Network_plus/PPI_HuMap_MAGMA_PoPS.txt")
tabb5 = read.table("/n/groups/price/kushal/GeneRank/Genes_by_X/Network_plus/PPI_STRING_MAGMA_PoPS.txt")

common_traits = colnames(tabb1)
common_genes = rownames(tabb1)

ccmat1 = c()
ccmat2 = c()
ccmat3 = c()
ccmat4 = c()
ccmat5 = c()
ccmat6 = c()
for(traitname in common_traits){
  dff = cbind(tabb1[, traitname], tabb2[, traitname], tabb3[, traitname], tabb4[, traitname], tabb5[, traitname])
  vecc = apply(dff, 1, max)
  vecc2 = apply(dff, 1, mean)
  vecc3 = apply(dff, 1, function(x) return(2/9*x[1] + 1/9*x[2]+1/9*x[3] + 1.5/9*x[4] + 3.5/9*x[4]))
  vecc4 = apply(dff[, c(1, 5)], 1, max)
  vecc5 = apply(dff[, c(1, 5)], 1, mean)
  vecc6 = apply(dff[, c(1, 5)], 1, function(x) return(0.2*x[1]+0.8*x[2]))

  qvecc = as.vector(normalize.quantiles.use.target(as.matrix(vecc), dff[,1], copy=T, subset=NULL))
  qvecc2 = as.vector(normalize.quantiles.use.target(as.matrix(vecc2), dff[,1], copy=T, subset=NULL))
  qvecc3 = as.vector(normalize.quantiles.use.target(as.matrix(vecc3), dff[,1], copy=T, subset=NULL))
  qvecc4 = as.vector(normalize.quantiles.use.target(as.matrix(vecc4), dff[,1], copy=T, subset=NULL))
  qvecc5 = as.vector(normalize.quantiles.use.target(as.matrix(vecc5), dff[,1], copy=T, subset=NULL))
  qvecc6 = as.vector(normalize.quantiles.use.target(as.matrix(vecc6), dff[,1], copy=T, subset=NULL))

  ccmat1 = cbind(ccmat1, qvecc)
  ccmat2 = cbind(ccmat2, qvecc2)
  ccmat3 = cbind(ccmat3, qvecc3)
  ccmat4 = cbind(ccmat4, qvecc4)
  ccmat5 = cbind(ccmat5, qvecc5)
  ccmat6 = cbind(ccmat6, qvecc6)

  cat("We are at trait:", traitname, "\n")
}

rownames(ccmat1) = common_genes
colnames(ccmat1) = common_traits

rownames(ccmat2) = common_genes
colnames(ccmat2) = common_traits

rownames(ccmat3) = common_genes
colnames(ccmat3) = common_traits

rownames(ccmat4) = common_genes
colnames(ccmat4) = common_traits

rownames(ccmat5) = common_genes
colnames(ccmat5) = common_traits

rownames(ccmat6) = common_genes
colnames(ccmat6) = common_traits

write.table(ccmat1, file = paste0("/n/groups/price/kushal/GeneRank/Genes_by_X/Network_plus/Combo_PPI_All_Max.txt"),
            row.names = T, col.names = T, sep = "\t", quote=F)
write.table(ccmat2, file = paste0("/n/groups/price/kushal/GeneRank/Genes_by_X/Network_plus/Combo_PPI_All_Mean.txt"),
            row.names = T, col.names = T, sep = "\t", quote=F)
write.table(ccmat3, file = paste0("/n/groups/price/kushal/GeneRank/Genes_by_X/Network_plus/Combo_PPI_All_WtdMean.txt"),
            row.names = T, col.names = T, sep = "\t", quote=F)
write.table(ccmat4, file = paste0("/n/groups/price/kushal/GeneRank/Genes_by_X/Network_plus/Combo_PPI_CoexpressDB_STRING_Max.txt"),
            row.names = T, col.names = T, sep = "\t", quote=F)
write.table(ccmat5, file = paste0("/n/groups/price/kushal/GeneRank/Genes_by_X/Network_plus/Combo_PPI_CoexpressDB_STRING_Mean.txt"),
            row.names = T, col.names = T, sep = "\t", quote=F)
write.table(ccmat6, file = paste0("/n/groups/price/kushal/GeneRank/Genes_by_X/Network_plus/Combo_PPI_CoexpressDB_STRING_WtdMean.txt"),
            row.names = T, col.names = T, sep = "\t", quote=F)






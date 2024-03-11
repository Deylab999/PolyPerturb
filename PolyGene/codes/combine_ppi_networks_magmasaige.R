tabb1 = read.table("/n/groups/price/kushal/GeneRank/Genes_by_X/Network_plus/PPI_CoExpressDB_MAGMA_SAIGE_pooled.txt")
tabb2 = read.table("/n/groups/price/kushal/GeneRank/Genes_by_X/Network_plus/PPI_CoExpressDB_MAGMA_SAIGE_indiv.txt")
tabb3 = read.table("/n/groups/price/kushal/GeneRank/Genes_by_X/Network_plus/PPI_STRING_MAGMA_SAIGE_pooled.txt")
tabb4 = read.table("/n/groups/price/kushal/GeneRank/Genes_by_X/Network_plus/PPI_STRING_MAGMA_SAIGE_indiv.txt")

common_traits = colnames(tabb1)
common_genes = rownames(tabb1)

ccmat1 = c()
ccmat2 = c()
ccmat3 = c()
ccmat4 = c()

for(traitname in common_traits){
  dff = cbind(tabb1[, traitname], tabb2[, traitname], tabb3[, traitname], tabb4[, traitname])
  vecc1 = apply(dff[, c(1, 3)], 1, function(x) return(0.2*x[1]+0.8*x[2]))
  vecc2 = apply(dff[, c(2, 4)], 1, function(x) return(0.2*x[1]+0.8*x[2]))

  qvecc1 = as.vector(normalize.quantiles.use.target(as.matrix(vecc1), dff[,1], copy=T, subset=NULL))
  qvecc2 = as.vector(normalize.quantiles.use.target(as.matrix(vecc2), dff[,1], copy=T, subset=NULL))

  ccmat1 = cbind(ccmat1, qvecc1)
  ccmat2 = cbind(ccmat2, qvecc2)
  cat("We are at trait:", traitname, "\n")
}

rownames(ccmat1) = common_genes
colnames(ccmat1) = common_traits

rownames(ccmat2) = common_genes
colnames(ccmat2) = common_traits

write.table(ccmat1, file = paste0("/n/groups/price/kushal/GeneRank/Genes_by_X/Network_plus/Net_MAGMA_SAIGE_STRING_CoexprressDB_pooled.txt"),
            row.names = T, col.names = T, sep = "\t", quote=F)
write.table(ccmat2, file = paste0("/n/groups/price/kushal/GeneRank/Genes_by_X/Network_plus/Net_MAGMA_SAIGE_STRING_CoexprressDB_indiv.txt"),
            row.names = T, col.names = T, sep = "\t", quote=F)


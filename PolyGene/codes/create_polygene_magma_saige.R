
# tabb0 = read.table("/n/groups/price/kushal/GeneRank/Genes_by_X/POPS-plus/POPS_0kb_qmatched_MAGMA.txt")
# tabb0[tabb0 < 0] = 0
# colnames(tabb0)[which(colnames(tabb1) == "PASS_Bipolar_Disorder")] = "PASS_BipolarDisorder_Ruderfer2018"
# colnames(tabb0)[which(colnames(tabb1) == "PASS_ChildOnsetAsthma_Ferreira2019")] = "UKB_460K.disease_ASTHMA_DIAGNOSED"


tabb1 = read.table("/n/groups/price/kushal/GeneRank/Genes_by_X/POPS-plus/POPs_MAGMA_SAIGE_eqwtd_0kb_qmatched_MAGMA.txt")
tabb1[tabb1 < 0] = 0
colnames(tabb1)[which(colnames(tabb1) == "PASS_Bipolar_Disorder")] = "PASS_BipolarDisorder_Ruderfer2018"
colnames(tabb1)[which(colnames(tabb1) == "PASS_ChildOnsetAsthma_Ferreira2019")] = "UKB_460K.disease_ASTHMA_DIAGNOSED"


tabb2 = read.table("/n/groups/price/kushal/GeneRank/Genes_by_X/Network_plus/Net_MAGMA_SAIGE_STRING_CoexprressDB_pooled.txt")
tabb3 = read.table("/n/groups/price/kushal/GeneRank/Genes_by_X/Network_plus/Net_MAGMA_SAIGE_STRING_CoexprressDB_indiv.txt")

common_traits = Reduce(intersect, list(colnames(tabb1), colnames(tabb2),
                                       colnames(tabb3)))
common_genes = Reduce(intersect, list(rownames(tabb1), rownames(tabb2),
                                      rownames(tabb3)))

ccmat1 = c()
ccmat2 = c()

for(traitname in common_traits){
  dff = cbind(tabb1[common_genes, traitname],
              tabb2[common_genes, traitname],
              tabb3[common_genes, traitname])
  dff[dff > 25] = 25
  qdff = normalize.quantiles.use.target(dff,dff[,1],copy=TRUE,subset=NULL)
  vecc = apply(qdff[, c(1, 2)], 1, mean)
  vecc2 = apply(qdff[, c(1, 3)], 1, mean)
  qvecc = as.vector(normalize.quantiles.use.target(as.matrix(vecc), qdff[,1], copy=T, subset=NULL))
  qvecc2 = as.vector(normalize.quantiles.use.target(as.matrix(vecc2), qdff[,1], copy=T, subset=NULL))
  ccmat1 = cbind(ccmat1, qvecc)
  ccmat2 = cbind(ccmat2, qvecc2)
  cat("We are at trait:", traitname, "\n")
}

rownames(ccmat1) = common_genes
colnames(ccmat1) = common_traits

rownames(ccmat2) = common_genes
colnames(ccmat2) = common_traits

write.table(ccmat1, file = paste0("/n/groups/price/kushal/GeneRank/Genes_by_X/PolyGene_plus/PolyGene_PPI_PoPS_justMAGMA_SAIGE_pooled.txt"),
            row.names = T, col.names = T, sep = "\t", quote=F)
write.table(ccmat2, file = paste0("/n/groups/price/kushal/GeneRank/Genes_by_X/PolyGene_plus/PolyGene_PPI_PoPS_justMAGMA_SAIGE_indiv.txt"),
            row.names = T, col.names = T, sep = "\t", quote=F)


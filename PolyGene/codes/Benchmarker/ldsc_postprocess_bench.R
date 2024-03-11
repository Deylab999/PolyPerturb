annot_cell="/n/groups/price/kushal/singlecellLDSC/data/ANNOTATIONS/Benchmarker/POPs_MAGMA_SAIGE_eqwtd_0kb_qmatched_MAGMA"
results_cell = "/n/groups/price/kushal/singlecellLDSC/data/LDSC_RESULTS/Benchmarker/POPs_MAGMA_SAIGE_eqwtd_0kb_qmatched_MAGMA/baselineLD_v2.1"
annot_modules = list.files(results_cell)
flag = 0
index_in_results = 1

library(data.table)

get_sd_annot = function(cell_path, annot_index = 1, flag=0){
  if(flag == 0){
    if(file.exists(paste0(cell_path, "/", "sd_annot_", annot_index, ".rda"))){
      sd_annot = get(load(paste0(cell_path, "/", "sd_annot_", annot_index, ".rda")))
      return(sd_annot)
    }else{
      flag = 1
    }}

  if(flag == 1){
    num = 0
    den = 0
    ll <- list.files(cell_path, pattern = ".annot.gz")
    for(m in 1:length(ll)){
      dat <- data.frame(fread(paste0("zcat ", cell_path, "/", ll[m])))
      num = num  + (nrow(dat)-1) * var(dat[,4+annot_index])
      den = den + (nrow(dat)-1)
      rm(dat)
    }
  }

  estd_sd_annot = sqrt(num/den)
  save(estd_sd_annot, file = paste0(cell_path, "/", "sd_annot_", annot_index, ".rda"))
  return(estd_sd_annot)
}


base_index = index_in_results
annot_names = list.dirs(paste0(results_cell, "/", annot_modules[1]), full.names = F)[-1]

tau_star_mat = matrix(0, length(annot_modules), length(annot_names))
ptau_star_mat = matrix(0, length(annot_modules), length(annot_names))

for(annot_id in 1:length(annot_modules)){
  final_df = c()
  annot_names = list.dirs(paste0(results_cell, "/", annot_modules[annot_id]), full.names = F)[-1]
  for(aa in 1:length(annot_names)){
    cell_path = paste0(annot_cell, "/", annot_modules[annot_id], "/", annot_names[aa])
    sd_annot1=get_sd_annot(cell_path, annot_index=base_index, flag = flag)
    result.file=paste0(results_cell, "/", annot_modules[annot_id], "/", annot_names[aa], "/",
                       annot_modules[annot_id], ".sumstats.part_delete")
    new_table=read.table(result.file,header=F)
    logfile = paste(results_cell, "/", annot_modules[annot_id], "/", annot_names[aa], "/",
                    annot_modules[annot_id],".sumstats.log", sep="")
    log = read.table(logfile,h=F,fill=T)
    h2g = as.numeric(as.character(log[which(log$V4=="h2:"),5]))
    Mref = 5961159
    coef1=sd_annot1*Mref/h2g
    sc=c()
    for(i in 1:dim(new_table)[1]){
      tau1=as.numeric(new_table[i,base_index])
      taus1=tau1*coef1
      sc=c(sc,taus1)
    }
    mean_sc=mean(sc)
    se_sc=sqrt(199**2/200*var(sc))
    p_sc=pnorm(abs(mean_sc/se_sc), 0, 1, lower.tail=F)*2
    tau_star_mat[annot_id, aa] = mean_sc
    ptau_star_mat[annot_id, aa] = p_sc
  }
  cat("We are at trait:", annot_modules[annot_id])
}

rownames(tau_star_mat) = annot_modules
colnames(tau_star_mat) = annot_names

rownames(ptau_star_mat) = annot_modules
colnames(ptau_star_mat) = annot_names



write.table(tau_star_mat, file = "/n/groups/price/kushal/GeneRank/output/POPs_MAGMA_SAIGE_eqwtd_0kb_qmatched_MAGMA_tau.txt",
            row.names = T, col.names = T, sep = "\t", quote=F)

write.table(ptau_star_mat, file = "/n/groups/price/kushal/GeneRank/output/POPs_MAGMA_SAIGE_eqwtd_0kb_qmatched_MAGMA_ptau.txt",
            row.names = T, col.names = T, sep = "\t", quote=F)


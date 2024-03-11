
ABC_Road_GI_bedgraph_calc <- function(scores,
                                 output_cell,
                                 tissuename = "BRN",
                                 output_bed = "temp.bed"){

  df_pre = data.frame(fread(paste0("/n/groups/price/kushal/extras/",
                                   "AllPredictions.AvgHiC.ABC0.015.minus150.withcolnames.ForABCPaper.txt.gz")))
  df_pre = df_pre[which(df_pre$class == "intergenic" | df_pre$clas == "genic"), ]
  if(tissuename != "ALL"){

  if(tissuename == "BRN"){
    tissuenames2 = c("neur", "Neur", "astro", "spinal", "Brain", "brain")
  }
  if(tissuename == "PANC"){
    tissuenames2 = c("pancrea")
  }
  if(tissuename == "BLD"){
    tissuenames2 = as.character(read.table("/n/groups/price/kushal/extras/ABC.listbloodQC.txt", header=F)[,1])
  }
  if(tissuename == "BLD_Bcell"){
    tissuenames2 = as.character(read.table("/n/groups/price/kushal/extras/ABC.listbloodQC_Bcell.txt", header=F)[,1])
  }
  if(tissuename == "BLD_Tcell"){
    tissuenames2 = as.character(read.table("/n/groups/price/kushal/extras/ABC.listbloodQC_Tcell.txt", header=F)[,1])
  }
  if(tissuename == "BLD_Monocyte"){
    tissuenames2 = as.character(read.table("/n/groups/price/kushal/extras/ABC.listbloodQC_Monocyte.txt", header=F)[,1])
  }
  if(tissuename == "BLD_DC"){
    tissuenames2 = as.character(read.table("/n/groups/price/kushal/extras/ABC.listbloodQC_DC.txt", header=F)[,1])
  }
  if(tissuename == "BLD_Megakaryocyte"){
    tissuenames2 = as.character(read.table("/n/groups/price/kushal/extras/ABC.listbloodQC_Megakaryocyte.txt", header=F)[,1])
  }
  if(tissuename == "BLD_RBC"){
    tissuenames2 = as.character(read.table("/n/groups/price/kushal/extras/ABC.listbloodQC_RBC.txt", header=F)[,1])
  }
  if(tissuename == "LNG"){
    tissuenames2 = c("lung")
  }
  if(tissuename == "GI"){
    tissuenames2 = c("stomach", "intestine", "colon", "esophagus", "Stomach", "Intestine", "Colon", "Rectal", "rectal",
                     "Duodenum", "Gastric")
  }
  if(tissuename == "HRT"){
    tissuenames2 = c("heart", "cardiac", "coronary")
  }
  if(tissuename == "FAT"){
    tissuenames2 = c("adipo", "Adipo")
  }
  if(tissuename == "LIV"){
    tissuenames2 = c("liver", "hepato", "HepG2")
  }
  if(tissuename == "KID"){
      tissuenames2 = c("kidney")
  }
  if(tissuename == "SKIN"){
      tissuenames2 = c("skin", "keratin")
  }
  tissue_ids = as.numeric(unlist(sapply(tissuenames2, function(x) return(grep(x, df_pre$CellType)))))
  }else if (tissuename == "ALL"){
    tissue_ids = 1:nrow(df_pre)
  }

  df = df_pre[tissue_ids, ]
  df2 = cbind.data.frame(df$chr, df$start, df$end, df$TargetGene)
  colnames(df2) = c("chr", "start", "end", "TargetGene")
  matched_ids = match(df2$TargetGene, names(scores))
  temp = as.numeric(scores)[matched_ids]
  temp[is.na(temp)] = 0
  final_bed1 = cbind(df2[,c(1:3)], temp)

  roadmap_meta = read.delim("/n/groups/price/kushal/extras/Roadmap_map_EID_names.txt",
                            header=F)

  if (tissuename != "ALL"){
     if (tissuename == "BLD_Bcell"){
     	Road_ids = unique(as.character(roadmap_meta[unlist(sapply("B ",
                                               function(x) return(grep(x, roadmap_meta[,3])))), 1]))
     } 
     else if (tissuename == "BLD_Tcell"){
     	       Road_ids = unique(as.character(roadmap_meta[unlist(sapply("T ",
                                               function(x) return(grep(x, roadmap_meta[,3])))), 1]))
     }
     else if (tissuename == "BLD_Monocyte"){
            Road_ids = unique(as.character(roadmap_meta[unlist(sapply("onocyte",
                                               function(x) return(grep(x, roadmap_meta[,3])))), 1]))
     }
     else if (tissuename == "BLD_DC"){
     	    Road_ids = unique(as.character(roadmap_meta[unlist(sapply("mononuclear",
                                               function(x) return(grep(x, roadmap_meta[,3])))), 1]))
     }
     else if (tissuename == "BLD_Megakaryocyte"){
            Road_ids = unique(as.character(roadmap_meta[unlist(sapply("mononuclear",
                                               function(x) return(grep(x, roadmap_meta[,3])))), 1]))
     }
     else if (tissuename == "BLD_RBC"){
            Road_ids = unique(as.character(roadmap_meta[unlist(sapply("mononuclear",
                                               function(x) return(grep(x, roadmap_meta[,3])))), 1]))
     } else{
	Road_ids =  unique(as.character(roadmap_meta[unlist(sapply(tissuename,
                                                               function(x) return(grep(x, roadmap_meta[,2])))), 1]))     	  
  }}else{
    Road_ids = unique(as.character(roadmap_meta[,1]))
  }

  Enhancer = c()
  for(ee in Road_ids){
    temp = read.table(paste0("/n/groups/price/kushal/extras/RoadmapLinks/",
                             "links_", ee, "_7_2.5.txt"))
    temp2 = read.table(paste0("/n/groups/price/kushal/extras/RoadmapLinks/",
                              "links_", ee, "_6_2.5.txt"))
    Enhancer = rbind(Enhancer, temp[,1:4], temp2[,1:4])
    cat("We processed file:", ee, "\n")
  }

  geneanno = read.csv("/n/groups/price/kushal/extras/gene_anno_unique_datefix.txt", sep = "\t")
  ff = geneanno$symbol[match(Enhancer[,4], geneanno$id)]
  tmp = cbind.data.frame(Enhancer, ff, 1)
  dff = tmp[,c(1:3, 5, 6)]

  matched_ids = match(dff[,4], names(scores))
  temp = as.numeric(scores)[matched_ids]
  temp[is.na(temp)] = 0

  final_bed2 = cbind(dff[,1:3], temp)
  colnames(final_bed1) = c("V1", "V2", "V3", "V4")
  colnames(final_bed2) = c("V1", "V2", "V3", "V4")
  final_bed = rbind(final_bed1, final_bed2)

  write.table(final_bed, paste0(output_cell, "/", output_bed),
              sep = "\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
}



ABC_Road_I_bedgraph_calc <- function(scores,
                                      output_cell,
                                      tissuename = "BRN",
                                      output_bed = "temp.bed"){

  df_pre = data.frame(fread(paste0("/n/groups/price/kushal/extras/",
                                   "AllPredictions.AvgHiC.ABC0.015.minus150.withcolnames.ForABCPaper.txt.gz")))
  df_pre = df_pre[which(df_pre$class == "intergenic"), ]
  if(tissuename != "ALL"){

    if(tissuename == "BRN"){
      tissuenames2 = c("neur", "Neur", "astro", "spinal", "Brain", "brain")
    }
    if(tissuename == "PANC"){
      tissuenames2 = c("pancrea")
    }
    if(tissuename == "BLD"){
      tissuenames2 = as.character(read.table("/n/groups/price/kushal/extras/ABC.listbloodQC.txt", header=F)[,1])
    }
    if(tissuename == "LNG"){
      tissuenames2 = c("lung")
    }
    if(tissuename == "GI"){
      tissuenames2 = c("stomach", "intestine", "colon", "esophagus", "Stomach", "Intestine", "Colon", "Rectal", "rectal",
                       "Duodenum", "Gastric")
    }
    if(tissuename == "HRT"){
      tissuenames2 = c("heart", "cardiac", "coronary")
    }
    if(tissuename == "FAT"){
      tissuenames2 = c("adipo", "Adipo")
    }
    if(tissuename == "LIV"){
      tissuenames2 = c("liver", "hepato", "HepG2")
    }
    if(tissuename == "KID"){
      tissuenames2 = c("kidney")
    }
    if(tissuename == "SKIN"){
      tissuenames2 = c("skin", "keratin")
  }

    tissue_ids = as.numeric(unlist(sapply(tissuenames2, function(x) return(grep(x, df_pre$CellType)))))
  }else if (tissuename == "ALL"){
    tissue_ids = 1:nrow(df_pre)
  }

  df = df_pre[tissue_ids, ]
  df2 = cbind.data.frame(df$chr, df$start, df$end, df$TargetGene)
  colnames(df2) = c("chr", "start", "end", "TargetGene")
  matched_ids = match(df2$TargetGene, names(scores))
  temp = as.numeric(scores)[matched_ids]
  temp[is.na(temp)] = 0
  final_bed1 = cbind(df2[,c(1:3)], temp)

  roadmap_meta = read.delim("/n/groups/price/kushal/extras/Roadmap_map_EID_names.txt",
                            header=F)

  if (tissuename != "ALL"){

    Road_ids =  unique(as.character(roadmap_meta[unlist(sapply(tissuename,
                                                               function(x) return(grep(x, roadmap_meta[,2])))), 1]))

    print(unique(as.character(roadmap_meta[unlist(sapply(tissuename,
                                                         function(x) return(grep(x, roadmap_meta[,2])))), 3])))
  }else{
    Road_ids = unique(as.character(roadmap_meta[,1]))
  }

  Enhancer = c()
  for(ee in Road_ids){
    temp = read.table(paste0("/n/groups/price/kushal/extras/RoadmapLinks/",
                             "links_", ee, "_7_2.5.txt"))
    Enhancer = rbind(Enhancer, temp[,1:4])
    cat("We processed file:", ee, "\n")
  }

  geneanno = read.csv("/n/groups/price/kushal/extras/gene_anno_unique_datefix.txt", sep = "\t")
  ff = geneanno$symbol[match(Enhancer[,4], geneanno$id)]
  tmp = cbind.data.frame(Enhancer, ff, 1)
  dff = tmp[,c(1:3, 5, 6)]

  matched_ids = match(dff[,4], names(scores))
  temp = as.numeric(scores)[matched_ids]
  temp[is.na(temp)] = 0

  final_bed2 = cbind(dff[,1:3], temp)
  colnames(final_bed1) = c("V1", "V2", "V3", "V4")
  colnames(final_bed2) = c("V1", "V2", "V3", "V4")
  final_bed = rbind(final_bed1, final_bed2)

  write.table(final_bed, paste0(output_cell, "/", output_bed),
              sep = "\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
}





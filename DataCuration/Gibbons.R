# Rename Gibbons data

data_dir ="/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/data/"

#list_zero_cols = list()
for(k in c(5,6,7,8)){
  feature_otu = read.csv(paste0(data_dir ,"CRC_k",k,"/kmer_table.txt"),
                   sep = "\t",stringsAsFactors = FALSE,header = TRUE)
  meta = read.csv(paste0(data_dir ,"CRC_k",k,"/metadata.txt"),
                         sep = "\t",stringsAsFactors = FALSE,header = TRUE)
  
  meta = meta[colnames(feature_otu)[which(colSums(feature_otu) != 0)],]
  feature_otu = feature_otu[,which(colSums(feature_otu) != 0)]
  
  print(dim(meta))
  print(dim(feature_otu))
  
  saveRDS(feature_otu,paste0(data_dir ,"Gibbonsr_max_k",k,"/feature_table.rds"))
  write.table(feature_otu ,paste0(data_dir,"Gibbonsr_max_k",k,"/feature_table.txt"),sep="\t",quote=FALSE)

  feature_otu = sweep(feature_otu, MARGIN = 2,STATS = colSums(feature_otu), "/")
  saveRDS(feature_otu,paste0(data_dir ,"Gibbonsr_max_k",k,"/feature_table_rel.rds"))
  write.table(feature_otu ,paste0(data_dir,"Gibbonsr_max_k",k,"/feature_table_rel.txt"),sep="\t",quote=FALSE)

  write.table(meta ,paste0(data_dir,"Gibbonsr_max_k",k,"/metadata.txt"),sep="\t",quote=FALSE)
  
  
  
}


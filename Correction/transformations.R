## example: run from main director (Microbiome_PCA_Correction)
# Rscript Correction/transformations.R ~/Documents/MicroBatch/microbatch_vc/data/Gibbonsr_max_k5 rds counts vst

#on cluster
# Rscript Correction/transformations.R /u/home/b/briscoel/project-halperin/MicroBatch/data/Gibbonsr_max_k5 rds counts vst

args = commandArgs(trailingOnly=TRUE)
data_dir = args[1] #~/Documents/MicroBatch/microbatch_vc/data/Gibbonsr_max_k5 
format = args[2] #rds
trans = args[3]
if(length(args)> 3){
  correction_method = args[4]
}else{
  correction_method = "standard"
}
packs = c("edgeR","compositions","DESeq2")
lapply(packs,library, character.only = TRUE)

require(compositions)
if(trans == "counts"){
  trans_string = ""
}else{
  trans_string = paste0("_",trans)
}

if(format == "rds"){
  
  feature_table =  readRDS(paste0(data_dir,"/feature_table",trans_string,".rds"))
  
}else{
  feature_table =  read.csv(paste0(data_dir,"/feature_table",trans_string,".txt"), sep = "\t",stringsAsFactors = FALSE,header=TRUE,row.names=1)
}
source(paste0("Correction/correction_source.R"),local=FALSE)

##Replace 0 with pseudocount
print("Starting CLR transformation of data")
pseudocount  = min(feature_table[feature_table!= 0])*0.65


num_zero = sum(feature_table == 0)
prop_zero = num_zero/(dim(feature_table)[1] * dim(feature_table)[2])

print("Adding pseudocount to zeroes")
feature_table[feature_table == 0] = pseudocount
print(pseudocount)
print("Number of zero values in data" )
print(num_zero)
print("Proportion of values in data that are zero" )
print(prop_zero)

if(correction_method  =="standard"){
  feature_table = t(clr(t(feature_table)))
  feature_table = data.frame(feature_table)
  feature_table = as.matrix(feature_table)
  
  
  saveRDS(feature_table,paste0(data_dir,"/feature_table",trans_string,"_clr.rds"))
  write.table(feature_table,paste0(data_dir,"/feature_table",trans_string,"_clr.txt"),sep="\t",quote=FALSE)
  
  feature_table = t(scale_custom(t(feature_table)))
  
  saveRDS(feature_table,paste0(data_dir,"/feature_table",trans_string,"_clr_scale.rds"))
  write.table(feature_table,paste0(data_dir,"/feature_table",trans_string,"_clr_scale.txt"),sep="\t",quote=FALSE)
  print(paste0("Transformation done and exported to ", data_dir))
}else{
  
  if(correction_method  =="logcpm"){
    y <- DGEList(counts=feature_table)
    keep <- filterByExpr(y)
    y <- y[keep,,keep.lib.sizes=FALSE]
    y <- calcNormFactors(y)
    
    logcpm <- cpm(y, log=TRUE)
    feature_table = as.matrix(logcpm)
  }else if(correction_method  =="vst"){
    feature_table <- as.matrix(varianceStabilizingTransformation(as.matrix(feature_table)))
    
  }
  saveRDS(feature_table,paste0(data_dir,"/feature_table",trans_string,"_", correction_method,".rds"))
  write.table(feature_table,paste0(data_dir,"/feature_table",trans_string,"_", correction_method,".txt"),sep="\t",quote=FALSE)
  
}


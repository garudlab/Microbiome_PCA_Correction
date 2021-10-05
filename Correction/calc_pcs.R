## example: run inside Corrections Folder
## example: run from main director (Microbiome_PCA_Correction)
# Rscript Correction/calc_pcs.R ~/Documents/MicroBatch/microbatch_vc/data/Gibbonsr_max_k5 rel_clr
# needed packages: bigstatsr

#install.packages("bigstatsr")
args = commandArgs(trailingOnly=TRUE)
data_dir = args[1] #~/Documents/MicroBatch/microbatch_vc/data/Gibbonsr_max_k5 
format = args[2]
transformation = args[3] #rel_clr

# FUNCTIONS
source(paste0("Correction/correction_source.R"))



##READ in data 
num_pcs_calc = 15

if(format == "rds"){
  feature_table =  readRDS(paste0(data_dir,"/feature_table_",transformation,".rds"))
  
}else{
  feature_table =  read.csv(paste0(data_dir,"/feature_table_",transformation,".txt"), sep = "\t",stringsAsFactors = FALSE,header=TRUE,row.names=1)
}

pca_score_results = pca_method(feature_table, num_pcs = num_pcs_calc)
pca_score = pca_score_results[[1]]
print("dim pca scores")
print(dim(pca_score))
saveRDS(pca_score,paste0(data_dir,"/pca_score_",transformation,".rds"))
write.table(pca_score_results[[2]],paste0(data_dir,"/eigenvalues_",transformation,".txt"))



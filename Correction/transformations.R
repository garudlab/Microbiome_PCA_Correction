## example: run from main director (Microbiome_PCA_Correction)
# Rscript Correction/transformations.R ~/Documents/MicroBatch/microbatch_vc/data/Gibbonsr_max_k5 rds

args = commandArgs(trailingOnly=TRUE)
data_dir = args[1] #~/Documents/MicroBatch/microbatch_vc/data/Gibbonsr_max_k5 
format = args[2] #rds
require(compositions)

if(format == "rds"){
  feature_table =  readRDS(paste0(data_dir,"/feature_table_rel.rds"))
  
}else{
  feature_table =  read.csv(paste0(data_dir,"/feature_table_rel.txt"), sep = "\t",stringsAsFactors = FALSE,header=TRUE,row.names=1)
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
feature_table = t(clr(t(feature_table)))
feature_table = data.frame(feature_table)
feature_table = as.matrix(feature_table)

saveRDS(feature_table,paste0(data_dir,"/feature_table_rel_clr.rds"))
write.table(feature_table,paste0(data_dir,"/feature_table_rel_clr.txt"),sep="\t",quote=FALSE)

feature_table = t(scale_custom(t(feature_table)))

saveRDS(feature_table,paste0(data_dir,"/feature_table_rel_clr_scale.rds"))
write.table(feature_table,paste0(data_dir,"/feature_table_rel_clr_scale.txt"),sep="\t",quote=FALSE)
print(paste0("Transformation done and exported to ", data_dir))

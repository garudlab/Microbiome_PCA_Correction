args = commandArgs(trailingOnly=TRUE)
print(args)

require(compositions)
main_dir ="/u/home/b/briscoel/project-halperin/MicroBatch/data/"
script_folder= "/u/home/b/briscoel/project-halperin/MicroBatch/RevisionSequence/"

folder =args[1] #"AGPr_max_k6" #  "AGPr_complete_otu"  #
data_dir = paste0(main_dir,folder,"/")
metadata_table = read.csv(paste0(data_dir,"metadata.txt"), sep = "\t",stringsAsFactors = FALSE,header=TRUE,row.names=1)
feature_table =  readRDS(paste0(data_dir,"feature_table_rel.rds"))
source(paste0(script_folder,"/correction_source.R"))

##Replace 0 with pseudocount
pseudocount  = min(feature_table[feature_table!= 0])*0.65
print("pseudo")
print(pseudocount)
print("sum zero and zero prop")
num_zero = sum(feature_table == 0)
prop_zero = num_zero/(dim(feature_table)[1] * dim(feature_table)[2])
print(num_zero)
print(prop_zero)
feature_table[feature_table == 0] = pseudocount
feature_table = t(clr(t(feature_table)))
feature_table = data.frame(feature_table)
feature_table = as.matrix(feature_table)

saveRDS(feature_table,paste0(data_dir,"/feature_table_rel_clr.rds"))
write.table(feature_table,paste0(data_dir,"/feature_table_rel_clr.txt"),sep="\t",quote=FALSE)

feature_table = t(scale_custom(t(feature_table)))

saveRDS(feature_table,paste0(data_dir,"/feature_table_rel_clr_scale.rds"))
write.table(feature_table,paste0(data_dir,"/feature_table_rel_clr_scale.txt"),sep="\t",quote=FALSE)


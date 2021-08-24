# /u/local/apps/submit_scripts/R_job_submitter.sh -n handlingAGP2.R -m 24 -t 24 -hp -v 3.6.0 


main_dir ="/u/home/b/briscoel/project-halperin/MicroBatch/data/"



# main_dir ="/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/data/"
# folder ="AGPr_complete_otu" 
# data_dir = paste0(main_dir,folder,"/")
# metadata = read.csv(paste0(data_dir,"metadata.txt"), sep = "\t",stringsAsFactors = FALSE,header=TRUE,row.names=1)
# table(metadata$bin_bowel_movement)
# 

folder ="AGPr_max_k8" # "AGPr_complete_otu" #
old_folder ="AGP_max_k8" #"AGP_complete_otu" # 
datype = "kmer" #"otu" 
data_dir = paste0(main_dir,folder,"/")
old_data_dir = paste0(main_dir,old_folder,"/")
metadata = read.csv(paste0(data_dir,"metadata.txt"), sep = "\t",stringsAsFactors = FALSE,header=TRUE,row.names=1)
data_table = readRDS(paste0(old_data_dir,datype,"_table.rds"))
print(dim(data_table))
print(colnames(data_table)[1:10])
print(row.names(metadata)[1:10])

data_table = data_table[,intersect(colnames(data_table),row.names(metadata))]
print(dim(data_table))
saveRDS(data_table,paste0(data_dir,"/feature_table.rds"))
write.table(data_table,paste0(data_dir,"/feature_table.txt"),sep="\t",quote=FALSE)

data_table = sweep(data_table, MARGIN = 2,STATS = colSums(data_table), "/")
print("TSS norm")
print(colSums(data_table)[1:10])
saveRDS(data_table,paste0(data_dir,"/feature_table_rel.rds"))
write.table(data_table,paste0(data_dir,"/feature_table_rel.txt"),sep="\t",quote=FALSE)

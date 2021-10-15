# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("edgeR")
# BiocManager::install("DESeq2")
# # 

args = commandArgs(trailingOnly=TRUE)
#args = c("Thomasr_complete_otu","edgeR")
print(args)

require(edgeR)
require(DESeq2)
main_dir ="/u/home/b/briscoel/project-halperin/MicroBatch/data/"
script_folder= "/u/home/b/briscoel/project-halperin/MicroBatch/RevisionSequence/"


# main_dir ="/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/data/"
# script_folder= "/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/RevisionSequence/"


folder =args[1]
data_dir = paste0(main_dir,folder,"/")
metadata_table = read.csv(paste0(data_dir,"metadata.txt"), sep = "\t",stringsAsFactors = FALSE,header=TRUE,row.names=1)

feature_table =  readRDS(paste0(data_dir,"feature_table.rds"))
source(paste0(script_folder,"/correction_source.R"))

##Replace 0 with pseudocount
pseudocount  = 1
print("pseudo")
print(pseudocount)
print("sum zero and zero prop")
num_zero = sum(feature_table == 0)
prop_zero = num_zero/(dim(feature_table)[1] * dim(feature_table)[2])
print(num_zero)
print(prop_zero)
feature_table[feature_table == 0] = pseudocount

## START logCPM or deseq2 transform
y <- DGEList(counts=feature_table)
keep <- filterByExpr(y)
y <- y[keep,,keep.lib.sizes=FALSE]
y <- calcNormFactors(y)

logcpm <- cpm(y, log=TRUE)
dim(logcpm)
feature_table_logcpm = as.matrix(logcpm)
range(feature_table)
feature_table_vsd <- as.matrix(varianceStabilizingTransformation(as.matrix(feature_table)))

# 
saveRDS(feature_table_logcpm ,paste0(data_dir,"/feature_table_logCPM.rds"))
write.table(feature_table_logcpm ,paste0(data_dir,"/feature_table_logCPM.txt"),sep="\t",quote=FALSE)
# 
# feature_table = t(scale_custom(t(feature_table)))
# 
saveRDS(feature_table_vsd,paste0(data_dir,"/feature_table_vsd.rds"))
write.table(feature_table_vsd,paste0(data_dir,"/feature_table_vsd.txt"),sep="\t",quote=FALSE)
# 

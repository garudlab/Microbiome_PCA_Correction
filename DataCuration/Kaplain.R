# test = read.csv("/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/data/AGPr_complete_otu/metadata.txt",
#          stringsAsFactors = FALSE,sep="\t" )
# table(test$)
# 
# 
data_dir ="/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/data/Kaplanr_complete_otu"
dir.create(data_dir)


meta_otu = readRDS("/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/data/Hispanic_otu/metadata.rds")
range(as.numeric(meta_otu$host_body_mass_index),na.rm = TRUE)
range(as.numeric(meta_otu$bmi_v2),na.rm = TRUE)
table(meta_otu$Assay.Type)
meta_otu = meta_otu[!grepl("BLANK",row.names(meta_otu)),]
sort(table(meta_otu$host_subject_id),decreasing = TRUE)[1:10]
meta_otu[meta_otu$host_subject_id == "G1321",]
meta_otu = meta_otu[(row.names(meta_otu)!="G1322A"),]
meta_otu = meta_otu[(row.names(meta_otu)!="G1085A"),]
dim(meta_otu)
feature_otu =  readRDS("/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/data/Hispanic_otu/otu_table.rds")

feature_otu = feature_otu[,row.names(meta_otu)]

saveRDS(feature_otu ,paste0(data_dir,"/feature_table.rds"))
write.table(feature_otu ,paste0(data_dir,"/feature_table.txt"),sep="\t",quote=FALSE)

feature_otu = sweep(feature_otu, MARGIN = 2,STATS = colSums(feature_otu), "/")

saveRDS(feature_otu ,paste0(data_dir,"/feature_table_rel.rds"))
write.table(feature_otu ,paste0(data_dir,"/feature_table_rel.txt"),sep="\t",quote=FALSE)
write.table(meta_otu ,paste0(data_dir,"/metadata.txt"),sep="\t",quote=FALSE)

table(meta_otu$extraction_robot..exp.)
#general_dir = "/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/data/"


general_dir = "/u/home/b/briscoel/project-halperin/MicroBatch/data/"
data_dir =paste0(general_dir,"Kaplanr_complete_otu/")
# 
feature_otu = readRDS(paste0(data_dir,"feature_table.rds"))


for(k in c(5,6,7,8)){
  print(k)
  input_dir = paste0(general_dir,"Hispanic_k",k,"/")
  meta_kmer = readRDS(paste0(input_dir,"metadata.rds"))
  data_dir =paste0(general_dir,"Kaplanr_max_k",k,"/")
  dir.create(data_dir)
  
  
  meta_kmer = meta_kmer[colnames(feature_otu),]
  
  print(meta_kmer[1:4,1:4])
  
  feature_kmer =  readRDS(paste0(input_dir,"kmer_table.rds"))
  feature_kmer = feature_kmer [,row.names(meta_kmer)]
  saveRDS(feature_kmer  ,paste0(data_dir,"feature_table.rds"))
  write.table(feature_kmer  ,paste0(data_dir,"feature_table.txt"),sep="\t",quote=FALSE)
  feature_kmer = sweep(feature_kmer, MARGIN = 2,STATS = colSums(feature_kmer), "/")
  
  saveRDS(feature_kmer ,paste0(data_dir,"feature_table_rel.rds"))
  write.table(feature_kmer ,paste0(data_dir,"feature_table_rel.txt"),sep="\t",quote=FALSE)
  write.table(meta_kmer ,paste0(data_dir,"metadata.txt"),sep="\t",quote=FALSE)
  
}
table(meta_kmer$extraction_robot..exp.,meta_kmer$yrsus_c2_v2)
table(meta_otu$extraction_robot..exp.,meta_otu$diabetes_self)
table(meta_otu$extraction_robot..exp.,meta_otu$diabetes_self_v2)
table(meta_otu$extraction_robot..exp.,meta_otu$diabetes_lab_v2.y)

table(meta_otu$extraction_robot..exp.,meta_otu$center)


hist(as.integer(meta_otu$yrsus_c2_v2))
hist(as.integer(meta_otu$us_born_v2))


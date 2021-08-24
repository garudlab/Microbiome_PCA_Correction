
require(lefser)
require(matrixStats)
require(UpSetR)
set.seed(0)
args = commandArgs(trailingOnly=TRUE)
local = TRUE

require(dplyr)
if(local){
  #args = c("Thomasr_complete_otu","rel","uncorrected")
  #args = c("Thomasr_complete_otu","rel_clr","pca4counts")
  #args = c("Kaplanr_complete_otu","rel","uncorrected")
  args = c("Kaplanr_complete_otu","rel_clr","pca2counts")
  
}

print(args)

main_dir ="/u/home/b/briscoel/project-halperin/MicroBatch/data/"
script_folder= "/u/home/b/briscoel/project-halperin/MicroBatch/RevisionSequence/"

if(local){
  main_dir ="/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/data/"
  script_folder= "/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/RevisionSequence/"
  
}

folder = args[1] # "AGPr_max_k5" #"AGPr_complete_otu" 
trans = args[2] #"rel"
correction= args[3] #"rel"

data_dir = paste0(main_dir,folder,"/")
metadata_table = read.csv(paste0(data_dir,"metadata.txt"), sep = "\t",stringsAsFactors = FALSE,header=TRUE,row.names=1)
if(correction == "uncorrected"){
  feature_table =  as.matrix(readRDS(paste0(data_dir,"feature_table_",trans,".rds")))
  
}else{
  feature_table =  readRDS(paste0(data_dir,"feature_table_",trans,"_", correction,".rds"))
}


require("compositions")

  #feature_table_counts = t(clrInv(z=t(feature_table)))
  #feature_table_counts = exp(feature_table)
  ####feature_table_counts = apply(feature_table_counts,2,as.numeric)
  #typeof(feature_table_counts)
  #range(feature_table_counts)
  #feature_table = as.matrix(feature_table_counts)  * 10000000
  # #typeof(feature_table)
feature_table = apply(feature_table,2,as.numeric)
feature_table = feature_table * 10000000


if(grepl("Thomas",folder)){
  study_names = unique(metadata_table$dataset_name)
  dataset_batch = "dataset_name"
  dataset_phenotype = "bin_crc_normal"
}else if(grepl("Kaplan",folder)){
  study_names = unique(metadata_table$extraction_robot..exp.)
  dataset_batch = "extraction_robot..exp."
  dataset_phenotype = "bmi_v2"
 
  
  new_pheno = sapply(metadata_table[,dataset_phenotype],function(x){
    if(is.na(x)){return(NA)}
    else if(x %in% c("not provided","not applicable","Not provided","Unspecified","")){return(NA)}
    else{
      #print(x)
      num_x = as.numeric(x)
      if((num_x < 25) & (num_x > 18.5)){
        return("Healthy")
      }else if((num_x > 30) & (num_x < 35)){
        return("Overweight")
      }else{
        return(NA)
      }
    }
    
  })
  table(new_pheno)
  metadata_table[,dataset_phenotype] = new_pheno
  
}


taxonomy = read.csv("/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/data/97_otu_taxonomy.txt",
                    sep = "\t",header=FALSE)

lefse_list = list()
#for( s in 1:3){
for( s in 1:length(study_names)){
 

  study_name = study_names[s]
  print(study_name )
  conditions = (metadata_table[, dataset_batch] == study_name & !is.na(metadata_table[,dataset_phenotype]))
  metadata_table_clean = metadata_table[conditions,]
  feature_table_clean = feature_table[,conditions]
  
  feature_table_clean = feature_table_clean[rowVars(feature_table_clean) != 0,]
  
  
  se0 <- SummarizedExperiment(assays=SimpleList(exprs=feature_table_clean),
                              colData=metadata_table_clean[,c(dataset_phenotype,dataset_batch)])
  
  #colData=metadata_table_clean[,c("bin_crc_normal","country")])
  
  res_group <- lefser(se0, groupCol = dataset_phenotype)
 
  if(grepl("Thomas",folder)){
    final_names = res_group$Names
  }else{
    # temp_names = res_group$Names
    # length(temp_names)
    # temp  = taxonomy %>% filter(V1 %in% temp_names)
    # final_names = as.character(temp$V2)
    final_names = res_group$Names
  }
  
  
  lefse_list[[study_name]] = final_names
}
dim(feature_table)
feature_table[1:4,1:4]
#install.packages("UpSetR")
pdf(paste0(data_dir ,"/upset_",trans,"_",correction,".pdf"))
upset(fromList(lefse_list),order.by="freq",nsets=length(study_names))

# upset(fromList(lefse_list),order.by="freq",nsets=length(study_names),
#       mainbar.y.max = 60)
dev.off()

saveRDS(lefse_list,paste0(data_dir ,"/upset_",trans,"_",correction,".rds"))


lefse_list = readRDS(paste0(data_dir ,"/upset_",trans,"_",correction,".rds"))
names(lefse_list) = gsub("_"," ",names(lefse_list))
pdf(paste0(data_dir ,"/upset_",trans,"_",correction,".pdf"))

upset(fromList(lefse_list),order.by="freq",nsets=7,text.scale=2,mainbar.y.max = 60)

# upset(fromList(lefse_list),order.by="freq",nsets=length(study_names),
#       mainbar.y.max = 60)
dev.off()
# trans = "rel"
# correction = "uncorrected"
# lefse_list = readRDS(paste0(data_dir ,"/upset_",trans,"_",correction,".rds"))
# pdf(paste0(data_dir ,"/upset_",trans,"_",correction,".pdf"))
# upset(fromList(lefse_list),order.by="freq",nsets=length(study_names),mainbar.y.max = 250)
# dev.off()

for( l in 1:length(lefse_list)){
  if( l == 1){
    get_intersection = lapply(lefse_list[[l]],function(x){return(1)})
    names(get_intersection) = lefse_list[[l]]
  }else{
    for(x in lefse_list[[l]]){
      if(x %in% names(get_intersection)){
        get_intersection[[x]] = get_intersection[[x]] +1
      }else{
        get_intersection[[x]] = 1
      }
    
    }
  }
  
  
}

sum(get_intersection > 1)
sum(get_intersection > 2)
sum(get_intersection > 3)

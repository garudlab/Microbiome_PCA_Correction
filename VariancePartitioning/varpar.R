args = commandArgs(trailingOnly=TRUE)
local = TRUE
if(local){
  #args = c("Thomasr_complete_otu","rel","study_condition") #study_condition
  #args = c("Gibbonsr_complete_otu","rel_clr","bin_crc")
  #args = c("Gibbonsr_complete_otu","rel","bin_crc")
  #args = c("Gibbonsr_max_k6","rel_clr","bin_crc_normal")
  args = c("AGPr_complete_otu","rel_clr","Instrument") #bin_antibiotic_last_year
}
print(args)


main_dir ="/u/home/b/briscoel/project-halperin/MicroBatch/data/"
script_folder= "/u/home/b/briscoel/project-halperin/MicroBatch/RevisionSequence/"

if(local){
  main_dir ="/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/data/"
  script_folder= "/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/RevisionSequence/"
  
}

folder = args[1] #"AGPr_max_k5" #"AGPr_complete_otu" 
trans = args[2] #"rel"
group_column = args[3] # "Instrument"
data_dir = paste0(main_dir,folder,"/")
metadata_table = read.csv(paste0(data_dir,"metadata.txt"), sep = "\t",stringsAsFactors = FALSE,header=TRUE,row.names=1)

# FUNCTIONS
source(paste0(script_folder,"/correction_source.R"))


### MAKE PCA COrrelation plot
library(corrplot)
library(variancePartition)
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("variancePartition")
if(folder == "Thomasr_complete_otu"){
  covariates = c('bin_crc_normal',"age","BMI","dataset_name","gender",'DNA_extraction_kit','country')
  main_phenotype = 'bin_crc_normal'
  
}

if(folder == "AGPr_complete_otu"){
  covariates = c("bin_antibiotic_last_year","bmi_corrected","age_corrected","race.x", "bin_alcohol_consumption",
                 "bin_omnivore_diet","bin_bowel_movement","Instrument","collection_year")
  main_phenotype = c("bin_antibiotic_last_year","bmi_corrected")
  
  
  new_pheno = sapply(metadata_table$bin_antibiotic_last_year,function(x){
    if(is.na(x)){return(NA)}
    if(x == "Yes"){return(1)}
    if(x == "No"){return(0)}
    
  })
  metadata_table$bin_antibiotic_last_year = new_pheno
}

if(grepl("Gibbonsr",folder)){
  
  covariates = c('bin_crc_normal',"age","sex","seq_meth","host_race","study","library_size")
  main_phenotype = 'bin_crc_normal'
  
  
  if(grepl("otu",folder)){
    new_pheno = sapply(metadata_table$bin_crc_normal,function(x){
      if(is.na(x)){return(NA)}
      if(x == "CRC"){return(1)}
      if(x == "H"){return(0)}
      
    })
    metadata_table$bin_crc_normal = unlist(new_pheno)
    
  }
  random_effects_tech = c("study","seq_meth") # "center_project_name","collection_days")#"Instrument",
  random_effects_bio = c("host_race","bin_crc_normal","sex") 
  fixed_effects_tech = c("library_size")
  fixed_effects_bio = c("age")#,"bmi_corrected")
  binary_vars = c("bin_crc_normal","sex")
  categorical_vars = c("study","seq_meth","host_race")
  numeric_vars = c("library_size","age") #"bmi_corrected",
  
  
  top_5 = c( "study","seq_meth","library_size","age","sex","host_race","bin_crc_normal")
  
}



## POLISHING METADATA
input_metadata_table = process_model_matrix(total_metadata = total_metadata,
                                          binary_vars=binary_vars,
                                          categorical_vars =categorical_vars,
                                          numeric_vars = numeric_vars)

technical_vars = c(random_effects_tech,fixed_effects_tech)
biological_vars = c(random_effects_bio,fixed_effects_bio)
random_effects_vars = c(random_effects_tech,random_effects_bio)
fixed_effects_vars = c(fixed_effects_tech,fixed_effects_bio)

for(n in numeric_vars){
  input_metadata_table[,n] = scale(input_metadata_table[,n])
}

for( r in random_effects_vars){
  input_metadata_table[,r] = as.character(input_metadata_table[,r])
}

##END  POLISHING METADATA



# VARPAR TEST
input_abundance_table = readRDS(paste0(data_dir,"feature_table_",trans,".rds"))
input_abundance_table  = as.matrix(input_abundance_table)
input_abundance_table = apply(input_abundance_table,2,as.numeric)



formula_random = paste0('~ (1| ',paste( random_effects_vars, collapse = ') + (1|'),")")
formula_fixed =  paste(fixed_effects_vars, collapse = ' + ')

formula_input = paste0(formula_random, " + ", formula_fixed)


result = fitExtractVarPartModel(formula = formula_input,
                                exprObj = input_abundance_table, data = data.frame(input_metadata_table))

 #"BMI"
result[,top_5]
p <- plotVarPart(result) + theme(text = element_text(size=16),axis.text.x= element_text(size=14),plot.title = element_text(size=17)) +
  scale_x_discrete(labels=top_5_pretty) + ggtitle(pretty_titles[t]) +
  xlab("Known technical variables") + ylab("Proportion of feature variance") 

ggsave(filename = paste0(plot_path,'/top5_plot_kmer_variance_partition_',varpar_types[t],'.pdf'),
       plot = p,width = 7, height = 6)
# END VARPAR


## PLOTTING:

vp5 = vp[,top_5] 

p <- plotVarPart(vp5) + theme(text = element_text(size=16),axis.text.x= element_text(size=14),plot.title = element_text(size=17)) +
  scale_x_discrete(labels=top_5_pretty) + ggtitle(pretty_titles[t]) +
  xlab("Known technical variables") + ylab("Proportion of feature variance") 

ggsave(filename = paste0(plot_path,'/top5_plot_kmer_variance_partition_',varpar_types[t],'.pdf'),
       plot = p,width = 7, height = 6)




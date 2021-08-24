rm(list = ls())
args = commandArgs(trailingOnly=TRUE)
print(args)
# HOFFMAN
#for method in DomainCorrect; do for trans in none; do for phen in bin_antibiotic_last_year; do /u/local/apps/submit_scripts/R_job_submitter.sh -n variance_partioning.R -m 10 -t 24 -v 3.6.0 -arg kmer -arg 7 -arg /u/home/b/briscoel/project-halperin/MicroBatch -arg AGP_max -arg "$method"filter_TRUE_trans_"$trans" -arg protect_"$phen" -arg BatchCorrected -arg 0 -arg 0 -arg 0; done; done; done

#args = c("otu", "WR_AD","~/Documents/MicroBatch/", "0-0.5","1-2","01/07/2016","DiseaseState","study")
# args = c("kmer", 5,'/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc',"AGP_max",
#          "minerva_first1filter_TRUE_trans_clr_scale","protect_bmi_corrected","BatchCorrected",
#          0,0,0)

# args = c("kmer", 6,'/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc',"Hispanic",
#          "minerva_first1filter_TRUE_trans_clr_scale","protect_diabetes3_v2",
#          "BatchCorrected",0,1,0)
# args = c("kmer", 7,'/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc',"Thomas",
#          "minerva_first2filter_TRUE_trans_clr_scale","protect_bin_crc_normal",
#          "BatchCorrected",1,0,0)

local = TRUE
# "bmc&DCC&combat&limma&clr&clr_pca1counts&clr_pca33counts"
# args = c("Thomasr_complete_otu","rel",
#          "clr_pca33counts",0,0,0)

args = c("Kaplanr_complete_otu","rel",
         "nocorrection&bmc&DCC&combat&limma&clr&clr_pca2&clr_pca33",0,0,0)
# args = c("otu", 6,'/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc',"CRC_thomas",
#          "minerva_first1filter_TRUE_trans_clr_scale","protect_bin_crc_normal",
#          "BatchCorrected",0,0,0)

# args = c("kmer", 7,'/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc',"CRC",
#          "minerva_first1filter_TRUE_trans_clr_scale","protect_bin_crc_normal",
#          "BatchCorrected",0,0,0)
# 
# args = c("otu", 7, "/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc", "AGP_otumatch_noabx",
#          "raw&bmc&ComBat&limma",'Instrument',"BatchCorrected",1)
#args[5] = "no_scale_clr&no_scale_no_clr"
# ============================================================================== #
# user input
#data_type = args[1]#"kmer"
#kmer_len = args[2]#6
#microbatch_folder = args[3]#'/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/'
#study_name = args[4]


folder = args[1]
trans = args[2]
methods_list = unlist(strsplit(args[3],"&"))#c("ComBat_with_batch2")#"pca_regress_out_scale","clr_pca_regress_out_no_scale","clr_pca_regress_out_scale") #)#,
use_quant_norm = as.logical(as.integer(args[4]))
use_std =  as.logical(as.integer(args[5]))
filter_low_counts =  as.logical(as.integer(args[6]))
apply_bootstrap = FALSE
bootstrap_prop = 0.80
# ============================================================================== #
# load packages and functions
require(varhandle)
#library(variancePartition)
require(matrixStats)
require(dplyr)
require(varhandle)
require(variancePartition)

main_dir ="/u/home/b/briscoel/project-halperin/MicroBatch/data/"
script_folder= "/u/home/b/briscoel/project-halperin/MicroBatch/RevisionSequence/"
data_dir = paste0(main_dir,folder ,"/")

if(local){
  main_dir ="/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/data/"
  script_folder= "/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/RevisionSequence/"
  data_dir = paste0(main_dir,folder ,"/")
  
}
source(paste0(script_folder,"correction_source.R"))
source(paste0(script_folder,"VariancePartitioning/varpar_source.R"))


#source(paste0(batch_script_folder,"/batch_correction_source.R"))
# ============================================================================== #
# define input folder
total_metadata  = read.csv(paste0(data_dir, "metadata.txt"), sep = "\t",header=TRUE)

 
# =========================================================================== #
#getting read depth

# ============================================================================== #
# read in data
batch_corrected_data = list()
batch_corrected_data_quant_norm = list()
batch_corrected_data_scale = list()
for(m in 1:length(methods_list)){
  if(methods_list[m] == "nocorrection"){
    feature_table = readRDS(paste0(data_dir, "feature_table_" , trans , ".rds" ))
    
  }else{
    feature_table = readRDS(paste0(data_dir, "feature_table_" , trans ,"_" , methods_list[m] ,".rds" ))
    
  }
  feature_table = apply(feature_table,2,as.numeric)
  batch_corrected_data[[methods_list[m]]] = feature_table
}
length(batch_corrected_data)
#table(total_metadata$antibiotics_current_use)
# ============================================================================== #
# quant norm data?
if(use_quant_norm){
  for(m in 1:length(methods_list)){
    print(Sys.time())
    print(m)
    batch_corrected_data_quant_norm[[methods_list[m]]] = quantile_norm(batch_corrected_data[[methods_list[m]]])
    print(Sys.time())
  }
}else if(use_std){
  for(m in 1:length(methods_list)){
    print(Sys.time())
    print(m)
    batch_corrected_data_scale[[methods_list[m]]] = t(scale(t(batch_corrected_data[[methods_list[m]]])))
    print(Sys.time())
  }
}

# ============================================================================== #
# make model matrix
if(grepl("AGP",folder)){
  replacement_year = total_metadata$collection_year
  replacement_year[replacement_year < 2012 | replacement_year > 2020] = NA
  total_metadata$collection_year = replacement_year
  
  #"collection_AM",
  binary_vars = c("bin_alcohol_consumption","bin_omnivore_diet","bin_antibiotic_last_year")#,"sex")
  categorical_vars = c("race.x","bin_bowel_movement",
                       "collection_year","Instrument")
  numeric_vars = c("bmi_corrected","age_corrected","librarysize")
  
  
}else if(grepl("Kaplan",folder)){
  collection_year = format(as.Date(total_metadata$collection_timestamp, format="%m/%d/%Y"), "%Y")
  total_metadata$collection_year = collection_year
  
  binary_vars = c("antibiotic","sex")
  categorical_vars = c("collection_year","hispanic_origin.x","frequency_bowel_movement.y","diabetes3_v2",
                       "mastermix_lot..exp.","processing_robot..exp.","extraction_robot..exp.",
                       "center","prep","tm300_8_tool..exp.","extractionkit_lot..exp.","plating..exp.","primer_date..exp.","primer_plate..exp.","run_prefix..exp.")
  
  
  # categorical_vars = c("hispanic_origin.x","frequency_bowel_movement.y","diabetes3_v2",
  #                      "mastermix_lot..exp.","processing_robot..exp.","extraction_robot..exp.",
  #                      "center","prep","tm300_8_tool..exp.","extractionkit_lot..exp.","plating..exp.","primer_date..exp.","primer_plate..exp.","run_prefix..exp.")
  # 
  
  numeric_vars = c("bmi_v2","age_v2.x","librarysize")

}else if(grepl("Thomas",folder)){
  binary_vars = c("bin_crc_normal","bin_crc_adenomaORnormal","gender")
  categorical_vars = c("dataset_name")#,"sequencing_platform")
  numeric_vars = c("number_reads","age") #"bmi_corrected",
  
  
}else if(grepl("CRC",folder)){
  binary_vars = c("bin_crc_normal","sex")
  categorical_vars = c("study","seq_meth","host_race")
  numeric_vars = c("library_size","age") #"bmi_corrected",
  
  
}else if(grepl("Thomas",folder)){
  # binary_vars = c("gender","LibraryLayout")
  # categorical_vars = c("study","Instrument",'multi_crc_adenoma_normal',
  #                      'CenterName','DNA_extraction_kit')
  # numeric_vars = c("LibrarySize")#,"age","BMI")
  
  binary_vars = c("gender","LibraryLayout","bin_crc_normal")
  categorical_vars = c("study","Instrument",'country',
                       'CenterName','DNA_extraction_kit')
  numeric_vars = c("LibrarySize")#,"age","BMI")
  
  # binary_vars = c("gender","LibraryLayout","bin_crc_normal")
  # categorical_vars = c("study","Instrument",'country',
  #                      'CenterName','DNA_extraction_kit',"CollectionYear")
  # numeric_vars = c("LibrarySize")#,"age","BMI")
  
  
  
  
}else if(grepl("T2D",folder)){
  binary_vars = c("sex","bin_t2d","seq_instrument")
  categorical_vars = c("study")
  numeric_vars = c("library_size","age")
  
}



total_metadata_mod = process_model_matrix(total_metadata = total_metadata,
                                          binary_vars=binary_vars,
                                          categorical_vars =categorical_vars,
                                          numeric_vars = numeric_vars)

# for(i in colnames(total_metadata_mod)[1:10]){
#   print(table(total_metadata_mod[,i]))
#   print(table(total_metadata[,i]))
# }

total_metadata_mod_formula = as.formula(paste0(" ~ ",paste(colnames(total_metadata_mod), collapse = " + ")))

# ============================================================================== #
# cleaning data


for(n in numeric_vars){
  total_metadata_mod[,n] = scale(total_metadata_mod[,n])
}
# if(grepl("AGP",study_name)){
#   
#   total_metadata_mod$librarysize = scale(total_metadata_mod$librarysize)
#   
#   total_metadata_mod$bmi_corrected =scale(total_metadata_mod$bmi_corrected)
#   
#   total_metadata_mod$age_corrected =scale(total_metadata_mod$age_corrected)
#   
#   
# }else if(grepl("Hispanic",study_name)){
#   total_metadata_mod$librarysize = scale(total_metadata_mod$librarysize)
#   
#   total_metadata_mod$bmi_v2 =scale(total_metadata_mod$bmi_v2)
#   
#   total_metadata_mod$age_v2.x =scale(total_metadata_mod$age_v2.x)
#   
#   
# }


# ============================================================================== #
# define fixed and random
if(grepl("AGP",folder)){
  # random_effects_tech = c("collection_year","Instrument") # "center_project_name","collection_days")#"Instrument",
  # random_effects_bio = c("race.x","bin_alcohol_consumption","bin_omnivore_diet","bin_antibiotic_last_year","bin_bowel_movement") #"diet_type.x","artificial_sweeteners"
  # 
  # fixed_effects_tech = c("librarysize")#,"collection_AM")
  # fixed_effects_bio = c("bmi_corrected","age_corrected") #"sex",
  
  random_effects_tech = c("collection_year","Instrument","race.x","bin_alcohol_consumption","bin_omnivore_diet","bin_bowel_movement") # "center_project_name","collection_days")#"Instrument",
  random_effects_bio = c("bin_antibiotic_last_year") #"diet_type.x","artificial_sweeteners"
  
  fixed_effects_tech = c("librarysize","age_corrected")#,"collection_AM")
  fixed_effects_bio = c("bmi_corrected") #"sex",
  
}else if(grepl("Kaplan",folder)){
  random_effects_tech = c("collection_year","mastermix_lot..exp.","processing_robot..exp.",
                          "extraction_robot..exp.", "center","prep") # "center_project_name","collection_days")#"Instrument",
  
  random_effects_bio = c("hispanic_origin.x","diabetes3_v2","antibiotic","frequency_bowel_movement.y","sex") 
  fixed_effects_tech = c("librarysize")
  fixed_effects_bio = c("bmi_v2","age_v2.x")

}else if(grepl("T2D",folder)){
  random_effects_tech = c("seq_instrument","study") # "center_project_name","collection_days")#"Instrument",
  
  random_effects_bio = c() 
  fixed_effects_tech = c("library_size")
  fixed_effects_bio = c("age","sex","bin_t2d")
  
  
}else if(grepl("Thomas",folder)){
  
  random_effects_tech = c("dataset_name")#,"sequencing_platform")
  random_effects_bio = c("bin_crc_normal","bin_crc_adenomaORnormal","gender") 
  fixed_effects_tech = c("number_reads")
  fixed_effects_bio = c("age")#"age","BMI")
  
}else if(grepl("Gibbons",folder)){
  
  random_effects_tech = c("study","seq_meth") # "center_project_name","collection_days")#"Instrument",
  random_effects_bio = c("host_race","bin_crc_normal","sex") 
  fixed_effects_tech = c("library_size")
  fixed_effects_bio = c("age")#,"bmi_corrected")
  
  
}





# ============================================================================== #
# make formula

technical_vars = c(random_effects_tech,fixed_effects_tech)
biological_vars = c(random_effects_bio,fixed_effects_bio)
random_effects_vars = c(random_effects_tech,random_effects_bio)
fixed_effects_vars = c(fixed_effects_tech,fixed_effects_bio)


for( r in random_effects_vars){
  total_metadata_mod[,r] = as.character(total_metadata_mod[,r])
}
  
  
formula_random = paste0('~ (1| ',paste( random_effects_vars, collapse = ') + (1|'),")")
formula_fixed =  paste(fixed_effects_vars, collapse = ' + ')

formula_input = paste0(formula_random, " + ", formula_fixed)

# ============================================================================== #
# var par
#length(batch_correocted_data_input)

#row.names(batch_corrected_data_input[[methods_list[4]]])
collect_var_pars_full_BC = list()
if(use_quant_norm){
  batch_corrected_data_input = batch_corrected_data_quant_norm
  
}else if(use_std){
  batch_corrected_data_input = batch_corrected_data_scale
}else{
  batch_corrected_data_input = batch_corrected_data
}


for(i in 1:length(batch_corrected_data_input)){
  colnames(batch_corrected_data_input[[methods_list[i]]]) = colnames(batch_corrected_data[[methods_list[i]]])
  print(methods_list[i])
  varPartMetaData = c()
  input_abundance_table = c()
  
  input_abundance_table = batch_corrected_data_input[[methods_list[i]]]
  input_metadata_table = total_metadata_mod
  
  
  common_samples = intersect(colnames(input_abundance_table ),row.names(input_metadata_table))
  print("common samples")
  print(length(common_samples))
  input_abundance_table  = input_abundance_table[,common_samples]
  input_metadata_table = input_metadata_table[common_samples,]
  
  print(dim(input_abundance_table))
  print(dim( input_metadata_table))
  
  if(apply_bootstrap){
    
    samples_picked = sample(1:ncol(input_abundance_table),as.integer(bootstrap_prop*ncol(input_abundance_table)))
    sample_names_picked = colnames(input_abundance_table)[samples_picked]
    
    # subsample
    input_abundance_table= input_abundance_table[,sample_names_picked]
    #sub_abundance_table_kmer = input_abundance_table_kmer[,samples_picked]
    input_metadata_table = input_metadata_table[sample_names_picked,]
  }
  
  print(dim(input_abundance_table))
  print(dim(input_metadata_table))
  # remove any features with 0 variance in uncorrected data
  yes_no_na = apply(input_metadata_table,1,function(x){
    any(is.na(x))
  })
  input_abundance_table = input_abundance_table[,!yes_no_na]
  input_metadata_table = input_metadata_table[!yes_no_na,]
  print(dim(input_abundance_table))
  print(dim(input_metadata_table))
  
  
  input_abundance_table = input_abundance_table[rowVars(as.matrix(input_abundance_table)) > 10e-9,]
  
  if(filter_low_counts){
    #rowvar = rowVars(input_abundance_table)
    
    #filter_at_least_two_samples_sub = (rowSums(input_abundance_table  > 0 ) > 2)
    filter_at_least_two_samples_sub = (rowVars(input_abundance_table) >2e-7 )
    
     #test = rowVars(input_abundance_table)
     #sum(test  > 2e-7)
     #hist(test,breaks=100)
    input_abundance_table = input_abundance_table[filter_at_least_two_samples_sub,]
    
  }

  #row.names(input_abundance_table) = paste0("OTU",1:nrow(input_abundance_table))

  



  # test = cor(input_abundance_table)
  # test2 <- (test == 1)
  # cols <- colSums(test2) 
  # sort(cols,decreasing=TRUE)[1:10]
  # 
  # 
  # fe <- apply(as.factor(input_metadata_table),2,as.numeric)
  # dim(fe)
  # 
  # test = cor(as.numeric(as.factor(input_metadata_table)))
  # test2 <- (test == 1)
  # cols <- colSums(test2) 
  # sort(cols,decreasing=TRUE)[1:10]
  # 
  # 
  # 
  # noncorrelated_samples = names(cols)[which(cols <= 1)]
  # print(paste0("uncorre samples", length(noncorrelated_samples)))
  # input_abundance_table <- input_abundance_table[,noncorrelated_samples]
  # input_metadata_table <- input_metadata_table[noncorrelated_samples,]
  # 

  #table(input_metadata_table[,c(4,5)])
  #table(input_metadata_table$)
  #table(total_metadata$sequencing_platform)
  #sum(is.na(input_abundance_table))
  t1 = Sys.time()
  varPartMetaData = fitExtractVarPartModel(formula = formula_input,
                                           exprObj = input_abundance_table, data = data.frame(input_metadata_table))
  print(Sys.time() - t1)
  collect_var_pars_full_BC[[methods_list[i]]] = varPartMetaData
  
  # write.table(as.matrix(varPartMetaData), paste0(input_folder ,"varpart",methods_list[i],".txt"),
  #             sep = "\t",quote = FALSE)

  saveRDS(varPartMetaData, paste0(data_dir ,"/varpart_quant",use_quant_norm ,"_",methods_list[i],"_filter_", filter_low_counts,".rds"))
  
}
head(varPartMetaData)

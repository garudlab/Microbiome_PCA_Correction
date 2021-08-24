args = commandArgs(trailingOnly=TRUE)
local = TRUE
if(local){
  #args = c("Thomasr_max_k6","rel_clr","dataset_name") #
  #args = c("Thomasr_complete_otu","logCPM","study_condition") #
  #args = c("Gibbonsr_complete_otu","rel","bin_crc")
  #args = c("Gibbonsr_complete_otu","rel_clr","study")
  args = c("Gibbonsr_max_k6","rel","study")
  #args = c("AGPr_complete_otu","rel_clr","Instrument") #bin_antibiotic_last_year
  #args = c("AGPr_max_k6","rel","Instrument") #bin_antibiotic_last_year
  #args = c("Kaplanr_complete_otu","rel","extraction_robot..exp.") #bin_antibiotic_last_year
  #args = c("Kaplanr_max_k6","rel_clr","extraction_robot..exp.") #bin_antibiotic_last_year
}

print(args)

main_dir ="/u/home/b/briscoel/project-halperin/MicroBatch/data/"
script_folder= "/u/home/b/briscoel/project-halperin/MicroBatch/RevisionSequence/"

if(local){
  main_dir ="/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/data/"
  script_folder= "/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/RevisionSequence/"
  
}

#dim(metadata_table)
folder = args[1] #"AGPr_max_k5" #"AGPr_complete_otu" 
trans = args[2] #"rel"
group_column = args[3] # "Instrument"
data_dir = paste0(main_dir,folder,"/")
metadata_table = read.csv(paste0(data_dir,"metadata.txt"), sep = "\t",stringsAsFactors = FALSE,header=TRUE,row.names=1)

# FUNCTIONS
source(paste0(script_folder,"/correction_source.R"))



##READ in data 
num_pcs_calc = 15

pca_score = readRDS(paste0(data_dir,"/pca_score_", trans,".rds"))
colnames(pca_score)  = paste0("PC",c(1:ncol(pca_score)))

#metadata_table = metadata_table[row.names(pca_score),]
print(dim(pca_score))
print(dim(metadata_table))

#row.names(metadata_table)
# MAKE PCA PLOT

## FILTER
#pca_score = pca_score[!is.na(metadata_table[,group_column]),]
#metadata_table = metadata_table[!is.na(metadata_table[,group_column]),]
dev.off()
metadata_table[,group_column] = gsub("_"," ",metadata_table[,group_column])
p <- pca_plot(pca_score,metadata_table, title = "PC 1 and 2",group_column=group_column,coord1=1,coord2=2)
#p <- p + ylim(-20,10)
#p <- p + xlim(-0.3,0)
p

ggsave(p,file=paste0(data_dir,"/pca_plot_",trans, "_", group_column,".pdf"),device ="pdf")

### MAKE PCA COrrelation plot
library(corrplot)
library(variancePartition)
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("variancePartition")
if(grepl("Thomasr",folder )){
  covariates = c('bin_crc_normal',"age","BMI","gender","dataset_name",'DNA_extraction_kit','country')
  main_phenotype = 'bin_crc_normal'
  
  if(grepl("max_k",folder)){
    new_pheno = sapply(metadata_table$bin_crc_normal,function(x){
      if(is.na(x)){return(NA)}
      if(x == "CRC"){return(1)}
      if(x == "H"){return(0)}
      
    })
    metadata_table$bin_crc_normal = unlist(new_pheno)
    
  }
  
}

if(grepl("Kaplanr",folder )){
  table(metadata_table$placeofbirth_group.x)
  covariates = c("bmi_v2","age_v2.x",'placeofbirth_group.x',"extraction_robot..exp.","center",'extractionkit_lot..exp.')
  main_phenotype = "bmi_v2"
  
  new_pheno = sapply(metadata_table$age_v2.x,function(x){
    if(is.na(x)){return(NA)}
    else if(x %in% c("not provided","not applicable","Not provided","Unspecified")){return(NA)}
    else{return(as.integer(x))}
    
  })
  metadata_table$age_v2.x = new_pheno
  
  new_pheno = sapply(metadata_table$bmi_v2,function(x){
    if(is.na(x)){return(NA)}
    else if(x %in% c("not provided","not applicable","Not provided","Unspecified")){return(NA)}
    else{return(as.numeric(x))}
    
  })
  metadata_table$bmi_v2 = new_pheno
  
  
}


if(grepl("AGPr",folder)){
  covariates = c("bin_antibiotic_last_year","bmi_corrected","age_corrected","race.x", "bin_alcohol_consumption",
                 "bin_omnivore_diet","bin_bowel_movement","Instrument","collection_year")
  main_phenotype = c("bin_antibiotic_last_year","bmi_corrected")
  

  new_pheno = sapply(metadata_table$bin_antibiotic_last_year,function(x){
    if(is.na(x)){return(NA)}
    if(x == "Yes"){return(1)}
    if(x == "No"){return(0)}
    
  })
  metadata_table$bin_antibiotic_last_year = new_pheno
  new_pheno = sapply(metadata_table$bmi_corrected,function(x){
    if(is.na(x)){return(NA)}
    else if(x %in% c("not provided","not applicable","Not provided","Unspecified")){return(NA)}
    else{return(as.numeric(x))}
    
  })
  metadata_table$bmi_corrected = new_pheno
}

if(grepl("Gibbonsr",folder)){

  covariates = c('bin_crc_normal',"age","sex","host_race","seq_meth","study","library_size")
  main_phenotype = 'bin_crc_normal'
  

  if(grepl("otu",folder)){
    new_pheno = sapply(metadata_table$bin_crc_normal,function(x){
      if(is.na(x)){return(NA)}
      if(x == "CRC"){return(1)}
      if(x == "H"){return(0)}
      
    })
    metadata_table$bin_crc_normal = unlist(new_pheno)
    
  }
  
}
#print(colnames(metadata_table))
metadata_pc = data.frame(metadata_table[,covariates],
                         pca_score)


pc_formula = as.formula(paste0(" ~ ",colnames(metadata_pc), collapse = " + "))

# metadata_pc = metadata_pc[c(1:50,10000:10200,15000:15042),]
# CanCorC = canCorPairs(formula = '~bin_antibiotic_last_year + ~bmi_corrected + ~age_corrected + 
#     ~race.x + ~bin_alcohol_consumption + ~bin_omnivore_diet + 
#         ~collection_year + 
#             ~PC1 + ~PC2 + ~PC3' , data = metadata_pc )
# num_pcs_calc = 3


CanCorC = canCorPairs(formula = pc_formula , data = metadata_pc )
copy_CanCorC = CanCorC
plotCorrMatrix( copy_CanCorC  )


pheno_confounding_column = CanCorC[1:(nrow(CanCorC)-num_pcs_calc),main_phenotype ]
CanCorC = CanCorC[1:(nrow(CanCorC)-num_pcs_calc),(nrow(CanCorC)-num_pcs_calc +1 ):ncol(CanCorC)]
CanCorC = cbind(CanCorC,pheno_confounding_column )
new_col = colnames(CanCorC)
new_col[(length(new_col) +1 - length(main_phenotype)): length(new_col)] = main_phenotype
colnames(CanCorC) = new_col

p_val_matrix = CanCorC
for(r in 1:nrow(CanCorC )){
  for(cl in 1:ncol(CanCorC )){
    #r = 1
    #cl =1 
    rname = row.names(CanCorC)[r]
    cname = colnames(CanCorC)[cl]
    
    #if(rname %in% categorical_covariates){
    r_input = model.matrix( as.formula(paste0(" ~ ",rname )),metadata_pc )
    
    all_cors = c()
    all_pvals = c()
    head(r_input)
    for( cl2 in 2:ncol(r_input)){
      
      #head(metadata_pc)
      
      cor = cor(r_input[,cl2],metadata_pc[row.names(r_input),cname],use = "pairwise.complete.obs", method = "spearman")
      if(is.na(cor)){
        cor = 0
      }
      cor_test = cor.test(r_input[,cl2],metadata_pc[row.names(r_input),cname],method = "spearman")
      if(is.na(cor_test$p.value)){
        cor_test$p.value = 1
      }
      all_pvals = c(all_pvals,cor_test$p.value)
      all_cors = c(all_cors, cor)
    }
    p_val_matrix[r,cl] = all_pvals[which.max(abs(all_cors))]

  }
}
max_val = 1
new_rn =  gsub("crc_normal","CRC status", row.names(CanCorC) )
new_rn = gsub("_"," ", new_rn)
new_rn = gsub("bin ", "", new_rn)
new_rn = gsub(" corrected","",new_rn)
new_rn = gsub("v2.x","",new_rn)
new_rn = gsub("v2","",new_rn)
new_rn = gsub("\\.."," ",new_rn)



new_cn =  gsub("crc_normal","CRC status", colnames(CanCorC) )
new_cn = gsub("_"," ", new_cn)
new_cn = gsub("bin ", "", new_cn)
new_cn = gsub(" corrected","",new_cn)
new_cn = gsub("v2.x","",new_cn)
new_cn = gsub("v2","",new_cn)
new_cn = gsub("\\.."," ",new_cn)

if(grepl("AGP",folder)){
  new_rn = gsub(".x","",new_rn)
  new_cn = gsub(".x","",new_cn)
}

new_rn = gsub("Instrument","sequencing instrument",new_rn)
row.names(CanCorC) = new_rn
row.names( p_val_matrix) = new_rn
colnames(CanCorC) = new_cn
colnames( p_val_matrix) = new_cn
CanCorC[is.na(CanCorC)] = 0
plot_object = list(CanCorC = CanCorC,p_val_matrix = p_val_matrix )

plot_object$CanCorC
saveRDS(plot_object,paste0(data_dir,"/","CanCorPlotObj_",trans, ".rds"))
#plot_object = readRDS(paste0(data_dir,"/","CanCorPlotObj_",trans, ".rds"))
write.table(plot_object$CanCorC,paste0(data_dir,"/","CanCor_",trans, ".txt"),sep="\t",
            quote = FALSE)



pdf(paste0(data_dir,"/","CanCor_",trans, ".pdf"))
dim(plot_object$CanCorC)
dim(plot_object$p_val_matrix)

require(RColorBrewer)
col2 = colorRampPalette(c("#FFFFFF", "#00569c")) #002b54  2B9EDE
#?corrplot
corrplot(plot_object$CanCorC,na.label = "X",
         na.label.col = "grey",
         p.mat = plot_object$p_val_matrix,insig = "label_sig",
         sig.level = c(.001, .01, .05), pch.cex = .9,pch.col = "black",
         tl.col="black",col.lim=c(0.0,1.0),is.corr=F,tl.cex = 1,
         cl.cex=1,cl.length = 6,cl.align.text = "c",cl.ratio=0.2,
         col=col2(20))
dim(plot_object$CanCorC)
dim(plot_object$p_val_matrix)
dev.off()

plot_object$p_val_matrix[1:4,1:4]

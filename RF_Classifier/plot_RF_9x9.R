args = commandArgs(trailingOnly=TRUE)
local = TRUE

require(dplyr)


margins_list = list()
margins_list[["Thomasr_max_k"]] =c(13,13)
margins_list[["Thomasr_complete_otu"]] =c(13,13)
margins_list[["Kaplanr_max_k"]] =c(13,13)
margins_list[["Kaplanr_complete_otu"]] =c(13,13)
margins_list[["Gibbonsr_max_k"]] =c(13,13)
margins_list[["Gibbonsr_complete_otu"]] =c(13,13)
margins_list[["AGPr_complete_otu"]] =c(23,10)
margins_list[["AGPr_max_k"]] =c(23,10)

notecex_list = list()
notecex_list[["Thomasr_max_k"]] = 0.8
notecex_list[["Thomasr_complete_otu"]] = 0.8
notecex_list[["Kaplanr_max_k"]] = 1
notecex_list[["Kaplanr_complete_otu"]] = 1
notecex_list[["Gibbonsr_max_k"]] = 1
notecex_list[["Gibbonsr_complete_otu"]] = 1
notecex_list[["AGPr_max_k"]] = 1
notecex_list[["AGPr_complete_otu"]] = 1

folders = c("AGPr_complete_otu","Thomasr_complete_otu","Gibbonsr_complete_otu","Kaplanr_complete_otu",
            "AGPr_max_k","Thomasr_max_k","Gibbonsr_max_k","Kaplanr_max_k")
transs = c("logcpm","vst", "rel_clr")
methods = c("combat","limma","bmc")
main_dir ="/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/data/"
script_folder= "/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/RevisionSequence/"
corrections_vec = c()
for(trans in transs){
  for(method in methods){
    corrections_vec = c(corrections_vec,paste0(trans,"_",method))
  }
  
}

folder_spec = list()

for(a in 1:length(folders)){ 
  print(folders[a])
#for(a in c(5)){ 
  #a = 7
  if(grepl("_k",folders[a])){
    
    if(grepl("Thomas",folders[a])){
      folders_a =  paste0(folders[a],c(6:7))
    }else if(grepl("Gibbons|AGP|Kaplan",folders[a])){
      folders_a =  paste0(folders[a],c(5:8))
    }
   
  }else{
    folders_a = folders[a]
  }
  
  corrections_val = list()
  
  for(fold_a in folders_a){
    args = c(fold_a,"True","AUC")
    folder = args[1] # "AGPr_max_k5" #"AGPr_complete_otu" 
    lodo = args[2]
    meas = args[3]
    if(grepl("Kaplan",folder)){
      meas = "pearson"
    }
    data_dir = paste0(main_dir,folder,"/")
    
    for(cori in 1:length(corrections_vec)){
      if(meas == "AUC"){
        pathy  = paste0(data_dir, "GRID_", meas,"_OUTPUT_" ,corrections_vec[cori] , "_lodo_" , lodo  , ".csv")
        metrics_pre =  read.csv(pathy, header=TRUE)
        
        
      }else{
        metrics_pre =  read.csv(paste0(data_dir, "PRED_OUTPUT_" ,corrections_vec[cori] , "_lodo_" , lodo  , ".csv"), header=TRUE)
        
        
      }

      
      if(meas == "AUC"){
        mean_metric = mean(metrics_pre$val_auc)
        metrics = data.frame(t(metrics_pre$test_auc))
        colnames(metrics) =paste0(meas,"_",metrics_pre$fold)
       
      }else{
        mean_metric  = mean(metrics_pre$val_pearson)
        metrics = data.frame(t(metrics_pre$test_pearson))
        colnames(metrics) =paste0(meas,"_",metrics_pre$test_dataset)
  
      }
      
      
      metrics$mean_metric = mean_metric 
      
      if(cori == 1){
        corrections_df = metrics
      }else{
        corrections_df = rbind(corrections_df, metrics)
      }
    }
    
    corrections_df$corrections = corrections_vec
    if(length(folders_a)> 1){
      corrections_val[[fold_a]] = corrections_df
    }
  }
  
  if(length(folders_a) > 1){
    val_list = list()
    for(f in 1:length(folders_a)){
      val_list[[f]] = corrections_val[[folders_a[f]]]$mean_metric
    }
    best_param = apply(do.call(cbind,val_list),1,function(x){which.max(x)})
    print(best_param)
    corrections_df = list()
    for(p in 1:length(best_param)){
      corrections_df[[p]] = corrections_val[[folders_a[best_param[p]]]][p,]
    }
    corrections_df = do.call( rbind, corrections_df)
  }
  require(gplots)
  nonnice_names =  corrections_vec
  nice_names =  gsub("_"," ",corrections_vec)
  nice_names =  gsub("logcpm","logCPM", nice_names )
  nice_names =  gsub("bmc","BMC", nice_names )
  nice_names =  gsub("rel clr","CLR", nice_names )
  nice_names =  gsub("combat","ComBat", nice_names )
  nice_names =  gsub("vst","VST", nice_names )
  
  
  AUC_results = corrections_df[,grepl(meas,colnames(corrections_df))]
  row.names(AUC_results) = corrections_vec
  row.names(AUC_results) = nice_names #"ComBat",
  AUC_results$Average = rowMeans(AUC_results)
  input = as.matrix(AUC_results)
  colnames(input) = gsub(paste0(meas,"_"), "", colnames(input))
  
  
  if(grepl("Thomas",folder )){
    input = input[,c("Average","FengQ_2015", "ThomasAM_2018b",  "ZellerG_2014", "YuJ_2015" ,  "ThomasAM_2018a" , "VogtmannE_2016"  ,"HanniganGD_2017" )]
  }
  
  colnames(input) = gsub("_", " ", colnames(input))
  colnames(input) = gsub("crc ", "", colnames(input))
  firstup <- function(x) {
    substr(x, 1, 1) <- toupper(substr(x, 1, 1))
    x
  }
  
  colnames(input) = sapply(colnames(input),firstup)
  
  input_str = apply(input,2, function(x){sprintf("%.2f",round(x,2))})
  
  # blue palette c("#FFFFFF", "#2B9EDE"
  
  
  data_dir = paste0(main_dir,folders_a[1],"/")
  
  pdf(paste0(data_dir,"/",meas,"LODO_Heatmap_",trans, ".pdf"))
  heatmap.2(t(input), trace="none", density="none", col=colorRampPalette(c("red","yellow")), cexRow=1.2, cexCol=1.2, 
            margins = margins_list[[folders[a]]],
            Rowv = FALSE, Colv =  "Rowv",cellnote=t(input_str),notecol="black",srtCol = 45,notecex=notecex_list[[folders[a]]])
  dev.off()
  
  custom_colors = c('#e32f27',"#C3FFCE",'#FF9300','#FFE800','#fdd0a2',"#9A33FF","#F133FF","#3341FF","#72C1FC")
  
  require("reshape2")
  require(ggplot2)
  library(ggpubr)
  to_plot = melt(input) 
  to_plot = to_plot %>% filter(Var2 != "Average")
  head(to_plot)
  p <- ggboxplot(to_plot, x = "Var1", y = "value",
                 fill = "Var1", palette = custom_colors) +xlab("Correction") + 
    ylab(paste0("Cross-validated ", meas) ) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1,hjust = 1),text = element_text(size=19),
          legend.position = "none",
          panel.grid.major.x = element_blank() ,
          # explicitly set the horizontal lines (or they will disappear too)
          panel.grid.major.y = element_line( size=.1, color="black" ))
  
  #
  p
  
  ggsave(plot=p,filename=paste0(data_dir,"/",meas,"_BOX_lodo",lodo,"_9x9_",folders_a[1], ".pdf"),width = 7,height = 5,units="in")
  
  saveRDS(p,paste0(data_dir,"/",meas,"_BOX_lodo",lodo,"_9x9_",folders_a[1], ".rds"))

  
  
  
}
  
  


  
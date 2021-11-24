
local = TRUE
main_dir ="/u/home/b/briscoel/project-halperin/MicroBatch/data/"
script_folder= "/u/home/b/briscoel/project-halperin/MicroBatch/RevisionSequence/"

if(local){
  main_dir ="/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/data/"
  script_folder= "/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/RevisionSequence/"
  
}

# folder = "Kaplanr_complete_otu"
# grouping = "batches"
# shortname = "HCHS"
# num_batches = 4
# trans_list = c(rep("rel",6), rep("counts",2), rep("rel_clr",2))
# correction_list = c("uncorrected","DCC","combatcounts","limmacounts","bmc",
#                     "clr","logcpm","vst","pca3counts","pca2counts")
# 

folder = "Thomasr_complete_otu"
shortname = "CRC-WGS"
grouping = "studies"
num_batches = 7
trans_list = c(rep("rel",6), rep("counts",2), rep("rel_clr",2))
correction_list = c("uncorrected","DCC","combatcounts","limmacounts","bmc",
                    "clr","logcpm","vst","pca3counts","pca4counts")
nice_names =  c("Uncorrected","DCC","ComBat","limma","BMC","CLR","logCPM","VST","Fixed PCA Correction","Tuned PCA Correction")

data_dir = paste0(main_dir,folder,"/")
lefse_stats = list()
for(i in 1:length(trans_list)){
  #i = 7
  
  lefse_list= readRDS(paste0(data_dir ,"/upset_",trans_list[i],"_",correction_list[i],".rds"))
  
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
  lefse_stats[[correction_list[i]]] = list()
  for( l in 1:num_batches){
    lefse_stats[[correction_list[i]]][[l]] = sum(get_intersection >= l)
  }

  
  
}

lefse_df = data.frame(do.call(rbind, lefse_stats))
row.names(lefse_df) = nice_names

lefse_df$methods = row.names(lefse_df)
library(tidyr)
biomarker_sharing <- lefse_df %>% gather("sharing","biomarker_num",1:num_batches)
biomarker_sharing$sharing = as.integer(gsub("X","",biomarker_sharing$sharing))
biomarker_sharing$biomarker_num = unlist(biomarker_sharing$biomarker_num)
require(ggplot2)
custom_colors = c('#e32f27',"#C3FFCE",'#FF9300','#FFE800','#fdd0a2',"#9A33FF","#F133FF","#3341FF","#72C1FC","#0093FF")

biomarker_sharing$methods = factor(biomarker_sharing$methods, levels = nice_names)
biomarker_sharing$log_biomarker_num = log(biomarker_sharing$biomarker_num)
g1 <- biomarker_sharing %>% ggplot( aes(x=sharing, y=biomarker_num, group=methods, color=methods)) +
  geom_line() +scale_color_manual(values=custom_colors)+
  ggtitle(paste0("Biomarker sharing: ",shortname)) +
  ylab("Number of LEFSe biomarkers") + xlab(paste0("Number of ",grouping," sharing biomarkers")) + 
  theme_bw() + theme(text = element_text(size=19))

g2 <- biomarker_sharing %>% ggplot( aes(x=sharing, y=log(biomarker_num), group=methods, color=methods)) +
  geom_line() +scale_color_manual(values=custom_colors)+
  ggtitle(paste0("Biomarker sharing: ",shortname)) +
  ylab("log(Number of LEFSe biomarkers)") + xlab(paste0("Number of ",grouping," sharing biomarkers")) + theme_bw() +
  theme(text = element_text(size=19))
        

ggsave(g1,filename = paste0(data_dir ,"/lefse_log", folder,".pdf"))

ggsave(g2,filename = paste0(data_dir ,"/lefse_reg", folder,".pdf"))





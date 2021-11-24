args = commandArgs(trailingOnly=TRUE)
local = TRUE

require(dplyr)
if(local){

  #Kaplan
  # args = c("Kaplanr_max_k5","rel","True","pearson")
  # # "bmc", "combat", "percentilenorm", "limma", "DCC", a
  # dtype_pca = "" # counts
  # pca_methods = c(paste0("clr_pca",c(c(2,33)),dtype_pca)) #,paste0("clr_pca",c(1:5)))

  #Kaplan otu
  args = c("Kaplanr_complete_otu","rel","True","pearson")
  dtype_pca = "" # counts
  pca_methods = c(paste0("clr_pca",c(c(1,2,3,4,5,33)),dtype_pca)) #,paste0("clr_pca",c(1:5)))



}

margins_list = list()
margins_list[["Gibbonsr_complete_otu"]] =c(16,17)
margins_list[["Thomasr_complete_otu"]] =c(13,15)
margins_list[["Thomasr_max_k7"]] =c(15,13)
margins_list[["AGPr_complete_otu"]] =c(5,13)
margins_list[["Kaplanr_complete_otu"]] =c(20,8)
margins_list[["Kaplanr_max_k6"]] =c(20,8)
margins_list[["Kaplanr_max_k5"]] =c(20,8)
margins_list[["Kaplanr_max_k7"]] =c(20,8)
margins_list[["Kaplanr_max_k8"]] =c(20,8)
folder_spec = list()
folder_spec[["AGPr_max_k5"]] = list(nocorrection = "AGPr_max_k5",
                                    combat =  "AGPr_max_k5",
                                    DCC =   "AGPr_max_k5",
                                    limma =  "AGPr_max_k5",
                                    bmc =   "AGPr_max_k5",
                                    clr = "AGPr_max_k8",
                                    clr_pca1counts = "AGPr_max_k5",
                                    clr_pca33counts = "AGPr_max_k6")
# lodo false
# folder_spec[["Kaplanr_max_k5"]] = list(nocorrection = "Kaplanr_max_k5",
#                                     combat =  "Kaplanr_max_k8",
#                                     DCC =   "Kaplanr_max_k5",
#                                     limma =  "Kaplanr_max_k8",
#                                     bmc =   "Kaplanr_max_k5",
#                                     clr =  "Kaplanr_max_k8",
#                                     clr_pca1 = "Kaplanr_max_k8",
#                                     clr_pca33 = "Kaplanr_max_k8")

# lodo true
folder_spec[["Kaplanr_max_k5"]] = list(nocorrection = "Kaplanr_max_k5",
                                    combat =  "Kaplanr_max_k8",
                                    DCC =   "Kaplanr_max_k5",
                                    limma =  "Kaplanr_max_k8",
                                    bmc =   "Kaplanr_max_k5",
                                    clr =  "Kaplanr_max_k8",
                                    clr_pca2 = "Kaplanr_max_k8",
                                    clr_pca33 = "Kaplanr_max_k8",
                                    logcpm =  "Kaplanr_max_k7",
                                    vst =  "Kaplanr_max_k8")

notecex_list= list()
notecex_list[["Gibbonsr_complete_otu"]] = 1
notecex_list[["Thomasr_complete_otu"]] = 1
notecex_list[["Thomasr_max_k7"]] = 1
notecex_list[["AGPr_complete_otu"]] = 1
notecex_list[["Kaplanr_complete_otu"]] = 1
notecex_list[["Kaplanr_max_k6"]] = 1
notecex_list[["Kaplanr_max_k5"]] = 1
notecex_list[["Kaplanr_max_k7"]] = 1
notecex_list[["Kaplanr_max_k8"]] = 1

notecex_text_list = notecex_list
notecex_text_list[["Thomasr_complete_otu"]] = 1.2

print(args)


main_dir ="/u/home/b/briscoel/project-halperin/MicroBatch/data/"
script_folder= "/u/home/b/briscoel/project-halperin/MicroBatch/RevisionSequence/"

if(local){
  main_dir ="/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/data/"
  script_folder= "/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/RevisionSequence/"
  
}

folder = args[1] # "AGPr_max_k5" #"AGPr_complete_otu" 
trans = args[2] #"rel"
lodo = args[3]
meas = args[4]
data_dir = paste0(main_dir,folder,"/")


#   c("clr_scale_pca",
# "clr_pca1roundcounts", "clr_pca1", "clr_pca2roundcounts", "clr_pca2", "clr_pca3roundcounts", "clr_pca3")
other_methods = c("nocorrection","DCC","combat","limma","bmc","clr","logcpm","vst" ) # "combat", 

corrections_vec = c(other_methods, pca_methods)
trans_methods = rep(trans,length(corrections_vec))

#c("nocorrection","bmc","clr_pca")#"clr_pcacounts","clr_scale_pca","clr_pca")
corrections_list = list()
for(cori in 1:length(corrections_vec)){
  #cori = 1
 if( folder %in% names(folder_spec)){
    metrics_pre =  read.csv(paste0(main_dir,folder_spec[[folder]][[corrections_vec[cori]]],"/", "PRED_OUTPUT_" ,trans_methods[cori], "_" , corrections_vec[cori] , "_lodo_" , lodo  , ".csv"), header=TRUE)
    
  }else{
    metrics_pre =  read.csv(paste0(data_dir, "PRED_OUTPUT_" ,trans_methods[cori], "_" , corrections_vec[cori] , "_lodo_" , lodo  , ".csv"), header=TRUE)
    
  }
  
  corrections_list[[corrections_vec[cori]]] = metrics_pre
  head(metrics_pre)
  if(meas == "score"){
    mean_metric = mean(metrics_pre$val_score)
    metrics = data.frame(t(metrics_pre$test_score))
    colnames(metrics) =paste0(meas,"_",metrics_pre$test_dataset)
    metrics$mean_metric = mean_metric 
    
  }else{
    mean_metric  = mean(metrics_pre$val_pearson)
    metrics = data.frame(t(metrics_pre$test_pearson))
    colnames(metrics) =paste0(meas,"_",metrics_pre$test_dataset)
    metrics$mean_metric = mean_metric 
    
  }
  
 
  if(cori == 1){
    corrections_df = metrics
  }else{
    corrections_df = rbind(corrections_df, metrics)
  }
}
corrections_df$corrections = corrections_vec


### INSPECT VAL AUC only keep pca result with highest val
if(any(grepl("pca", corrections_vec))){
  
  grepl("pca",corrections_df$corrections)
  corrections_temp = corrections_df %>% filter(grepl("pca",corrections) & !grepl("33",corrections))
  pca_method_best = corrections_temp$corrections[which.max(unlist(corrections_temp %>% select(mean_metric )))]
  
  
}


corrections_vec = c(other_methods,  paste0("clr_pca33",dtype_pca), pca_method_best)
row.names(corrections_df) = corrections_df$corrections
#corrections_df = corrections_dzf %>% filter(corrections %in% corrections_vec)
corrections_df = corrections_df[corrections_vec,]
#install.packages("gplots")
require(gplots)

nonnice_names =  c("nocorrection","DCC","combat","limma","bmc", "clr", "logcpm","vst",paste0("clr_pca33",dtype_pca),pca_method_best) 
nice_names =  c("Uncorrected","DCC","ComBat","limma","BMC","CLR","logCPM","VST","Fixed PCA Correction","Tuned PCA Correction")
presence_index = which(nonnice_names %in% corrections_vec)
nice_names  = nice_names[presence_index ]

print(corrections_df)


if(lodo == "True"){
  AUC_results = corrections_df[,grepl(meas,colnames(corrections_df))]
  row.names(AUC_results) = corrections_vec
  row.names(AUC_results) = nice_names#"ComBat",
  AUC_results$Average = rowMeans(AUC_results)
  input = as.matrix(AUC_results)
  colnames(input) = gsub(paste0(meas,"_"), "", colnames(input))
  
  
  if(grepl("Thomas",folder )){
    input = input[,c("Average" ,"FengQ_2015", "ThomasAM_2018b",  "ZellerG_2014", "YuJ_2015" ,  "ThomasAM_2018a" , "VogtmannE_2016"  ,"HanniganGD_2017" )]
  }
  if(grepl("Kaplan",folder )){
    input = input[,c("Average","HOWE_KF1",   "HOWE_KF2",   "HOWE_KF3",   "HOWE_KF4"    )]
  }
  
  colnames(input) = gsub("_", " ", colnames(input))
  colnames(input) = gsub("crc ", "", colnames(input))
  firstup <- function(x) {
    substr(x, 1, 1) <- toupper(substr(x, 1, 1))
    x
  }
  
  colnames(input) = sapply(colnames(input),firstup)
  
  input_str = apply(input,2, function(x){sprintf("%.2f",round(x,2))})
  
  ##2B9EDE
  pdf(paste0(data_dir,"/",meas,"LODO_Heatmap_",trans, ".pdf"))
  heatmap.2(t(input), trace="none", density="none", col=colorRampPalette(c("red", "yellow")), cexRow=1.4, cexCol= 1.6, 
            margins = margins_list[[folder]],
            Rowv = FALSE, Colv =  "Rowv",cellnote=t(input_str),notecol="black",srtCol = 45,notecex=notecex_list[[folder]])
  dev.off()
 
 
  # require("ComplexHeatmap")
  # Heatmap(input, name = "mat", rect_gp = gpar(col = "white", lwd = 2),
  #         column_title = "set cell borders")
}

if(lodo == "False"){
  AUC_results = corrections_df[,grepl(meas,colnames(corrections_df))]
  row.names(AUC_results) = corrections_vec
  row.names(AUC_results) = nice_names
  AUC_results$Average = rowMeans(AUC_results)
  input = as.matrix(AUC_results)
  colnames(input) = gsub(meas, "", colnames(input))
  input_str = apply(input,2, function(x){sprintf("%.2f",round(x,2))})
  pdf(paste0(data_dir,"/",meas,"_Heatmap_",trans, ".pdf"))
  heatmap.2(input, trace="none", density="none", col=colorRampPalette(c("white", "blue")), cexRow=1, cexCol=1, margins = c(5,13),
            Rowv = FALSE, Colv =  "Rowv",cellnote=input_str,notecol="black")
  dev.off()
}




# built comparisons
my_comparisons_dcc  = list()
my_comparisons_dcc[[1]] = c("DCC","Fixed PCA Correction")
my_comparisons_uncorr = list()
comparison_candidates = row.names(AUC_results)[!(row.names(AUC_results) %in% c("Fixed PCA Correction","DCC"))]

for(cv in 1:length(comparison_candidates)){
  my_comparisons_uncorr[[cv]] = c("Uncorrected" ,comparison_candidates[cv])
}



require("reshape2")
to_plot = melt(input) %>% filter(Var2 != "Average")
require(ggplot2)
library(ggpubr)
#install.packages("ggpubr")
custom_colors = c('#e32f27',"#C3FFCE",'#FF9300','#FFE800','#fdd0a2',"#9A33FF","#F133FF","#3341FF","#72C1FC","#0093FF")[presence_index ]
#custom_colors = c('#e32f27',"#C3FFCE",'#FF9300','#FFE800','#fdd0a2',"#9A33FF","#F133FF","#3341FF","#72C1FC","#0093FF")
# 
# stat_compare_means(comparisons = my_comparisons_dcc ,bracket.size = 0,ref.group = "DCC",method = "t.test",label = "p.signif",paired=TRUE) + 
#   stat_compare_means(comparisons = my_comparisonas_uncorr, bracket.size = 0,ref.group = "Uncorrected",method = "t.test",label = "p.signif",paired=TRUE) + 
#   
#palette = "jco"
p <- ggboxplot(to_plot, x = "Var1", y = "value",
               fill = "Var1", palette = custom_colors) +xlab("Correction") + 
  
  #stat_compare_means(ref.group = "DCC",method = "t.test",label = "p.signif",paired=TRUE,col = "#808080",vjust=1,method.args = list(alternative = "less")) +
  stat_compare_means(ref.group = "DCC",method = "t.test",label = "p.signif",paired=TRUE,col = "#86B78F",vjust=1,method.args = list(alternative = "greater")) +
  stat_compare_means(ref.group = "DCC",method = "t.test",label = "p.signif",paired=TRUE,col = "#808080",vjust=1,method.args = list(alternative = "less"))
p
  
p <- ggboxplot(to_plot, x = "Var1", y = "value",
               fill = "Var1", palette = custom_colors) +xlab("Correction") + 
  
  #stat_compare_means(ref.group = "DCC",method = "t.test",label = "p.signif",paired=TRUE,col = "#808080",vjust=1,method.args = list(alternative = "less")) +
  stat_compare_means(ref.group = "Uncorrected",method = "t.test",label = "p.signif",paired=TRUE,col = "#e32f27",vjust=1,
                     method.args = list(alternative = "greater"),hide.ns=TRUE,size=10) +
  stat_compare_means(ref.group = "Uncorrected",method = "t.test",label = "p.signif",paired=TRUE,col = "#808080",vjust=1,
                     method.args = list(alternative = "less"),hide.ns=TRUE) +
  #stat_compare_means(ref.group = "DCC",method = "t.test",label = "p.signif",paired=TRUE,col = "#86B78F",vjust=1,method.args = list(alternative = "greater")) +
  # stat_compare_means(ref.group = "DCC",method = "t.test",label = "p.signif",paired=TRUE,col = "#808080",method.args = list(alternative = "greater"))+

  theme(text = element_text(size=13))+
   ylab(paste0("Cross-validated " , meas) ) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1,hjust = 1),text = element_text(size=18),legend.position = "none",
        panel.grid.major.x = element_blank() ,
          # explicitly set the horizontal lines (or they will disappear too)
          panel.grid.major.y = element_line( size=.1, color="black" )) 
#range(to_plot$value)
if(grepl("Kaplanr",folder) & lodo == "True"){
  p<- p + ylim(-0.1,0.3) 
} 
if(grepl("Kaplanr",folder) & lodo == "False"){
    p<- p + ylim(-0.15,0.3) 
}  

if(grepl("AGPr",folder) & lodo == "False"){
  p<- p + ylim(-0.1,0.2) 
}  
   

ggsave(plot=p,filename=paste0(data_dir,"/",meas,"_BOX_lodo",lodo,"_trans_",trans, ".pdf"),width = 7,height = 5,units="in")
saveRDS(p,paste0(data_dir,"/",folder,"_",meas,"_BOX_",trans, "_lodo_",lodo,  ".rds"))




if(lodo == "False"){
  saveRDS(p,paste0(data_dir,"/",meas,"_BOX_",trans, ".rds"))
}  

arrange_time = FALSE
if(arrange_time){
  main_dir ="/u/home/b/briscoel/project-halperin/MicroBatch/data/"
  
  
  folder_list = c("AGPr_complete_otu","Thomasr_complete_otu")
  plot_list = list()
  for(f in 1:length(folder_list)){
    
    data_dir = paste0(main_dir,folder_list[f],"/")
    if(f > 7){
      plot_list[[f]] = readRDS(paste0(data_dir,"/","AUC_BOX_",trans, ".rds")) 
      
      
    }else{
      plot_list[[f]] = readRDS(paste0(data_dir,"/","AUC_BOX_",trans, ".rds")) + theme(axis.title.x=element_blank(),
                                                                                      axis.text.x=element_blank(),
                                                                                      axis.ticks.x=element_blank())
      
      
    }
    
    
  }
  
  ggarrange(bxp, dp, bp + rremove("x.text"), 
            labels = c("A", "B", "C"),
            ncol = 2, nrow = 2)
  
}

# df_means = aggregate(value ~ Var1, data=to_plot,FUN = "mean")
# colnames(df_means) = c("Method","mean")
# df_sd = aggregate(value ~ Var1, data=to_plot,FUN = "sd")
# colnames(df_sd) = c("Method","sd")
# to_plot2 <- merge(df_means,df_sd,by=c("Method"))
# 
# row.names(to_plot2 ) = to_plot2$Method
# to_plot2 = to_plot2[row.names(AUC_results),]
# 
# 
# p <- ggplot(to_plot2,aes(x=Method,y = color = Method))+
#   geom_boxplot(aes(lower=mean-sd,upper=mean+sd,middle=mean,ymin=mean-3*sd,ymax=mean+3*sd),stat="identity")+
#   scale_fill_manual( c('#e32f27',"#C3FFCE",'#FF9300','#FFE800','#fdd0a2',"#72C1FC","#0093FF"))+ 
#   stat_compare_means(comparisons = my_comparisons,method = "wilcox.test",label = "p.signif") + 
#   theme(axis.text.x = element_text(angle = 45, vjust = 1,hjust = 1),legend.position = "none")
# p
#   
print(corrections_df)

print(cbind(rowMeans(corrections_df[,1:(ncol(corrections_df)-1)]),corrections_df[,ncol(corrections_df)]))

range(to_plot$value)


args = commandArgs(trailingOnly=TRUE)
local = TRUE

require(dplyr)
if(local){
  args = c("AGPr_complete_otu","rel","False","AUC")
  
}
dtype_pca = "counts"
extra_trans = "clr_" #"clr_"
margins_list = list()
margins_list[["Gibbonsr_complete_otu"]] =c(16,17)
margins_list[["Gibbonsr_max_k5"]] =c(16,17)
margins_list[["Gibbonsr_max_k7"]] =c(16,17)
margins_list[["Gibbonsr_max_k8"]] =c(16,17)
margins_list[["Thomasr_complete_otu"]] =c(15,13)
margins_list[["Thomasr_max_k7"]] =c(15,13)
margins_list[["AGPr_complete_otu"]] =c(5,13)
margins_list[["AGPr_max_k8"]] =c(5,13)

margins_list[["Kaplanr_complete_otu"]] =c(5,13)
folder_spec = list()
folder_spec[["AGPr_max_k5"]] = list(nocorrection = "AGPr_max_k5",
                                    combat =  "AGPr_max_k5",
                                    DCC =   "AGPr_max_k5",
                                    limma =  "AGPr_max_k5",
                                    bmc =   "AGPr_max_k7",
                                    clr =   "AGPr_max_k5",
                                    clr_pca1counts = "AGPr_max_k7",
                                    clr_pca33counts = "AGPr_max_k7")
folder_spec[["Gibbonsr_max_k5"]] = list(nocorrection = "Gibbonsr_max_k8",
                                    combat =  "Gibbonsr_max_k8",
                                    DCC =   "Gibbonsr_max_k8",
                                    limma =  "Gibbonsr_max_k8",
                                    bmc =   "Gibbonsr_max_k8",
                                    clr = "Gibbonsr_max_k8",
                                    clr_pca1 = "Gibbonsr_max_k8",
                                    clr_pca33 = "Gibbonsr_max_k8")
folder_spec[["Thomasr_max_k6"]] = list(nocorrection = "Thomasr_max_k7",
                                        combat = "Thomasr_max_k7",
                                        DCC =   "Thomasr_max_k7",
                                        limma =  "Thomasr_max_k7",
                                        bmc =   "Thomasr_max_k7",
                                        clr = "Thomasr_max_k7",
                                        clr_pca3counts = "Thomasr_max_k7",
                                        clr_pca33counts = "Thomasr_max_k7")
#folder_spec[["Kaplanr_max_k5"]] = list(test =c(1,2), test2 = c(3,4))
notecex_list= list()
notecex_list[["Gibbonsr_complete_otu"]] = 1
notecex_list[["Gibbonsr_max_otu"]] = 1
margins_list[["Gibbonsr_max_k7"]] =1
margins_list[["Gibbonsr_max_k8"]] =1
notecex_list[["Thomasr_complete_otu"]] = 1
notecex_list[["Thomasr_max_k7"]] = 1
notecex_list[["AGPr_complete_otu"]] = 1
notecex_list[["Kaplanr_complete_otu"]] = 1
notecex_list[["AGPr_max_k8"]] = 1
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



# "bmc", "combat", "percentilenorm", "limma", "DCC", a
#pca_methods = c(paste0("clr_pca",c(1,2,3,4,5,33),dtype_pca)) #,paste0("clr_pca",c(1:5)))

pca_methods = c(paste0(extra_trans,"pca",c(1,33),dtype_pca)) #,paste0("clr_pca",c(1:5)))
#   c("clr_scale_pca",
# "clr_pca1roundcounts", "clr_pca1", "clr_pca2roundcounts", "clr_pca2", "clr_pca3roundcounts", "clr_pca3")
other_methods = c("nocorrection","DCC","combat","limma","bmc","clr") #,"logCPM","vsd" ) # , 
corrections_vec = c(other_methods, pca_methods)
trans_methods = rep(trans,length(corrections_vec))
#trans_methods[(length(trans_methods)-5): length(trans_methods)] = "vsd"
#c("nocorrection","bmc","clr_pca")#"clr_pcacounts","clr_scale_pca","clr_pca")
corrections_list = list()
for(cori in 1:length(corrections_vec)){
  if( folder %in% names(folder_spec)){
    pathy =  paste0(main_dir,folder_spec[[folder]][[corrections_vec[cori]]],"/", "GRID_", meas,"_OUTPUT_" ,trans_methods[cori], "_" , corrections_vec[cori] , "_lodo_" , lodo  , ".csv")
    
  }else{
    pathy  = paste0(data_dir, "GRID_", meas,"_OUTPUT_" ,trans_methods[cori], "_" , corrections_vec[cori] , "_lodo_" , lodo  , ".csv")
    
    
  }
  metrics_pre =  read.csv(pathy, header=TRUE)
  print(pathy)
  corrections_list[[corrections_vec[cori]]] = metrics_pre
  if(meas == "AUC"){
    mean_val_auc = mean(metrics_pre$val_auc)
    
  }else{
    mean_val_auc = mean(metrics_pre$val_f1)
    
  }
  metrics = data.frame(t(metrics_pre$test_auc))
  colnames(metrics) =paste0(meas,"_",metrics_pre$fold)
  metrics$mean_val_auc = mean_val_auc
 
  if(cori == 1){
    corrections_df = metrics
  }else{
    corrections_df = rbind(corrections_df, metrics)
  }
}
corrections_df$corrections = corrections_vec


### INSPECT VAL AUC only keep pca result with highest val
if(any(grepl("pca", corrections_vec))){
  corrections_temp = corrections_df %>% filter(grepl("pca",corrections))
  pca_method_best = corrections_temp$corrections[which.max(unlist(corrections_temp %>% select(mean_val_auc)))]
  
  
}
print(pca_method_best)

corrections_vec = c(other_methods,  paste0(extra_trans,"pca33",dtype_pca), pca_method_best)
row.names(corrections_df) = corrections_df$corrections
corrections_df = corrections_df[corrections_vec,]


#install.packages("gplots")
require(gplots)
nonnice_names =  c("nocorrection","DCC","combat","limma","bmc","clr","logCPM","vsd", paste0(extra_trans,"pca33",dtype_pca),pca_method_best) 
nice_names =  c("Uncorrected","DCC","ComBat","limma","BMC","CLR","logCPM","VSD","Fixed PCA Correction","Tuned PCA Correction")

presence_index = which(nonnice_names %in% corrections_vec)
nice_names  = nice_names[presence_index ]




if(lodo == "True"){
  AUC_results = corrections_df[,grepl(meas,colnames(corrections_df))]
  row.names(AUC_results) = corrections_vec
  row.names(AUC_results) = nice_names #"ComBat",
  AUC_results$Average = rowMeans(AUC_results)
  input = as.matrix(AUC_results)
  colnames(input) = gsub(paste0(meas,"_"), "", colnames(input))
  
  
  if(grepl("Thomas",folder )){
    input = input[,c("FengQ_2015", "ThomasAM_2018b",  "ZellerG_2014", "YuJ_2015" ,  "ThomasAM_2018a" , "VogtmannE_2016"  ,"HanniganGD_2017","Average"  )]
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
  pdf(paste0(data_dir,"/",meas,"LODO_Heatmap_",trans, ".pdf"))
  heatmap.2(input, trace="none", density="none", col=colorRampPalette(c("red","yellow")), cexRow=1.2, cexCol=1.2, 
            margins = margins_list[[folder]],
            Rowv = FALSE, Colv =  "Rowv",cellnote=input_str,notecol="black",srtCol = 45,notecex=notecex_list[[folder]])
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
  heatmap.2(input, trace="none", density="none", col=colorRampPalette(c("red", "yellow")), cexRow=1, cexCol=1, margins = c(5,13),
            Rowv = FALSE, Colv =  "Rowv",cellnote=input_str,notecol="black")
  dev.off()
}




# built comparisons
my_comparisons  = list()
for(cv in 1:(length(row.names(AUC_results))-1)){
  my_comparisons[[cv]] = c(row.names(AUC_results)[1],row.names(AUC_results)[cv+1])
}
require("reshape2")
to_plot = melt(input) %>% filter(Var2 != "Average")
require(ggplot2)
library(ggpubr)
#install.packages("ggpubr")
custom_colors = c('#e32f27',"#C3FFCE",'#FF9300','#FFE800','#fdd0a2',"#72C1FC","#72C1FC","#72C1FC","#72C1FC","#0093FF")[presence_index]


#palette = "jco"
p <- ggboxplot(to_plot, x = "Var1", y = "value",
               fill = "Var1", palette = custom_colors) +xlab("Correction") + 
  theme(text = element_text(size=15))+
   ylab(paste0("Cross-validated ", meas) ) + 
  stat_compare_means(ref.group = "DCC",method = "t.test",label = "p.signif",paired=TRUE,col = "#86B78F",vjust=1) + 
  stat_compare_means(ref.group = "Uncorrected",method = "t.test",label = "p.signif",paired=TRUE,col = "#e32f27") + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1,hjust = 1),legend.position = "none",
        panel.grid.major.x = element_blank() ,
          # explicitly set the horizontal lines (or they will disappear too)
          panel.grid.major.y = element_line( size=.1, color="black" )) 


if(grepl("Gibbonsr",folder) & lodo == "False"){
  p<- p + ylim(0.55,0.93) 
}  
if(grepl("Thomasr",folder) & lodo == "False"){
  p<- p + ylim(0.55,0.90) 
}  
if(grepl("AGPr",folder) & lodo == "False"){
  p<- p + ylim(0.47,0.72) 
}  

ggsave(plot=p,filename=paste0(data_dir,"/",meas,"_BOX_lodo",lodo,"_trans_",trans, ".pdf"),width = 7,height = 5,units="in")

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
corrections_df

print(cbind(rowMeans(corrections_df[,1:(ncol(corrections_df)-1)]),corrections_df[,ncol(corrections_df)]))



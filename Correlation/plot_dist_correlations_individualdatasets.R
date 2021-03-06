# folder_temp = "Gibbonsr_complete_otu"
# trans_temp = "rel_clr"
# data_dir = paste0(main_dir,folder_temp,"/")
# test = readRDS(paste0(data_dir,"CanCorPlotObj_",trans_temp, ".rds"))
# test$CanCorC
require(dplyr)
TOTAL_PCS=15
appendage ="_max_k6" #   "_complete_otu" #     "
folder = c(rep(paste0("Gibbonsr",appendage),2),
           rep(paste0("Thomasr",appendage),2),
           rep(paste0("AGPr",appendage),2),
           rep(paste0("Kaplanr",appendage),2))
pretty_folders = c("CRC-16S","CRC-WGS","AGP","HCHS")
# folder = c(rep("Gibbonsr_complete_otu",2),rep("Kaplanr_complete_otu",2),
#            rep("AGPr_complete_otu",2),rep("Thomasr_complete_otu",2))
trans = rep(c("rel","rel_clr"),4)
trans_pretty =rep(c("None","CLR"),4)
trans_pretty_ = c("None","CLR")
only_confounders = FALSE
only_biology = TRUE
phenotype_range =  list(1,1,1,1,c(1,3),c(1,3),1,1)
confounder_range = list(c(5,6,7),c(5,6,7),
                        c(5,6,7),c(5,6,7),
                        c(8,9),c(8,9),
                        c(4,5,6),c(4,5,6)) # technical confounders only not biological



local = TRUE
main_dir ="/u/home/b/briscoel/project-halperin/MicroBatch/data/"
if(local){
  main_dir ="/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/data/"
  
}
require(reshape2)
require(ggplot2)
for(b in 1:length(folder)){

  data_dir = paste0(main_dir,folder[b],"/")
  plot_object = readRDS(paste0(data_dir,"CanCorPlotObj_",trans[b], ".rds"))
  if(only_confounders){
    title_time = "Technical variables only"
    wavy_melt =melt( plot_object$CanCorC[confounder_range[[b]],1:TOTAL_PCS,drop=FALSE])
  }else if(only_biology){
    title_time = "Phenotype variables only"
    wavy_melt =melt( plot_object$CanCorC[phenotype_range[[b]],1:TOTAL_PCS,drop=FALSE])
  }else{
    title_time = "All variables"
    wavy_melt = melt(plot_object$CanCorC) 
  }
  mean(wavy_melt$value)
  
  colnames(wavy_melt) = c("Variable","PC","Correlation")
  wavy_melt$Transformation = trans_pretty[b]
  wavy_melt$study =  folder[b]

  if(b == 1){
    data_to_plot = wavy_melt
  }else{
    data_to_plot = rbind(data_to_plot,wavy_melt)
  }
  
 
  
}

data_to_plot$Transformation = factor(data_to_plot$Transformation,levels =trans_pretty_)

median(unlist(data_to_plot %>% 
       filter(study == folder[5],Transformation == "CLR") %>% 
       select(Correlation)))

head(data_to_plot)
unique_folder = unique(folder)
plot_list = list()
for( b in 1:length(unique_folder)){
  temp_data_to_plot = data_to_plot %>% filter(study == unique_folder[b])
  x_test = temp_data_to_plot %>% filter(Transformation == "None") %>% select(Correlation)
  y_test = temp_data_to_plot %>% filter(Transformation == "CLR") %>% select(Correlation)
  
  #ks = ks.test(x = as.numeric(x_test$Correlation),y=as.numeric(y_test$Correlation),alternative = "greater")
  head(x_test )
  head(y_test)
  ks = wilcox.test(x = as.numeric(x_test$Correlation),y=as.numeric(y_test$Correlation), paired = TRUE, 
                   alternative = "less")
  
  
  mu = data.frame(mean =c(mean(x_test$Correlation ),mean(y_test$Correlation ) ),Transformation = trans_pretty_)
  #dev.off()
  g <- ggplot(temp_data_to_plot , aes(x=Correlation, fill=Transformation,color=Transformation)) +
    geom_histogram(position="identity", alpha=0.5,binwidth=0.01)+
    geom_vline(data=mu, aes(xintercept=mean, color=Transformation),
               linetype="dashed") + theme_bw() +
    theme(aspect.ratio=1,text = element_text(size=10),axis.text.x = element_text(size=15),
          axis.text.y = element_text(size=15)) +
    annotate(geom = 'text', label = paste0('p-value =\n',round(ks$p.value,4)), 
             x = Inf, y = Inf, hjust = 1.2, vjust = 1.2,size=6) + 
    ggtitle(paste0(unique_folder))
  
  if(only_confounders){
    ggsave(paste0(main_dir,unique_folder[b],appendage,"CorrDist_Confounders.pdf"),plot=g)
  }else if(only_biology){
    ggsave(paste0(main_dir,unique_folder[b],appendage,"CorrDist_Biology.pdf"),plot=g)
  }else{
    ggsave(paste0(main_dir,unique_folder[b],appendage,"CorrDist_AllVars.pdf"),plot=g)
  }
  
  g <- ggplot(temp_data_to_plot , aes(x=Correlation, fill=Transformation,color=Transformation)) +
    geom_histogram(position="identity", alpha=0.5,binwidth=0.01)+
    geom_vline(data=mu, aes(xintercept=mean, color=Transformation),
               linetype="dashed") + theme_bw() +
    theme(aspect.ratio=1,text = element_text(size=5),axis.text.x = element_text(size=5),
          axis.text.y = element_text(size=5),legend.position = "none") +
    annotate(geom = 'text', label = paste0('p-value =\n',round(ks$p.value,4)), 
             x = Inf, y = Inf, hjust = 1.2, vjust = 1.2,size=2.5) + 
    ylab("Number of correlation pairs") + 
    xlab(paste0("PC Correlation in ", pretty_folders[b]))
  
  plot_list[[b]] = g
  
 
  
  
}


# colnames(wavy_melt) = c("Variable","PC","Correlation","Transformation")
# wavy_melt$Transformation = factor(wavy_melt$Transformation,levels = c("None","CLR"))
# p <- ggplot(data_to_plot , aes(x=Correlation, fill=Transformation,color=Transformation)) +
#   geom_density(alpha=0.5,position = "identity") + theme_bw() +
#   theme(aspect.ratio=1,text = element_text(size=18),axis.text.x = element_text(size=15),
#         axis.text.y = element_text(size=15))
# p
#dev.off()

require(gridExtra)
require(ggplot2)
require(cowplot)

coords_x = rep(0.5,4)
coords_y = c(-0.12,0.13,0.38,0.63)
height_ = 0.5
scale1 = 0.5


g <- ggdraw() + draw_plot(plot_list[[1]],coords_x[1], coords_y[1], 0.5, 0.5,scale=scale1) +
  draw_plot(plot_list[[2]],coords_x[2], coords_y[2], 0.5, 0.5,scale = scale1)  + 
  draw_plot(plot_list[[3]],coords_x[3], coords_y[3], 0.5, 0.5,scale = scale1) + 
  draw_plot(plot_list[[4]],coords_x[4], coords_y[4], 0.5, 0.5,scale = scale1)
data_dir = "/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/data/"
ggsave(g,file=paste0(data_dir,"/",appendage,"only_confounders",
                     only_confounders,"only_biology",only_biology,"corr_plot.pdf"),device ="pdf")




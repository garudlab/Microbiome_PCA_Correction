# folder_temp = "Gibbonsr_complete_otu"
# trans_temp = "rel_clr"
# data_dir = paste0(main_dir,folder_temp,"/")
# test = readRDS(paste0(data_dir,"CanCorPlotObj_",trans_temp, ".rds"))
# test$CanCorC
require(dplyr)
folder = c(rep("Gibbonsr_complete_otu",2),rep("Kaplanr_complete_otu",2),
           rep("AGPr_complete_otu",2),rep("Thomasr_complete_otu",2))
trans = rep(c("rel","rel_clr"),4)
trans_pretty =rep(c("None","CLR"),4)
trans_pretty_ = c("None","CLR")
only_confounders = FALSE
only_biology = FALSE
phenotype_range =  list(c(1,1),c(1,1),c(1,1),c(1,1),c(1,3),c(1,3),c(1,1),c(1,1))
confounder_range = list(c(5,6,7),c(5,6,7),c(4,5,6),c(4,5,6),c(8,9),c(8,9),c(5,6,7),c(5,6,7)) # technical confounders only not biological
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
    wavy_melt =melt( plot_object$CanCorC[confounder_range[[b]],1:15,drop=FALSE])
  }else if(only_biology){
    title_time = "Phenotype variables only"
    wavy_melt =melt( plot_object$CanCorC[phenotype_range[[b]],1:15,drop=FALSE])
  }else{
    title_time = "All variables"
    wavy_melt = melt(plot_object$CanCorC) 
  }
  
  colnames(wavy_melt) = c("Variable","PC","Correlation")
  wavy_melt$Transformation = trans_pretty[b]

  if(b == 1){
    data_to_plot = wavy_melt
  }else{
    data_to_plot = rbind(data_to_plot,wavy_melt)
  }
  
 
  
}

data_to_plot$Transformation = factor(data_to_plot$Transformation,levels =trans_pretty_)
x_test = data_to_plot %>% filter(Transformation == "None") %>% select(Correlation)
y_test = data_to_plot %>% filter(Transformation == "CLR") %>% select(Correlation)

ks = ks.test(x = as.numeric(x_test$Correlation),y=as.numeric(y_test$Correlation),alternative = "greater")

# colnames(wavy_melt) = c("Variable","PC","Correlation","Transformation")
# wavy_melt$Transformation = factor(wavy_melt$Transformation,levels = c("None","CLR"))
# p <- ggplot(data_to_plot , aes(x=Correlation, fill=Transformation,color=Transformation)) +
#   geom_density(alpha=0.5,position = "identity") + theme_bw() +
#   theme(aspect.ratio=1,text = element_text(size=18),axis.text.x = element_text(size=15),
#         axis.text.y = element_text(size=15))
# p
dev.off()
mu = data.frame(mean =c(mean(x_test$Correlation ),mean(y_test$Correlation ) ),Transformation = trans_pretty_)
g <- ggplot(data_to_plot , aes(x=Correlation, fill=Transformation,color=Transformation)) +
  geom_histogram(position="identity", alpha=0.5,binwidth=0.01)+
  geom_vline(data=mu, aes(xintercept=mean, color=Transformation),
             linetype="dashed") + theme_bw() +
  theme(aspect.ratio=1,text = element_text(size=18),axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15)) +
  annotate(geom = 'text', label = paste0('KS p-value =\n',round(ks$p.value,4)), 
           x = Inf, y = Inf, hjust = 1.2, vjust = 1.2,size=6) + 
  ggtitle(title_time)

g

if(only_confounders){
  ggsave(paste0(main_dir,"CorrDist_Confounders.pdf"),plot=g)
}else if(only_biology){
  ggsave(paste0(main_dir,"CorrDist_Biology.pdf"),plot=g)
}else{
  ggsave(paste0(main_dir,"CorrDist_AllVars.pdf"),plot=g)
}

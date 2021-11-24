#DiseaseState
gen_path ="/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/data/" 

paths = c(   paste0(gen_path,"Gibbonsr_complete_otu/pca_plot_rel_study.rds"),
             paste0(gen_path,"Gibbonsr_complete_otu/pca_plot_rel_clr_study.rds"),
             paste0(gen_path,"Gibbonsr_max_k6/pca_plot_rel_study.rds"),
             paste0(gen_path,"Gibbonsr_max_k6/pca_plot_rel_clr_study.rds"),
             paste0(gen_path,"Thomasr_complete_otu/pca_plot_rel_dataset_name.rds"),
             paste0(gen_path,"Thomasr_complete_otu/pca_plot_rel_clr_dataset_name.rds"),
             paste0(gen_path,"Thomasr_max_k6/pca_plot_rel_dataset_name.rds"),
             paste0(gen_path,"Thomasr_max_k6/pca_plot_rel_clr_dataset_name.rds"),
             paste0(gen_path,"Kaplanr_complete_otu/pca_plot_rel_extraction_robot..exp..rds"),
             paste0(gen_path,"Kaplanr_complete_otu/pca_plot_rel_clr_extraction_robot..exp..rds"),
             paste0(gen_path,"Kaplanr_max_k6/pca_plot_rel_extraction_robot..exp..rds"),
             paste0(gen_path,"Kaplanr_max_k6/pca_plot_rel_clr_extraction_robot..exp..rds"),
             paste0(gen_path,"AGPr_complete_otu/pca_plot_rel_Instrument.rds"),
             paste0(gen_path,"AGPr_complete_otu/pca_plot_rel_clr_Instrument.rds"),
             paste0(gen_path,"AGPr_max_k6/pca_plot_rel_Instrument.rds"),
             paste0(gen_path,"AGPr_max_k6/pca_plot_rel_clr_Instrument.rds"))

# paths = c("/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/data/Thomasr_max_k6/pca_plot_rel_dataset_name.rds",
#           "/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/data/Thomasr_max_k6/pca_plot_rel_clr_dataset_name.rds")

require(gridExtra)
require(ggplot2)
require(cowplot)


plot_list = list()
for(p_num in c(1:length(paths))){
  plot_list[[p_num]] = readRDS(paths[p_num]) 
}

require("ggplotify")
width_list = list()
for(p_num in 1:length(paths)){
  print(p_num)
  #p_num=1
  p_temp = plot_list[[p_num]] + theme(text = element_text(size=5),plot.title = element_blank(),
                                         axis.text.x = element_text(size=5),
                                         axis.text.y = element_text(size=5),
                                      legend.text=element_text(size=1))
  
  q <- ggplot_build(p_temp)
  
  q$data[[1]]$size = 0.2
  
  q <- ggplot_gtable(q)
  #plot(q)
  p_temp = as.ggplot(q)
  
  if(grepl("k6",paths[p_num]) & grepl("rel_clr",paths[p_num])){
    print('nuthin')
    width_list <- c(width_list,1)
  }else{
    p_temp <- p_temp  + theme(legend.position = "none")  
    width_list <- c(width_list,0.5)
  }
  #p_temp <- p_temp + geom_point(size = 0.2) 
  assign(paste0("p",p_num),p_temp)
  #p1 + theme(legend.position = "none")  
  #p_temp
  
}

coords_x = rep(c(-0.07,0.18,0.43,0.66),4)
coords_y = c(rep(-0.10,4),rep(0.15,4),rep(0.40,4),rep(0.63,4))
height_ = 0.5
#p1 
# p1 <- p1 + theme(legend.position = "none",plot.margin = margin(2,.8,2,.8, "cm"))
# p2 <- p2 + theme(plot.margin = margin(2,.8,2,.8, "cm"))
# grid.arrange(p1, p2, nrow = 1)
scale1 = 0.7
scale2 = 0.73
g <- ggdraw() + draw_plot(p1,coords_x[1], coords_y[1], 0.5, 0.5,scale=scale1) +
  draw_plot(p2,coords_x[2], coords_y[2], 0.5, 0.5,scale = scale1)  + 
  draw_plot(p3,coords_x[3], coords_y[3], 0.5, 0.5,scale = scale1) + 
  draw_plot(p4,coords_x[4], coords_y[4], 0.5, 0.5,scale = scale1)+ 
  draw_plot(p5,coords_x[5], coords_y[5], 0.5, 0.5,scale = scale2)  + 
  draw_plot(p6,coords_x[6], coords_y[6], 0.5, 0.5,scale = scale2)  + 
  draw_plot(p7,coords_x[7], coords_y[7], 0.5, 0.5,scale = scale2)  + 
  draw_plot(p8,coords_x[8], coords_y[8], 0.5, 0.5,scale = scale2) + 
  draw_plot(p9,coords_x[9], coords_y[9], 0.5, 0.5,scale = scale1) +
  draw_plot(p10,coords_x[10], coords_y[10], 0.5, 0.5,scale = scale1)  + 
  draw_plot(p11,coords_x[11], coords_y[11], 0.5, 0.5,scale = scale1)  + 
  draw_plot(p12,coords_x[12], coords_y[12], 0.5, 0.5,scale = scale1) +
  draw_plot(p13,coords_x[13], coords_y[13], 0.5, 0.5,scale = scale1)  + 
  draw_plot(p14,coords_x[14], coords_y[14], 0.5, 0.5,scale = scale1)  + 
  draw_plot(p15,coords_x[15], coords_y[15], 0.5, 0.5,scale = scale1)  + 
  draw_plot(p16,coords_x[16], coords_y[16], 0.5, 0.5,scale = scale1) 
data_dir = "/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/data/"
ggsave(g,file=paste0(data_dir,"/pca_plot.pdf"),device ="pdf")

  
  

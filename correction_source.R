scale_custom <- function(mat){
  mat_means = apply(mat,2,mean)
  mat_sd = apply(mat,2,sd)
  
  mat_sub_mean = sweep(mat,2,mat_means,FUN = "-" )
  mat_sub_mean_div_sd = sweep(mat_sub_mean,2,mat_sd,FUN = "/" )
  
  return(mat_sub_mean_div_sd)
  
}


pca_method <- function(input,num_pcs){
  require("bigstatsr")
  colnames_orig_input = colnames(input)
  myFBM = as_FBM(t(input), type = c("double"))
  
  
  t1 = Sys.time()
  
  svd_result = big_SVD(myFBM,k=num_pcs)

  print(Sys.time()-t1)
  
  pca_score <-  svd_result$u %*% diag(svd_result$d)
  print("eigen values")
  print(svd_result$d)
  eigenvalues = svd_result$d
  
  row.names(pca_score) = colnames_orig_input
  return(list(pca_score,eigenvalues))
}

regress_out <- function(pc_score,data,index_pcs_correct){
  
  model_residuals<-lm(as.matrix(data) ~ pc_score[,index_pcs_correct] ) 
  
  extracted_residuals <- residuals(model_residuals)
  return(t(extracted_residuals))
  
  
  #pc_score = pca_score
  #data = t(feature_table)
  #index_pcs_correct = c(1:num_pcs_regress)

  #test = cbind(1,pc_score[,index_pcs_correct]) %*% model_residuals$coefficients  + extracted_residuals 
  #dim(test)
  #dim(data)
  #test[1:4,1:4]
  #data[1:4,1:4]
  
}

pca_plot <- function(df,meta,title,group_column,coord1,coord2){
  
  # df = pca_score
  # meta = metadata_table
  # title = "PC 1 and 2"
  # group_column=group_column
  # coord1=1
  # coord2=2
  
  require(ggplot2)
  to_plot = data.frame(plotx = df[,coord1], ploty = df[,coord2], group = meta[,group_column])
  
  p<-ggplot(to_plot,aes(x=plotx,y=ploty,color=group)) + ggtitle(title) 
  p<-p + geom_point(size = 2.5) + theme_bw()  + xlab(paste0("PC",coord1)) + ylab(paste0("PC",coord2))
  p <-p + coord_fixed(ratio=1) + 
    theme(aspect.ratio=1,text = element_text(size=15),axis.text.x = element_text(size=15),
          axis.text.y = element_text(size=15),legend.text=element_text(size=15))#,legend.title = element_text("Study"))
  return(p)
    #dev.off()
}


percentile_norm <- function(df_otu,df_meta,replace_zeroes,case_class,control_class){
  require(dplyr)
  # df_otu = subset_df_otu_rel_ab
  # df_meta = subset_df_meta
  # replace_zeroes = FALSE
  # 
  
  #df_otu = df_otu_rel_ab
  #df_meta = df_meta
  #replace_zeroes=FALSE
  
  
  set.seed(3)
  dat = df_otu
  meta = df_meta
  
  #replace_zeroes=FALSE
  # replace zeros with runif
  if(replace_zeroes){
    df_zero_replace = apply(dat,2,function(x){
      replace(x,which(x == 0),runif(length(which(x == 0)),0,10^-9))
    })
    dat = df_zero_replace
  }else{
    dat  = df_otu
  }
  
  studies <- meta %>% pull(study) %>% unique
  
  dat_percentile_collection = list()
  for (s in studies) {
    #s= "crc_baxter"
    #s = "CCIS"
    
    dat_in_study <- dat[, as.character(meta %>% filter(study==s) %>% pull(Sample_ID)),drop=FALSE]
    controls <- dat[, as.character(meta %>% filter(study==s) %>% 
                                     filter(DiseaseState==control_class) %>% pull(Sample_ID)),drop=FALSE]
    n_control = ncol(controls)
    dat_percentile = dat_in_study
    
    for(tax in 1:nrow(dat_in_study )){
      #tax = 1
      dat_percentile[tax,] = unlist(sapply(dat_in_study[tax,],function(sample){
        #sample = dat_in_study[tax,1]
        (sum(controls[tax,,drop=FALSE] < sample)/n_control + (1-sum(controls[tax,,drop=FALSE] > sample)/n_control))/2
      }))
      
    }
    dat_percentile_collection[[s]] = dat_percentile
    print(paste0(s," percentile norm complete"))
    
    
  }
  concat_dat = do.call(cbind,  dat_percentile_collection)
  colnames(concat_dat) <- colnames(dat)
  return(concat_dat*100)
}


correct_limma <- function(mat,batch_labels,batch_labels2 = NULL){
  require(limma)
  #?removeBatchEffect
  input = removeBatchEffect( x=mat , batch= batch_labels,batch2 = batch_labels2)
  return(input)
}
correct_bmc <- function(mat,batch_labels){
  
  
  require(dplyr)
  corrected_mat = mat
  unique_batches= unique(batch_labels)
  for( b in 1:length(unique_batches)){
    samples = colnames(mat)[batch_labels == unique_batches[b]]
    batch_mat = mat[,samples]
    corrected_mat[,samples] = sweep(mat[,samples],MARGIN = 1, rowMeans(batch_mat))
  }
  
  return(corrected_mat)
  
}
correct_ComBat <- function(mat, batch_labels,model_matrix=NULL){
  require(sva)
  
  
  #mat <- data$df_otu_corrected
  #range(mat)
  # make continuous
  #log(mat + 1) #mat #
  input = ComBat( dat=mat, batch = batch_labels,mod = model_matrix)
  return(input)
}

correct_DCC <- function(mat, batch_labels){
  batch_labels_factor = factor(batch_labels)
  batch_mat = model.matrix(~batch_labels_factor )
  return(t(resid(lm(t(mat) ~ batch_mat))))
  
}


#output_folder = paste0(output_folder, "_",methods_list[m])
#dir.create(output_folder)
#saveRDS(new_metadata,paste0(output_folder,"/metadata.rds"))
#write.table(new_metadata,paste0(output_folder,"/metadata.txt"),sep="\t",quote=FALSE)



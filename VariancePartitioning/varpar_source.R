process_model_matrix <- function(total_metadata =NULL,binary_vars=NULL,categorical_vars = NULL,numeric_vars = NULL,integer_vars = NULL,
                                 label_pos_or_neg=NULL,target_label=NULL){
  #,
  for(b_v in binary_vars){
    #b_v = "antibiotic"
    #target_label = "1"
    #label_pos_or_neg = 1
    
    data_na_included = as.character(total_metadata[,b_v])
    data_na_included[data_na_included == "Other" | data_na_included == "Not provided" | data_na_included == "other" 
                     | data_na_included == '' | data_na_included == 'not applicable' | data_na_included == 'not provided'] = NA
    #temp = to.dummy(as.factor(data_na_included),paste0(b_v,"_"))[,-1,drop=FALSE]
    
    if(length(target_label) > 0){
      if( label_pos_or_neg == 1){
        print("Positive")
        data_na_included = sapply(data_na_included,function(x){
          if(is.na(x)){return(x)}
          else if(x %in% target_label){return(1)}
          else{return(0)}
        })
      }else{
        data_na_included = sapply(data_na_included,function(x){
          if(is.na(x)){return(x)}
          else if(x %in% target_label){return(0)}
          else{return(1)}
        })
      }
      
      
    }
    #might keep
    #temp = as.integer(as.factor(data_na_included))
    temp = as.character(as.integer(as.factor(data_na_included))-1)
    #print(colnames(temp))
    assign(b_v ,temp)
  }
  for(c_v in categorical_vars){
    #print(table(total_metadata[,c_v]))
    data_na_included = as.character(total_metadata[,c_v])
    data_na_included[data_na_included == "Other" | data_na_included == "Not provided" | data_na_included == "other" 
                     | data_na_included == '' | data_na_included == 'not applicable' | data_na_included == 'not provided'
                     | data_na_included == 'Not applicable'| data_na_included == 'Unspecified'] = NA
    #temp = to.dummy(as.factor(data_na_included),paste0(c_v,"_"))[,-1,drop=FALSE]
    #print(colnames(temp))
    assign(c_v ,data_na_included)
    
  }
  
  
  for(n_v in numeric_vars){
    #print(sort(table(as.numeric(total_metadata[,n_v])),decreasing=TRUE)[1:3])
    assign(n_v ,as.numeric(total_metadata[,n_v]))
  }
  
  for(i_v in integer_vars){
    #print(sort(table(as.numeric(total_metadata[,n_v])),decreasing=TRUE)[1:3])
    assign(i_v ,as.integer(total_metadata[,i_v]))
  }
  
  list_vars = mget(c(binary_vars,categorical_vars,numeric_vars,integer_vars))
  df_vars =data.frame(list_vars)
  row.names(df_vars) = row.names(total_metadata)
  return(df_vars)
}
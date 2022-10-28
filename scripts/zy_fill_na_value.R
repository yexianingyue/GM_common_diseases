library(dplyr)

# 针对表中的数值类型，删除NA值大于cutoff的列，且使用中位数填充NA
zy_fill_na_value <- function(in_matrix=NA, cutoff=1/4, del_na=TRUE, FUN = median){
  # 只处理数值类型
  del_count = 0
  total_num = nrow(in_matrix)
  for(i in 1:ncol(in_matrix)){
    na_num = sum(is.na(in_matrix[,i]))
    fill_value = FUN(in_matrix[,i], na.rm=T)
    if (na_num/total_num < cutoff){
      in_matrix[,i][is.na(in_matrix[,i])] = fill_value
    }else{
      del_count = del_count + 1
    }
  }
  # 储存将要删除的列
  if (isTRUE(del_na)){
    del_col = c()
    for(i in 1:ncol(in_matrix)){
      if(is.numeric(in_matrix[,i])){
        if(is.na(sum(in_matrix[,i]))){
          del_col = c(del_col, i)
        }
      }
    }
    message("delete ",del_count," items")
    return(as.data.frame(in_matrix[,-del_col], check.names=F))
  }else{
    return(as.data.frame(in_matrix, check.names=F))
  }
}

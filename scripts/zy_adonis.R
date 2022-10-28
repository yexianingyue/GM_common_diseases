get_adjusted_r2 <- function(adonis_object) {
    # n_observations <- ncol(adonis_object$coef.sites) # 用于adonis
    n_observations = tail(adonis_object$Df,1)+1 # 用于adonis2
    d_freedom <- adonis_object$Df[1]
    r2 <- adonis_object$R2[1]
    adjusted_r2 <- RsquareAdj(r2, n_observations, d_freedom)
    adjusted_r2
}
zy_adonis <- function(query=NA, target=NA, method="bray"){
  # 对于target每一列一个样本
  # 对于query，每行一个样本
  names_ = colnames(query)
  nrow_q = nrow(query)
  result <- matrix(NA,ncol=4,nrow=length(names_),
                   dimnames = list(names_, c("name","r2","pvalue","adjust.R2")))
  for(c_ in names_){
    message(c_)
    ngs = length(query[,c_])
    if(ngs == 1 | ngs==nrow_q){
        next
    }
    if(is.character(query[,c_])){
        query[,c_] = as.factor(query[,c_])
    }
    x = na.omit(query[,c_])
    rm_index = attr(x,"na.action")
    # 针对query， 删除有NA值的样本
    if( !is.null(rm_index)){
      ado = adonis2(t(target[,-rm_index])~query[-rm_index,c_], method = method)
    }else{
      ado = adonis2(t(target)~query[,c_], method = method)
    }
    r2 = ado$R2[1]
    p = ado$`Pr(>F)`[1]
    q2 = get_adjusted_r2(ado)
    result[c_,] = c(c_, r2,p, q2)
  }

  result[,2:4] = apply(result[,2:4],2,as.numeric)
  result = as.data.frame(result)
  result
}

zy_parallel_adonis <- function(query=NA, target=NA, method="bray", p=2){
  # 对于target每一列一个样本
  # 对于query，每行一个样本
  library(parallel)
  names_ = colnames(query)
  result <- matrix(NA,ncol=4,nrow=length(names_),
                   dimnames = list(names_, c("name","r2","pvalue","adjust.R2")))
  for(c_ in names_){
    message(c_)
    ngs = length(query[,c_])
    if(ngs == 1 | ngs==nrow_q){
      next
    }
    if(is.character(query[,c_])){
        query[,c_] = as.factor(query[,c_])
    }
    x = na.omit(query[,c_])
    rm_index = attr(x,"na.action")
    # 针对query， 删除有NA值的样本
    if( !is.null(rm_index)){
      ado = adonis2(t(target[,-rm_index])~query[-rm_index,c_], method = method)
    }
    else{
      ado = adonis2(t(target)~query[,c_], method = method)
    }
    r2 = ado$R2[1]
    p = ado$`Pr(>F)`[1]
    q2 = get_adjusted_r2(ado)
    result[c_,] = c(c_, r2,p, q2)
  }
  
  result = as.data.frame(result)
  result[,2:4] = apply(result[,2:4], 2, as.numeric)
  result
}

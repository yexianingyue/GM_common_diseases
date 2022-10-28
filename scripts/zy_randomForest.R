library(randomForest, quietly = T, warn.conflicts=T)
library(tidyselect, quietly = T)
library(dplyr, quietly = T)


zy_Leave_One_Out_RandomFroestClassification <- function(rf_dt=NA, rf_map=NA,
                                          zy_sample="zy_RF_temp_ID", group=NA,
                                          nspecies=NA,seed=123,
                                          ntree=999){

  # 是否给定物种数量
  if(! is.na(nspecies)){
    rf_dt = rf_dt[1:nspecies,] %>% filter(rowSums(rf_dt) !=0)
  }
  set.seed(seed)
  # 用于二分类
    rf_map[,group] = as.factor(rf_map[,group])
    rf_dt = rf_dt[,rf_map[,zy_sample]]
    # rf_map$zy_temp_Group = ifelse(rf_map[,group] == Control,"Control","Disease")
    temp_result = list()
    x = randomForest(x= t(rf_dt), y = rf_map[,group], ntree=ntree,importance=TRUE, proximity=TRUE)
    for(i in 1:ncol(rf_dt)){
        message("sample: ", i)
        test_sample = rf_map[i,]
        test_dt = rf_dt[,test_sample[,zy_sample]]
        train_sample = rf_map[-i,]
        train_dt = rf_dt[,train_sample[,zy_sample]]
        message("train")
        #fit = randomForest(train_sample[,group]~.,data=t(train_dt), ntree=ntree,importance=TRUE, proximity=TRUE)
        fit = randomForest(x= t(train_dt), y = train_sample[,group], ntree=ntree,importance=TRUE, proximity=TRUE)
        message("test\n")
        pred = as.data.frame(predict(fit, t(test_dt), type='prob'))
        pred[,zy_sample] = test_sample[,zy_sample]
        pred[,group] = test_sample[,group]
        temp_result = append(temp_result, list(pred))
    }
    temp_result <- do.call("rbind", temp_result)
    list(pred=temp_result, imporv = x$importance)
}


zy_format_class_name <- function(rf_dt=NA, rf_map=NA, zy_sample=NA){
  row.names(rf_dt) = make.names(row.names(rf_dt))
  colnames(rf_dt) = make.names(colnames(rf_dt))
  rf_map$zy_RF_temp_ID = make.names(rf_map[,zy_sample])
  rf_map[,zy_sample] = make.names(rf_map[,zy_sample])
  rf_dt = rf_dt[,rf_map$zy_RF_temp_ID]
  return(list(rf_dt=rf_dt, rf_map=rf_map))
}

zy_RF <- function(train_dt=NA, train_map=NA, test_dt=NA, test_map=NA, 
                  ntree=999,seed=123,
                  zy_sample="zy_RF_temp_ID", group=NA,
                  test_sample=NA, test_group=NA,
                  train_sample=NA, train_group=NA){
  if(!is.na(zy_sample)){
    test_sample = train_sample = zy_sample
    test_group = train_group = group
  }
  group_fact = unique(test_map[,test_group])
  test_map[,test_group] = factor(test_map[,test_group], levels=group_fact)
  train_map[,train_group] = factor(train_map[,train_group], levels=group_fact)
  train_dt = train_dt[, train_map[,train_sample]]
  test_dt  = test_dt[, test_map[, test_sample]]
  set.seed(seed)
  rf <- randomForest(train_map[,train_group]~., data=t(train_dt), importance=T, proximity=T,ntree=ntree )
  pred <- as.data.frame(predict(rf,t(test_dt),type = "prob"))
  result = merge(test_map[,c(test_sample,test_group)], pred, by.x=test_sample, by.y="row.names")
  result
}


zy_RF_two_class <- function(rf_dt=NA, rf_map=NA, 
                            zy_sample="zy_RF_temp_ID", group=NA, 
                            ntree=999, cross_n = 10,
                            nspecies = NA,
                            seed=123){
  # 返回每个样本的预测值
  #
  # 是否给定物种数量
  if(! is.na(nspecies)){
    rf_dt = rf_dt[1:nspecies,] %>% filter(rowSums(.) !=0)
  }
  set.seed(seed)
  gs = rf_map %>%
    dplyr::group_by(across({{group}})) %>%
    summarise(value=n()) %>%
    as.data.frame()
  
  g1 <- rf_map %>%
    filter(across({{group}})==gs[1,1]) %>%
    mutate(rf_temp_cross_n=rep(sample(1:cross_n), gs[1,2]/cross_n+1)[1:gs[1,2]])
  
  g2 <- rf_map %>%
    filter(across({{group}})==gs[2,1]) %>%
    mutate(rf_temp_cross_n = rep(sample(1:cross_n), gs[2,2]/cross_n+1)[1:gs[2,2]])
  
  rf_map = rbind(g1,g2)
  rf_map[,group] = as.factor(rf_map[,group])
  
  predict_result = list() # 储存预测的结果

  for(i in 1:cross_n){
    test_sample = rf_map[rf_map$rf_temp_cross_n == i,]
    test_dt = rf_dt[,test_sample[,zy_sample]]
    
    train_sample = rf_map[rf_map$rf_temp_cross_n != i,]
    train_dt = rf_dt[,train_sample[,zy_sample]]
    fit = randomForest(train_sample[,group]~.,data=t(train_dt), ntree=ntree,importance=TRUE, proximity=TRUE)
    pred = as.data.frame(predict(fit, t(test_dt), type='prob'))
    predict_result = append(predict_result, list(pred)) # 储存每个样本预测到的结果
  }
  predict_result <- do.call("rbind", predict_result)
  predict_result = merge(rf_map[,c(zy_sample, group)], predict_result, by.x=zy_sample, by.y="row.names")
  predict_result
}

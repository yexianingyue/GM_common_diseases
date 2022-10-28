library(ggplot2)
library(pROC)
library(dplyr)
library(data.table) # SetDT
map_name <- function(roc.list){
  # 返回映射后的新名字
  oc = c() # old name
  nc = c() # new name
  for(rb in names(roc.list)){
    # b = signif(ci(roc.list[[rb]], of="auc")*100, digits=3)
    b = sprintf("%0.1f",ci(roc.list[[rb]], of="auc")*100)
    c = paste(rb, " (",b[2], "%)\t95% CI: " , b[1],"%-",b[3],"%", sep="")
    oc = c(oc,rb)
    nc = c(nc,c)
  }
  names(nc) = oc
  nc
}

calc_auc <- function(dt, pred=NA, true=NA, group=NA,acc=F,
                     boot_n=2000){

  roc.list = list()
  if(is.na(group)){
      roc.list['AUC'] = list( roc(dt[,true], dt[,pred]))
      grps = "AUC"
  }else{
      grps = unique(dt[,group])
      if(length(grps) == 1){
        roc.list['AUC'] = list(roc(dt[,true], dt[,pred]))
      }else{
        for(g in grps){
          temp_dt = dt[dt[,group]==g,]
          roc.list[as.character(g)] = list(  roc(temp_dt[,true], temp_dt[,pred]))
        }
      }
  }
  names_ = names(roc.list)
  result_auc = matrix(NA, ncol=6,nrow=length(grps),
                      dimnames=list(names_, c("low","auc","high","low_acc","acc","high_acc")))
  for(i in names_){
    b = sprintf("%0.4f",ci(roc.list[[i]], of="auc")*100)
    result_auc[i,c(1,2,3)] = as.numeric(b)
    if(isTRUE(acc)){
      ac = ci.coords(roc.list[[i]], x="best", ret="accuracy", transpose=F)
      ac = sprintf("%0.4f",unlist(ac)*100)
      ac[2] = sprintf("%0.4f",coords(roc.list[[i]], x="best", ret="accuracy", transpose=F)*100)
      result_auc[i,c(4,5,6)] = ac
    }
  }
  list(table=as.data.frame(result_auc))
}

plot_roc <- function(dt, pred=NA, true=NA, group=NA,
                     fill=FALSE,
                     cols = NA, conf_level=0.95, boot_n=2000){
  # 为了防止曾经修改过小数位等，这边先备份，最后再恢复
  old_scipen = getOption("scipen")
  old_digits = getOption("digits")
  options(scipen=0)
  options(digits=7)
  
  roc.list = list()
  if(is.na(group)){
    if(is.na(cols)){cols = "darkblue"}
    roc.list['AUC'] = list(roc(dt[,true], dt[,pred]))
  }else{
    grps = unique(dt[,group])
    if(is.na(cols)){
      cols=c(1:length(grps))
    }
    for(g in grps){
      temp_dt = dt[dt[,group]==g,]
      roc.list[as.character(g)] = list(  roc(temp_dt[,true], temp_dt[,pred]))
    }
  }
  new_name_map <- map_name(roc.list)
  p <- ggroc(roc.list)+
    theme_bw()+
    geom_abline(slope=1, intercept=1,
                linetype="dashed",color="gray",alpha=0.7)+
    scale_color_manual(values=cols, labels=new_name_map)
  
  if(isTRUE(fill)){
    ci.list <- lapply(roc.list, function(rocobj)
      setDT(
        data.frame(
          ci.se(rocobj, specificities=seq(0, 1, 0.1)), check.names=F)
        ,keep.rownames = T)
    )
    data_ci <- bind_rows(ci.list, .id="plot_group")
    data_ci$rn = as.numeric(data_ci$rn)
    p <- p+
      geom_ribbon(data=data_ci,aes(x=rn, ymin=`2.5%`, ymax=`97.5%`, fill=plot_group),
                  alpha=.3,
                  inherit.aes = F)+ # 必须有参数inherit.aes
      scale_fill_manual(values=cols, labels=new_name_map)
  }
  
  # 设置为原来的小数位等
  options(scipen = old_scipen)
  options(digits = old_digits)
  list(plot=p, ROC=roc.list, labels=new_name_map)
}


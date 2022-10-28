# 已经算好的fold-change去计算，因为计算qvalue时会有过滤
zy_table_to_volcano = function(){
  
}

taxo_sample_to_grps <- function(sample_map=NA, group=NA,ID=NA){
  grps = unique(sample_map[,group])
  grps_list = list()
  for(grp in grps){
    gg <- sample_map[which(sample_map[,group] == grp), ID]
    grps_list[[grp]] = gg
  }
  return(grps_list)
}

# 直接用院士表计算，使用pvalue表示
zy_raw_profile_to_volcano = function(dt=NA, sample_map=NA, group=NA, ID=NA, cutoff=10){
  # ID -> ID columns name
  # gorup -> how to group data
  # dt -> profile
  # sample_map -> mapping file
  # cutoff -> fold-change 
  dt = dt[, sample_map[,ID]]
  grps <- taxo_sample_to_grps(sample_map, "Group", "Sample")
  com = t(combn(names(grps),2))
  nspecies = nrow(dt)
  names = rownames(dt)
  result = matrix(NA,nrow = nrow(com)*nspecies, ncol = 8,
                  dimnames = list(NULL,c("name","g1","g2","m1","m2","enriched","fold_change","pvalue")))
  nr = 1
  for (n in 1:nspecies){
    temp_dt = dt[n,]
    for(c in 1:nrow(com)){
      g1 = com[c,1]
      g2 = com[c,2]
      dt1 = as.matrix(temp_dt[,grps[[g1]]])
      dt2 = as.matrix(temp_dt[,grps[[g2]]])
      m1 = mean(dt1)
      m2 = mean(dt2)
      enrich = ifelse(m1>m2,g1,g2)
      fold = max(m1/m2,m2/m1)
      p = wilcox.test(dt1,dt2)$p.value
      result[nr,] = c(names[n], g1, g2, m1, m2,enrich,fold,p)
      nr = nr+1
    }
  }
  result = as.data.frame(result)
  for(x in c("m1","m2","fold_change","pvalue")){
    result[,x] = as.numeric(result[,x])
  }
  result
}

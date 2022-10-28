library(dplyr)

zy_pvalue = function(dt=NA, sample_map=NA, group=NA, ID=NA){
    # ID -> ID columns name
    # gorup -> how to group data
    # dt -> profile
    # sample_map -> mapping file
    dt = dt[, sample_map[,ID]]
    grps = unique(sample_map[,group])
    com = t(combn(grps,2))
    nspecies = nrow(dt)
    names_ = rownames(dt)
    # Avg -> 平均数
    # Avg.weighted.g1 -> 这个分组的加权平均数
    result = data.frame(matrix(NA,nrow = nrow(com)*nspecies, ncol = 16,
                    dimnames = list(NULL,c("name","g1","g2","Avg.g1","Avg.g2","fold_change","enriched",
                                           "Avg.weighted.g1","Avg.weighted.g2","all.avg","all.var","pvalue",
                                           "count1","count2",
                                           "rank1.avg", "rank2.avg"))))
    nr = 1
    for (n in 1:nspecies){
        temp_dt = dt[n,]
        for(c in 1:nrow(com)){
            g1 = com[c,1]
            g2 = com[c,2]
            g1s = sample_map[which(sample_map[,group] == g1), ID]
            g2s = sample_map[which(sample_map[,group] == g2), ID]
            dt1 = as.matrix(temp_dt[,g1s])
            dt2 = as.matrix(temp_dt[,g2s])
            c1 = sum(dt1 != 0 )
            c2 = sum(dt2 != 0)
            m1 = mean(dt1)
            m2 = mean(dt2)
            ag1 = sum(dt1)/c1
            ag2 = sum(dt2)/c2
            am = mean(c(dt1,dt2))
            a_var=var(c(dt1,dt2))
            p = wilcox.test(dt1,dt2)$p.value
            fold_change = ifelse(m1>m2, m1/m2, m2/m1)
            enriched = ifelse(m1>m2, g1,g2)
            
            m = sample_map[which(sample_map[,group] %in% c(g1,g2)), ID]
            all_rank = rank(temp_dt[,m])
            rank1 = all_rank[colnames(dt1)]
            rank2 = all_rank[colnames(dt2)]
            rank1.avg = mean(rank1)
            rank2.avg = mean(rank2)

            result[nr,] = c(names_[n], g1, g2, m1, m2, fold_change, enriched, ag1, ag2, am,a_var, p, c1, c2,rank1.avg, rank2.avg)
            #data.frame(name = names[n],g1=g1, g2=g2,mean1 = m1, mean2=m2,pvalue=p, count1=c1, count2= c2)
            nr = nr+1
        }
    }
    result
}

zy_qvalue = function(dt=NA, sample_map=NA, group=NA, ID=NA,
                     method="BH",min_count=0, min_avg=0,min_fd=0){
  # min_count 至少有一个分组有这么多样本
  # avg 总体的平均含量阈值
  # fd fold-change阈值
  result <- as.data.frame(zy_pvalue(dt, sample_map, group=group, ID=ID))
  result[,c(4:6,8:15)] = lapply(result[,c(4:6,8:15)], as.numeric)
  result$qvalue = NA
  # 筛选出要计算qvalue的数据
  result1 <- result %>%
    dplyr::filter((count1 >= min_count | count2 >= min_count) & fold_change >= min_fd & all.avg >= min_avg)
  result1$qvalue = p.adjust(result1$pvalue, method=method)
  # 不计算qvalue的数据与计算完的qvalue数据合并
  result2 <- result %>%
    dplyr::filter((count1 < min_count & count2 < min_count) | fold_change < min_fd | all.avg < min_avg) %>%
    rbind(result1)
  result2
}

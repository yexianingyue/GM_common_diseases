# 健康人富集的菌求和
# 病人富集的菌求和
# 相减/除

dt = read.table("../00.data/enriched.merged.profile", sep="\t", header=T, check.names = F, row.names=1)
marker = read.table("../00.data/enriched.merged.group", sep="\t", header=T)
sample_map = read.table("../00.data/sample.info.20220103.tsv",sep="\t", header=T, check.names=F)

sample_map <- sample_map %>% filter(Project_2 == "Discover") %>% dplyr::select(Sample,Group,Project)

#dt = log10(dt+1e-6)
dt = dt+1e-6
dt1 = dt[marker$name, sample_map$Sample]

con_list = marker %>% filter(Group=="Control", Level == "species") %>% dplyr::select(name)
dis_list = marker %>% filter(Group=="Disease", Level == "species") %>% dplyr::select(name)

con_dt = data.frame(con = colSums(dt1[con_list$name,]))
dis_dt = data.frame(dis = colSums(dt1[dis_list$name,]))

data = data.frame(merge(con_dt, dis_dt, by='row.names'), row.names=1)
data$sub = data$dis-data$con

ps = unique(sample_map$Project)

result = rbind()
for(i in ps){
 x <- sample_map %>% 
   filter(Project == i) %>% 
   merge(data, by.x='Sample', by.y='row.names')
 x_summ <- x %>% filter(Group == "Control") %>% summarise(var=var(sub), sd=sd(sub), mean=mean(sub))
 
 # f0
 y_con <- x %>% filter(Group == "Control") %>% mutate(score=(sub-x_summ$mean)/x_summ$sd)
 y_dis <- x %>% filter(Group == "Disease") %>% mutate(score=(sub-x_summ$mean)/x_summ$sd)
 
 # f1
 # y_con <- x %>% filter(Group == "Control") %>% mutate(score=scale(sub))
 # y_dis <- x %>% filter(Group == "Disease") %>% mutate(score=scale(sub))
 
 # f2
 # y_con <- x %>% filter(Group == "Control") %>% mutate(score=(sub-x_summ$mean)/x_summ$var)
 # y_dis = x %>% filter(Group=="Disease") %>% mutate(score=(sub-x_summ$mean)/x_summ$var)
 result = rbind(result, y_con, y_dis)

}


#roc(result$Group~result$sub)

project_order = read.table("../00.data/project_order.list", sep="\t", header=T)
result$Project = factor(result$Project, levels=unique(project_order$project))
ggplot(result,aes(x=Project, y=con, fill=Group))+
  geom_boxplot()+
  theme_bw()+
  theme(axis.text.x = element_text(angle=90,hjust=1))+
  scale_fill_manual(values=structure(c("blue", "red"), names=c("Control","Disease")))

stack_result = rbind()
for(i in unique(result$Project)){
  temp_dt = result %>% filter(Project == i)
  dc = temp_dt %>% filter(Group == "Control")
  dd = temp_dt %>% filter(Group == "Disease")
  Pcon = wilcox.test(dc$con,dd$con)$p.value
  Ccon = mean(dc$con)
  Dcon = mean(dd$con)
  
  Pdis = wilcox.test(dc$dis,dd$dis)$p.value
  Cdis = mean(dc$dis)
  Ddis = mean(dd$dis)
  temp_result = data.frame(name=i,
                           "Mean(Control-con)"=Ccon, "Mean(Disease-con)"=Dcon,"Pval-con"=Pcon,
                           "Mean(Control-dis)"=Cdis, "Mean(Disease-dis)"=Ddis,"Pval-dis"=Pdis,
                           check.names=F
                           )
  stack_result = rbind(stack_result, temp_result)
  
}
#xx = result %>% group_by(Project, Group) %>% summarise(mean=mean(score), sd=sd(score), var=var(score))

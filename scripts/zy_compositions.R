library(ggplot2)
library(reshape2)
library(vegan)

zy_compositions <- function(dt=NA, top_N = 10, 
                               title="Composition",taxo.color = NA,width=0.9,
                               order_func = "order", order_n = 1, label_order=NA){
  # order_n -> 按照第几个丰度排序
  # order_func -> 按照什么对样本排序
  total_color = c("#ed1299","#09f9f5","#246b93","#cc8e12","#d561dd","#c93f00","#ddd53e","#4aef7b","#e86502","#9ed84e","#39ba30","#6ad157","#8249aa","#99db27","#e07233","#ff523f","#ce2523","#f7aa5d","#cebb10","#03827f","#931635","#373bbf","#a1ce4c","#ef3bb6","#d66551","#1a918f","#ff66fc","#2927c4","#7149af","#57e559","#8e3af4","#f9a270","#22547f","#db5e92","#edd05e","#6f25e8","#0dbc21","#280f7a","#6373ed","#5b910f","#7b34c1","#0cf29a","#d80fc1","#dd27ce","#07a301","#167275","#391c82","#2baeb5","#925bea","#63ff4f")
  if(is.na(taxo.color)){
    taxo.color = total_color[1:top_N]
  }
  
  dt = dt[rowSums(dt)!=0,] # 删除所有样本都没有的物种
  dt = dt[order(rowMeans(dt), decreasing=T),]
  
  dt$taxo_temp=rownames(dt)
  if(nrow(dt)>=top_N){dt$taxo_temp[top_N:nrow(dt)] = 'other'} # 如果物种数大于top_N
  data = data.frame(aggregate(. ~ taxo_temp, data=dt, sum), row.names=1, check.names=F)
  data = data[order(rowMeans(data), decreasing=T),]
  dl = reshape2::melt(as.matrix(data))
  # 样本排序方式
  if(order_func %in% c("order","cluster","specific")){
    if (order_func == "order"){
      # 按照含量第order_n个物种对样本进行排序
      label_order = data[order_n,] %>% t() %>% as.data.frame() %>% arrange_all() %>% rownames()
      #label_order = col_names[order(data[order_n,], decreasing=T )]
    }else if(order_func == "cluster"){
      # 按照样本聚类
      otu.dist=vegdist(t(data), method="bray")
      hc = hclust(otu.dist)
      label_order = hc$labels[hc$order]
    }else if(order_func == "specific"){
      label_order = label_order
    }
  }
  
  # 分类排序
  tax_ord = rownames(data)
  if(order_n != 1){
    tax_ord = levels(forcats::fct_relevel(tax_ord, tax_ord[order_n]))
  }
  dl$Var1 = factor(dl$Var1, level=rev(tax_ord))
  dl$Var2 = factor(dl$Var2, level=label_order)
  
  p <- ggplot(dl, aes(x=Var2, y=value, fill=Var1))+
    geom_bar(stat='identity',width=width)+
    theme_bw()+
    theme(panel.grid = element_blank(), axis.text.x = element_text(angle=45, hjust=1))+
    ggtitle(label=title)+
    scale_fill_manual(values=taxo.color)+
    scale_y_continuous(expand = c(0,0))
  
  p
}

zy_group_compositions <- function(dt=NA, sample_map=NA, ID=NA, group=NA, top_N = 10, 
                         title="Composition",taxo.color = NA,width=0.9,label_order=NA,
                         order_func = "order", order_n = 1){
    # order_n -> 按照第几个丰度排序
    # order_func -> 按照什么对样本排序
    total_color = c("#ed1299","#09f9f5","#246b93","#cc8e12","#d561dd","#c93f00","#ddd53e","#4aef7b","#e86502","#9ed84e","#39ba30","#6ad157","#8249aa","#99db27","#e07233","#ff523f","#ce2523","#f7aa5d","#cebb10","#03827f","#931635","#373bbf","#a1ce4c","#ef3bb6","#d66551","#1a918f","#ff66fc","#2927c4","#7149af","#57e559","#8e3af4","#f9a270","#22547f","#db5e92","#edd05e","#6f25e8","#0dbc21","#280f7a","#6373ed","#5b910f","#7b34c1","#0cf29a","#d80fc1","#dd27ce","#07a301","#167275","#391c82","#2baeb5","#925bea","#63ff4f")
    if(is.na(taxo.color)){
      taxo.color = total_color[1:top_N]
    }
    
    dt = dt[, sample_map[,ID] ]
    dt = dt[rowSums(dt)!=0,] # 删除所有样本都没有的物种
    dt = dt[order(rowMeans(dt), decreasing=T),]

    dt$taxo_temp=rownames(dt)
    if(nrow(dt)>=top_N){dt$taxo_temp[top_N:nrow(dt)] = 'other'} # 如果物种数大于top_N
    data = data.frame(aggregate(. ~ taxo_temp, data=dt, sum), row.names=1, check.names=F)
    data = data[order(rowMeans(data), decreasing=T),]
    dl = melt(as.matrix(data))
    
    # 样本排序方式
    if(order_func %in% c("order","cluster","specific")){
      if (order_func == "order"){
        # 按照含量第order_n个物种对样本进行排序
        label_order = data[order_n,] %>% t() %>% as.data.frame() %>% arrange_all() %>% rownames()
      }else if(order_func == "cluster"){
        # 按照样本聚类
        otu.dist=vegdist(t(data), method="bray")
        hc = hclust(otu.dist)
        label_order = hc$labels[hc$order]
      }else if(order_func == "specific"){
        label_order = label_order
      }
    }
    
    # 分类排序
    tax_ord = rownames(data)
    if(order_n != 1){
      tax_ord = levels(forcats::fct_relevel(tax_ord, tax_ord[order_n]))
    }
    
    dm = merge(dl, sample_map, by.x='Var2', by.y=ID)
    dm$Var1 = factor(dm$Var1, level=rev(tax_ord))
    dm$Var2 = factor(dm$Var2, level=label_order)
  
    p <- ggplot(dm, aes(x=Var2, y=value, fill=Var1))+
        geom_bar(stat='identity',width=width)+
        theme_bw()+
        theme(panel.grid = element_blank(),
              strip.placement = "outside",
              axis.text.x = element_text(angle=90, hjust=1))+
        facet_grid(.~get(`group`), scale='free',space = 'free_x'
                   ,switch = "both" # 标签在下
                   )+
        ggtitle(label=title)+
        scale_fill_manual(values=taxo.color)+
        scale_y_continuous(expand = c(0,0))
  
    p
}

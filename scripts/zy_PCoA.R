library(vegan)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(dplyr)

zy_pcoa <- function(dt=NA, sample_map=NA, group=NA, ID=NA, sample.color=NA,
                    ado_method="bray", pca_method="bray",
                    levels=0.95,star_plot=F,ellipse_plot=T,
                    title="PCoA", x="V1", y="V2"){
  
    # 对齐profile和分组的样本名称
    dt = dt[,sample_map[,ID]]
    
    ## colors 
    if (is.na(sample.color) || is.nan(sample.color)){
        sample.color = c(1:length(unique(sample_map[,group])))
    }
    # 统计每个分组各有多少,作为新的图例
    group_summ <- sample_map %>% dplyr::select(all_of(group)) %>% dplyr::group_by(across({{group}})) %>% summarise(count=n()) %>% mutate(new_label=paste(get(`group`), " (", count, ")", sep=""))
    new_label <- structure(group_summ$new_label,names=unlist(group_summ[,group]))
    
    message(paste(length(sample.color), "of groups to plot"))
    ## adonis
    ado = adonis(t(dt)~sample_map[,group], method = ado_method)
    ado_r2 = round(ado$aov.tab$R2[1], digits = 4)
    ado_p = ado$aov.tab$`Pr(>F)`[1]

    ## PCoA
    otu.dist = vegdist(t(dt), method = pca_method)
    pcoa = cmdscale(otu.dist, k=10, eig=T)
    eigs = signif(pcoa$eig/sum(pcoa$eig), 4)*100
    point = pcoa$points

    xlab = paste(x, " (",eigs[1],"%)", sep="")
    ylab = paste(y, " (",eigs[2],"%)", sep="")
    title = paste(title, "\nR2=",ado_r2,"\npvalue=", ado_p, sep="")
    
    dm = merge(point, sample_map, by.x='row.names', by.y=ID)
    p1 <- ggscatter(data=dm, x=x,y=y,
                    color=group,
                    star.plot = star_plot,
                    ellipse.level = levels,ellipse = ellipse_plot,
                    xlab=xlab,
                    ylab=ylab,
                    title=title
                    )+
    theme_bw()+
    geom_vline(xintercept=0)+
    geom_hline(yintercept=0)+
    theme(panel.grid = element_blank(),
          text = element_text(color="black"),
          axis.text = element_text(color="black"),
          axis.ticks = element_line(color="black", size=0.25),
          panel.border = element_rect(colour="black", size=0.25))+
    scale_fill_manual(values=sample.color, guide="none")+
    scale_color_manual(values=sample.color, labels=new_label)
    list(plot=p1, new_label=new_label)
}

zy_pcoa_with_arrow <- function(dt=NA, sample_map=NA, group=NA, ID=NA, sample.color=NA,
                    ado_method="bray", pca_method="bray",
                    levels=0.95,
                    title="PCoA", x="V1", y="V2"){
  
  # 对齐profile和分组的样本名称
  dt = dt[,sample_map[,ID]]
  
  ## colors 
  if (is.na(sample.color)){
    sample.color = c(1:length(unique(sample_map[,group])))
  }
  ## adonis
  ado = adonis(t(dt)~sample_map[,group], method = ado_method)
  ado_r2 = round(ado$aov.tab$R2[1], digits = 4)
  ado_p = ado$aov.tab$`Pr(>F)`[1]
  
  #-------------------------------------------------
  # pcoa
  otu.dist=vegdist(t(dt), method=pca_method)
  dt = dt[,colnames(as.matrix(otu.dist))]
  pca <- cmdscale(otu.dist, eig=TRUE, k = 10)
  #-----------
  # get species fit
  fit <- envfit(pca, t(dt), permutations = 3)
  fit_val <- scores(fit, display = c("vectors"))
  fit_val <- fit_val*vegan::ordiArrowMul(fit_val, fill = 2)
  fit_val <- fit_val[head(order(sqrt((fit_val[,1])^2+(fit_val[,2])^2),decreasing = T),n = 10),]
  #fit_val <- fit_val[dt %>% arrange(desc(rowMeans((.)))) %>% head(10) %>% rownames(),]
  
  eigs <-round(pca$eig/sum(pca$eig)*100,digits = 2)
  sample_axis <- as.data.frame(pca$points)
  sample_axis <- sample_axis[sample_map[,ID],]
  sample_axis$Group <- sample_map[,group]
  xlab = paste(x, " (",eigs[1],"%)", sep="")
  ylab = paste(y, " (",eigs[2],"%)", sep="")
  title = paste(title, "\nR2=",ado_r2,"\npvalue=", ado_p, sep="")
  
  vector_fd<-function(tab,sam){
    ms<-max(sam[,1]^2,sam[,2]^2)
    mt<-max(tab[,1]^2,tab[,2]^2)
    sqrt(ms)/sqrt(mt)
  }
  Change_axis <- vector_fd(fit_val,sample_axis)*0.5
  
  
  ggplot(sample_axis, aes(x = V1, y = V2, color = .data[[group]])) +
    stat_ellipse(aes(x = V1, y = V2, fill = .data[[group]]), geom = "polygon", alpha = 0.2, level = levels) +
    geom_point(aes(fill=.data[[group]]),size = 1, alpha = 0.3, shape=21) +
    scale_color_manual(values = sample.color) +
    scale_fill_manual(values = sample.color) +
    geom_segment(data=data.frame(fit_val), 
                 aes(x=0,y=0,xend=Dim1*Change_axis, yend=Dim2*Change_axis), 
                 arrow=arrow(length=unit(0.2,"cm")), 
                 color='black',alpha=1)  + 
    geom_text_repel(data=data.frame(fit_val), aes(Dim1*Change_axis, Dim2*Change_axis, label=rownames(fit_val)),
                    color='black',alpha=1)+
    geom_hline(yintercept = 0)+
    geom_vline(xintercept = 0)+
    theme_bw()+
    theme(panel.grid = element_blank(),
          text = element_text(color="black"),
          axis.text = element_text(color="black"),
          axis.ticks = element_line(color="black", size=0.25),
          panel.border = element_rect(colour="black", size=0.25))+
    labs(x = xlab,
         y = ylab,
         title = title) 
  
}


zy_dbrda <- function(dt=NA, sample_map=NA, group=NA, ID=NA, sample.color=NA,
                    ado_method="bray", pca_method="bray",
                    levels=0.95,ellipse_plot=F,star_plot=F,
                    title="dbRDA", x="CAP1", y="MDS1"){
  
  # 对齐profile和分组的样本名称
  dt = dt[,sample_map[,ID]]
  
  ## colors 
  if (is.na(sample.color) || is.nan(sample.color)){
    sample.color = c(1:length(unique(sample_map[,group])))
  }
  # 统计每个分组各有多少,作为新的图例
  group_summ <- sample_map %>% dplyr::select(all_of(group)) %>% dplyr::group_by(across({{group}})) %>% summarise(count=n()) %>% mutate(new_label=paste(get(`group`), " (", count, ")", sep=""))
  new_label <- structure(group_summ$new_label,names=unlist(group_summ[,group]))
  
  message(paste(length(sample.color), "of groups to plot"))
  ## adonis
  ado = adonis2(t(dt)~sample_map[,group], method = ado_method)
  ado_r2 = round(ado$R2[1], digits = 4)
  ado_p = ado$`Pr(>F)`[1]
  
  ## dbrda
  otu.dist = vegdist(t(dt), method = pca_method)
  db_rda = capscale(otu.dist ~ sample_map[,group])
  db_rda_score = scores(db_rda, choices = 1:10)
  point = db_rda_score$sites
  
  eigs = round(summary(eigenvals(db_rda))[2,]*100, digits=4)
  
  xlab = paste(x, " (", eigs[1], "% )", sep="")
  ylab = paste(y, " (", eigs[2], "% )", sep="")
  title = paste(title, "\nR2=",ado_r2,"\npvalue=", ado_p, sep="")
  
  dm = merge(point, sample_map, by.x='row.names', by.y=ID)
  p1 <- ggscatter(data=dm, x=x,y=y,
                  color=group,
                  ellipse.level = levels,ellipse = ellipse_plot,
                  star.plot = star_plot,
                  xlab=xlab,
                  ylab=ylab,
                  title=title
  )+
    theme_bw()+
    geom_vline(xintercept=0)+
    geom_hline(yintercept=0)+
    theme(panel.grid = element_blank(),
          text = element_text(color="black"),
          axis.text = element_text(color="black"),
          axis.ticks = element_line(color="black", size=0.25),
          panel.border = element_rect(colour="black", size=0.25))+
    scale_fill_manual(values=sample.color, guide="none")+
    scale_color_manual(values=sample.color, labels=new_label)
  list(plot=p1, new_label=new_label)
}

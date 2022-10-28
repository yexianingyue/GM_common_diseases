library(vegan)
library(ggpubr)

# 物种累计曲线
get_pan_data=function(dd){
    dd.curve=specaccum(t(dd), method = "random",permutations=1000)
    dd.curve.data=data.frame(Sites=dd.curve$sites, Richness=dd.curve$richness, SD=dd.curve$sd)
    #dd.curve.data$label=rep(nm, nrow(dd.curve.data))
    dd.curve.data
}


zy_nspecies = function(dt=NA, sample_map = NA
                       ,group="Group", ID="Sample"
                       ,sample.color=NA
                       ,title="Rarefaction curve analysis"
                       ){

    ## colors 
    if (is.na(sample.color) || is.nan(sample.color)){
        sample.color = c(1:length(unique(sample_map[,group])))
    }
    message(paste(length(sample.color), "of groups to plot"))

    dt = dt[,sample_map[,ID]]

    # 循环对每个分组执行函数
    result = rbind()
    for(g in unique(sample_map[,group])){
        temp_map = sample_map[which(sample_map[,group]==g),]
        temp_dt = dt[,temp_map[,ID]]
        temp_result = get_pan_data(temp_dt)
        temp_result[,group] = g
        result = rbind(result, temp_result)
    }


    p = ggplot(data=result, aes(x=Sites, y=Richness,color=.data[[group]]))+
        geom_line(size=1)+
        geom_errorbar(aes(ymax=Richness+SD, ymin=Richness-SD), width=.25)+
        theme_bw()+
        theme(panel.grid = element_blank())+
        scale_color_manual(values=sample.color)+
        xlab("Number of samples")+
        ylab("Number of species")+
        ggtitle(title)
    p
}


sigFunc = function(x){
    if(x < 0.001){"***"} 
    else if(x < 0.01){"**"}
    else if(x < 0.05){"*"}
    else{NA}}

zy_alpha = function(dt=NA, sample_map=NA, group="Group", ID="Sample", # 必须参数
                    index="shannon", # 计算参数
                    sample.color=NA, # 美化参数
                    box_width=0.5, # 箱式图宽度
                    title="alpha diversity", # 文字参数,
                    violin = F
                    ){
    # pvalue给的是非精确计算exact=F
    ## colors 
    if (is.na(sample.color) || is.nan(sample.color)){
        sample.color = c(1:length(unique(sample_map[,group])))
    }
    message(paste(length(sample.color), "of groups to plot"))

    ## align dt and group
    dt = dt[,sample_map[,ID]]
    dt = dt[rowSums(dt)!=0,]

    #alpha
    alpha = data.frame(alpha = diversity(t(dt),index=index))
    dm = merge(alpha,sample_map, by.x='row.names', by.y=ID)
    comp = combn(as.character(unique(dm[,group])),2,list)
    
    p = ggplot(dm, aes(x=.data[[group]], y=alpha,fill=.data[[group]]))
    if(isTRUE(violin)){
        p <- p+
            geom_violin()+
            geom_boxplot(width=box_width, fill="white", position = position_dodge2(preserve = 'single')
                         ,outlier.shape = 21,outlier.fill=NA, outlier.color="#c1c1c1")
    }else{
        p <- p+ 
            geom_boxplot(position = position_dodge2(preserve = 'single')
                     ,outlier.shape = 21,outlier.fill=NA, outlier.color="#c1c1c1")
    }
    
    p <- p+
        theme_bw()+
        theme(panel.grid = element_blank())+
        scale_fill_manual(values=sample.color)+
        #geom_signif(comparisons =comp,test='wilcox.test',test.args=list(exact=F),step_increase = 0.1,map_signif_level=sigFunc)+
        geom_signif(comparisons =comp,test='wilcox.test',test.args=list(exact=F),step_increase = 0.1)+
        ggtitle(title)

    p
}

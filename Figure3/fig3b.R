library(ComplexHeatmap)
library(dplyr)
library(tidyselect)
library(reshape2)
library(circlize)

rm(list=ls())

# define function for heatmap
sigFun <- function(x){
  if(is.na(x)){NA}
  if(x<0.01){"+"
  }else if(x<0.05){
    "*"
  }
}

load("colors.map")
project_map = read.table("project.group",sep="\t", header=T, row.names=1, check.names=F)
marker = read.table("enriched.species.group",sep="\t", check.names=F, header=T)
dt = read.table("wilcox.tsv",sep="\t", check.names=F)
dt$log2_fd = ifelse(dt$enriched == "Disease", 0-log2(dt$fold_change), log2(dt$fold_change))
pq <- data.frame(dcast(dt, name~Project, value.var="qval"), row.names=1, check.names=F)
pv <- data.frame(dcast(dt, name~Project, value.var="pvalue"), row.names=1, check.names=F)
fd <- data.frame(dcast(dt, name~Project,value.var='log2_fd'), row.names=1, check.names=F)

#filter
marker = marker %>% filter(known != "unknown")
rownames(marker) = marker$genomes
dt <- dt %>% 
  filter(name %in% marker$genomes)

fd = fd[marker$genomes,]
pq = pq[marker$genomes,colnames(fd)]



# Legend
col_legend = colorRamp2(c(-1,-0.5,0, 0.5,1),
                        c("#8E0152","#D5589D", "#ffffff","#6EAE36","#276419")
)
lgd = Legend(col_fun = col_legend, title = "foo")


# 聚类
rod = hclust(dist(fd,method = "canb"),method = "complete")$order
cod = hclust(dist(t(fd),method = "maximum"),method = "complete")$order
fd = fd[rod, cod]
pq = pq[rod, cod]
marker = marker[rownames(fd),]

# significant barplot
barplot_dt = matrix(0,ncol=2, nrow=nrow(fd),
                    dimnames = list(
                      rownames(fd),
                      c("Control","Disease")
                    ))

for( i in rownames(pq)){
  cc = sum((pq[i,fd[i,]>0]<0.05)+0) # count -> Control
  cd = sum((pq[i,fd[i,]<0]<0.05)+0) # count -> Disease
  barplot_dt[i,"Control"] = cc
  barplot_dt[i,"Disease"] = cd
}

# barplot annotation
row_ha = rowAnnotation(
  bar=anno_barplot(barplot_dt)
)



# 分割
rsp = data.frame(enriched=marker$Group)
rsp = factor(rsp$enriched, levels=c("Control","Disease"))

# heatmap
cm = colnames(fd)
rownames(fd) = marker$name
Heatmap(as.matrix(fd),
        row_split = rsp, border=T,
        cluster_rows = F, cluster_columns = F,
        show_row_names = T, show_column_names = T,
        row_names_gp = gpar(fongsize=1),
        column_names_gp = gpar(col = colors[project_map[cm,'group']]),
        use_raster=F,
        col = col_legend,
        show_heatmap_legend = F,
        right_annotation = row_ha,
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%s", sigFun(pq[i, j])), x, y=ifelse(pq[i, j]<0.01,y,y+unit(-0.01,"npc")), 
                    gp = gpar(fontsize = ifelse(pq[i, j]<0.01,15,20))
          )
        }
        )
# draw legend
draw(lgd, x=unit(0.9,"npc"), y=unit(0.1, "npc"))
# 再加一层菌株分类

source("../scipts/zy_plot_ROC.R")
in_g = "sample.info.20220103.tsv"
low_list = c("LiR_2021.CA", "ZhongH_2019.pT2D", "ZhongH_2019.tnT2D","YeZ_2018.BD",
             "LiuP_2021.MG", "YeZ_2020.VKH", "LuW_2018.HIV", "WangM_2019.ASD",
             "WanY_2021.ASD", "QianY_2020.PD", "MaoL_2021.PD", "ZhuF_2020.SCZ")
sample_map = read.table(in_g, sep="\t", header=T, check.names=F)
sample_map$Sample = make.names(sample_map$Sample)

# fig 5a
dt = read.table("all_predict.tsv", sep="\t")
temp_map <- sample_map %>%
    filter(Project_2 == "Discover") %>%
    mutate(group=ifelse(Project %in% low_list, "low", "high")) %>%
    select(Sample, Group, group) %>%
    unique()
temp_map1 <- sample_map %>%
    filter(Project_2 == "Discover") %>%
    mutate(group="all") %>%
    select(Sample, Group, group) %>%
    unique() %>%
    rbind(temp_map)
dt1 <- dt %>%
    filter(V1==260, V2 %in% c(3751, 9516, 1004))
    data = merge(temp_map1, dt1, by.x="Sample", by.y="V3")
    data1 = data %>%
    filter((group=="all" & V2==3751)|(group=="high" & V2==9516)|(group=="low" & V2==1004))

x = plot_roc(data1 , pred="V4", true="Group", group="group") # plot
y = calc_auc(data1 , pred="V4", true="Group", group="group", acc=T) # auc table

# fig 5b
tax_color = c("#c51b7d","#762a83","#7f0000","#542788","#b35806","#b2182b","#2166ac","#b8e186","#fdbb84","#fee8c8","#fdd49e","#f1b6da","#9970ab","#4393c3","#ef6548","#fc8d59","#8073ac","#d6604d","#f4a582","#d7301f","#b30000","#4d9221","#fdb863","white")
names(tax_color) = c("f__Enterococcaceae","f__Erysipelatoclostridiaceae","f__Lachnospiraceae","f__Lactobacillaceae","f__Mycoplasmoidaceae","f__Oscillospiraceae","f__Ruminococcaceae","f__Streptococcaceae","g__Acetatifactor","g__Agathobacter","g__Coprococcus","g__Enterococcus","g__Erysipelatoclostridium","g__Faecalibacterium","g__Fusicatenibacter","g__Lachnospira","g__Lactobacillus","g__Lawsonibacter","g__Oscillibacter","g__Roseburia","g__Sellimonas","g__Streptococcus","g__Ureaplasma","other")
dt = read.table("./impvar_species_scaled.tsv", sep="\t", header=T)
dt$name = factor(dt$name, levels=rev(dt$name))
dt$Enriched = as.factor(dt$Enriched)
table(dt$taxo)
p <- ggplot(dt,aes(y=name))+
    theme_bw()+
    geom_bar(aes(x=MeanDecreaseAccuracy), stat="identity", fill="#75aadb", color="black")+
    theme(axis.text.y=element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank())
# fig 5c
dt = read.table("boxplot.data",sep="\t", check.names=F)
sample_map = read.table("sample.info.20220103.tsv", sep="\t", header=T, check.names=F)

sample_map <- sample_map %>%
    filter(Project_2 == "Discover") %>%
    dplyr::select(Sample, Group, Project) %>%
    unique()

sample_map$Sample = make.names(sample_map$Sample)

data = merge(dt, sample_map, by.x='row.names', by.y='Sample')

# density

density_data <- data %>%
    dplyr::select(Row.names, Group, Control)%>%
    unique()

density_summ <- density_data %>%
    group_by(Group) %>%
    summarise(grp.mean=mean(Control), grp.median=median(Control))

ggplot(density_data, aes(x=Control, color=Group, fill=Group))+
    geom_density(alpha=0.6)+
    geom_vline(data=density_summ, aes(xintercept=grp.mean, color=Group))+
    theme_bw()


# fig 5d
load("marker_scaled.modle.RData")
load("marker.modle.RData")
test_dt = read.table("validation.marker_OTU.relative.profile", sep="\t", header=T,row.names = 1,check.names=F)
test_map = read.table("validation.group.20220721", sep="\t", header=T)

result = rbind()
for(g in unique(test_map$group1)){
    temp_test_map <- test_map %>% 
    filter(group1 == g) %>%
    select(Sample, Group) %>%
    unique()
    temp_test_dt = test_dt[,temp_test_map$Sample]
    y <- zy_format_class_name(temp_test_dt, temp_test_map, "Sample")
    pred <- as.data.frame(predict(fit,t(y$rf_dt),type = "prob"))
    temp_result <- merge(y$rf_map[,c("Sample","Group")], pred, by.x="Sample", by.y="row.names")
    temp_result$group1 = g
    result = rbind(result, temp_result)
}
x = c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6")
plot_roc(result, pred="Disease", true="Group", group="group1", cols = x)

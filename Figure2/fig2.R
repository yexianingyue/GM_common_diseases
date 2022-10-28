load("ggplot_observed_foldchange.RData") #  fig 2a
load("ggplot_shannon_foldchange.RData") #  fig 2b
load("ggplot_adonis.RData") # fig 2c
load("ggplot_auc.RData") # fig 2d

remove_y <- theme(
  axis.text.y = element_blank(),
  axis.ticks.y = element_blank(),
  axis.title.y = element_blank()
)

fig1_labels <- adonis_plot+theme(aspect.ratio = 1000)

figs = list(
  fig1_labels,
  observed_foldchange+remove_y,
  shannon_foldchange+remove_y,
  adonis_plot+remove_y,
  auc_plot+remove_y
)

x = ggarrange(plotlist = figs,
              nrow = 1,
              legend='bottom',
              common.legend = TRUE)
x

ggsave(x, "fig2.pdf", width=20, hight=8)

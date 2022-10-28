library(ggplot2)
library(ggpubr)
load("ggplot_study_count_plot.RData") # fig 1a
load("ggplot_composition_phylum.RData") # fig 1b
load("ggplot_composition_genus.RData") # fig 1c
load("ggplot_composition_species.RData") # fig 1d

remove_y <- theme(
                  axis.text.y = element_blank(),
                  axis.ticks.y = element_blank(),
                  axis.title.y = element_blank()
                  )

fig1_labels <- study_count_plot+theme(aspect.ratio = 1000)

figs = list(
            fig1_labels,
            study_count_plot+remove_y,
            composition_phylum+remove_y,
            composition_genus+remove_y,
            composition_species+remove_y
            )

x = ggarrange(plotlist = figs,
              nrow = 1,
              legend='bottom',
              common.legend = TRUE)
x

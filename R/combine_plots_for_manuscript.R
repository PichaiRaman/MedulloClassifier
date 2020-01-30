# script to combine plots for Figure 2 and Supp Figure 1
library(ggpubr)

load('results/Fig_S1A_S1B.RData')
load('results/Fig_S2A.RData')
load('results/Fig_S2B.RData')
load('results/Fig_2B.RData')
load('results/Fig_2C.RData')

ggexport(ggarrange(s1a, s1b, labels = c("A", "B"), 
                   align = "h"), width = 10, height = 5, 
         filename = "results/plots/SuppFigure1.pdf")
ggexport(ggarrange(s2a, s2b, labels = c("A", "B"), 
                   align = "h"), width = 10, height = 5,
         filename = "results/plots/SuppFigure2.pdf")

ggexport(ggarrange(ggarrange(fig2b[[4]], NULL, nrow = 2, heights = c(8, 2)), 
                   ggarrange(fig2c, NULL, nrow = 2, heights = c(4, 6)), 
                   nrow = 1, 
                   widths = c(10, 5), labels = c("A", "B"), common.legend = T),
         filename = "results/plots/Figure2.pdf", width = 10)
# ggexport(ggarrange(ggarrange(fig1b[[4]], heights = 10), 
#                    ggarrange(fig1c, NULL, nrow = 2, heights = c(4, 6)), 
#                    nrow = 1, 
#                    widths = c(10, 5), labels = c("A", "B"), common.legend = T),
#          filename = "results/plots/Figure2.pdf", width = 10)

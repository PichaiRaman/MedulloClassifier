# script to combine plots for Figure 2 and Supp Figure 1
library(ggpubr)

load('results/Fig_S1A_S1B_2B_S1C.RData')
load('results/Fig_S1D.RData')
load('results/Fig_2C.RData')

ggexport(ggarrange(s1a, s1b, s1c, s1d, labels = c("A", "B", "C", "D"), 
                   align = "h"), 
         filename = "results/plots/SuppFigure1.pdf")
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

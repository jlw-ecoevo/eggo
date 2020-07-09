

## JLW 2020 - Comparison of two gut MAG datasets

# Load Packages ----------------------------------------------------------------


library(dplyr)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)

# Load Data --------------------------------------------------------------------

setwd("~/eggo/Data/eggo-data/")
load("CodonStatistics_Pasolli.RData")
mag_df1 <- mag_df %>% subset(nHE>=10) %>% 
  subset(!is.na(d)) %>% subset(Body.Site=="Stool")
load("CodonStatistics_Nayfach.RData")
mag_df2 <- mag_df %>% subset(nHE>=10)

# Plot -------------------------------------------------------------------------

setwd("~/eggo/Figs")
pdf("Pasolli_vs_Nayfach.pdf",width=7,height=5)
ggplot(mag_df1,aes(x=d,fill="Pasolli et al.")) + geom_density(alpha=0.5) + 
  geom_density(data=mag_df2,aes(fill="Nayfach et al."),alpha=0.5) + 
  scale_x_log10(limit=c(0.05,100)) + theme_bw() + labs(fill="MAGs Set") + 
  xlab("Predicted Minimal Doubling Time (Hours)")
dev.off()

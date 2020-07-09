

library(dplyr)
library(data.table)
library(ggplot2)
library(ggpubr)


setwd("~/eggo/Data")
load("CodonStatistics_Okie.RData")

water_df <- metagenome_df %>% subset(Final=="X" & Source=="W")
wilcox.test(d.adj ~ Fertilized, data=water_df) 
wilcox.test(d ~ Fertilized, data=water_df) 

p1 <- ggplot(water_df,aes(x=Fertilized,group=Fertilized,y=d)) + 
  scale_y_log10(limits=c(0.5,2.5)) + 
  geom_boxplot() + geom_point(color="red",size=3) +theme_bw() +
  ylab("Median Minimal Doubling Time (d)") + 
  ggtitle("Unweighted") 

p2 <- ggplot(water_df,aes(x=Fertilized,group=Fertilized,y=d.adj)) + 
  scale_y_log10(limits=c(0.5,2.5)) + 
  geom_boxplot() + geom_point(color="blue",size=3) +theme_bw() +
  ylab("Median Minimal Doubling Time (d)") + 
  ggtitle("Weighted by Relative Abundance")

p3 <- ggplot(water_df,aes(x=d,y=d.adj)) + 
  geom_point() + 
  geom_abline(intercept = 0, slope = 1, lty = 2) + 
  theme_bw() + 
  xlab("Unweighted Minimal Doubling Time") + 
  ylab("Abundance-Weighted Minimal Doubling Time") + 
  scale_x_log10() + 
  scale_y_log10()

setwd("~/eggo/Figs")
pdf("Okie_water.pdf",width=7,height=8)
ggarrange(ggarrange(p1,p2,ncol=2,labels=c("(a)","(b)")),
          p3,
          nrow = 2,
          labels=c("","(c)"),
          hjust=0,
          vjust=1)
dev.off()


# Okie et al. Samples Analyzed in Paper (excluding low read samples) -----------

paper_df <- metagenome_df %>% subset(Paper=="X")

wilcox.test(d.adj ~ Fertilized, data=paper_df) 

p1 <- ggplot(paper_df,aes(x=Fertilized,group=Fertilized,y=d)) + 
  scale_y_log10(limits=c(0.5,2.5)) + 
  geom_boxplot() + geom_point(color="red",size=3) +theme_bw() +
  ylab("Median Minimal Doubling Time (d)") + 
  ggtitle("Unweighted") 

p2 <- ggplot(paper_df,aes(x=Fertilized,group=Fertilized,y=d.adj)) + 
  scale_y_log10(limits=c(0.5,2.5)) + 
  geom_boxplot() + geom_point(color="blue",size=3) +theme_bw() +
  ylab("Median Minimal Doubling Time (d)") + 
  ggtitle("Weighted by Relative Abundance")

p3 <- ggplot(paper_df,aes(x=d,y=d.adj)) + 
  geom_point() + 
  geom_abline(intercept = 0, slope = 1, lty = 2) + 
  theme_bw() + 
  xlab("Unweighted Minimal Doubling Time") + 
  ylab("Abundance-Weighted Minimal Doubling Time") + 
  scale_x_log10() + 
  scale_y_log10()

setwd("~/eggo/Figs")
pdf("Okie_paper.pdf",width=7,height=8)
ggarrange(ggarrange(p1,p2,ncol=2,labels=c("(a)","(b)")),
          p3,
          nrow = 2,
          labels=c("","(c)"),
          hjust=0,
          vjust=1)
dev.off()


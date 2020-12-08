

library(dplyr)
library(data.table)
library(ggplot2)
library(ggpubr)


setwd("~/eggo/Data")
load("CodonStatistics_RedSea.RData")

metagenome_df$Nutrient <- gsub("addition","",metagenome_df$Nutrient)
metagenome_df$Nutrient[grep("-C",metagenome_df$Treatment)] <- 
  paste(metagenome_df$Nutrient[grep("-C",metagenome_df$Treatment)],"(continuous addition)")
metagenome_df$Nutrient[grep("-I",metagenome_df$Treatment)] <- 
  paste(metagenome_df$Nutrient[grep("-I",metagenome_df$Treatment)],"(initial addition)")

mean_c <- mean(metagenome_df$d.adj[metagenome_df$Nutrient=="control mesocosm"])
p1 <- ggplot(metagenome_df %>% 
               group_by(Time,Nutrient) %>% 
               summarise(d.adj=mean(d.adj)), 
             aes(x=Time, y=d.adj)) +
  geom_line(alpha=0.5,lwd=2) +  scale_y_log10() +scale_x_date(date_labels = "%Y %b %d") + 
  facet_wrap(vars(Nutrient),nrow=1) + geom_hline(yintercept=mean_c,lty=2,color="red") +
  geom_point(data = metagenome_df) +
  theme_bw() + ylab("Median Minimal Doubling Time (d)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
  ggtitle("Weighted") + 
  geom_vline(xintercept = as.Date("2013-02-03"), lty=2, color="blue")


mean_c <- mean(metagenome_df$d[metagenome_df$Nutrient=="control mesocosm"])
p2 <- ggplot(metagenome_df %>% 
               group_by(Time,Nutrient) %>% 
               summarise(d=mean(d)), 
             aes(x=Time, y=d)) +
  geom_line(alpha=0.5,lwd=2) +  scale_y_log10() +scale_x_date(date_labels = "%Y %b %d") + 
  facet_wrap(vars(Nutrient),nrow=1) + geom_hline(yintercept=mean_c,lty=2,color="red") +
  geom_point(data = metagenome_df) +
  theme_bw() + ylab("Median Minimal Doubling Time (d)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
  ggtitle("Unweighted") + 
  geom_vline(xintercept = as.Date("2013-02-03"), lty=2, color="blue")

ggarrange(p1,p2,nrow=2)


setwd("~/eggo/Figs")
pdf("RedSea.pdf",width=8,height=8)
ggarrange(p1,p2,
          labels=c("(a)","(b)"),
          nrow=2,
          hjust=0,
          vjust=1)
dev.off()

p1 <- ggplot(metagenome_df,aes(x=d,y=d.adj)) + 
  geom_point() + 
  geom_abline(intercept = 0, slope = 1, lty = 2) + 
  theme_bw() + 
  xlab("Raw Predicted Minimal Doubling Time") + 
  ylab("Abundance-Weighted Predicted Minimal Doubling Time") + 
  scale_x_log10() + 
  scale_y_log10()
p2 <- ggplot(metagenome_df,aes(x=1,y=d.adj-d)) + 
  geom_violin(fill="gray") + geom_boxplot(width=0.1) +
  theme_bw() + geom_hline(yintercept = 0,lty=2) + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  xlab("") + ylab(expression("Difference in Estimate (" * d["Weighted"] - d["Unweighted"] * ")"))
setwd("~/eggo/Figs")
pdf("RedSea_diff.pdf",width=7,height=5)
ggarrange(p1,p2,ncol=2,widths = c(3,1))
dev.off()

t.test(metagenome_df$d, metagenome_df$d.adj, paired = TRUE)

########### Poster

mean_c <- mean(metagenome_df$d.adj[metagenome_df$Nutrient=="control mesocosm"])
init_metagenomes <- metagenome_df %>% 
  subset(Nutrient %in% c("control mesocosm",
                         "N,P  (initial addition)",
                         "N,P,Si  (initial addition)") & nHE>100)
init_metagenomes$Nutrient <- init_metagenomes$Nutrient %>%
  gsub(pattern = "[(]initial addition[)]", replace = "Addition") %>%
  gsub(pattern = ",", replace = "+")
p3 <- ggplot(init_metagenomes %>% 
               group_by(Time,Nutrient) %>% 
               summarise(d.adj=mean(d.adj)), 
             aes(x=Time, y=d.adj)) +
  geom_line(alpha=0.5,lwd=2) +  scale_y_log10() +scale_x_date(date_labels = "%Y %b %d") + 
  facet_wrap(vars(Nutrient),nrow=1) + geom_hline(yintercept=mean_c,lty=2,color="red") +
  #geom_point(data = init_metagenomes) +
  theme_pubclean() + ylab("Median Minimal Doubling Time (d)") +
  theme(axis.text.x = element_text(angle = 25, hjust = 1, vjust = 1),
        text = element_text(size=18)) + 
  geom_vline(xintercept = as.Date("2013-02-03"), lty=2, color="blue") + 
  ggtitle("Nutrient Enrichment Experiment in Marine Mesocosms [2]")


p1 <- ggplot(metagenome_df,aes(x=d,y=d.adj)) + 
  geom_point() + 
  geom_abline(intercept = 0, slope = 1, lty = 2) + 
  theme_pubclean() + 
  xlab("Unweighted Predictions (Hours)") + 
  ylab("Abundance-Weighted Predictions") + 
  scale_x_log10() + 
  scale_y_log10() +
  theme(text = element_text(size=18))
p2 <- ggplot(metagenome_df,aes(x=1,y=d.adj-d)) + 
  geom_violin(fill="gray") + geom_boxplot(width=0.1) +
  theme_pubclean() + geom_hline(yintercept = 0,lty=2) + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        text = element_text(size=18)) +
  xlab("") + ylab(expression("Difference in Estimate (" * d["Weighted"] - d["Unweighted"] * ")"))

setwd("~/eggo/Figs")
pdf("RedSea_poster.pdf",width=15,height=5)
ggarrange(p3,p1,ncol=2,widths = c(5,3))
dev.off()
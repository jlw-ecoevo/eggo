
library(ggplot2)
library(dplyr)
library(data.table)
library(ggpubr)

merge.easy <- function(df1,df2,key){
  df1 <- data.table(df1,key=key)
  df2 <- data.table(df2,key=key)
  return(as.data.frame(unique(merge(df1,df2,all.x=TRUE,by=.EACHI,allow.cartesian=TRUE)),stringsAsFactors=F))
}

x <- read.csv("~/eggo/Data/condensed_species_NCBI.csv")
x$Vol_lo <- pi*((x$d1_lo/2)^2)*x$d2_lo

y <- x %>% subset(!is.na(Vol_lo)) %>% subset(select=c(species,genus,Vol_lo,d1_lo,d2_lo))
names(y) <- c("Species","Genus","Vol","Dim1","Dim2")

setwd("~/eggo/Data/eggo-data")
load("CodonStatistics_RefSeq.RData")

y_spp <- y %>% group_by(Species) %>% summarise(Vol=mean(Vol),Dim1=mean(Dim1),Dim2=mean(Dim2))
spp_d <- assembly_df %>% group_by(Species) %>% summarize(d=mean(d))
spp_d <- merge.easy(spp_d,y_spp,key="Species")

ggplot(spp_d,aes(x=Vol,y=d)) + geom_point() + 
  scale_y_log10(limits=c(0.01,100)) + scale_x_log10() + 
  geom_smooth(method="lm") + theme_pubclean() + 
  geom_hline(yintercept=5,color="red",lty=2)

lm_dv <- lm(log10(d)~log10(Vol),data=spp_d)
summary(lm_dv)

pv <- ggplot(spp_d,aes(x=Vol,y=d)) + geom_point(alpha=0.25) +
  scale_y_log10(limits=c(0.01,100)) + scale_x_log10() +
  geom_smooth() + theme_pubclean()+ 
  geom_hline(yintercept=5,color="red",lty=2) +
  xlab("Volume (assuming cylinder)") + 
  ylab("Max. Doubling Time (Hours)")

pd1 <- ggplot(spp_d,aes(x=Dim1,y=d)) + geom_point(alpha=0.25) +
  scale_y_log10(limits=c(0.01,100)) + scale_x_log10() +
  geom_smooth() + theme_pubclean()+ 
  geom_hline(yintercept=5,color="red",lty=2)+
  xlab("Shorter Dimension (um)") + 
  ylab("Max. Doubling Time (Hours)")

pd2 <- ggplot(spp_d,aes(x=Dim2,y=d)) + geom_point(alpha=0.25) +
  scale_y_log10(limits=c(0.01,100)) + scale_x_log10() +
  geom_smooth() + theme_pubclean()+ 
  geom_hline(yintercept=5,color="red",lty=2)+
  xlab("Longer Dimension (um)") + 
  ylab("Max. Doubling Time (Hours)")

pd1d2 <- ggplot(spp_d,aes(x=Dim1,y=Dim2)) + geom_point(alpha=0.25) +
  scale_y_log10(limits=c(0.01,100)) + scale_x_log10() +
  geom_smooth() + theme_pubclean()+ 
  geom_hline(yintercept=5,color="red",lty=2)+
  xlab("Shorter Dimension (um)") + 
  ylab("Longer Dimension (um)")

ggarrange(pd1,pd2,pv,pd1d2,ncol=2,nrow=2)

# Not significant
fisher.test(table(spp_d$Dim1>quantile(spp_d$Dim1,0.99,na.rm=T),
                    spp_d$d<quantile(spp_d$d,0.01,na.rm=T)))

fisher.test(table(spp_d$Dim2>quantile(spp_d$Dim2,0.99,na.rm=T),
                  spp_d$d<quantile(spp_d$d,0.01,na.rm=T)))

fisher.test(table(spp_d$Vol>quantile(spp_d$Vol,0.99,na.rm=T),
                  spp_d$d<quantile(spp_d$d,0.01,na.rm=T)))

y_genus <- y %>% group_by(Genus) %>% summarise(Vol=mean(Vol),Dim1=mean(Dim1),Dim2=mean(Dim2))
genus_d <- assembly_df %>% group_by(Genus) %>% summarize(d=mean(d))
genus_d <- merge.easy(genus_d,y_genus,key="Genus")

pv <- ggplot(genus_d,aes(x=Vol,y=d)) + geom_point(alpha=0.25) +
  scale_y_log10(limits=c(0.01,100)) + scale_x_log10() +
  geom_smooth() + theme_pubclean()+ 
  geom_hline(yintercept=5,color="red",lty=2) +
  xlab("Volume (assuming cylinder)") + 
  ylab("Max. Doubling Time (Hours)")

pd1 <- ggplot(genus_d,aes(x=Dim1,y=d)) + geom_point(alpha=0.25) +
  scale_y_log10(limits=c(0.01,100)) + scale_x_log10() +
  geom_smooth() + theme_pubclean()+ 
  geom_hline(yintercept=5,color="red",lty=2)+
  xlab("Shorter Dimension (um)") + 
  ylab("Max. Doubling Time (Hours)")

pd2 <- ggplot(genus_d,aes(x=Dim2,y=d)) + geom_point(alpha=0.25) +
  scale_y_log10(limits=c(0.01,100)) + scale_x_log10() +
  geom_smooth() + theme_pubclean()+ 
  geom_hline(yintercept=5,color="red",lty=2)+
  xlab("Longer Dimension (um)") + 
  ylab("Max. Doubling Time (Hours)")

pd1d2 <- ggplot(genus_d,aes(x=Dim1,y=Dim2)) + geom_point(alpha=0.25) +
  scale_y_log10(limits=c(0.01,100)) + scale_x_log10() +
  geom_smooth() + theme_pubclean()+ 
  geom_hline(yintercept=5,color="red",lty=2)+
  xlab("Shorter Dimension (um)") + 
  ylab("Longer Dimension (um)")

setwd("~/eggo/Figs")
pdf("Madin_Size.pdf",width=8,height=7)
ggarrange(pd1,pd2,pv,pd1d2,ncol=2,nrow=2)
dev.off()

# Not significant
fisher.test(table(genus_d$Dim1>quantile(genus_d$Dim1,0.99,na.rm=T),
                  genus_d$d<quantile(genus_d$d,0.01,na.rm=T)))

fisher.test(table(genus_d$Dim2>quantile(genus_d$Dim2,0.99,na.rm=T),
                  genus_d$d<quantile(genus_d$d,0.01,na.rm=T)))

fisher.test(table(genus_d$Vol>quantile(genus_d$Vol,0.99,na.rm=T),
                  genus_d$d<quantile(genus_d$d,0.01,na.rm=T)))

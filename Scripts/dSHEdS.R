
# dSHe/dS summary

# Load  packages ---------------------------------------------------------------

library(dplyr)
library(reshape2)
library(ggplot2)
library(data.table)
library(ggpubr)
library(ggpointdensity)

merge.easy <- function(df1,df2,key){
  df1 <- data.table(df1,key=key)
  df2 <- data.table(df2,key=key)
  return(as.data.frame(unique(merge(df1,df2,all.x=TRUE,by=.EACHI,allow.cartesian=TRUE)),stringsAsFactors=F))
}

# Load data --------------------------------------------------------------------

dnds <- read.delim("~/eggo/Data/summary_dNdS_grouped.tbl",
                   head=F,stringsAsFactors = F,sep=" ")
names(dnds) <- c("Pair","dNdS")
r <- read.delim("~/eggo/Data/summary_dS_ribo_grouped.tbl",
                head=F,stringsAsFactors = F,sep=" ")
names(r) <- c("Pair","dSHE")
nr <- read.delim("~/eggo/Data/summary_dS_noribo_grouped.tbl",
                 head=F,stringsAsFactors = F,sep=" ")
names(nr) <- c("Pair","dS")
atgc <- read.delim("~/eggo/Data/atgcgenomelist.csv",sep="\t",
                   head=F,stringsAsFactors = F)
load("~/eggo/Data/EGGO.RData")
EGGO$Assembly <- substr(EGGO$Assembly,1,15) 

load("~/eggo/Data/gRodon-files/GrowthRates.rda")
load("~/eggo/Data/gRodon-files/Accession2Species.rda")
names(spp_acc) <- c("Assembly","Species")
acc_d <- merge.easy(spp_acc,d,key="Species")
d <- EGGO %>% subset(select=c(Assembly,d))

# summarise --------------------------------------------------------------------

# Limit to pairs w/ 0.01 < dS < 1
pairs_rnr <- merge.easy(r,nr,key="Pair")
x <- pairs_rnr %>% subset(dS > 0.01 & dS < 1)

atgc_genomes <- data.frame(ATGC=character(),
                           Genome=character(),
                           stringsAsFactors = F)
for(i in 1:nrow(atgc)){
  genomes <- strsplit(atgc$V2[i],split=",") %>% 
    unlist() %>% 
    gsub(pattern=".*[.]GCF_",replace="GCF_")
  if(length(genomes)==0){genomes <- atgc$V2[i]}
  print(i)
  print(genomes)
  atgc_genomes <- rbind(atgc_genomes,
                        data.frame(ATGC=atgc$V1[i],
                                   Assembly=genomes,
                                   stringsAsFactors = F))
}
atgc_genomes <- merge.easy(atgc_genomes,d,key="Assembly") %>% subset(!is.na(d))
atgc <- atgc_genomes %>% group_by(ATGC) %>% summarise(d=mean(d)) %>% as.data.frame(stringsAsFactors=F)
atgc$Copiotroph <- atgc$d < 5
atgc_genomes$d <- NULL
x$ATGC <- x$Pair %>% gsub(pattern="[.].*",replace="")
x <- merge.easy(x,atgc,key="ATGC") %>% subset(!is.na(Copiotroph))
y <- x
y$d <- NULL
y$Copiotroph <- NULL
atgc_genomes2 <- merge.easy(atgc_genomes,acc_d,key="Assembly")
atgc <- atgc_genomes2 %>% group_by(ATGC) %>% summarise(d=mean(d,na.rm=T)) %>% 
  as.data.frame(stringsAsFactors=F) %>% subset(!is.na(d))
atgc$Copiotroph <- atgc$d < 5
y <- merge.easy(y,atgc,key="ATGC") %>% subset(!is.na(d))

# Plot -------------------------------------------------------------------------

p1 <- ggplot(y,aes(x=dSHE/dS)) + 
  geom_density(aes(fill=factor(Copiotroph)),alpha=0.75,lwd=0.5) + 
  scale_x_log10(limits=c(0.1,10)) + 
  theme_pubclean() + 
  geom_vline(xintercept=1,col="black",lty=2,lwd=1.2) + 
  #geom_vline(xintercept=mean(y$dSHE[y$Copiotroph==T]/y$dS[y$Copiotroph==T],na.rm=T),col="red",lty=1) + 
  #geom_vline(xintercept=mean(y$dSHE[y$Copiotroph==F]/y$dS[y$Copiotroph==F],na.rm=T),col="blue",lty=1) + 
  xlab(expression(dS[HE]/dS)) +
  scale_fill_brewer(palette="Set1",labels=c("Oligotroph","Copiotroph"),direction=-1) +
  labs(fill="") + 
  theme(legend.position = c(0.8,0.8))


p2 <- ggplot(y,aes(x=factor(Copiotroph),y=dSHE/dS,fill=Copiotroph)) + 
  geom_boxplot() + 
  scale_y_log10() + 
  geom_hline(yintercept=1,lty=2) +
  theme_pubclean() + 
  xlab("") + 
  ylab(expression(dS[HE]/dS))+
  scale_x_discrete(labels=c("FALSE"="Oligotroph","TRUE"="Copiotroph")) + 
  scale_fill_brewer(palette="Set1",labels=c("Oligotroph","Copiotroph"),direction=-1) + 
  theme(text = element_text(size=12),
        axis.text.x = element_text(angle=30, hjust=1),
        legend.position = "none")

# Against Ne and s proxies -----------------------------------------------------

z <- merge.easy(y,dnds,key="Pair") %>% subset(dNdS < 10 & dNdS > 0.001)
hist(log10(z$dNdS),breaks=100)

p3 <- ggplot(z,aes(x=(1/d)*(1/dNdS),y=dSHE/dS,fill=Copiotroph)) +
  geom_point(alpha=0.1,pch=21) + 
  scale_y_log10() + scale_x_log10() + 
  geom_hline(yintercept=1,lty=2,col="black") +
  geom_smooth(aes(color=Copiotroph),method="lm") + 
  theme_classic2() +
  labs(fill="",color="") + 
  theme(legend.position = c(0.2,0.9)) + 
  xlab(expression("(" * 1/d * ")" * "*" * "(" * dS/dN * ")" * " (Proxy for " * sN[e] * ")")) + 
  ylab(expression(dS[HE]/dS))  + 
  scale_fill_brewer(palette="Set1",labels=c("Oligotroph","Copiotroph"),direction=-1) +
  scale_color_brewer(palette="Set1",labels=c("Oligotroph","Copiotroph"),direction=-1)

p4 <- ggplot(z,aes(x=(1/dNdS),y=dSHE/dS,fill=Copiotroph)) +
  geom_point(alpha=0.1,pch=21) + 
  scale_y_log10() + scale_x_log10() + 
  geom_hline(yintercept=1,lty=2,col="black") +
  geom_smooth(aes(color=Copiotroph),method="lm") + 
  theme_classic2() +
  labs(fill="",color="") + 
  theme(legend.position = "none") + 
  xlab(expression("(" * dS/dN * ")" * " (Proxy for " * N[e] * ")")) + 
  ylab(expression(dS[HE]/dS))  + 
  scale_fill_brewer(palette="Set1",labels=c("Oligotroph","Copiotroph"),direction=-1) +
  scale_color_brewer(palette="Set1",labels=c("Oligotroph","Copiotroph"),direction=-1)

# Save -------------------------------------------------------------------------

setwd("~/eggo/Figs")
png("dSHEdSNe.png",width=800,height=1000)
ggarrange(ggarrange(p1,p2,ncol=2,labels=c("(a)","(b)"),widths = c(3,1)),
          ggarrange(p3,p4,ncol=2,labels=c("(c)","(d)"),widths=c(2,1)),
          nrow=2)
dev.off()

# Test significance ------------------------------------------------------------

hist((z$dNdS))
hist(log10(z$dNdS))

hist(z$dSHE/z$dS)
hist(log10(z$dSHE/z$dS))

summary(lm(log10(dSHE/dS)~log10((1/d)*(1/dNdS)),data=z))
summary(lm(log10(dSHE/dS)~log10((1/d)*(1/dNdS)),data=z[z$Copiotroph==T,]))

wilcox.test((dSHE/dS)~Copiotroph,data=y)
t.test((dSHE/dS)~Copiotroph,data=y)

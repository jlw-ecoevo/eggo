


## JLW 2020 

set.seed(123456)

# Load Packages ----------------------------------------------------------------

library(dplyr)
library(data.table)
library(ggplot2)
library(ggpubr)
library(MASS)

# Helper Functions -------------------------------------------------------------

merge.easy <- function(df1,df2,key){
  df1 <- data.table(df1,key=key)
  df2 <- data.table(df2,key=key)
  return(unique(merge(df1,df2,all.x=TRUE,by=.EACHI,allow.cartesian=TRUE)))
}

rgrep <- function(big,small_vec){
  small_vec[lapply(small_vec,grepl,x=big) %>% unlist()]
}

boxcoxTransform <- function(x, lambda, back_transform = F) {
  if (back_transform == TRUE) {
    (x*lambda +1)^(1/lambda)  %>% return()
  } else {
    (((x^lambda) - 1) / lambda) %>% return()
  }
}


# Compile Dataset --------------------------------------------------------------

setwd("~/eggo/Data/eggo-data/")
load("CodonStatistics_RefSeq.RData")

setwd("~/eggo/Data/gRodon-files/")
load("Accession2Species.rda")
load("GrowthRates.rda")
load("CodonStatistics.rda")
cu <- cu %>% mutate_all(unlist)

# Merge growth data
rownames(spp_acc) <- spp_acc$V1 %>% gsub(pattern="[.].*",replace="")
cu$Accession <- cu$File %>% gsub(pattern="[.].*",replace="")
cu$Species <- spp_acc[cu$Accession,"V2"]
cu$Species <- lapply(cu$Species,rgrep,small_vec=d$Species) %>%
  lapply("[",1) %>% unlist()
cu <- merge.easy(cu,d,key="Species")

# Merge genera
genus_df <- assembly_df %>% subset(select=c(Assembly,Genus,SppName)) %>% 
  mutate(Accession = gsub("[.].*","",Assembly))
cs_df <- merge.easy(cu,genus_df,key="Accession") %>% subset(!is.na(d))


# Actual d ---------------------------------------------------------------------

sppsub_df <- cs_df %>% subset(!is.na(Species)) %>% group_by(Species) %>% sample_n(1)
genus_multi <- names(table(sppsub_df$Genus))[table(sppsub_df$Genus)>1]
pair_df <- sppsub_df %>% subset(Genus %in% genus_multi) %>% 
  group_by(Genus) %>% sample_n(2) %>% group_by(Genus) %>% 
  summarise(d1=d[1],d2=d[2]) %>% subset(Genus!="")
pair_df$Ratio <- pair_df$d1/pair_df$d2
pair_df$Diff <- pair_df$d1 - pair_df$d2

p1 <- ggplot(pair_df,aes(x=d1,y=d2)) + geom_point(alpha=0.5) + 
  scale_x_log10() + scale_y_log10() + theme_bw() + geom_smooth() + 
  geom_abline(intercept=0,slope=1,lty=2) + 
  xlab("Minimal Doubling Time of Organism 1") + 
  ylab("Minimal Doubling Time of Organism 2")
bc_milc <- boxcox(d1~d2,data=pair_df)
lambda <- bc_milc$x[which.max(bc_milc$y)]
lm(boxcoxTransform(d1, lambda) ~ d2,data=pair_df) %>% summary()


# Inferred d Pairs -------------------------------------------------------------

#Remove any large outliers (d>100)
sppsub_df <- assembly_df %>% subset(d<100) %>% group_by(SppName) %>% sample_n(1)
genus_multi <- names(table(sppsub_df$Genus))[table(sppsub_df$Genus)>1]
pair_df <- sppsub_df %>% subset(Genus %in% genus_multi) %>% 
  group_by(Genus) %>% sample_n(2) %>% group_by(Genus) %>% 
  summarise(d1=d[1],d2=d[2]) %>% subset(Genus!="")
pair_df$Ratio <- pair_df$d1/pair_df$d2
pair_df$Diff <- pair_df$d1 - pair_df$d2

p2 <- ggplot(pair_df,aes(x=d1,y=d2)) + geom_point(alpha=0.5) + 
  scale_x_log10() + scale_y_log10() + theme_bw() + geom_smooth() + 
  geom_abline(intercept=0,slope=1,lty=2) + 
  xlab("Estimated Minimal Doubling Time of Organism 1") + 
  ylab("Estimated Minimal Doubling Time of Organism 2")
bc_milc <- boxcox(d1~d2,data=pair_df)
lambda <- bc_milc$x[which.max(bc_milc$y)]
lmsum <- lm(boxcoxTransform(d1, lambda) ~ d2,data=pair_df) %>% summary()
lmsum$coefficients
lmsum$adj.r.squared

# Put it together --------------------------------------------------------------

setwd("~/eggo/Figs")
pdf("Pair_Prediction.pdf",width=10,height=5)
ggarrange(p1,p2,ncol=2,labels=c("(a)","(b)"))
dev.off()

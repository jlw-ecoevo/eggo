

## JLW 2020 - Does streamlining impact our model?

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



# CUB, Consistency, CPB
setwd("~/eggo/Data/gRodon-files/")
load("CodonStatistics.rda")
load("Accession2Species.rda")
load("GrowthRates.rda")
cu <- cu %>% mutate_all(unlist)
# d <- d %>%  subset(Extremophile!="X")

setwd("~/eggo/Data/")

#genome lengths
genome_len <- read.delim("genome_lengths_vs.txt", sep = "\t", header = F)
names(genome_len) <- c("Name","Length")
genome_len$Accession <- genome_len$Name %>% gsub(pattern="[.].*",replace="")
genome_len$Name <- NULL

#length coding and GC
codgc <- read.delim("gc_and_length_vs.tsv", sep = "\t", header = F)
names(codgc) <- c("Name","LengthCoding","GC")
codgc$Accession <- codgc$Name %>% gsub(pattern="[.].*",replace="")
codgc$Name <- NULL

gclen <- merge.easy(genome_len,codgc,key="Accession")
gclen$PercentCoding <- gclen$LengthCoding/gclen$Length

cu$Accession <- cu$File %>% gsub(pattern="[.].*",replace="")
cs_df <- merge.easy(gclen,cu,key="Accession")

#Match genome accessions with species for dataset matching
rownames(spp_acc) <- spp_acc$V1
cs_df$Species <- spp_acc[cs_df$Accession,"V2"]
cs_df$Species <- lapply(cs_df$Species,rgrep,small_vec=d$Species) %>% 
  lapply("[",1) %>% unlist()

# Merge datasets
cs_df <- merge.easy(cs_df,d,key="Species")

# Average CUB estimates over species, add phylum info
vs_avg <- cs_df %>% group_by(Species) %>% 
  summarise_all(mean,na.rm=T) %>% 
  subset(!is.na(Species)) 

#Get Resulduals
bc_milc <- boxcox(d~CUBHE+ConsistencyHE+CPB,data=vs_avg)
lambda_milc <- bc_milc$x[which.max(bc_milc$y)]
gRodon_model_base <- 
  lm(boxcoxTransform(d, lambda_milc) ~ CUBHE+ConsistencyHE+CPB,data=vs_avg)
vs_avg$Residual <- gRodon_model_base$residuals


p1 <- ggplot(vs_avg,aes(x=GC,y=d)) + 
  geom_point() + stat_smooth() + scale_y_log10() + 
  theme_bw() + ylab("Doubling Time (Hours)")
p2 <- ggplot(vs_avg,aes(x=PercentCoding,y=d)) + 
  geom_point() + stat_smooth() + scale_y_log10() + 
  theme_bw() + ylab("Doubling Time (Hours)") + 
  geom_hline(yintercept = 0,lty=2) + 
  xlab("Percent Coding")

p3 <- ggplot(vs_avg,aes(x=GC,y=Residual)) + 
  geom_point() + stat_smooth() + theme_bw() + 
  geom_hline(yintercept = 0,lty=2)
p4 <- ggplot(vs_avg,aes(x=PercentCoding,y=Residual)) + 
  geom_point() + stat_smooth() + theme_bw() + 
  geom_hline(yintercept = 0,lty=2) + 
  xlab("Percent Coding")

setwd("~/eggo/Figs")
pdf("gRodon_streamlining.pdf",width=7,height=7)
ggarrange(ggarrange(p1,p2,nrow=1),
          ggarrange(p3,p4,nrow=1),
          nrow=2)
dev.off()



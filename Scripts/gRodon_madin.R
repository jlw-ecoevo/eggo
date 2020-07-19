
library(dplyr)
library(data.table)
library(ggplot2)
library(ggpubr)
library(ggpointdensity)

merge.easy <- function(df1,df2,key){
  df1 <- data.table(df1,key=key)
  df2 <- data.table(df2,key=key)
  return(unique(merge(df1,df2,all.x=TRUE,by=.EACHI,allow.cartesian=TRUE)))
}

rgrep <- function(big,small_vec){
  small_vec[lapply(small_vec,grepl,x=big) %>% unlist()]
}


setwd("~/eggo/Data")
load("gRodon-files/GrowthRates_Madin.rda")
dm_all <- d
load("gRodon-files/GrowthRates.rda")
d_all <- unique(d)
load("eggo-data/CodonStatistics_Training.RData")
load("gRodon-files/Accession2Species_Madin.rda")
spp_acc_m <- spp_acc
load("gRodon-files/Accession2Species.rda")
load("EGGO.RData")


d <- d_all %>% subset(Extremophile=="")
dm <- dm_all %>% subset(Extremophile==FALSE)

rownames(spp_acc) <- spp_acc$V1 %>% gsub(pattern="[.].*",replace="")
vs_df$Accession <- vs_df$Assembly %>% gsub(pattern="[.].*",replace="")
vs_df$Spp <- spp_acc[vs_df$Accession,"V2"]
vs_df$Species <- lapply(vs_df$Spp,rgrep,small_vec=d$Species) %>%
  lapply("[",1) %>% unlist()
vs_df$Species[vs_df$Spp %in% d$Species] <- vs_df$Spp[vs_df$Spp %in% d$Species]
names(vs_df)[names(vs_df)=="d"] <- "d.vs"
vs_df <- merge.easy(vs_df,d,key="Species") %>% subset(!is.na(Species))
stat_data <- vs_df %>%
  group_by(Species) %>%
  summarise_all(mean,na.rm=T) %>%
  subset(!is.na(Species))
plot(stat_data$d,stat_data$d.madin,log="xy")

rownames(spp_acc_m) <- spp_acc_m$V1 %>% gsub(pattern="[.].*",replace="")
madin_df$Accession <- madin_df$Assembly %>% gsub(pattern="[.].*",replace="")
madin_df$Spp <- spp_acc_m[madin_df$Accession,"V2"]
dm$Species <- dm$species
madin_df$Species <- lapply(madin_df$Spp,rgrep,small_vec=dm$Species) %>%
  lapply("[",1) %>% unlist()
madin_df$Species[madin_df$Spp %in% dm$Species] <- madin_df$Spp[madin_df$Spp %in% dm$Species]
names(madin_df)[names(madin_df)=="d"] <- "d.vs"
madin_df <- merge.easy(madin_df,dm,key="Species") %>% subset(!is.na(Species))
stat_data_madin <- madin_df %>%
  group_by(Species) %>%
  summarise_all(mean,na.rm=T) %>%
  subset(!is.na(Species))
plot(stat_data_madin$d,stat_data_madin$d.vs,log="xy")

pMM <- ggplot(stat_data_madin,aes(x=d,y=d.madin)) + geom_point(alpha=0.5) + 
  scale_x_log10() + scale_y_log10() + theme_pubclean() + 
  geom_smooth(color="darkgrey") + xlab("Empirical Min. Doubling Time (Madin et al.)") + 
  ylab("Predicted Min. Doubling Time (Trained on Madin et al.)") + 
  geom_abline(slope = 1,intercept = 0,lty=2) + 
  geom_vline(xintercept = 5,lty=2,color="red")

pVSVS <- ggplot(stat_data,aes(x=d,y=d.vs)) + geom_point(alpha=0.5) + 
  scale_x_log10() + scale_y_log10() + theme_pubclean() + 
  geom_smooth(color="darkgrey") + xlab("Empirical Min. Doubling Time (Vieira-Silva et al.)") + 
  ylab("Predicted Min. Doubling Time (Trained on Vieira-Silva et al.)") + 
  geom_abline(slope = 1,intercept = 0,lty=2) + 
  geom_vline(xintercept = 5,lty=2,color="red")

pMVS <- ggplot(stat_data,aes(x=d,y=d.madin)) + geom_point(alpha=0.5) + 
  scale_x_log10() + scale_y_log10() + theme_pubclean() + 
  geom_smooth(color="darkgrey") + xlab("Empirical Min. Doubling Time (Vieira-Silva et al.)") + 
  ylab("Predicted Min. Doubling Time (Trained on Madin et al.)") + 
  geom_abline(slope = 1,intercept = 0,lty=2) + 
  geom_vline(xintercept = 5,lty=2,color="red")

pVSM <- ggplot(stat_data_madin,aes(x=d,y=d.vs)) + geom_point(alpha=0.5) + 
  scale_x_log10() + scale_y_log10() + theme_pubclean() + 
  geom_smooth(color="darkgrey") + xlab("Empirical Min. Doubling Time (Madin et al.)") + 
  ylab("Predicted Min. Doubling Time (Trained on Vieira-Silva et al.)") + 
  geom_abline(slope = 1,intercept = 0,lty=2) + 
  geom_vline(xintercept = 5,lty=2,color="red")

setwd("~/eggo/Figs")
pdf("Madin_vs_VieiraSilva_fits.pdf",width=10,height=10)
ggarrange(pVSVS,pVSM,
          pMVS,pMM,
          nrow=2,
          ncol=2,
          labels=c("(a)",
                   "(b)",
                   "(c)",
                   "(d)"), 
          hjust = -1, 
          vjust = 1)
dev.off()

# Actual Comparison ------------------------------------------------------------

d <- d_all
dm <- dm_all

rownames(d) <- d$Species
dm$Species <- NA
dm$Species <- d[dm$species,"Species"]
dm$SpeciesInd <- NA
dm$SpeciesInd[is.na(dm$Species)] <- 
  lapply(dm$species[is.na(dm$Species)],grep,x=d$Species) %>% lapply("[",1) %>% unlist()
dm$Species[!is.na(dm$SpeciesInd)] <- d$Species[dm$SpeciesInd[!is.na(dm$SpeciesInd)]]
dmd <- merge.easy(dm,d,key="Species")

p1 <- ggplot(dmd, aes(x = d.y, y = d.x)) +
  geom_point() + 
  scale_y_log10(limits=c(0.1,500)) + scale_x_log10(limits=c(0.1,500)) + 
  geom_abline(slope=1, intercept = 0, lty = 2) + 
  geom_smooth(color = "grey") + 
  xlab("Empirical Min. Doubling Time (Vieira-Silva et al.)") + 
  ylab("Empirical Min. Doubling Time (Madin et al.)") + 
  theme_classic2()

# ----------------------------------------------------------------------------
  
p2 <- ggplot(EGGO, aes(x = d, y = d.madin)) +
  geom_pointdensity(adjust=1) + 
  scale_colour_gradient(low="white",high="black",trans = "log") +
  scale_y_log10(limits=c(0.1,500)) + scale_x_log10(limits=c(0.1,500)) + 
  geom_abline(slope=1, intercept = 0, lty = 2) + 
  geom_smooth(color = "grey") + 
  theme_classic2() +
  theme(legend.position = c(0.3,0.7)) + 
  xlab("Predicted Min. Doubling Time (Trained on Vieira-Silva et al.)") + 
  ylab("Predicted Min. Doubling Time (Trained on Madin et al.)")  + 
  labs(color = "Neighbor Density")


setwd("~/eggo/Figs")
#pdf("Madin_vs_VieiraSilva.pdf",width=12,height=6)
tiff(file = "Madin_vs_VieiraSilva.tiff", 
     width = 12, 
     height = 6, 
     units = "in", 
     res = 800,
     compression="lzw")
ggarrange(p1,p2,ncol=2,labels=c("(a)","(b)"), hjust = 0, vjust = 1)
dev.off()

cor.test(EGGO$d,EGGO$d.madin)
cor.test(dmd$d.x,dmd$d.y)

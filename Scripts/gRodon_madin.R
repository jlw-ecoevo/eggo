
## JLW 2020 - Comparison Against Growrthpred Results on Original Viera-Silva et 
## al. Dataset

# Load Packages ----------------------------------------------------------------

library(dplyr)
library(data.table)
library(MASS)
library(ggplot2)
library(ggpubr)
library(car)
library(lmtest)
library(reshape2)
library(psych)

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

setwd("~/eggo/Data")

# gRodon fitting data
load("gRodon-files/GrowthRates_Madin.rda")
load("gRodon-files/CodonStatistics_Madin.rda")
load("gRodon-files/Accession2Species_Madin.rda")
cu <- cu %>% mutate_all(unlist)
names(d)[1] <- "Species"
d <- d %>% as.data.frame(stringsAsFactors=F)

# Merge datasets
rownames(spp_acc) <- spp_acc$V1 %>% gsub(pattern="[.].*",replace="")
cu$Accession <- cu$File %>% gsub(pattern="[.].*",replace="")
cu$Spp <- spp_acc[cu$Accession,"V2"]
cu$Species <- lapply(cu$Spp,rgrep,small_vec=d$Species) %>%
  lapply("[",1) %>% unlist()
cu$Species[cu$Spp %in% d$Species] <- cu$Spp[cu$Spp %in% d$Species]
cu <- merge.easy(cu,d,key="Species") %>% subset(!is.na(Species))

# Average CUB estimates over species
stat_data <- cu %>%
  subset(Extremophile == FALSE) %>%
  group_by(Species) %>%
  summarise_all(mean,na.rm=T) %>%
  subset(!is.na(Species))


# Box Cox Transformation -------------------------------------------------------

#linear models
m_milc <- lm(d~CUBHE+ConsistencyHE+CPB,data=stat_data)
#transformation
bc_milc <- boxcox(d~CUBHE+ConsistencyHE+CPB,data=stat_data)
lambda_milc <- bc_milc$x[which.max(bc_milc$y)]

# re-run with transformation
mnew_milc <- 
  lm(boxcoxTransform(d, lambda_milc) ~ CUBHE+ConsistencyHE+CPB,data=stat_data)

#look at residuals
ggqqplot(m_milc$residuals)
ggqqplot(mnew_milc$residuals)

smm <- summary(mnew_milc)
smm$coefficients
smm$r.squared

# Likelihood Ratio Test --------------------------------------------------------

# Nested Models
milc1 <- 
  lm(boxcoxTransform(d, lambda_milc) ~ CUBHE,data=stat_data)
milc2 <- 
  lm(boxcoxTransform(d, lambda_milc) ~ CUBHE+ConsistencyHE,data=stat_data)
milc3 <- 
  lm(boxcoxTransform(d, lambda_milc) ~ CUBHE+ConsistencyHE+CPB,data=stat_data)

# LR Test
m.test1 <- lrtest(milc1, milc2)
m.test2 <- lrtest(milc1, milc3)
m.test3 <- lrtest(milc2, milc3)

m.test1
m.test2
m.test3

# Model Comparison/Performance -------------------------------------------------

summary(mnew_milc)
cor(boxcoxTransform(stat_data$d, lambda_milc),mnew_milc$fitted.values)

plot(mnew_milc$fitted.values,boxcoxTransform(stat_data$d, lambda_milc))

#Full gRodon
gRodon_model_base <- 
  lm(boxcoxTransform(d, lambda_milc) ~ CUBHE+ConsistencyHE+CPB,data=stat_data)

# Partial genome mode
gRodon_model_partial <- 
  lm(boxcoxTransform(d, lambda_milc) ~ CUBHE+ConsistencyHE,data=stat_data)

# Metagenome mode
gRodon_model_meta <- 
  lm(boxcoxTransform(d, lambda_milc) ~ CUBHE,data=stat_data)

# Save compiled dataset
stat_data$dGR <- boxcoxTransform(gRodon_model_base$fitted.values,
                              lambda_milc,
                              back_transform = TRUE)
stat_data$dGRP <- boxcoxTransform(gRodon_model_partial$fitted.values,
                              lambda_milc,
                              back_transform = TRUE)
stat_data$dGRM <- boxcoxTransform(gRodon_model_meta$fitted.values,
                               lambda_milc,
                               back_transform = TRUE)
stat_data$Residual <- gRodon_model_base$residuals
stat_data$ResidualM <- gRodon_model_meta$residuals
save(stat_data,file="stat_data_madin.RData")

# Plot Performance Nicely ------------------------------------------------------


p1 <- ggplot(stat_data,aes(x=d,y=dGR)) + geom_point(alpha=0.5) + 
  scale_x_log10() + scale_y_log10() + theme_bw() + 
  geom_smooth(color="darkgrey") + xlab("Empirical Minimal Doubling Time (Hours)") + 
  ylab("Predicted Minimal Doubling Time (Hours)") + 
  geom_abline(slope = 1,intercept = 0,lty=2) + 
  geom_vline(xintercept = 5,lty=2,color="red")


p2 <- ggplot(stat_data,aes(x=d,y=CUBHE)) + geom_point(alpha=0.5) + 
  scale_x_log10()  + theme_bw() + 
  geom_smooth(color="darkgrey") + xlab("Empirical Minimal Doubling Time (Hours)") + 
  ylab("Codon Usage Bias (Ribosomal Proteins)") + 
  geom_vline(xintercept = 5,lty=2,color="red")

p3 <- ggplot() +theme_minimal()

setwd("~/eggo/Figs")
pdf("gRodon_performance_Madin.pdf",width=8.75,height=3.75)
ggarrange(p3,
          ggarrange(p1,p2,nrow=1,
                    labels = c("(a)","(b)"),
                    hjust=-1,
                    vjust=-0.5),
          nrow=2,heights = c(1,9))
dev.off()


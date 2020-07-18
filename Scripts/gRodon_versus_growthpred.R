
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
load("gRodon-files/GrowthRates.rda")
load("gRodon-files/CodonStatistics.rda")
load("gRodon-files/Accession2Species.rda")
cu <- cu %>% mutate_all(unlist)

# Growthpred on VS genomes
vsgp <-  read.delim("growthpred-predictions/growthpred_VS.tbl",
                    head = F,stringsAsFactors = F,sep="\t")
vsgpm <-  read.delim("growthpred-predictions/growthpred_VS_metagenome.tbl",
                     head = F,stringsAsFactors = F,sep="\t")

# Phylum info
spp_phy <- read.delim("phyla.tbl",stringsAsFactors = F,header = F)
names(spp_phy) <- c("Spp","Phylum")

# Merge growth data
rownames(spp_acc) <- spp_acc$V1 %>% gsub(pattern="[.].*",replace="")
cu$Accession <- cu$File %>% gsub(pattern="[.].*",replace="")
cu$Species <- spp_acc[cu$Accession,"V2"]
cu$Species <- lapply(cu$Species,rgrep,small_vec=d$Species) %>%
  lapply("[",1) %>% unlist()
cu <- merge.easy(cu,d,key="Species")

# Merge growthpred predictions
names(vsgp) <- c("Accession","dGP")
vsgp$Accession <- vsgp$Accession %>% gsub(pattern="[.].*",replace="")
cu <- merge.easy(cu,vsgp,key="Accession")
names(vsgpm) <- c("Accession","dGPM")
vsgpm$Accession <- vsgpm$Accession %>% gsub(pattern="[.].*",replace="")
cu <- merge.easy(cu,vsgpm,key="Accession")

# Average CUB estimates over species
stat_data <- cu %>%
  subset(Extremophile == "") %>%
  group_by(Species) %>%
  summarise_all(mean,na.rm=T) %>%
  subset(!is.na(Species))

#Merge phylum
stat_data$Spp <- stat_data$Species %>% gsub(pattern=" ",replace="")
stat_data <- merge(stat_data,spp_phy,key="Spp")

# Box Cox Transformation -------------------------------------------------------

#linear models
m_milc <- lm(d~CUBHE+ConsistencyHE+CPB,data=stat_data)
m_gp <- lm(d~dGP,data=stat_data)

#transformation
bc_milc <- boxcox(d~CUBHE+ConsistencyHE+CPB,data=stat_data)
lambda_milc <- bc_milc$x[which.max(bc_milc$y)]
bc_gp <- boxcox(d~dGP,data=stat_data)
lambda_gp <- bc_gp$x[which.max(bc_gp$y)]

# re-run with transformation
mnew_milc <- 
  lm(boxcoxTransform(d, lambda_milc) ~ CUBHE+ConsistencyHE+CPB,data=stat_data)
mnew_gp <- 
  lm(boxcoxTransform(d, lambda_gp) ~ boxcoxTransform(dGP, lambda_gp),data=stat_data)

#look at residuals
ggqqplot(m_milc$residuals)
ggqqplot(mnew_milc$residuals)
ggqqplot(m_gp$residuals)
ggqqplot(mnew_gp$residuals)

smm <- summary(mnew_milc)
smm$coefficients
smm$r.squared

# ANOVA COmparison of Models ---------------------------------------------------

anova(lm(boxcoxTransform(d, lambda_milc) ~ boxcoxTransform(dGP, lambda_gp),data=stat_data),
      lm(boxcoxTransform(d, lambda_milc) ~ CUBHE+ConsistencyHE+CPB,data=stat_data))

anova(lm(boxcoxTransform(d, lambda_gp) ~ boxcoxTransform(dGP, lambda_gp),data=stat_data),
      lm(boxcoxTransform(d, lambda_gp) ~ CUBHE+ConsistencyHE+CPB,data=stat_data))

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

summary(mnew_gp)
cor(boxcoxTransform(stat_data$d, lambda_gp),mnew_gp$fitted.values)

plot(mnew_milc$fitted.values,boxcoxTransform(stat_data$d, lambda_milc))
plot(mnew_gp$fitted.values,boxcoxTransform(stat_data$d, lambda_gp))


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
save(stat_data,file="stat_data.RData")

# Plot Performance Nicely ------------------------------------------------------


p1 <- ggplot(stat_data,aes(x=d,y=dGR)) + geom_point(alpha=0.5) + 
  scale_x_log10() + 
  scale_y_log10() + 
  theme_pubclean() + 
  geom_smooth(color="darkgrey") + xlab("Empirical Minimal Doubling Time (Hours)") + 
  ylab("Predicted Minimal Doubling Time (Hours)") + 
  geom_abline(slope = 1,intercept = 0,lty=2) + 
  geom_vline(xintercept = 5,lty=2,color="red")


p2 <- ggplot(stat_data,aes(x=d,y=CUBHE)) + geom_point(alpha=0.5) + 
  scale_x_log10()  + theme_pubclean() + 
  geom_smooth(color="darkgrey") + xlab("Empirical Minimal Doubling Time (Hours)") + 
  ylab("Codon Usage Bias (Ribosomal Proteins)") + 
  geom_vline(xintercept = 5,lty=2,color="red")

p3 <- ggplot() +theme_minimal()

setwd("~/eggo/Figs")
pdf("gRodon_performance.pdf",width=8.75,height=3.75)
ggarrange(p3,
          ggarrange(p1,p2,nrow=1,
                    labels = c("(a)","(b)"),
                    hjust=0,
                    vjust=-1.5),
          nrow=2,heights = c(1,9))
dev.off()


pS1 <- ggplot(stat_data,aes(x=d,y=dGR)) + 
  geom_point(aes(y = dGR, col = "gRodon"),alpha=0.5) + 
  geom_smooth(aes(y = dGR, col = "gRodon"), se = F) +
  geom_point(aes(y = dGP, col = "growthpred"),alpha=0.5) + 
  geom_smooth(aes(y = dGP, col = "growthpred"), se = F) +
  #geom_smooth(aes(y = dGPM, col = "growthpred M"), se = F) +
  scale_x_log10() + scale_y_log10() + theme_bw() + 
  geom_abline(slope = 1,intercept = 0,lty=2) + 
  xlab("Actual Minimal Doubling Time (Hours)") + 
  ylab("Predicted Minimal Doubling Time (Hours)") + 
  labs(color = "Prediction Method") 



# Plain Old Cross Validation ---------------------------------------------------


doCV <- function(stat_data,lambda_gp,n_fold=10){
  stat_data$Fold <- sample(1:n_fold,nrow(stat_data),replace = T)
  err <- list()
  errp <- list()
  errm <- list()
  err_gp <- list()
  err_gpm <- list()
  for(f in 1:n_fold){
    #Train and test set
    train_df <- stat_data %>% subset(Fold != f)
    test_df <- stat_data %>% subset(Fold == f)
    
    # Box cox
    bc_p <- boxcox(d ~ CUBHE + ConsistencyHE + CPB, data = train_df)
    lambda_p <- bc_p$x[which.max(bc_p$y)]
    
    # Predict fold
    p_lm <- lm(boxcoxTransform(d, lambda_p) ~ 
                 CUBHE + ConsistencyHE + CPB, data = train_df)
    pred <- predict.lm(p_lm, newdata = test_df)
    
    p_lm <- lm(boxcoxTransform(d, lambda_p) ~ 
                 CUBHE + ConsistencyHE, data = train_df)
    predp <- predict.lm(p_lm, newdata = test_df)
    
    p_lm <- lm(boxcoxTransform(d, lambda_p) ~ 
                 CUBHE, data = train_df)
    predm <- predict.lm(p_lm, newdata = test_df)
    
    #Calculate errors
    err[[f]] <- mean((pred - boxcoxTransform(test_df$d, lambda_p))^2)
    errm[[f]] <- mean((predm - boxcoxTransform(test_df$d, lambda_p))^2)
    errp[[f]] <- mean((predp - boxcoxTransform(test_df$d, lambda_p))^2)
    err_gp[[f]] <- mean((boxcoxTransform(test_df$dGP, lambda_gp) -
                           boxcoxTransform(test_df$d, lambda_gp))^2)
    err_gpm[[f]] <- mean((boxcoxTransform(test_df$dGPM, lambda_gp) -
                           boxcoxTransform(test_df$d, lambda_gp))^2)
  }
  return(list(gRodon=mean(unlist(err)),
              gRodon_partial=mean(unlist(errp)),
              gRodon_meta=mean(unlist(errm)),
              growthpred=mean(unlist(err_gp)),
              growthpred_meta=mean(unlist(err_gpm))))
}

cv_rep <- replicate(100,doCV(stat_data,lambda_gp))
err_df <- data.frame(gRodon=cv_rep["gRodon",]%>%unlist(),
                     gRodonP=cv_rep["gRodon_partial",]%>%unlist(),
                     gRodonM=cv_rep["gRodon_meta",]%>%unlist(),
                     GP=cv_rep["growthpred",]%>%unlist(),
                     GPM=cv_rep["growthpred_meta",]%>%unlist(),
                     stringsAsFactors = FALSE)

x <- melt(err_df)

p1 <- ggplot(x,aes(x=variable,group=variable,y=value)) + 
  geom_boxplot() + theme_pubclean() + 
  ylab("MSE from Cross Validation") + xlab("") + 
  scale_x_discrete(labels = c('gRodon',
                              'gRodon partial mode',
                              'gRodon metagenome mode',
                              'growthpred',
                              'growthpred metagenome mode')) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))


# Blocked Cross Validation with Phyla ------------------------------------------

phyla <- unique(stat_data$Phylum)

err <- list()
err_gp <- list()
err_gr <- list()
n <- numeric()
for(p in phyla){
  #Train and test set
  train_df <- stat_data %>% subset(Phylum != p)
  test_df <- stat_data %>% subset(Phylum == p)
  
  # Box cox
  bc_p <- boxcox(d ~ CUBHE + ConsistencyHE + CPB, data = train_df)
  lambda_p <- bc_p$x[which.max(bc_p$y)]
  
  # Predict fold
  p_lm <- lm(boxcoxTransform(d, lambda_p) ~ 
               CUBHE + ConsistencyHE + CPB, data = train_df)
  pred <- predict.lm(p_lm, newdata = test_df)
  
  #Calculate errors
  err[[p]] <- mean((pred - boxcoxTransform(test_df$d, lambda_p))^2)
  err_gp[[p]] <- mean((boxcoxTransform(test_df$dGP, lambda_gp) -
                         boxcoxTransform(test_df$d, lambda_gp))^2)
  err_gr[[p]] <- mean((boxcoxTransform(test_df$dGR, lambda_milc) -
                         boxcoxTransform(test_df$d, lambda_milc))^2)
  n[p] <- nrow(test_df)
}

err_df_blocked <- data.frame(GP=unlist(err_gp),
                             GR=unlist(err_gr),
                             gRodon=unlist(err),
                             n=n,
                             stringsAsFactors = F)

p2 <- ggplot(err_df_blocked,aes(y=GP,x=gRodon,size=n)) + 
  geom_point() + geom_abline(slope=1,intercept=0,lty=2) + 
  theme_classic2() + ylab("MSE growthpred") + 
  xlim(0,2) + ylim(0,2) + box() + 
  xlab("Blocked CV MSE gRodon (By Phylum)") + 
  theme(legend.position = c(0.85,0.5),
        legend.box.background = element_rect(color="black", size=1)) 

# Performance Panel Figure -----------------------------------------------------

setwd("~/eggo/Figs")
pdf("gRodon_vs_growthpred_loess.pdf",width=7,height=3.5)
pS1
dev.off()


setwd("~/eggo/Figs")
pdf("gRodon_vs_growthpred_cv.pdf",width=8,height=4)
ggarrange(p1,
                    p2,
                    nrow=1,
                    labels = c("(a)","(b)"),
                    hjust=0,
                    vjust=1)
dev.off()

# Classifier -------------------------------------------------------------------

p1 <- ggplot(stat_data,aes(x=dGR)) + 
  scale_x_log10() + 
  geom_density(aes(x=dGR,fill="Predictions (gRodon)"),alpha=0.5,lwd=1) + 
  geom_density(aes(x=d,fill="Actual Values"),alpha=0.5,lwd=1) +
  theme_bw() + geom_vline(xintercept=5,lty=2,color="red") + 
  xlab("Minimal Doubling Time (Hours)") + 
  theme(legend.title = element_blank()) + 
  scale_fill_brewer(palette = "Accent")

getKappa <- function(thresh){
  cohen.kappa(data.frame(dGR=log10(stat_data$dGR)>thresh,
                         d=log10(stat_data$d)>thresh))
}

x <- seq(-0.6,2,0.1)
y <- lapply(x,getKappa) %>% lapply("[[","kappa") %>% unlist()
plot(x,y,type="l")
abline(v=log10(5))


p2 <- ggplot(data.frame(x=10^x,y=y),aes(x=x,y=y)) + geom_line() + 
  theme_bw() + ylab(expression("Cohen's " * kappa)) + 
  xlab("Threshold Doubling Time (Hours)") + scale_x_log10()  + 
  geom_vline(xintercept=5,lty=2,color="red")

setwd("~/eggo/Figs")
pdf("Classifier.pdf",width=7,height=6)
ggarrange(p1,p2,nrow=2,
          labels = c("(a)","(b)"))
dev.off()


library(dplyr)
library(data.table)
library(ape)
library(ggplot2)
library(ggpubr)
library(castor)
library(MASS)
library(reshape2)
library(parallel)
library(psych)

set.seed(123456)

merge.easy <- function(df1,df2,key){
  df1 <- data.table(df1,key=key)
  df2 <- data.table(df2,key=key)
  return(unique(merge(df1,df2,all.x=TRUE,by=.EACHI,allow.cartesian=TRUE)) %>%
           as.data.frame())
}

boxcoxTransform <- function(x, lambda, back_transform = F) {
  if (back_transform == TRUE) {
    (x*lambda +1)^(1/lambda)  %>% return()
  } else {
    (((x^lambda) - 1) / lambda) %>% return()
  }
}

sqErrBoxCox <- function(d1,d2,lambda){
  (boxcoxTransform(d1,lambda = lambda) -
     boxcoxTransform(d2,lambda = lambda))^2
}

getNeighborhoodPrediction <- function(i, tree, values, n = 5, threshold = 1e-8){
  if(i%%100==0){print(i)}
  to_drop <- i
  nt <- numeric(n)
  dnt <- numeric(n)
  for(neighbor in 1:n){
    nt[neighbor] <- find_nearest_tips(tree, 
                                      target_tips = tree$tip.label[-to_drop])$nearest_tip_per_tip[i]
    dnt[neighbor] <- find_nearest_tips(tree, 
                                       target_tips = tree$tip.label[-to_drop])$nearest_distance_per_tip[i]
    to_drop <- c(to_drop, nt[neighbor])
  }
  
  tip_weights <- 1/(dnt+threshold)
  vals <- values[nt]
  return(c(Actual = values[i], 
           Predicted = weightedGeometricMean(vals, tip_weights), 
           MeanDistance = mean(dnt)))
}

weightedGeometricMean <- function(x, w, na.rm = TRUE, zero.propagate = FALSE){
  if(any(x < 0, na.rm = TRUE)){
    return(NaN)
    
  }
  if(zero.propagate){
    if(any(x == 0, na.rm = TRUE)){
      return(0)
      
    }
    return(exp(weighted.mean(log(x), w, na.rm=na.rm)))
    
  } else {
    nonzero_ind <- x>0 & w>0
    w <- w[nonzero_ind]
    x <- x[nonzero_ind]
    return(exp(weighted.mean(log(x), w, na.rm=na.rm)))
    
  }
}

setwd("~/eggo/Data")
load("gRodon-files/sysdata.rda")
load("EGGO.RData")
tip_names <- read.delim("16S_seqnames.txt", sep = "\t", header = F, stringsAsFactors = F)
names(tip_names) <- c("File","Tip")
tip_names$ID <- tip_names$File %>% 
  gsub(pattern = "_rna_from_.*", replace = "") %>% 
  gsub(pattern = ".fa.ffn.16S*", replace = "") %>% 
  gsub(pattern = "_prokka.ffn.16S*", replace = "")
EGGO <- EGGO %>% mutate(ID = Assembly) %>% 
  mutate_all(unlist) %>%
  subset(select=c(ID,d,LowerCI,UpperCI))
tip_meta <- merge.easy(tip_names,EGGO,key="ID")
rownames(tip_meta) <- tip_meta$Tip
x <- read.tree("16S.cut80.tree")


#Sample One 16S Per Organism
tip_meta_sub <- tip_meta %>% 
  subset(Tip %in% x$tip.label & !is.na(d)) %>% 
  group_by(File) %>% 
  sample_n(1) %>%
  as.data.frame(stringsAsFactors=FALSE)
rownames(tip_meta_sub) <- tip_meta_sub$Tip
x_sub <- drop.tip(x, which(!(x$tip.label %in% tip_meta_sub$Tip)))

# #Predict - very slow
# neighbor_pred_list <- mclapply(1:Ntip(x_sub), 
#                                getNeighborhoodPrediction, 
#                                tree=x_sub, 
#                                values=tip_meta_sub[x_sub$tip.label,"d"],
#                                mc.cores = 6)
# neighbor_pred <- do.call("rbind",neighbor_pred_list) %>%
#   as.data.frame(stringsAsFactors = F)
# setwd("~/eggo/Data")
# save(neighbor_pred, file = "neighbors.RData")
# 
# nt <- numeric(length(x$tip.label))
# dnt <- numeric(length(x$tip.label))
# for(i in 1:length(x$tip.label)){
#   if(i%%100==0){print(i)}
#   nt[i] <- find_nearest_tips(x, target_tips = x$tip.label[-i])$nearest_tip_per_tip[i]
#   dnt[i] <- find_nearest_tips(x, target_tips = x$tip.label[-i])$nearest_distance_per_tip[i]
# }
# 
# x_dist <- data.frame(Name = x$tip.label,
#                      Nearest = nt,
#                      NearestDist = dnt)
# x_dist$Name <- as.character(x_dist$Name)
# x_dist$NearestName <- x$tip.label[x_dist$Nearest]
# x_dist$d1 <- tip_meta[x_dist$Name,"d"]
# x_dist$d2 <- tip_meta[x_dist$NearestName,"d"]
# x_dist$Org1 <- tip_meta[x_dist$Name,"ID"]
# x_dist$Org2 <- tip_meta[x_dist$NearestName,"ID"]
# setwd("~/eggo/Data")
# save(x_dist,file="nearest_neighbors.RData")

setwd("~/eggo/Data")
load("neighbors.RData")


#Sample 10000 tips
n <- 10000
x_subsub <- drop.tip(x_sub, sample(1:Ntip(x_sub), Ntip(x_sub) - n, replace = FALSE))

#Distances
x_dist_matrix <- cophenetic.phylo(x_subsub)
x_dist_matrix[upper.tri(x_dist_matrix)] <- NA
x_dist <- melt(x_dist_matrix, 
               varnames = c('Tip1', 'Tip2'), 
               na.rm = TRUE) %>%
  mutate(Tip1 = as.character(Tip1)) %>%
  mutate(Tip2 = as.character(Tip2))
x_dist$d1 <- tip_meta[x_dist$Tip1,"d"]
x_dist$d2 <- tip_meta[x_dist$Tip2,"d"]
x_dist$Org1 <- tip_meta[x_dist$Tip1,"ID"]
x_dist$Org2 <- tip_meta[x_dist$Tip2,"ID"]
x_dist$dErr <- (log10(x_dist$d1)-log10(x_dist$d2))^2
x_dist$dErr2 <- abs((x_dist$d1)-(x_dist$d2))


p1 <- ggplot(neighbor_pred,aes(x=Actual,y=Predicted)) + 
  geom_hex(aes(fill = stat(log10(count)),alpha=stat(log10(count))), bins=100)+ 
  geom_abline(slope = 1, intercept = 0, lty = 2, alpha=0.5) +
  scale_x_log10() + scale_y_log10() + theme_bw() + 
  xlab("Doubling Time (d)") + 
  ylab(expression("Predicted Doubling Time (" * d["predicted"] * ")")) + 
  theme(legend.position = "bottom")

p2 <- ggplot(x_dist,aes(x = value + 1e-8, y = dErr2)) + 
  geom_hex(aes(fill = stat(log10(count)),alpha=stat(log10(count)))) + 
  scale_y_log10(limit=c(1e-1,1e2)) +scale_x_log10() + geom_smooth(color="black") + theme_bw() + 
  geom_hline(yintercept = 0.5, lty = 2, color = "red") + 
  geom_vline(xintercept = 0.2, lty = 2, color = "red") + 
  ylab(expression("Absolute Error (|" * d[1] - d[2] * "|)")) + 
  xlab("Patristic Distance") + 
  theme(legend.position = "bottom")

setwd("~/eggo/Figs")
pdf("16S_Neighbors_and_Pairs.pdf",width=14,height=5)
ggarrange(p1,p2,ncol=2,labels=c("(a)","(b)"),hjust = 0)
dev.off()

setwd("~/eggo/Figs")
pdf("16S_PairDistances.pdf",width=7,height=5)
p2
dev.off()

setwd("~/eggo/Data")
load("nearest_neighbors.RData")

xd <- x_dist %>% subset(Org1!=Org2) %>% subset(d1<100 & d2<100)
xd$ld1 <- log10(xd$d1)
xd$ld2 <- log10(xd$d2)
xd$dErr <- (xd$ld1-xd$ld2)^2

p1 <- ggplot(xd,aes(x=d1,y=d2)) + 
  geom_hex(aes(fill = stat(log10(count)),alpha=stat(log10(count))), bins=100)+ 
  geom_abline(slope = 1, intercept = 0, lty = 2, alpha=0.5) +
  scale_x_log10() + scale_y_log10() + theme_bw() + 
  xlab("Doubling Time (d)") + 
  ylab(expression("Predicted Doubling Time (" * d["Closest Relative"] * ")")) + 
  theme(legend.position = "bottom")

setwd("~/eggo/Figs")
pdf("16S_Closest.pdf",width=7,height=5)
p1
dev.off()
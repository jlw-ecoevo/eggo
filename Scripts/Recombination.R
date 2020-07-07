
## JLW 2020 - CUB Assoc. W/ Recomb?

# Load Packages ----------------------------------------------------------------

library(dplyr)
library(data.table)
library(ggplot2)
library(ggpubr)
library(reshape2)

# Helper Functions -------------------------------------------------------------

merge.easy <- function(df1,df2,key){
  df1 <- data.table(df1,key=key)
  df2 <- data.table(df2,key=key)
  return(unique(merge(df1,df2,all.x=TRUE,by=.EACHI,allow.cartesian=TRUE)))
}

# Load/Organize Data -----------------------------------------------------------

setwd("~/eggo/Data")
load("PhiRecombTest.RData") # Generated in Weissman et al. 2019, PLoS Genet.
load("ATGCCUB.RData")
phi_df$ATGC.COG <- paste(phi_df$ATGC,phi_df$COG,sep=".")
phi_cub <- merge.easy(phi_df,atgc_cub,key="ATGC.COG")

ribo_prot <- readLines("ribosomal_protein_COGs.txt")
phi_cub$Ribo <- phi_cub$COG %in% ribo_prot

phi_cub$R <- phi_cub$PvalCorrectedBH<5e-2

cub_RvsNR <- phi_cub %>% group_by(ATGC,R) %>% summarise(MILC=mean(MILC)) %>% 
  as.data.frame(stringsAsFactors=F)
MILC_RvsNR <- reshape(cub_RvsNR,idvar = "ATGC",timevar = "R",direction = "wide")

# Plot -------------------------------------------------------------------------

p1 <- ggplot(MILC_RvsNR,aes(x=MILC.FALSE,y=MILC.TRUE)) + geom_point() + 
  geom_abline(slope=1,intercept = 0,lty=2) + theme_bw() + 
  xlim(0.45,0.75) + ylim(0.45,0.75) +
  xlab("CUB Non-Recombining Genes") + ylab("CUB Recombining Genes")
p2 <- ggplot(MILC_RvsNR,aes(x=1,y=MILC.TRUE-MILC.FALSE)) + 
  geom_violin(fill="gray") + geom_boxplot(width=0.1) +
  theme_bw() + geom_hline(yintercept = 0,lty=2) + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  xlab("") + ylab("(CUB Recombining Genes)-(CUB Non-Recombining Genes)")
setwd("~/eggo/Figs")
pdf("CUB_Recomb_ByGene.pdf",width=7,height=5)
ggarrange(p1,p2,ncol=2,widths = c(3,1))
dev.off()

# Just Ribosomal Proteins ------------------------------------------------------

table(ribo=phi_cub$Ribo,recomb=phi_cub$R) %>% chisq.test()

ribo_cub <- phi_cub %>% subset(Ribo==T)
cub_RvsNR <- ribo_cub %>% group_by(ATGC,R) %>% summarise(MILC=mean(MILC)) %>% as.data.frame(stringsAsFactors=F)
MILC_RvsNR <- reshape(cub_RvsNR,idvar = "ATGC",timevar = "R",direction = "wide")
p1 <- ggplot(MILC_RvsNR,aes(x=MILC.FALSE,y=MILC.TRUE)) + geom_point() + 
  geom_abline(slope=1,intercept = 0,lty=2) + theme_bw() + 
  xlim(0.45,0.75) + ylim(0.45,0.75) +
  xlab("CUB Non-Recombining Genes") + ylab("CUB Recombining Genes") +
  ggtitle("Ribosomal Proteins Only")
p2 <- ggplot(MILC_RvsNR,aes(x=1,y=MILC.TRUE-MILC.FALSE)) + 
  geom_violin(fill="gray") + geom_boxplot(width=0.1) +
  theme_bw() + geom_hline(yintercept = 0,lty=2) + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  xlab("") + ylab("(CUB Recombining Genes)-(CUB Non-Recombining Genes)")
setwd("~/eggo/Figs")
pdf("CUB_Recomb_ByGene_ribosomal.pdf",width=7,height=5)
ggarrange(p1,p2,ncol=2,widths = c(3,1))
dev.off()

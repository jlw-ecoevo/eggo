


## JLW 2020 - Do big cells grow fast? 

# Load Packages ----------------------------------------------------------------

library(dplyr)
library(data.table)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)

# Helper Function --------------------------------------------------------------

merge.easy <- function(df1,df2,key){
  df1 <- data.table(df1,key=key)
  df2 <- data.table(df2,key=key)
  return(unique(merge(df1,df2,all.x=TRUE,by=.EACHI,allow.cartesian=TRUE)))
}

# Load Data --------------------------------------------------------------------

setwd("~/eggo/Data/eggo-data")
load("CodonStatistics_GORG.RData")
sag_df$CUBgenomic <- sag_df$CUBgenomic %>% lapply(median) %>% unlist()
sag_df_all <- sag_df %>% mutate_all(unlist)
sag_df <- sag_df_all %>% subset(nHE>=10)

setwd("~/eggo/Data/")
load("GORGsize.RData")

sag_df <- merge.easy(sag_df,gorg_sizes,key="ID")
sdf <- sag_df %>% subset(!is.na(CellDiameter)) %>% 
  subset(CellDiameter!="<0.5")%>% 
  mutate(CellDiameter=as.numeric(CellDiameter))
sdf <- sdf %>% mutate(CellDiameter=as.numeric(CellDiameter))

# Stats ------------------------------------------------------------------------

table(sdf$CellDiameter < quantile((sdf$CellDiameter),c(0.99)),
      sdf$d > quantile((sdf$d),c(0.01))) %>% 
  fisher.test()

# Plot -------------------------------------------------------------------------

p1 <- ggplot(sdf,aes(y=CellDiameter,x=d)) + 
  geom_point(alpha=0.5) + scale_y_log10() + scale_x_log10() + geom_smooth() + 
  theme_bw() + xlab(expression("Cell Diameter (" * mu * "m)")) +
  xlab("Predicted Minimal Doubling Time (Hours)") + 
  geom_hline(yintercept=quantile((sdf$CellDiameter),c(0.99)),lty=2,col="red") +
  geom_vline(xintercept=quantile((sdf$d),c(0.01)),lty=2,col="red")

p2 <- ggplot(sdf,aes(x=d,fill=(CellDiameter>quantile((sdf$CellDiameter),c(0.99))))) + 
  geom_density(alpha=0.5) + theme_bw() +
  scale_x_log10() + labs(fill="Top 1% Cell Size?") + 
  xlab("Predicted Minimal Doubling Time (Hours)") +
  theme(legend.position = "top") +
  scale_fill_manual(values = brewer.pal(4,"Set1")[1:2])

p3 <- ggplot(sdf %>% subset(CellDiameter>0),
       aes(x=CellDiameter,fill=d<quantile((sdf$d),c(0.01)))) + 
  geom_density(alpha=0.5) + theme_bw() + 
  labs(fill="Top 1% Growth Rate?") + 
  theme(legend.position = "top") +
  xlab(expression("Cell Diameter (" * mu * "m)")) +
  scale_fill_manual(values = brewer.pal(4,"Set1")[4:3])

setwd("~/eggo/Figs")
pdf("GORGSize.pdf",width=7,height=10)
ggarrange(p1,
          ggarrange(p2,
                    p3,
                    ncol=2,
                    labels=c("(b)","(c)")),
          nrow=2,
          labels=c("(a)",""))
dev.off()


# nHE Relationship -------------------------------------------------------------

setwd("~/eggo/Figs")
pdf("GORG_cutoff.pdf", width = 7, height = 5)
ggplot(sag_df_all,aes(x = nHE, y = d)) + geom_point() + 
  scale_y_log10() + geom_smooth() + theme_bw() + 
  geom_vline(xintercept = 10, color = "red", lty = 5) + 
  ylab("Predicted Minimal Doubling Time (Hours)") + 
  xlab("Number of Ribosomal Proteins")
dev.off()
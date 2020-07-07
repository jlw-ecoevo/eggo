
## JLW 2020 - Growth/Error Assoc. W/ Ne?

# Load Packages ----------------------------------------------------------------

library(dplyr)
library(data.table)
library(ggplot2)
library(ggpubr)

# Helper Functions -------------------------------------------------------------

merge.easy <- function(df1,df2,key){
  df1 <- data.table(df1,key=key)
  df2 <- data.table(df2,key=key)
  return(unique(merge(df1,df2,all.x=TRUE,by=.EACHI,allow.cartesian=TRUE)))
}

# Load Data --------------------------------------------------------------------

setwd("~/eggo/Data")
load("stat_data.RData")
ne <- read.csv("Ne_Bobay2018.csv", head = T, stringsAsFactors = F)
ne$Species[lapply(ne$Species,grepl,x=stat_data$Species) %>% 
             lapply(sum) %>% unlist() %>% as.logical()] <- 
  stat_data$Species[lapply(ne$Species,grep,x=stat_data$Species) %>% unlist()]
stat_data <- merge.easy(stat_data, ne, key = "Species")

# Stats ------------------------------------------------------------------------

summary(lm(d ~ Ne, data=stat_data))
summary(lm(Residual ~ Ne, data=stat_data))

# Nice Plots -------------------------------------------------------------------

p1 <- ggplot(stat_data,aes(x=Ne,y=d)) + geom_point() + 
  geom_smooth(color="darkgray") + 
  ylab("Doubling Time (Hours)") + scale_y_log10() +
  theme_bw()
p2 <- ggplot(stat_data,aes(x=Ne,y=Residual)) + geom_point() + 
  geom_smooth(color="darkgray") + 
  ylab("Residual (gRodon)") + 
  geom_hline(yintercept = 0, lty =2) + 
  theme_bw()
setwd("~/eggo/Figs")
pdf("Ne.pdf",width=7,height=3.5)
ggarrange(p1,
          p2,
          nrow=1,
          labels = c("(a)","(b)"),
          hjust=0)
dev.off()


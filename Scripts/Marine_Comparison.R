


## JLW 2020 - EGGO Marine SAGs, MAGs, and Isolate Genomes

# Load Packages ----------------------------------------------------------------


library(dplyr)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)

# Load Data --------------------------------------------------------------------

setwd("~/eggo/Data/eggo-data")
load("CodonStatistics_GORG.RData")
load("CodonStatistics_MarRef.RData")
isolate_df <- assembly_df
load("CodonStatistics_Delmont.RData")
mag_df_delmont <- mag_df
load("CodonStatistics_Tully.RData")
sag_df$CUBgenomic <- sag_df$CUBgenomic %>% lapply(median) %>% unlist()
sag_df <- sag_df %>% mutate_all(unlist)
sag_df <- sag_df %>% subset(nHE>=10)
mag_df <- mag_df %>% subset(!is.na(d))
mag_df <- mag_df %>% subset(nHE>=10)
mag_df_delmont <- mag_df_delmont %>% subset(nHE>=10)

# Plot -------------------------------------------------------------------------

p1 <- ggplot(NULL,aes(x=d)) + 
  geom_density(data=isolate_df,aes(fill="MarRef Genomes"),
               alpha=0.75,color="black") + 
  geom_density(data=mag_df,aes(fill="Tully et al. MAGs"),
               alpha=0.75,color="black") + 
  geom_density(data=mag_df_delmont,aes(fill="Delmont et al. MAGs"),
               alpha=0.75,color="black") + 
  geom_density(data=sag_df,aes(fill="GORG-tropics SAGs"),
               alpha=0.75,color="black") + 
  scale_x_log10(limit=c(1,50)) + theme_pubclean() +
  scale_fill_brewer(palette = "Dark2") + 
  theme(legend.title = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(size=8),
        legend.key.size = unit(1.3, "lines")) + 
  geom_vline(xintercept = 5, lty = 2, color = "red") + 
  xlab("Predicted Minimal Doubling Time (Hours)") 

p1l <- ggplot(NULL,aes(x=d)) + 
  geom_density(data=isolate_df,aes(fill="MarRef Genomes"),
               alpha=0.5,color="white") + 
  geom_density(data=mag_df,aes(fill="Tully et al. MAGs"),
               alpha=0.5,color="white") + 
  geom_density(data=mag_df_delmont,aes(fill="Delmont et al. MAGs"),
               alpha=0.5,color="white") + 
  geom_density(data=sag_df,aes(fill="GORG-tropics SAGs"),
               alpha=0.5,color="white") + 
  theme_bw() + xlim(0,30) +
  scale_fill_brewer(palette = "Dark2") + 
  theme(legend.title = element_blank(),legend.position = "bottom") + 
  geom_vline(xintercept = 5, lty = 2, color = "red") + 
  xlab("Predicted Minimal Doubling Time (Hours)") 

x <- data.frame(Fast=rep(c(F,T),4),
           Data=c("Isolate","Isolate",
                  "MAG (Tully)","MAG (Tully)",
                  "MAG (Delmont)","MAG (Delmont)",
                  "SAG","SAG"),
           Value=c(table(isolate_df$d<5)/nrow(isolate_df),
                   table(mag_df$d<5)/nrow(mag_df),
                   table(mag_df_delmont$d<5)/nrow(mag_df_delmont),
                   table(sag_df$d<5)/nrow(sag_df)))
p2 <- ggplot(x%>%subset(Fast==T),aes(x=Data,y=Value,fill=Data)) + 
  geom_bar(position="stack", stat="identity",alpha=0.75, color = "black") +
  theme_pubclean() + xlab("") + 
  ylab(expression("Proportion Doubling Times  <5 Hours")) +
  # theme(axis.text.x = element_text(angle = 90,hjust=1),legend.position = "none") +
  # scale_fill_brewer(palette = "Dark2",direction = -1)+
  theme(legend.position = "none") +
  scale_fill_manual(values = brewer.pal(4,"Dark2")[c(3,1,4,2)]) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

setwd("~/eggo/Figs")
pdf("MarineComparison.pdf",width=8,height=5)
ggarrange(p1,p2,ncol=2,widths=c(4,1),
          labels=c("(a)","(b)"),
          hjust=0,
          vjust=1)
dev.off()

setwd("~/eggo/Figs")
pdf("MarineComparison_nolog.pdf",width=7,height=5)
ggarrange(p1l,p2,ncol=2,widths=c(5,1))
dev.off()

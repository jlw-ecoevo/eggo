


## JLW 2020 - EGGO Gut MAGs and Isolate Genomes

# Load Packages ----------------------------------------------------------------


library(dplyr)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)

# Load Data --------------------------------------------------------------------


setwd("~/eggo/Data/eggo-data")
load("CodonStatistics_Poyet.RData")
isolate_df1 <- assembly_df
load("CodonStatistics_Zou.RData")
isolate_df2 <- assembly_df
load("CodonStatistics_Pasolli.RData")
mag_df <- mag_df %>% subset(nHE>=10)
all_mag <- mag_df
mag_df <- mag_df %>% subset(!is.na(d)) %>% subset(Body.Site=="Stool")

# Plot Comparison --------------------------------------------------------------

p1 <- ggplot(NULL,aes(x=d)) + 
  geom_density(data=isolate_df1,aes(fill="Poyet et al Isolates"),
               alpha=0.5,color="white") + 
  geom_density(data=isolate_df2,aes(fill="Zou et al Isolates"),
               alpha=0.5,color="white") + 
  geom_density(data=mag_df,aes(fill="Pasolli et al MAGs"),
               alpha=0.5,color="white") + 
  scale_x_log10(limit=c(0.09,15)) + theme_bw() +
  scale_fill_brewer(palette = "Dark2") + 
  theme(legend.title = element_blank(),legend.position = "bottom") + 
  geom_vline(xintercept = 5, lty = 2, color = "red") + 
  xlab("Predicted Minimal Doubling Time (Hours)") 

x <- data.frame(Fast=rep(c(F,T),3),
           Data=c("Isolate (Poyet)","Isolate (Poyet)","Isolate (Zou)","Isolate (Zou)","MAG","MAG"),
           Value=c(table(isolate_df1$d>5)/nrow(isolate_df1),
                   table(isolate_df2$d>5)/nrow(isolate_df2),
                   table(mag_df$d>5)/nrow(mag_df)))
p2 <- ggplot(x%>%subset(Fast==T),aes(x=Data,y=Value,fill=Data)) + 
  geom_bar(position="stack", stat="identity",alpha=0.75) +
  theme_bw() + xlab("") + 
  ylab(expression("Proportion Doubling Times  >5 Hours")) +
  # theme(axis.text.x = element_text(angle = 90,hjust=1),legend.position = "none") +
  # scale_fill_brewer(palette = "Dark2",direction = -1)+
  theme(legend.position = "none") +
  scale_fill_manual(values = brewer.pal(3,"Dark2")[c(2,3,1)]) + 
  ggpubr::rotate_x_text()

# Comparison of microbes assoc. w/ non-westernized microbiomes -----------------

mag_df$NonWesternized <- "All"
mag_df$NonWesternized[mag_df$Non.westernized.enriched == "Yes"] <- 
  "Non-Westernized"
table(mag_df$NonWesternized)


t.test(log10(d)~NonWesternized,data=mag_df)

p3 <- ggplot(mag_df,aes(y=d,x=NonWesternized,group=NonWesternized)) + 
  geom_violin(fill="gray") + geom_boxplot(width=0.5) + 
  scale_y_log10() + theme_bw() + 
  ylab("Predicted Minimal Doubling Time (Hours)") +
  xlab("") + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1))

# Put it together --------------------------------------------------------------

setwd("~/eggo/Figs")
pdf("GutComparison.pdf",width=7,height=4)
ggarrange(p1,p2,p3,ncol=3,widths=c(4,1,1),
          labels=c("(a)","(b)","(c)"),
          hjust=0,
          vjust=1)
dev.off()

# Plot by body site ------------------------------------------------------------

table(all_mag$Body.Site)


bs_mag <- all_mag %>% subset(Body.Site %in% c("Stool",
                                              "Skin"))
p1 <- ggplot(bs_mag,aes(x=d,fill=Body.Site)) + 
  geom_density(alpha=0.5) + scale_x_log10() + theme_bw() +
  scale_fill_manual(values = brewer.pal(5,"Set1")[1:2]) +
  labs(fill="") + theme(legend.position = "bottom") + 
  xlab("Predicted Minimal Doubling Time (Hours)") + 
  geom_vline(xintercept = 5, color="red", lty=2)

bs_mag <- all_mag %>% subset(Body.Site %in% c("Oral cavity",
                                              "Airways",
                                              "Vagina"))
p2 <- ggplot(bs_mag,aes(x=d,fill=Body.Site)) + 
  geom_density(alpha=0.5) + scale_x_log10() + theme_bw() +
  scale_fill_manual(values = brewer.pal(5,"Set1")[3:5]) +
  labs(fill="") + theme(legend.position = "bottom") + 
  xlab("Predicted Minimal Doubling Time (Hours)") + 
  geom_vline(xintercept = 5, color="red", lty=2)

setwd("~/eggo/Figs")
pdf("GutComparison_BodySite.pdf",width=7,height=4)
ggarrange(p1,p2,nrow=1,labels = c("(a)","(b)"),
          hjust=0,vjust=1)
dev.off()


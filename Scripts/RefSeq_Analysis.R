

## JLW 2020 - EGGO RefSeq Analysis

# Load Packages ----------------------------------------------------------------


library(dplyr)
library(ggplot2)
library(mclust)
library(ggpubr)
library(RColorBrewer)

set.seed("123456")

# Load Data --------------------------------------------------------------------

setwd("~/eggo/Data/eggo-data")
load("CodonStatistics_Parks.RData")
mag_df <- mag_df %>% subset(nHE>=10)
load("CodonStatistics_RefSeq.RData")

#Remove Genomes w/ Unrealistic #'s of Ribosomal Proteins
quantile(assembly_df$nHE,c(0.005,0.995))
assembly_df <- assembly_df %>% subset(nHE>50 & nHE<70)

# Plot Overall Distribution ----------------------------------------------------

setwd("~/eggo/Figs")
pdf("RefSeq_Growth_All.pdf",width=7,height=5)
ggplot(assembly_df,aes(x=d)) + geom_density() + 
  scale_x_log10(limit=c(0.05,100)) + theme_bw() +
  geom_vline(xintercept = 5, color = "red", lty = 2) + 
  xlab("Predicted Minimal Doubling Time (Hours)")
dev.off()

# Subsample by genus -----------------------------------------------------------


# genus_df <- assembly_df %>% group_by(Genus) %>% sample_n(1)
genus_df <- assembly_df %>% group_by(Genus,Phylum) %>% summarise_all(mean)

# Cluster growth rates ---------------------------------------------------------

x <- log10(genus_df$d)
m <- Mclust(x)
10^m$parameters$mean


# Plot Clusters ----------------------------------------------------------------

p <- (colSums(m$z)/sum(m$z))

cl_df <- data.frame(x=10^seq(-2,3,0.01),
                    cl1=dnorm(seq(-2,3,0.01),
                              mean=m$parameters$mean[1],
                              sd=sqrt(m$p$variance$sigmasq[1]))*p[1],
                    cl2=dnorm(seq(-2,3,0.01),
                              mean=m$parameters$mean[2],
                              sd=sqrt(m$p$variance$sigmasq[2]))*p[2])
p1 <- ggplot() + geom_density(data=genus_df,aes(x=d),lwd=2) + 
  scale_x_log10(limit=c(0.05,100)) +
  # geom_line(data=cl_df,aes(x=x,y=cl1,color="Cluster 1"),lwd=1) + 
  geom_polygon(data=cl_df,aes(x=x, y=cl1, fill="Cluster 1"),alpha=0.75, color="black") +
  # geom_line(data=cl_df,aes(x=x,y=cl2,color="Cluster 2"),lwd=1) + 
  geom_polygon(data=cl_df,aes(x=x, y=cl2, fill="Cluster 2"), alpha=0.75, color="black") +
  theme_pubclean() + scale_color_brewer(palette = "Set1")  +
  geom_vline(xintercept = 5, color = "red", lty = 2) + 
  xlab("Predicted Minimal Doubling Time (Hours)")  + 
  theme(legend.position = "right") + 
  labs(fill="")


setwd("~/eggo/Figs")
pdf("RefSeq_vs_Parks.pdf",width=10,height=5)
ggplot() + 
  geom_density(data=genus_df,aes(x=d,fill="RefSeq Assemblies"),alpha=0.5,color="white") + 
  geom_density(data=mag_df,aes(x=d,fill="Parks et al. MAGs"),alpha=0.5,color="white") +
  scale_x_log10(limit=c(0.05,100)) +
  theme_bw() + scale_color_brewer(palette = "Set1")  +
  geom_vline(xintercept = 5, color = "red", lty = 2) + 
  xlab("Predicted Minimal Doubling Time (Hours)")  + 
  labs(fill="")
dev.off()


# Plot by phylum ---------------------------------------------------------------

table(genus_df$Phylum) %>% sort()
pdf1_genus <- genus_df %>% subset(Phylum %in% c("Proteobacteria",
                                                "Bacteroidetes",
                                               "Actinobacteria",
                                               "Firmicutes" ))
pdf2_genus <- genus_df %>% subset(Phylum %in% c("Chloroflexi",
                                               "Planctomycetes",
                                               "Cyanobacteria"))
p2 <- ggplot(pdf1_genus,aes(x=d,fill=Phylum)) + 
  geom_density(alpha=0.6) + scale_x_log10(limit=c(0.05,100))  +
  geom_vline(xintercept = 5, color = "red", lty = 2) + 
  theme_pubclean() + theme(legend.position = "bottom",
                     legend.text = element_text(size=7),
                     legend.key.size = unit(0.7, "lines")) + 
  labs(fill="") + 
  xlab("Predicted Minimal Doubling Time (Hours)") + 
  scale_fill_manual(values = brewer.pal(7,"Set1")[1:4])

p3 <- ggplot(pdf2_genus,aes(x=d,fill=Phylum)) + 
  geom_density(alpha=0.6) + scale_x_log10(limit=c(0.05,100))  +
  geom_vline(xintercept = 5, color = "red", lty = 2) + 
  theme_pubclean() + theme(legend.position = "bottom",
                     legend.text = element_text(size=7),
                     legend.key.size = unit(0.8, "lines")) + 
  labs(fill="") + 
  xlab("Predicted Minimal Doubling Time (Hours)") + 
  scale_fill_manual(values = brewer.pal(7,"Set1")[5:7])


# Put it together --------------------------------------------------------------

setwd("~/eggo/Figs")
pdf("RefSeq_Growth.pdf",width=8,height=6)
ggarrange(p1,
          ggarrange(p2,
                    p3,
                    nrow=1,
                    labels=c("(b)","(c)"),
                    hjust=0,
                    vjust=1),
          nrow=2,
          labels=c("(a)",""),
          hjust=0,
          vjust=1)
dev.off()


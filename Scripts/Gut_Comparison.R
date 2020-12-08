


## JLW 2020 - EGGO Gut MAGs and Isolate Genomes

# Load Packages ----------------------------------------------------------------


library(dplyr)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(data.table)
library(phytools)

merge.easy <- function(df1,df2,key){
  df1 <- data.table(df1,key=key)
  df2 <- data.table(df2,key=key)
  return(as.data.frame(unique(merge(df1,df2,all.x=TRUE,by=.EACHI,allow.cartesian=TRUE)),stringsAsFactors=F))
}
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

#Taxonomy
setwd("~/eggo/Data")
ar <- read.delim("gut.gtdbtk.ar122.summary.tsv", 
                 sep = "\t", 
                 head = T, 
                 stringsAsFactors = F) %>% 
  subset(select = c(user_genome, classification))
bac <- read.delim("gut.gtdbtk.bac120.summary.tsv", 
                  sep = "\t", 
                  head = T, 
                  stringsAsFactors = F) %>% 
  subset(select = c(user_genome, classification))
tax <- rbind(ar,bac)
tax$Phylum <- tax$classification %>% 
  gsub(pattern=".*;p__",replace="") %>% 
  gsub(pattern=";.*",replace="")
tax$Class <- tax$classification %>% 
  gsub(pattern=".*;c__",replace="") %>% 
  gsub(pattern=";.*",replace="")
tax$Order <- tax$classification %>% 
  gsub(pattern=".*;o__",replace="") %>% 
  gsub(pattern=";.*",replace="")
tax$Family <- tax$classification %>% 
  gsub(pattern=".*;f__",replace="") %>% 
  gsub(pattern=";.*",replace="")

#Keep full lineage for subset
tax$Genus <- tax$classification %>% 
  gsub(pattern=";s__.*",replace="") 


# Plot Comparison --------------------------------------------------------------

p1 <- ggplot(NULL,aes(x=d)) + 
  geom_density(data=isolate_df1,aes(fill="Poyet et al Isolates"),
               alpha=0.75,color="black") + 
  geom_density(data=isolate_df2,aes(fill="Zou et al Isolates"),
               alpha=0.75,color="black") + 
  geom_density(data=mag_df,aes(fill="Pasolli et al MAGs"),
               alpha=0.75,color="black") + 
  #scale_x_log10(limit=c(0.09,15)) +
  xlim(0,15) + 
  theme_pubclean() + 
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
  geom_bar(position="stack", stat="identity",alpha=0.75, color = "black") +
  xlab("") + theme_pubclean() + # + theme_bw() + 
  ylab(expression("Proportion Oligotrophs")) +
  # theme(axis.text.x = element_text(angle = 90,hjust=1),legend.position = "none") +
  # scale_fill_brewer(palette = "Dark2",direction = -1)+
  theme(legend.position = "none") +
  scale_fill_manual(values = brewer.pal(3,"Dark2")[c(2,3,1)])+ 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.title.y = element_text(size=9))

# Comparison of microbes assoc. w/ non-westernized microbiomes -----------------

mag_df$NonWesternized <- "All"
mag_df$NonWesternized[mag_df$Non.westernized.enriched == "Yes"] <- 
  "Non-Westernized"
table(mag_df$NonWesternized)


t.test(log10(d)~NonWesternized,data=mag_df)

p3 <- ggplot(mag_df,aes(y=d,x=NonWesternized,group=NonWesternized)) + 
  geom_violin(fill="gray",alpha=0.75) + geom_boxplot(width=0.3) + 
  scale_y_log10() + theme_pubclean() +#theme_bw() + 
  ylab("Predicted Minimal Doubling Time (Hours)") +
  xlab("") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

# Put it together --------------------------------------------------------------


setwd("~/eggo/Figs")
pdf("GutComparison.pdf",width=8,height=5)
ggarrange(p1,
          p2,
          p3,
          ncol=3,widths=c(4,1,1),
          labels=c("(a)","(b)","(c)"),
          hjust=0,
          vjust=1,
          align="hv")
dev.off()

# Plot by body site ------------------------------------------------------------

table(all_mag$Body.Site)


bs_mag <- all_mag %>% subset(Body.Site %in% c("Stool",
                                              "Skin"))
p1b <- ggplot(bs_mag,aes(x=d,fill=Body.Site)) + 
  geom_density(alpha=0.5) + scale_x_log10() + theme_bw() +
  scale_fill_manual(values = brewer.pal(5,"Set1")[1:2]) +
  labs(fill="") + theme(legend.position = "bottom") + 
  xlab("Predicted Minimal Doubling Time (Hours)") + 
  geom_vline(xintercept = 5, color="red", lty=2)

bs_mag <- all_mag %>% subset(Body.Site %in% c("Oral cavity",
                                              "Airways",
                                              "Vagina"))
p2b <- ggplot(bs_mag,aes(x=d,fill=Body.Site)) + 
  geom_density(alpha=0.5) + scale_x_log10() + theme_bw() +
  scale_fill_manual(values = brewer.pal(5,"Set1")[3:5]) +
  labs(fill="") + theme(legend.position = "bottom") + 
  xlab("Predicted Minimal Doubling Time (Hours)") + 
  geom_vline(xintercept = 5, color="red", lty=2)

setwd("~/eggo/Figs")
pdf("GutComparison_BodySite.pdf",width=7,height=4)
ggarrange(p1b,p2b,nrow=1,labels = c("(a)","(b)"),
          hjust=0,vjust=1)
dev.off()

# By Taxonomy ------------------------------------------------------------------

mag_df <- mag_df %>% subset(select=names(isolate_df1))
isolate_df1$user_genome <- paste0(isolate_df1$Assembly,"_genomic.fna")
isolate_df2$user_genome <- paste0(isolate_df2$Assembly,"_genomic.fna")
mag_df$user_genome <- paste0(mag_df$Assembly,".fa")

isolate_df1$Source <- "Poyet et al."
isolate_df2$Source <- "Zou et al."
mag_df$Source <- "Pasolli et al."

genomes <- rbind(isolate_df1,
                 isolate_df2,
                 mag_df)

genomes <- merge.easy(genomes,tax,key="user_genome")

library(ape)
library(ggtree)
library(phylolm)
library(sensiPhy)
library(dispRity)

tree <- read.tree("/home/jake/eggo/Data/gut.gtdbtk.bac120.classify.tree")
genomes <- genomes %>% subset(d != "NaN")
tree.bac <- drop.tip(tree,which(!tree$tip.label %in% genomes$user_genome))
meta <- genomes %>% subset(user_genome %in% tree.bac$tip.label)

#Test assoc w/ phylogenetic logistic regression
rownames(meta) <- meta$user_genome
meta$Isolate <- (meta$Source != "Pasolli et al.")*1
tree.test <- remove.zero.brlen(tree.bac)
isolate_mod <- phyloglm(Isolate~(d<5),
                        data=meta,
                        phy=tree.test,
                        btol=20,
                        log.alpha.bound = 5,
                        method="logistic_MPLE")
summary(isolate_mod)
table(meta$d<5)/nrow(meta)
table(meta$d<5,meta$Isolate)

tree.fam <- drop.tip(tree.bac, which(!tree.bac$tip.label %in% meta$user_genome))
d <- data.frame(id=meta$user_genome,Phylum=meta$Phylum)
d$Phylum <- d$Phylum %>% gsub(pattern="_.*",replace="")
d$Phylum[d$Phylum %in% names(table(d$Phylum)[table(d$Phylum)<25])] <- "Other"
p <- ggtree(tree.fam,layout="rectangular",size=0.001)
nColor <- 12
p <- p %<+% d + geom_tippoint(aes(color=Phylum),alpha=0.75,size=1) + 
  scale_color_manual(values=carto_pal(nColor, "Safe"))

heat_data <- meta %>% subset(select=c(d,Isolate)) %>% as.data.frame()
rownames(heat_data) <- meta$user_genome
heat_data$d <- factor(heat_data$d<5)
heat_data$Isolate <- factor(heat_data$Isolate==1)

colnames(heat_data) <- c("Copiotroph","Isolate")
ph <- gheatmap(p+ guides(colour = guide_legend(override.aes = list(size=4,alpha=1))),heat_data,
               low="lightgray",
               high=brewer.pal(4,"Dark2")[3],
               width=0.5,
               offset=-0.25,
               colnames=T,
               color=NA,
               font.size=3,
               hjust=0,
               colnames_angle = 20,
               colnames_offset_y=80,
               colnames_position = "top") + 
  scale_fill_manual(values=c("lightgray",brewer.pal(4,"Dark2")[3])) + 
  theme_tree() + 
  labs(fill="") + 
  theme(legend.text=element_text(size=8), 
        legend.key.height=unit(.5, "cm"),
        legend.key.width=unit(.4, "cm"))
ph

setwd("~/eggo/Figs")
pdf("growth_vs_isolate_phy_gut.pdf",height=8,width=8)
ggarrange(ggarrange(p1,p2,ncol=2,widths=c(5,2),
                    labels=c("(a)","(b)"),
                    hjust=0,
                    vjust=1),
          ggarrange(ph,
                    p3,
                    ncol=2,
                    widths=c(3,1),
                    labels=c("(c)","(d)"),
                    hjust=0),
          nrow = 2,
          heights = c(2,4))
dev.off()

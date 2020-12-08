


## JLW 2020 - EGGO Marine SAGs, MAGs, and Isolate Genomes

# Load Packages ----------------------------------------------------------------


library(dplyr)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(data.table)
library(phytools)
library(rcartocolor)

merge.easy <- function(df1,df2,key){
  df1 <- data.table(df1,key=key)
  df2 <- data.table(df2,key=key)
  return(as.data.frame(unique(merge(df1,df2,all.x=TRUE,by=.EACHI,allow.cartesian=TRUE)),stringsAsFactors=F))
}

give.n <- function(x){
  return(c(y = median(x)*1.05, label = length(x))) 
  # experiment with the multiplier to find the perfect position
}

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

#Taxonomy
setwd("~/eggo/Data")
ar <- read.delim("gtdbtk.ar122.summary.tsv", 
                 sep = "\t", 
                 head = T, 
                 stringsAsFactors = F) %>% 
  subset(select = c(user_genome, classification))
bac <- read.delim("gtdbtk.bac120.summary.tsv", 
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
  scale_x_log10(limit=c(1,30)) + theme_pubclean() +
  scale_fill_brewer(palette = "Dark2") + 
  theme(legend.title = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(size=8),
        legend.key.size = unit(1.3, "lines")) + 
  geom_vline(xintercept = 5, lty = 2, color = "red") + 
  xlab("Predicted Minimal Doubling Time (Hours)") 

p1l <- ggplot(NULL,aes(x=d)) + 
  geom_density(data=isolate_df,aes(fill="MarRef Genomes"),
               alpha=0.7,color="black") + 
  geom_density(data=mag_df,aes(fill="Tully et al. MAGs"),
               alpha=0.7,color="black") + 
  geom_density(data=mag_df_delmont,aes(fill="Delmont et al. MAGs"),
               alpha=0.7,color="black") + 
  geom_density(data=sag_df,aes(fill="GORG-tropics SAGs"),
               alpha=0.7,color="black") + 
  theme_pubclean() +
  xlim(0,25) +
  scale_fill_brewer(palette = "Dark2") + 
  theme(legend.title = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(size=8),
        legend.key.size = unit(0.8, "lines")) + 
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
  ylab(expression("Proportion Copiotrophs")) +
  # theme(axis.text.x = element_text(angle = 90,hjust=1),legend.position = "none") +
  # scale_fill_brewer(palette = "Dark2",direction = -1)+
  theme(legend.position = "none") +
  scale_fill_manual(values = brewer.pal(4,"Dark2")[c(3,1,4,2)]) + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1,size=7),
        axis.title.y = element_text(size=9))

setwd("~/eggo/Figs")
pdf("MarineComparison.pdf",width=8,height=5)
ggarrange(p1,p2,ncol=2,widths=c(4,1),
          labels=c("(a)","(b)"),
          hjust=0,
          vjust=1)
dev.off()

setwd("~/eggo/Figs")
pdf("MarineComparison_density.pdf",width=6,height=6)
p1
dev.off()

setwd("~/eggo/Figs")
pdf("MarineComparison_nolog.pdf",width=7,height=5)
ggarrange(p1l,p2,ncol=2,widths=c(5,1))
dev.off()

# By Taxonomy ------------------------------------------------------------------

isolate_df$user_genome <- paste0(isolate_df$Assembly,"_cds_from_genomic.fna.fna")
mag_df$user_genome <- paste0(mag_df$Assembly,"_cds_from_genomic.fna.fna")
mag_df_delmont$user_genome <- paste0(mag_df_delmont$Assembly,".fa.fna")
sag_df$user_genome <- sag_df$ID %>%
  paste0(.,"_contigs.fasta.fna")
sag_df$ID <- NULL

isolate_df$Source <- "MarRef"
mag_df$Source <- "Tully et al."
mag_df_delmont$Source <- "Delmont et al."
sag_df$Source <- "GORG-Tropics"

genomes <- rbind(isolate_df,
                 mag_df,
                 mag_df_delmont,
                 sag_df)

genomes <- merge.easy(genomes,tax,key="user_genome")

library(ape)
library(ggtree)
library(phylolm)
library(sensiPhy)
library(dispRity)

setwd("~/eggo/Data/")
tree <- read.tree("gtdbtk.bac120.classify.tree")
boot_trees <-  read.tree("RAxML_bootstrap.marine_gtdb_raxml_fewboot")

tree.bac <- drop.tip(tree,which(!tree$tip.label %in% genomes$user_genome))
subphy <- function(phy,tip_names){
  drop.tip(phy,which(!phy$tip.label %in% tip_names))
}
boot.bac <- lapply(boot_trees,subphy,tip_names=genomes$user_genome)
class(boot.bac) <- "multiPhylo"
meta <- genomes %>% subset(user_genome %in% tree.bac$tip.label)

#Test assoc w/ phylogenetic logistic regression
rownames(meta) <- meta$user_genome
meta$Isolate <- (meta$Source == "MarRef")*1
tree.test <- remove.zero.brlen(tree.bac)
isolate_mod <- phyloglm(Isolate~(d<5),
                        data=meta,
                        phy=tree.test,
                        btol=20,
                        method="logistic_MPLE")
summary(isolate_mod)


# Sensitivity of phylogenetic regression (slow)
influ_sense <- influ_phyglm(Isolate~(d<5), phy = tree.test,
                     data = meta, track=T)
samp_sense <- samp_phyglm(Isolate~(d<5), phy = tree.test,
            data = meta, track=T)
clade_sense <- clade_phyglm(Isolate~(d<5), phy = tree.test,
                     data = meta, clade.col = "Phylum", track=T)
tree_sense <- tree_phyglm(Isolate~(d<5), phy = boot.bac,
                   data = meta, n.tree = 10, track=T)
tree_influ_sense <- tree_influ_phyglm(Isolate~(d<5), phy = boot.bac,
                                      data = meta, n.tree = 10, track=T)
tree_clade_sense <- tree_clade_phyglm(Isolate~(d<5), phy = boot.bac,
                                      data = meta, n.tree = 10, track=T, clade.col = "Phylum")
tree_samp_sense <- tree_samp_phyglm(Isolate~(d<5), phy = boot.bac,
                                    data = meta, n.tree = 10, track=T)
setwd("~/eggo/Data")
save(influ_sense,
     samp_sense,
     clade_sense,
     tree_sense,
     tree_influ_sense,
     tree_clade_sense,
     tree_samp_sense,
     file="sentivitity_phyloglm.RData")

setwd("~/eggo/Data")
load("sentivitity_phyloglm.RData")

tree_influ_sense$sensi.estimates[tree_influ_sense$sensi.estimates[,"species"]=="AG-333-P16_contigs.fasta.fna",]

summary(clade_sense)

p_is <- ggplot(influ_sense$sensi.estimates,aes(x=estimate)) + 
  geom_histogram() + 
  theme_pubclean() +
  xlab("Estimated Relationship") + 
  theme(axis.title = element_text(size=9))

p_ss <- ggplot(samp_sense$sensi.estimates,aes(x=n.percent,group=n.percent,y=estimate)) + 
  geom_boxplot() + 
  theme_pubclean() +
  ylab("Estimated Relationship") + 
  xlab("Percent Data Dropped") + 
  theme(axis.title = element_text(size=9))

p_cs <- ggplot(clade_sense$sensi.estimates,aes(x=estimate)) + 
  geom_histogram() + 
  theme_pubclean() +
  xlab("Estimated Relationship") + 
  theme(axis.title = element_text(size=9))

p_tr <- ggplot(tree_sense$sensi.estimates,aes(x=estimate)) + 
  geom_histogram(binwidth=0.025) + 
  theme_pubclean() +
  xlab("Estimated Relationship") + 
  theme(axis.title = element_text(size=9)) + 
  xlim(0,1.5)

p_tr_is <- ggplot(tree_influ_sense$sensi.estimates,aes(x=estimate,fill=iteration)) + 
  geom_histogram(alpha=0.5) + 
  theme_pubclean() +
  xlab("Estimated Relationship") + 
  theme(axis.title = element_text(size=9)) +
  labs(fill="Bootstrap Replicate") + 
  xlim(0,1.5)

p_tr_cs <- ggplot(tree_clade_sense$sensi.estimates,aes(x=estimate,fill=iteration)) + 
  geom_histogram(alpha=0.5) + 
  theme_pubclean() +
  xlab("Estimated Relationship") + 
  theme(axis.title = element_text(size=9)) +
  labs(fill="Bootstrap Replicate") + 
  xlim(0,1.5)

p_tr_ss <- ggplot(tree_samp_sense$sensi.estimates,aes(x=factor(n.percent),y=estimate,fill=iteration)) + 
  geom_boxplot() + 
  theme_pubclean() +
  ylab("Estimated Relationship") + 
  xlab("Percent Data Dropped") + 
  theme(axis.title = element_text(size=9)) +
  labs(fill="Bootstrap Replicate") + 
  ylim(0,1.5)

meta_genus <- meta %>% 
  group_by(Genus) %>%
  summarise(d=mean(d),
            Isolate=max(Isolate),
            user_genome=sample(user_genome,1),
            Phylum=unique(Phylum))
tree.fam <- drop.tip(tree.bac, which(!tree.bac$tip.label %in% meta_genus$user_genome))
d <- data.frame(id=meta_genus$user_genome,Phylum=meta_genus$Phylum)
d$Phylum <- d$Phylum %>% gsub(pattern="_.*",replace="")
d$Phylum[d$Phylum %in% names(table(d$Phylum)[table(d$Phylum)<25])] <- "Other"
p <- ggtree(tree.fam,layout="rectangular",size=0.001)
nColor <- 12
p <- p %<+% d + geom_tippoint(aes(color=Phylum),alpha=0.75,size=1) + 
  scale_color_manual(values=carto_pal(nColor, "Safe"))


heat_data <- meta_genus %>% subset(select=c(d,Isolate)) %>% as.data.frame()
rownames(heat_data) <- meta_genus$user_genome
heat_data$d <- factor(heat_data$d<5)
heat_data$Isolate <- factor(heat_data$Isolate==1)

colnames(heat_data) <- c("Copiotroph","MarRef")
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
               colnames_offset_y=15,
               colnames_position = "top") + 
  # scale_fill_brewer(palette="Set1") + 
  scale_fill_manual(values=c("lightgray",brewer.pal(4,"Dark2")[3])) + 
  theme_tree() + 
  labs(fill="") + 
  theme(legend.text=element_text(size=8), 
            legend.key.height=unit(.5, "cm"),
            legend.key.width=unit(.4, "cm"))
            #legend.position=c(.13, y=.945))
ph
        

setwd("~/eggo/Figs")
pdf("growth_vs_isolate_phy.pdf",height=8,width=8)
ggarrange(ggarrange(p1l,p2,ncol=2,widths=c(5,2),
                    labels=c("(a)","(b)"),
                    hjust=0,
                    vjust=1),
          ggarrange(ph,
                    ggarrange(p_is,
                              p_cs,
                              p_ss,
                              nrow=3,
                              labels=c("(d)","(e)","(f)"),
                              hjust=0.8),
                    ncol=2,
                    widths=c(3,2),
                    labels=c("(c)",""),
                    hjust=0),
          nrow = 2,
          heights = c(2,4))
dev.off()

setwd("~/eggo/Figs")
pdf("tree_boot.pdf",height=8,width=10)
ggarrange(p_tr,
          p_tr_is,
          p_tr_cs,
          p_tr_ss,
          ncol=2,
          nrow=2,
          labels = c("(a)","(b)","(c)","(d)"))
dev.off()



setwd("~/eggo/Figs")
pdf("marine_slides.pdf",height=8,width=12)
ggarrange(ggarrange(p1l + theme(axis.title.y = element_text(size=12)),
                    p2 + theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1, size=10),
                             axis.title.y = element_text(size=12)),
                    nrow=2),
          ph,
          ncol=2,
          widths=c(1,1))
dev.off()

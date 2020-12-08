
## JLW 2020 - Functional Enrichment of Copiotrophs/Oligotrophs (CAZymes)

# Load Packages ----------------------------------------------------------------

library(RColorBrewer)
library(dplyr)
library(Hmisc)
library(reshape2)
library(ggplot2)
library(data.table)
library(ggpubr)

merge.easy <- function(df1,df2,key){
  df1 <- data.table(df1,key=key)
  df2 <- data.table(df2,key=key)
  return(as.data.frame(unique(merge(df1,df2,all.x=TRUE,by=.EACHI,allow.cartesian=TRUE)),stringsAsFactors=F))
}

# Load Data --------------------------------------------------------------------

setwd("~/eggo/Data")
tblout <- read.delim("refseq_filt.caz",sep="",head=F,stringsAsFactors = F) %>% subset(V1 != "#")
names(tblout) <- c("Gene","Hit","pval")
gene_genome <- read.delim("genes_genomes_hits.txt",sep="",head=F,stringsAsFactors = F)
names(gene_genome) <- c("Genome","Gene")
tblout <- merge.easy(tblout,gene_genome,key="Gene")
tblout$Accession <- tblout$Genome %>% substr(1,15)
proteins <- read.delim("n_proteins.tbl",head=F,stringsAsFactors = F,sep="")
names(proteins) <- c("Genomes","nProt")
proteins$Accession <- proteins$Genomes %>% substr(1,15)
annot <- read.delim("CAZyDB.07302020.fam-activities.txt",sep="\t",stringsAsFactors = F,skip=1,head=F)

# 
genomes <- readLines("genus_genomes.txt") %>% substr(1,15)
load("~/eggo/Data/EGGO.RData")
EGGO$Accession <- EGGO$Assembly %>% substr(1,15)
EGGO <- EGGO %>% subset(Accession %in% genomes)
copiotrophs <- EGGO$Accession[EGGO$d<5] 
oligotrophs <- EGGO$Accession[EGGO$d>=5] 

#Organize data
x <- tblout %>% subset(select=c(Hit,Accession)) %>% unique()
x$Copiotroph <- x$Accession %in% copiotrophs
x$Oligotroph <- x$Accession %in% oligotrophs
table(x$Copiotroph,x$Oligotroph)
x <- x %>% subset(Oligotroph==T | Copiotroph==T)
y <- x %>% group_by(Hit, Copiotroph) %>% count()
z <- reshape(as.data.frame(y),direction = "wide",idvar="Hit",timevar="Copiotroph")
z$n.FALSE[is.na(z$n.FALSE)] <- 0
z$n.TRUE[is.na(z$n.TRUE)] <- 0
z$nGenomes.FALSE <- length(oligotrophs)
z$nGenomes.TRUE <- length(copiotrophs)

# Fisher exact test ------------------------------------------------------------

for(i in 1:nrow(z)){
  contingency_table <- matrix(c(z[i,"n.TRUE"],
                                z[i,"n.FALSE"],
                                z[i,"nGenomes.TRUE"]-z[i,"n.TRUE"],
                                z[i,"nGenomes.FALSE"]-z[i,"n.FALSE"]),
                              ncol=2)
  z$pval[i] <- fisher.test(contingency_table)$p.value
  z$P.TRUE[i] <- z[i,"n.TRUE"]/z[i,"nGenomes.TRUE"]
  z$P.FALSE[i] <- z[i,"n.FALSE"]/z[i,"nGenomes.FALSE"]
}
z$pval.corrected <- p.adjust(z$pval,method="BH")

# Plot -------------------------------------------------------------------------

p1 <- ggplot(z,aes(x=P.FALSE,y=P.TRUE)) + 
  geom_point(aes(fill=(pval.corrected<0.01),shape=(pval.corrected<0.01)),alpha=0.75) + 
  geom_abline(slope=1,intercept=,lty=2) + 
  theme_pubclean() + 
  #geom_smooth(color="darkgrey") +
  #scale_fill_brewer(palette="Set3",direction=-1) + 
  xlab("Prevalence in Oligotroph Genomes") + 
  ylab("Prevalence in Copiotroph Genomes") + 
  scale_shape_manual(values=c(23,21),
                     name = expression(p<0.01),
                     labels = c("FALSE","TRUE")) +
  scale_fill_manual(values = c(brewer.pal(3,"Set3")[2],brewer.pal(3,"Set3")[1]),
                    name = expression(p<0.01),
                    labels = c("FALSE","TRUE")) +
  theme(legend.position = c(0.75,0.2),
        legend.background = element_rect(size=0.1,color="black"))

p2 <- ggplot(z,aes(y=P.TRUE/P.FALSE,x=1)) + 
  geom_violin(fill="grey") + geom_boxplot(width=0.5) + 
  scale_y_log10() + 
  geom_hline(yintercept=1,lty=2) + 
  theme_pubclean() + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y = element_text(size=10)) + 
  ylab("(Prevalence in Copiotrophs)/(Prevalence in Oligotrophs)")

#Not a number of proteins thing
boxplot(proteins$nProt[proteins$Accession %in% copiotrophs],
        proteins$nProt[proteins$Accession %in% oligotrophs])
t.test(proteins$nProt[proteins$Accession %in% copiotrophs],
       proteins$nProt[proteins$Accession %in% oligotrophs])

z$Hit <- z$Hit %>% gsub(pattern=".hmm",replace="")
names(annot) <- c("Hit", "Definition")
z <- merge.easy(z,annot,key="Hit")
z$Effect <- z$P.TRUE/z$P.FALSE

z$Class <- substr(z$Hit,1,2)
z$Class[z$Class=="CB"] <- "CBM"

z$Cat <- z$Effect>1
z_sig <- z %>% 
  subset(pval.corrected<0.01) %>% 
  group_by(Class,Cat) %>% 
  count() %>% 
  subset(Class %in% c("GH",
                      "GT",
                      "PL",
                      "CE",
                      "AA",
                      "CBM"))
z_sig$Cat[z_sig$Cat==TRUE] <- "Copiotroph"
z_sig$Cat[z_sig$Cat==FALSE] <- "Oligotroph"

z_sig$Class[z_sig$Class=="GH"] <- "Glycoside hydrolase"
z_sig$Class[z_sig$Class=="GT"] <- "Glycosyltranferase"
z_sig$Class[z_sig$Class=="PL"] <- "Polysaccharide lyase"
z_sig$Class[z_sig$Class=="CE"] <- "Carbohydrate esterase"
z_sig$Class[z_sig$Class=="AA"] <- "Auxilary activities"
z_sig$Class[z_sig$Class=="CBM"] <- "Carbohydrate-Binding Modules"

z_sig <- rbind(as.data.frame(z_sig,stringsAsFacotrs=F),
      data.frame(Class="Carbohydrate-Binding Modules",Cat="Oligotroph",n=0))

pC <- ggplot(z_sig[z_sig$Cat=="Copiotroph",], aes(x = Class, Cat)) +
  geom_tile(aes(fill = n)) +
  geom_text(aes(label = round(n, 3))) +
  # As our values are continuous, we'll use scale_fill_continuous instead of scale_fill_manual
  scale_fill_continuous(low = "gray", high = "aquamarine", trans="log") + 
  ylab("Copiotroph") + 
  theme_pubclean() +
  theme(legend.position = "none",
        axis.text.y = element_text(size = 0),
        axis.text.x = element_text(size = 0),
        axis.ticks.x = element_blank(),
        plot.margin = unit(c(.3,1,0,1), "cm"),
        axis.title.y = element_text(size=10)) +
  xlab("") + 
  scale_x_discrete(expand=c(0,0)) + 
  scale_y_discrete(expand=c(0,0))

pO <- ggplot(z_sig[z_sig$Cat=="Oligotroph",], aes(x = Class, Cat)) +
  geom_tile(aes(fill = n)) +
  geom_text(aes(label = round(n, 3))) +
  # As our values are continuous, we'll use scale_fill_continuous instead of scale_fill_manual
  scale_fill_continuous(low = "gray", high = "aquamarine", trans="log") + 
  ylab("Oligotroph") + 
  theme_pubclean() +
  theme(legend.position = "none",
        axis.text.y = element_text(size = 0),
        plot.margin = unit(c(0,1,0,1), "cm"),
        axis.text.x = element_text(angle = 30,hjust=1,size=10),
        axis.title.y = element_text(size=10)) +
  xlab("")+ 
  scale_x_discrete(expand=c(0,0)) + 
  scale_y_discrete(expand=c(0,0))


setwd("~/eggo/Figs")
pdf("CAZyme_enrichment.pdf",width=7,height=7)
ggarrange(ggarrange(p1,p2,ncol=2,widths=c(3,1), labels=c("(a)","(b)"),hjust=0),
          ggarrange(pC,pO,nrow=2,heights=c(1,2)),
          nrow=2,
          heights=c(3,2),
          labels=c("","(c)"),
          hjust=0)
dev.off()

table(substr(annot$Hit,1,2))

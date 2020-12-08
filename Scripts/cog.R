
## JLW 2020 - Functional Enrichment of Copiotrophs/Oligotrophs (eggnogmapper)

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

# annotations
tblout <- read.delim("genusgenomes.cogs",sep="",head=F,stringsAsFactors = F) #%>% subset(substr(V2,1,2)!="ar")
names(tblout) <- c("Genome","Hit")
tblout$Accession <- tblout$Genome %>% substr(1,15)

# number of proteins
proteins <- read.delim("n_proteins.tbl",head=F,stringsAsFactors = F,sep="")
names(proteins) <- c("Genomes","nProt")
proteins$Accession <- proteins$Genomes %>% substr(1,15)

#COG classes
arcog_def <- read.delim("ar14.arCOGdef.tab",sep="\t",head=T,stringsAsFactors = F) %>% 
  subset(select=c(CLU,FUNC_ONE))
cog_def <- read.delim("cog-20.def.tab",sep="",head=F,stringsAsFactors = F) %>% 
  subset(select=c(V1,V2))
names(arcog_def) <- c("Hit","FuncClass")
names(cog_def) <- c("Hit","FuncClass")
cog_def <- rbind(cog_def,arcog_def)

# functional class definitions
cog_cat_def <- read.delim("fun-20.tab",head = F, sep="\t",stringsAsFactors = F)
names(cog_cat_def) <- c("FuncClass","X","Def")

# genome list
genomes <- readLines("genus_genomes.txt") %>% substr(1,15)

# Enrichment
sig_df <- read.csv("function_enrichment_propmixed_significant.csv")
names(cog_cat_def)[1] <- "COG"
sdf <- merge.easy(sig_df,cog_cat_def,key="COG")
sdf$Name <- paste0(sdf$COG,": ",sdf$Def)

#Growth Data
load("~/eggo/Data/EGGO.RData")
EGGO$Accession <- EGGO$Assembly %>% substr(1,15)
EGGO <- EGGO %>% subset(Accession %in% genomes)
copiotrophs <- EGGO$Accession[EGGO$d<5] 
oligotrophs <- EGGO$Accession[EGGO$d>=5] 



# Prevalence Significance Test -------------------------------------------------

#organize data
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

# Fisher's exact test for significance
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

#Not a number of proteins thing
boxplot(proteins$nProt[proteins$Accession %in% copiotrophs],
        proteins$nProt[proteins$Accession %in% oligotrophs])
t.test(proteins$nProt[proteins$Accession %in% copiotrophs],
       proteins$nProt[proteins$Accession %in% oligotrophs])

# Appendix Figure (full) -------------------------------------------------------

p1a <- ggplot(z,aes(x=P.FALSE,y=P.TRUE)) + 
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
  theme(legend.position = c(0.9,0.3),
        legend.background = element_rect(size=0.1,color="black"),
        axis.title.y = element_text(size = 10))

z_sig <- merge.easy(z,cog_def,key="Hit") %>% subset(pval.corrected<0.01) %>% 
  subset(FuncClass %in% setdiff(c(LETTERS[1:22],"X"),c("R","S")))
x <- table(z_sig$FuncClass,z_sig$P.TRUE>z_sig$P.FALSE)
z_sig$Copio <- factor(z_sig$P.TRUE>z_sig$P.FALSE)
z_sig$FuncClass <- factor(z_sig$FuncClass)
z_cat <- z_sig %>% group_by(Copio,FuncClass,.drop=FALSE) %>%
  summarise(n=length(pval.corrected)) %>% as.data.frame(stringsAsFactors=F)
z_cat_diff <- z_cat %>% group_by(FuncClass) %>% summarise(n=(n[Copio==F]-n[Copio==T]))
z_cat_diff$Copio <- NA
names(cog_cat_def)[1] <- "FuncClass"
z_plot <- merge.easy(rbind(z_cat,z_cat_diff),cog_cat_def,key="FuncClass")
z_plot$Name <- paste0(z_plot$FuncClass,": ",z_plot$Def)

p3a <- ggplot(z_plot, aes(x = Copio, Name)) +
  geom_tile(aes(fill = n)) +
  geom_text(aes(label = round(n, 3))) +
  # As our values are continuous, we'll use scale_fill_continuous instead of scale_fill_manual
  #scale_fill_continuous(low = "violetred", high = "aquamarine") +
  xlab("") + 
  ylab("") + 
  theme_pubclean() +
  theme(legend.position = "none",
        axis.text.y = element_text(size = 10, angle = 0),
        axis.text.x = element_text(angle=30,hjust=1)) +
  scale_fill_gradientn(colours = c("violetred", "gray", "aquamarine"),
                       values = scales::rescale(c(-100 , 0, 150))) +
  scale_x_discrete(expand=c(0,0), labels=c("Overrepresented in Oligotrophs","Overrepresented in Copiotrophs","Difference"))


z_sig <- merge.easy(z,cog_def,key="Hit") %>% subset(pval.corrected<0.01) 
z_sig$P <- (z_sig$n.FALSE + z_sig$n.TRUE)/(z_sig$nGenomes.FALSE + z_sig$nGenomes.TRUE)
z_sig$Copio <- factor(z_sig$P.TRUE>z_sig$P.FALSE)

p2a <- ggplot(z_sig,aes(x=Copio,group=Copio,y=P)) + 
  geom_violin(fill="gray") +
  geom_boxplot(width=0.03) + 
  theme_pubclean() + 
  ylab("Overall Prevalence of Significant Families") + 
  scale_x_discrete(labels=c("Overrepresented in Oligotrophs","Overrepresented in Copiotrophs")) + 
  xlab("")

setwd("~/eggo/Figs")
pdf("Function_enrichment_presenceabsence_eggnogmapper.pdf",width=16,height=8)
ggarrange(ggarrange(p1a,
                    p2a,
                    nrow=2,
                    labels=c("(a)","(b)"),
                    hjust=0),
          p3a,
          ncol=2,
          labels=c("","(c)"),
          widths=c(2,3))
dev.off()

# In Text Figure ---------------------------------------------------------------

x <- z_plot %>% subset(!is.na(Copio)) %>% group_by(FuncClass) %>% summarise(n=sum(n))
many_hits <- x$FuncClass[x$n>=150]
z_plot2 <- z_plot %>% 
  subset(FuncClass %in% many_hits) %>%
  subset(!is.na(Copio))
z_plot2$NameReorder <- factor(z_plot2$Name,
                                 levels=c("G: Carbohydrate transport and metabolism",
                                "E: Amino acid transport and metabolism",
                                "L: Replication, recombination and repair",
                                "J: Translation, ribosomal structure and biogenesis",
                                "F: Nucleotide transport and metabolism",
                                "C: Energy production and conversion"))

p3 <- ggplot(z_plot2,aes(x = Copio, NameReorder)) +
  geom_tile(aes(fill = n)) +
  geom_text(aes(label = round(n, 3))) +
  # As our values are continuous, we'll use scale_fill_continuous instead of scale_fill_manual
  #scale_fill_continuous(low = "violetred", high = "aquamarine") +
  xlab("") + 
  ylab("") + 
  theme_pubclean() +
  theme(legend.position = "none",
        axis.text.y = element_text(angle = 0,size=8,face="bold"),
        axis.text.x = element_text(angle=30,hjust=1,size=8,face="bold")) +
  scale_fill_gradientn(colours = c("gray", "aquamarine")) +
  scale_x_discrete(expand=c(0,0), labels=c("Overrepresented in Oligotrophs","Overrepresented in Copiotrophs"))
p3


# p1 <- ggplot(sdf, aes(x= reorder(Name, -abs(Diff)),y=Diff)) + 
#   geom_bar(stat="identity") + 
#   theme_pubclean() +
#   ylab(expression(E[Copiotrophs]-E[Oligotrophs])) + 
#   xlab("") + 
#   theme(axis.text.x = element_text(angle=70,hjust=1,vjust=1),
#         axis.title.y = element_text(size=18))

p1 <- ggplot(sdf, aes(x= reorder(Name, -abs(Diff)),y=Diff, fill = abs(Diff)>0.005)) + 
  geom_bar(stat="identity") + 
  theme_pubclean() +
  ylab(expression(E[Copiotrophs]-E[Oligotrophs])) + 
  xlab("") + 
  theme(axis.text.x = element_text(angle=70,hjust=1,vjust=1,face="bold"),
        axis.title.y = element_text(size=14),
        legend.position = "none") +
  scale_fill_manual(values=c("grey35","darkred")) 

z_def <- merge.easy(z,cog_def,key="Hit") 
sig <- c("C","K","G","L")
z_def <- z_def %>% subset(FuncClass %in% sig)
p2 <- ggplot(z_def,aes(x=P.TRUE-P.FALSE,y=-log10(pval.corrected),fill=FuncClass)) + 
  geom_point(pch=21,size=2,alpha=0.9) + 
  theme_pubclean() + 
  geom_hline(yintercept=-log10(0.01),lty=2) + 
  xlab(expression(P[Copiotrophs]-P[Oligotrophs])) + 
  ylab(expression(-log[10](p[adjusted]))) + 
  scale_fill_manual(values=c("#1b9e77",
                             "#d95f02",
                             "#7570b3",
                             "#e7298a")) + 
  scale_alpha_manual(values=c(0.1,0.7)) +
  theme(legend.position = c(0.89,0.3),
        legend.background = element_rect(size=0.1,color="black"),
        axis.title = element_text(size = 14)) +
  labs(fill="Function")

setwd("~/eggo/Figs")
pdf("Funtion_COGs.pdf",height=8,width=8)
ggarrange(p1,
          ggarrange(p2,p3,nrow=2,labels=c("(b)","(c)")),
          ncol=2,labels=c("(a)",""))
dev.off()


# For slides

p1 <- ggplot(sdf, aes(x= reorder(Name, -abs(Diff)),y=Diff)) + 
  geom_bar(stat="identity") + 
  theme_pubclean() +
  ylab("Enrichment Copiotrophs vs. Oligotrophs (% Genes)") + 
  xlab("") + 
  theme(axis.text.x = element_text(angle=60,hjust=1,vjust=1),
        axis.title.y = element_text(size=12))

p3 <- ggplot(z_plot2,aes(x = Copio, NameReorder)) +
  geom_tile(aes(fill = n)) +
  geom_text(aes(label = round(n, 3))) +
  # As our values are continuous, we'll use scale_fill_continuous instead of scale_fill_manual
  #scale_fill_continuous(low = "violetred", high = "aquamarine") +
  xlab("") + 
  ylab("") + 
  theme_pubclean() +
  theme(legend.position = "none",
        axis.text.y = element_text(angle = 0,size=10),
        axis.text.x = element_text(angle=20,hjust=1,size=12)) +
  scale_fill_gradientn(colours = c("gray", "aquamarine")) +
  scale_x_discrete(expand=c(0,0), labels=c("Overrepresented in Oligotrophs","Overrepresented in Copiotrophs"))


setwd("~/eggo/Figs")
pdf("Funtion_COGs_slides.pdf",height=8,width=12)
ggarrange(p1,
          p3,
          ncol=2,#labels=c("(a)","(b)"),
          widths=c(5,4),
          hjust=0,
          vjust=1)
dev.off()
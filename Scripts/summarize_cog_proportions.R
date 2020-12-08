

## JLW 2020 - Compile COG Class Proportions

# Load Packages ----------------------------------------------------------------

library(dplyr)
library(ggplot2)
library(ggpubr)
library(boot)
library(data.table)
library(ggrepel)

samp_mean <- function(x, i) {
  mean(x[i],na.rm=T)
}

bootCI <- function(x,upper=T,R=100){
  if(upper==T){
    return(boot.ci(boot(x,statistic=samp_mean,R=R),type="norm")$normal[2])
  } else {
    return(boot.ci(boot(x,statistic=samp_mean,R=R),type="norm")$normal[3])
  }
}

merge.easy <- function(df1,df2,key){
  df1 <- data.table(df1,key=key)
  df2 <- data.table(df2,key=key)
  return(as.data.frame(unique(merge(df1,df2,all.x=TRUE,by=.EACHI,allow.cartesian=TRUE)),stringsAsFactors=F))
}


# Load Data --------------------------------------------------------------------


load("~/eggo/Data/EGGO.RData")
EGGO$Accession <- EGGO$Assembly %>% substr(1,15)
copiotrophs <- EGGO$Accession[EGGO$d<5] 
oligotrophs <- EGGO$Accession[EGGO$d>=5] 

# Process Annotations from eggnogmapper (slow, skip for preloaded) -------------

setwd("~/annotations/")
annot_files <- list.files(pattern="*.out")

prop_df <- data.frame(Accession=character(),
           COG=character(),
           prop=numeric(),
           stringsAsFactors = F)
for(f in annot_files){
  annot_i <- read.delim(f,sep="\t",head=F,stringsAsFactors=F) %>%
    group_by(V3) %>%
    count()
  prop_df <- rbind(prop_df,
                   data.frame(Accession=substr(f,1,15),
                              COG=annot_i$V3,
                              prop=annot_i$n/sum(annot_i$n) ))
}

setwd("~/eggo/Data")
save(prop_df,file="COG_proportions.RData")

# Loaed Processed Annotations --------------------------------------------------

setwd("~/eggo/Data")
load("COG_proportions.RData")

# Deal w/ mixed categories -----------------------------------------------------

prop_df <- prop_df %>% 
  subset(COG %in% setdiff(c(LETTERS[1:22],"X"),c("R","S"))) %>%
  subset(Accession %in% c(copiotrophs,oligotrophs))
prop_df$Copiotroph <- prop_df$Accession %in% copiotrophs  
prop_mix_df <- data.frame(Accession=character(),
                      COG=character(),
                      prop=numeric(),
                      stringsAsFactors = F)
for(acc_i in unique(prop_df$Accession)){
  for(COG_i in setdiff(c(LETTERS[1:22],"X"),c("R","S"))){
    prop_mix_df <- rbind(prop_mix_df,
                         data.frame(Accession=acc_i,
                                    COG=COG_i,
                                    prop=sum(prop_df$prop[prop_df$Accession==acc_i & grepl(COG_i,prop_df$COG)])))
  }
}

prop_mix_df <- prop_mix_df %>% 
  subset(COG %in% setdiff(c(LETTERS[1:22],"X"),c("R","S"))) %>%
  subset(Accession %in% c(copiotrophs,oligotrophs))
prop_mix_df$Copiotroph <- prop_mix_df$Accession %in% copiotrophs 

# Plot proportions -------------------------------------------------------------

prop_summ <- prop_mix_df %>% group_by(Copiotroph,COG) %>% summarize(mean.prop=mean(prop),
                                                                upper=bootCI(prop,upper=T),
                                                                lower=bootCI(prop,upper=F))  %>%
  as.data.frame(stringsAsFactors=F)
psw <- reshape(prop_summ,
               direction = "wide",
               idvar = "COG", 
               timevar = "Copiotroph")


cog_cat_def <- read.delim("~/oligotrophy/Data/fun-20.tab",head = F, sep="\t",stringsAsFactors = F)
names(cog_cat_def) <- c("COG","X","Def")
psw <- merge.easy(psw,cog_cat_def,key="COG")
psw$Name <- paste0(psw$COG,": ",psw$Def)

setwd("~/eggo/Figs")
pdf("Function_enrichment_proportionmixed_eggnogmapper.pdf",width=8,height=8)
ggplot(psw,aes(x=mean.prop.FALSE,y=mean.prop.TRUE,label=Name,color=Name)) + 
  geom_point() +
  geom_errorbar(aes(ymin = lower.TRUE,ymax = upper.TRUE)) + 
  geom_errorbarh(aes(xmin = lower.FALSE,xmax = upper.FALSE)) + 
  geom_abline(slope=1,intercept=0,lty=2) + 
  theme_pubclean()  + 
  geom_text_repel(force=5,
                  size=2.5,
                  arrow = arrow(length=unit(0.01, "inches"),type="closed")) + 
  xlab("Mean Proportion of Genes Assigned to Functional Class in Oligotroph Genomes") + 
  ylab("Mean Proportion of Genes Assigned to Functional Class in Copiotroph Genomes") + 
  theme(legend.position="none") 
dev.off()

# Test for significance --------------------------------------------------------

test_df <- data.frame(COG=character(),
                      pval=numeric(),
                      Diff=numeric())
for(COG_i in unique(prop_mix_df$COG)){
  df <- prop_mix_df %>% subset(COG==COG_i)
  test_df <- rbind(test_df,
                   data.frame(COG=COG_i,
                              pval=wilcox.test(df$prop[df$Copiotroph],df$prop[!df$Copiotroph])$p.value,
                              Diff=mean(df$prop[df$Copiotroph])-mean(df$prop[!df$Copiotroph])))
}
test_df$p.adj <- p.adjust(test_df$pval)

sig_df <- test_df %>% subset(p.adj<0.01)
setwd("~/eggo/Data")
write.csv(sig_df,file="function_enrichment_propmixed_significant.csv")




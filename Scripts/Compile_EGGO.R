
# JLW 2020 - Put together EGGO datasets

# Load Packages ----------------------------------------------------------------

library(dplyr)

# Load Datasets/ Add Metadata --------------------------------------------------

setwd("~/eggo/Data/eggo-data")

# GORG
load("CodonStatistics_GORG.RData")
gorg <- sag_df %>% subset(select=c(nHE,
                                   FilteredSequences, 
                                   d, 
                                   LowerCI, 
                                   UpperCI, 
                                   d.madin, 
                                   LowerCI.madin, 
                                   UpperCI.madin, 
                                   ID)) %>% mutate_all(unlist)%>% 
  mutate(Source = "GORG-tropics",
         Mode = "Partial",
         Type = "SAG",
         Environment = "Marine Surface")
names(gorg)[names(gorg)=="ID"] <- "Assembly"

# MarRef
load("CodonStatistics_MarRef.RData")
marref <- assembly_df %>% subset(select=c(nHE,
                                          FilteredSequences, 
                                          d, 
                                          LowerCI, 
                                          UpperCI, 
                                          d.madin, 
                                          LowerCI.madin, 
                                          UpperCI.madin, 
                                          Assembly)) %>% 
  mutate(Source = "MarRef",
         Mode = "Full",
         Type = "Isolate",
         Environment = "Marine")

#Nayfach
load("CodonStatistics_Nayfach.RData")
nayfach <- mag_df %>% subset(select=c(nHE,
                                      FilteredSequences, 
                                      d, 
                                      LowerCI, 
                                      UpperCI, 
                                      d.madin, 
                                      LowerCI.madin, 
                                      UpperCI.madin, 
                                      Assembly)) %>% 
  mutate(Source = "Nayfach et al.",
         Mode = "Partial",
         Type = "MAG",
         Environment = "Human Gut")

#Pasolli
load("CodonStatistics_Pasolli.RData")
pasolli <- mag_df %>% subset(select=c(nHE,
                                      FilteredSequences, 
                                      d, 
                                      LowerCI, 
                                      UpperCI, 
                                      d.madin, 
                                      LowerCI.madin, 
                                      UpperCI.madin, 
                                      Assembly)) %>% 
  mutate(Source = "Pasolli et al.",
         Mode = "Partial",
         Type = "MAG",
         Environment = "Human Microbiome")

#Poyet
load("CodonStatistics_Poyet.RData")
poyet <- assembly_df %>% subset(select=c(nHE,
                                         FilteredSequences, 
                                         d, 
                                         LowerCI, 
                                         UpperCI, 
                                         d.madin, 
                                         LowerCI.madin, 
                                         UpperCI.madin, 
                                         Assembly)) %>% 
  mutate(Source = "Poyet et al.",
         Mode = "Full",
         Type = "Isolate",
         Environment = "Human Gut")

#refseq
load("CodonStatistics_RefSeq.RData")
refseq <- assembly_df %>% subset(select=c(nHE,
                                          FilteredSequences, 
                                          d, 
                                          LowerCI, 
                                          UpperCI, 
                                          d.madin, 
                                          LowerCI.madin, 
                                          UpperCI.madin, 
                                          Assembly)) %>% 
  mutate(Source = "RefSeq Assemblies",
         Mode = "Full",
         Type = "Isolate",
         Environment = "") %>%
  subset(nHE>50 & nHE<70)

#tully
load("CodonStatistics_Tully.RData")
tully <- mag_df %>% subset(select=c(nHE,
                                    FilteredSequences, 
                                    d, 
                                    LowerCI, 
                                    UpperCI, 
                                    d.madin, 
                                    LowerCI.madin, 
                                    UpperCI.madin, 
                                    Assembly)) %>% 
  mutate(Source = "Tully et al.",
         Mode = "Partial",
         Type = "MAG",
         Environment = "Marine")

#Delmont
load("CodonStatistics_Delmont.RData")
delmont <- mag_df %>% subset(select=c(nHE,
                                      FilteredSequences, 
                                      d, 
                                      LowerCI, 
                                      UpperCI, 
                                      d.madin, 
                                      LowerCI.madin, 
                                      UpperCI.madin, 
                                      Assembly)) %>% 
  mutate(Source = "Delmont et al.",
         Mode = "Partial",
         Type = "MAG",
         Environment = "Marine")

#Parks
load("CodonStatistics_Parks.RData")
parks <- mag_df %>% subset(select=c(nHE,
                                    FilteredSequences, 
                                    d, 
                                    LowerCI, 
                                    UpperCI, 
                                    d.madin, 
                                    LowerCI.madin, 
                                    UpperCI.madin, 
                                    Assembly)) %>% 
  mutate(Source = "Parks et al.",
         Mode = "Partial",
         Type = "MAG",
         Environment = "")

#zou
load("CodonStatistics_Zou.RData")
zou <- assembly_df %>% subset(select=c(nHE,
                                       FilteredSequences, 
                                       d, 
                                       LowerCI, 
                                       UpperCI, 
                                       d.madin, 
                                       LowerCI.madin, 
                                       UpperCI.madin, 
                                       Assembly)) %>% 
  mutate(Source = "Zou et al.",
         Mode = "Full",
         Type = "Isolate",
         Environment = "Human Gut")

rm("mag_df","sag_df","assembly_df")

# Put Together -----------------------------------------------------------------

EGGO <- rbind(refseq,
              parks,
              gorg,
              tully,
              delmont,
              marref,
              nayfach,
              pasolli,
              poyet,
              zou) %>% 
  subset(nHE>=10)

table(EGGO$Environment)
table(EGGO$Type)
table(EGGO$Source)
sum(refseq$d>100)

setwd("~/eggo/Data")
save(EGGO,file="EGGO.RData")
write.csv(EGGO,file="EGGO.csv")


# Get genome list for eggnog mapping

seed(123456)

library(dplyr)

setwd("~/eggo/Data/eggo-data")
load("CodonStatistics_RefSeq.RData")

sub_genus <- assembly_df %>% group_by(Genus) %>% sample_n(1)
#sub_family <- assembly_df %>% group_by(Family) %>% sample_n(1)

setwd("~/eggo/Data/")
write(substr(sub_genus$Assembly,1,15),file="subsampled_assemblies.txt")

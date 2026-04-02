{
  library(ggplot2)
  library(dplyr)
  library(stringr)
  library(ggpubr)
  library(tidyverse)
  library(geomtextpath)
  library(scone)
  library(tidyr)
  library(changepoint)
  library(reshape2)
  setwd(dir = "E:/r_file")
} #lib
#ZW gametologs pair dS 
synZW <- read.table(file = "zwgame.identity.alignPerc") # with gametolog pairs
dndsZW <- read.table('psi.zw.game.dnds.info',header = F,row.names = 1)
geneloci <- read.table(file="psi.geneID.finalinfo.txt")
geneloci <- geneloci %>% filter(V1 == "chrZ" | V1 == "chrW")
synZW <- synZW %>% mutate(
  zID = geneloci$V5[match(synZW$V5,geneloci$V2)],
  dN = dndsZW$V2[match(zID,rownames(dndsZW))],
  dS = dndsZW$V3[match(zID,rownames(dndsZW))],
  w = dndsZW$V4[match(zID,rownames(dndsZW))]
) %>% mutate(
  strata = case_when(V6 <= 1517136 ~ "PAR",
                     V5 >= 1517136 & V6 <= 3178082 ~ "S3",
                     V5 >= 3178082 & V6 <= 10.4e6 ~ "S1",
                     V5 >= 10.4e6 & V6 <= 12.5e6 ~ "S2",
                     V5 >= 12.5e6 ~ "S0"
  )) 
ggplot(synZW, aes(x = V5/1e6, y = dS))+
  geom_point(aes(color = strata),size =3) + theme_test()+ 
  scale_color_manual(values=c("#a7d49e","#ed7d61","#fede8e","#8cb7df")) +
  geom_vline(aes(xintercept=1.51),linetype=2)+ guides(color = F)+
  geom_vline(aes(xintercept=10.4),linetype=2) +
  geom_vline(aes(xintercept=3.17),linetype=2) + geom_vline(aes(xintercept=12.5),linetype=2) + ylim(0,1)
ggplot(synZW,
       aes(x = strata, y = dS)) +
  geom_violin(aes(fill = strata), trim = FALSE) + geom_point() +
  stat_summary(fun = mean, geom = "crossbar", 
               width = 0.2, color = "black", 
               fatten = 2) +
  stat_compare_means(method = "wilcox.test", 
                     comparisons = list(c("S0", "S1"),
                                        c("S1", "S2"),
                                        c("S2", "S3")), 
                     method.args = list(alternative = "greater")) + 
  theme_test() + ylim(0,0.5)

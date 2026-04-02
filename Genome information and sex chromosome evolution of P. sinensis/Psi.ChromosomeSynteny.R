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
#Chromosome synteny   
psi2GgaSC <- read.table(file = "psisex-chick.delta.filt.coords", skip = 4)  
psi2GgaZ <- psi2GgaSC %>% filter(V13 == "chrZ" &V5 >= 500)
p1 <- ggplot(data = psi2GgaZ) +
  geom_segment(aes(x=V3/1000000,xend=V1/1000000,y=1,yend=2)) + theme_test() + xlim(0,40) +
  guides(color="none") 
psizw <- read.table("psizw.nucmer.coords", skip = 4)
psizw <- psizw %>% filter(
  V7 > 60 & V7 < 90 & V5 >= 1000
) %>% mutate(
  strata = case_when(V2 <= 1517136 ~ "PAR",
                     V1 >= 1517136 & V2 <= 3178082 ~ "S3",
                     V1 >= 3178082 & V2 <= 10.4e6 ~ "S1",
                     V1 >= 10.4e6 & V2 <= 12.5e6 ~ "S2",
                     V1 >= 12.5e6 ~ "S0"
  ))
p2 <- ggplot(data = psizw) +
  geom_segment(aes(x=V4/1000000,xend=V2/1000000,y=1,yend=2,color=strata)) + theme_test() + xlim(0,40)+ 
  geom_vline(aes(xintercept=1.51),linetype=2)+
  geom_vline(aes(xintercept=3.17),linetype=2)+
  geom_vline(aes(xintercept=10.4),linetype=2)+
  geom_vline(aes(xintercept=12.5),linetype=2) +
  guides(color="none") 
ggarrange(p1,p2, ncol =1) 

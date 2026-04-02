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
#Faster Z
dndsAll <- read.table("psi2tse.dnds.pair.txt") # pair-wise w against T. s. elegans (1-to-1 RBH orthologs)
geneloci <- read.table("trans.loci")
zsub <- dndsAll %>% filter(V1 %in% subset(geneloci, V1 == "chrZ")$V4) %>% 
  filter(V4 >= 0.01 & V4 <= 2) %>% mutate(
    chr = case_when(
      V7 <= 1.5e6 ~ "Z-PAR",
      V6 > 1.5e6 & V7 < 3.17e6 ~ "Z-S3",
      V6 > 3.17e6 & V7 < 10.4e6 ~ "Z-S1",
      V6 > 10.4e6 & V7 < 12.5e6 ~ "Z-S2",
      V6 > 12.5e6 ~ "Z-S0"
    ))
wsub <- dndsAll %>% filter(V1 %in% subset(geneloci, V1 == "chrW")$V4) %>% 
  filter(V4 >= 0.01 & V4 <= 2) %>%
  mutate(chr = "chrW") 
asub <- dndsAll %>% 
  filter(V1 %in% subset(geneloci, V1 != "chrZ" & V1 != "chrW" & !grepl("Hic", V1))$V4) %>%
  filter(V4 >= 0.01 & V4 <= 2) %>% 
  mutate(chr = gsub("chr","", V5),
         chr = ifelse(chr < 9, "micro", "macro")) 
dfsub <- rbind(asub,zsub,wsub) %>% mutate(
  V2 = ifelse(V2 < 0 ,0, V2)
)
ggplot(dfsub,aes(x = chr, y=V2,fill=chr))+
  geom_violin(trim = FALSE) + theme_bw() + 
  stat_summary(fun = mean, geom = "crossbar", 
               width = 0.2, color = "white", 
               fatten = 2)  +
  stat_compare_means(method = "wilcox.test", 
                     comparisons = list(c("Z-S0", "macro"),
                                        c("Z-S1", "macro"),
                                        c("chrW", "macro")), 
                     method.args = list(alternative = "greater"))+
  theme_test() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        axis.text.x= element_text(size=10),
        axis.text.y= element_text(size=15))+
  labs(x="",y="dN/dS") + ylim(-0.5,1)

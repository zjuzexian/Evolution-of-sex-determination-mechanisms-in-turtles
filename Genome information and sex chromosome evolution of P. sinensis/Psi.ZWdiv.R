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
#ZW similarity
zw.simi <- read.table(file = "zw.simi.lastz.final100k")
colnames(zw.simi)<-c("Chr","Pos","Mis","All","Similarity")
zw.simi <- zw.simi %>%
    filter(All>5000) %>% mutate(strata = 
                                  case_when(
                                    Pos <= 1.51e6 ~ "PAR",
                                    Pos > 1.51e6 & Pos < 3.17e6 ~ "S3",
                                    Pos > 3.17e6 & Pos < 10.4e6 ~ "S1",
                                    Pos > 10.4e6 & Pos < 12.5e6 ~ "S2",
                                    Pos > 12.5e6 ~ "S0"
                                  )) %>% mutate(
                                    color_group = ifelse(Similarity > 89.8, "high",
                                                         ifelse(Similarity >= 87 & Similarity <= 89.8, "mid", "low"))
                                  )
ggplot(zw.simi,aes(x=Pos/1e6,y=Similarity)) + 
    geom_point(aes(size = All,color= strata)) + 
    geom_smooth(method = "loess", span = .15, se =F, color = "#dcdcdc") +
    labs(x="Position(Mbp)",y="Z-W Similarity")+ theme_bw() + theme(legend.position = "none") +
  scale_color_manual(values = c("PAR"="#add69d","S3"="#dc6a67","S2"="#f3d68a",
                                "S1"="#8ac4eb","S0"="#8ac4eb")) + theme_test()+ ylim(78,94)+
  geom_vline(xintercept = c(1.51,3.17,10.4,13.1),linetype =2)+
  scale_size(range = c(.5, 4))
ggplot(zw.simi,aes(x=strata,y=Similarity)) + 
  geom_boxplot(aes(color= strata)) +
  stat_summary(fun = mean, geom = "crossbar", 
               width = 0.2, color = "black", 
               fatten = 2) +
  stat_compare_means(method = "wilcox.test", 
                     comparisons = list(c("S0", "S1"),
                                        c("S1", "S2"),
                                        c("S2", "S3")), 
                     method.args = list(alternative = "less"))
cptDf <- subset(zw.simi %>% filter(Similarity >= 78 & Similarity <= 94)) %>% select(x = 2, y = 5) 
cpt <- cpt.mean(cptDf$y, method = "PELT") 
cpts(cpt)
zw.simi %>% filter(Similarity >= 78 & Similarity <= 94) %>% filter(row_number(.) %in% cpts(cpt))

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
#DNA polymorphism (pai)
pai <- read.table(file = "pai500k-bialle.windowed.pi",header=T)
cptDf <- subset(pai,pai$CHROM=="chrZ"&pai$N_VARIANTS>200) %>% select(2, 3, 5) %>% mutate(
  x = (BIN_START+BIN_END)/2e6,
  y = PI
) %>% select(4,5)
cpt <- cpt.meanvar(cptDf$y, method = "PELT") 
pai %>% filter(CHROM == "chrZ") %>% filter(row_number(.) %in% cpts(cpt))
subpai <- subset(pai,pai$CHROM=="chrZ"&pai$N_VARIANTS>0) %>% mutate(
  strata = case_when(BIN_END <= 1.5e6 ~ "PAR",
                     BIN_START >= 1.5e6 & BIN_END <= 3.0e6 ~ "S3",
                     BIN_START >= 3.0e6 & BIN_END <= 10.5e6 ~ "S1",
                     BIN_START >= 10.5e6 & BIN_END <= 12.5e6 ~ "S2",
                     BIN_START >= 12.5e6 ~ "S0"
  )) 
ggplot(subpai,aes(x=(BIN_START)/1e6,y=PI*1e4))+ 
  geom_point(aes(color = strata))  + guides(size="none") +
  scale_color_manual(values=c("#a7d49e","#ed7d61","#fede8e","#8cb7df")) +
  geom_smooth(method = "loess", span = .4)+ theme_bw()+xlim(0,40)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        axis.text.x= element_text(size=15),
        axis.text.y= element_text(size=15)) + 
  geom_vline(aes(xintercept=1.51),linetype=2) +
  geom_vline(aes(xintercept=3.17),linetype=2) +
  geom_vline(aes(xintercept=10.4),linetype=2) +
  geom_vline(aes(xintercept=12.5),linetype=2) +
  geom_hline(aes(yintercept = mean(subset(subpai, strata == "S1")$PI*1e4)))+
  geom_hline(aes(yintercept = mean(subset(subpai, strata == "S0")$PI*1e4)))+
  guides(color="none") + theme_test() +
  labs(x="",y="") + ylim(0, 4)
ggplot(subpai %>% filter(strata %in% c("S0", "S1", "S2")),
       aes(x = strata, y = PI * 1e4)) +
  geom_violin(aes(fill = strata), trim = FALSE) +
  stat_summary(fun = mean, geom = "crossbar", 
               width = 0.2, color = "black", 
               fatten = 2) +
  stat_compare_means(method = "wilcox.test", 
                     comparisons = list(c("S0", "S1"),
                                        c("S1", "S2"),
                                        c("S2", "S3")), 
                     method.args = list(alternative = "greater")) +
  theme_test()

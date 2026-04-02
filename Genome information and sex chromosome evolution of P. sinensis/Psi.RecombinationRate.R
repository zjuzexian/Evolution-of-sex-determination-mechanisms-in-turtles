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
#RC rate
rc <- read.table(file = "chrZ-rcrate-200k",header=T)
rc1 <- subset(rc,rc$CIL>=0) %>%
  mutate(strata = case_when(
    End/1e6 <= 1.51 ~ "par",
    End/1e6 > 1.51 & End/1e6 <= 3.17 ~ "s3",
    End/1e6 > 3.17 & End/1e6 <= 10.4 ~ "s1",
    End/1e6 > 10.4 & End/1e6 <= 12.5 ~ "s2",
    End/1e6 > 12.5 ~ "s0"
  ))
rc1 <- rc1 %>%
  mutate(pos = (Start + End) / 2 / 1e6)  
ggplot(rc1, aes(x = Start/1e6)) +
  geom_ribbon(aes(ymin = CIL, ymax = CIR), fill = "#d0d9e8", alpha = 0.6) +
  geom_step(aes(y = Rho), color = "#2a6fba", size = 1) +
  labs(
    x = "Position on Chromosome (Mb)",
    y = "Recombination rate (ρ)"
  ) +
  theme_test() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  )+ geom_vline(aes(xintercept=1.51),linetype=2)+ geom_vline(aes(xintercept=10.4),linetype=2)+
  geom_vline(aes(xintercept=3.17),linetype=2) + geom_vline(aes(xintercept=12.5),linetype=2) 
ggplot(data = rc1, aes(x = strata, y = Rho, color = strata),shape =21) +
  geom_boxplot()  +
  theme_bw() +
  stat_compare_means(
    comparisons = list(c("s0", "s1"),c("s1","s2")),
    method = "wilcox.test",
    method.args = list(alternative = "greater")
  ) + 
  scale_color_manual(
    values = c("par" = "#a7d49e",
               "s3" = "#ed7d61",
               "s2" = "#fede8e",
               "s1" = "#8cb7df",
               "s0" = "#8cb7df"))
cptDf <- rc1 %>% select(x=7, y=3)
cpt <- cpt.mean(cptDf$y, method = "PELT") 
rc1 %>% filter(row_number(.) %in% cpts(cpt))

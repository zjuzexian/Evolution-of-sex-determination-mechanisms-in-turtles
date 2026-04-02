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
depthall <- read.table(file = "psi.dna.fm.100k.noN.txt") #e.g. 100k sliding window
cptDf <- subset(depthall,depthall$V1=="chrZ") %>% select(2, 3, 4, 6) %>% mutate(
  x = (V2+V3)/2e6,
  y = V6/V4
) %>% select(5,6)
cpt <- cpt.meanvar(cptDf$y, method = "PELT") 
depthall %>% filter(V1 == "chrZ") %>% filter(row_number(.) %in% cpts(cpt))
zdepth <- depthall %>% filter(V1 == "chrZ") %>% mutate(
  strata = case_when(V3 <= 1517136 ~ "PAR",
                     V2 >= 1517136 & V3 <= 3178082 ~ "S3",
                     V2 >= 3178082 & V3 <= 10.4e6 ~ "S1",
                     V2 >= 10.4e6 & V3 <= 12.5e6 ~ "S2",
                     V2 >= 12.5e6 ~ "S0"
  ))
strata_means <- zdepth %>%
  mutate(pos = (V2 + V3) / 2e6,
         ratio = V6 / V4) %>%
  group_by(strata) %>%
  summarise(
    mean_ratio = mean(ratio, na.rm = TRUE),
    x_pos = mean(pos, na.rm = TRUE)
  )
ggplot(zdepth, aes((V2 + V3) / 2e6, V6 / V4)) + 
  geom_point(aes(color = strata), size = 1) +
  stat_smooth(method = "loess", alpha = 0.3, span = 0.08) +
  labs(x = "", y = "DNA coverage") +
  theme_bw() +
  scale_color_manual(values = c("#a7d49e", "#ed7d61", "#fede8e", "#8cb7df", "#8cb7df")) +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 15)
  ) +
  guides(color = "none") +
  ylim(0, 2.5) +
  geom_hline(data = strata_means, aes(yintercept = mean_ratio, color = strata), linetype = "dashed")

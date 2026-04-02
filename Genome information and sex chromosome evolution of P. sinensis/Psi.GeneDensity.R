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
#Gene density
geneDens <- read.table(file = "psi.geneDens.100k.txt")
geneDensZ <- geneDens %>% filter(V1 == "chrZ") %>% mutate(
  strata = case_when(V3 <= 1517136 ~ "PAR",
                     V2 >= 1517136 & V3 <= 3178082 ~ "S3",
                     V2 >= 3178082 & V3 <= 10.4e6 ~ "S1",
                     V2 >= 10.4e6 & V3 <= 12.5e6 ~ "S2",
                     V2 >= 12.5e6 ~ "S0"
  ))
ggplot(geneDensZ, aes(x = (V2+V3)/2e6, y = V5)) +
  geom_bar(stat = "identity", aes(fill = strata)) + theme_test()

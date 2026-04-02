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
#Dosage compensation
exp <- read.table(file = "psi.tpmmean.NewFilt.txt",header = T,sep = ',',row.names = 1)
geneloci <- read.table(file = "psi.geneID.finalinfo.txt")
targetExp <- exp %>%
  rownames_to_column(var = "gene") %>%
  pivot_longer(cols = -gene, names_to = "sample", values_to = "expression") %>%
  mutate(
    tissue = str_sub(sample, 1, 2),
    sample_num = as.numeric(str_extract(sample, "\\d{2}$"))
  ) %>%
  filter(!is.na(sample_num) & sample_num >= 16) %>%
  group_by(gene, tissue) %>%
  summarise(mean_expression = mean(expression, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    chr = geneloci$V1[match(gene, geneloci$V4)],
    st = geneloci$V2[match(gene, geneloci$V4)],
    ed = geneloci$V3[match(gene, geneloci$V4)]
  ) %>%
  filter(chr != "chrW") %>%
  mutate(
    chr = if_else(chr == "chrZ", "chrZ", "Autosome"),
    strata = case_when(
      chr == "chrZ" & ed <= 1517136 ~ "PAR",
      chr == "chrZ" & st > 1517136 & ed <= 3178082 ~ "S3",
      chr == "chrZ" & st >= 3178082 & ed <= 10.4e6 ~ "S1",
      chr == "chrZ" & st >= 10.4e6 & ed <= 12.5e6 ~ "S2",  
      chr == "chrZ" & st >= 12.5e6 ~ "S0",
      chr == "Autosome" ~ "NA",
      TRUE ~ NA_character_
    )
  ) 
dfToshow <- targetExp %>% filter(tissue %in% c("FB", "MB")) %>% filter(mean_expression > 5)
ggplot(data = dfToshow, aes(x = strata, y = log2(mean_expression))) +
  geom_violin(trim = FALSE, aes(fill = tissue)) +
  stat_summary(
    aes(group = tissue),
    fun = mean,
    geom = "crossbar",
    width = 0.5,
    position = position_dodge(width = 0.9),
    color = "white",
    fatten = 2
  ) +
  stat_compare_means(aes(fill = tissue),
                     method = "wilcox.test",
                     method.args = list(alternative = "greater")
  ) +
  theme_test()

dfToshowMerge <- targetExp %>% filter(tissue %in% c("FB", "MB")) %>%
  filter(mean_expression > 5) %>% filter(strata != "PAR") 
ggplot(data = dfToshowMerge, aes(x = chr, y = log2(mean_expression))) +
  geom_violin(trim = FALSE, aes(fill = tissue)) +
  stat_summary(
    aes(group = tissue),
    fun = mean,
    geom = "crossbar",
    width = 0.5,
    position = position_dodge(width = 0.9),
    color = "white",
    fatten = 2
  ) +
  stat_compare_means(aes(fill = tissue),
                     method = "wilcox.test",
                     method.args = list(alternative = "less")
  ) +
  theme_test()



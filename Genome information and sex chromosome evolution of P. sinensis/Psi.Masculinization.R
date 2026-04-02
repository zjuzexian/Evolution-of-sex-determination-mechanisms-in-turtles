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
#Masculinization of Z-linked genes   
psiDEgene <- read.table(file = "psi.gonad.DEgene.txt")
psigeneloci <- read.table(file="psi.geneID.finalinfo.txt")
psiDEgeneSub <- psiDEgene %>%
  mutate(
    type = case_when(
      log2FoldChange >= 1 & padj <= 0.05 ~ "Testis-biased", 
      log2FoldChange <= -1 & padj <= 0.05 ~ "Ovary-biased",
      TRUE ~ "Non-biased"
    ),
    chr = psigeneloci$V1[match(rownames(.), psigeneloci$V4)],
    st = psigeneloci$V2[match(rownames(.), psigeneloci$V4)],
    ed = psigeneloci$V3[match(rownames(.), psigeneloci$V4)]
  ) %>%
  mutate(
    chr = case_when(
      chr %in% paste0("chr", 1:8) ~ "macro",
      chr %in% paste0("chr", 9:32) ~ "micro",
      TRUE ~ chr
    )
  ) %>%
  filter(chr != "chrW") %>% filter(!grepl("HiC", chr)) # get DE genes by chr and types
testisGprop <- psiDEgeneSub %>% group_by(type, chr) %>%
  summarise(freq = n(), .groups = "drop") %>%
  group_by(chr) %>%
  mutate(perc = freq / sum(freq) * 100)
ggplot(data = testisGprop, aes(x = chr, y = perc, fill = type)) +
  geom_bar(stat = "identity") + labs(x = "", y = "% Genes") +
  theme_test() + scale_fill_manual(values=c("#a7d49e","#ed7d61","#8cb7df"))
{
M_freq <- dfForFreq %>%
  filter(type == "M") %>%
  group_by(chr_group = ifelse(chr == "chrZ", "chrZ", "other")) %>%
  summarise(M = sum(freq))
Total_freq <- dfForFreq %>%
  group_by(chr_group = ifelse(chr == "chrZ", "chrZ", "other")) %>%
  summarise(total = sum(freq))
M_vs_nonM <- left_join(M_freq, Total_freq, by = "chr_group") %>%
  mutate(nonM = total - M) %>%
  select(-total)
mat <- as.matrix(M_vs_nonM %>% column_to_rownames("chr_group"))
chisq.test(mat)
} # Significances

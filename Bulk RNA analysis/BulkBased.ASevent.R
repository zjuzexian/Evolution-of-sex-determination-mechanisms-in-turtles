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
#AS event in turtles
tseAS <- read.table("tse.AS.filt.count.txt",header = F,fill = T) 
psiAS <- read.table("psi.AS.filt.count.txt",header = F,fill = T) 
orthlogs <- read.table("psi-red.orth-1to1.50")
tseASOrth <- tseAS %>% filter(V1 %in% orthlogs$V2) %>%
  group_by(V3, V4) %>% mutate(
    V5 = as.numeric(V2)
  ) %>%
  summarise(V2_sum = sum(V5), .groups = "drop") %>% mutate(
    species = "tse"
  ) 
psiASOrth <- psiAS %>% filter(V1 %in% orthlogs$V1) %>%
  group_by(V3, V4) %>% mutate(
    V5 = as.numeric(V2)
  ) %>%
  summarise(V2_sum = sum(V5), .groups = "drop") %>% mutate(
    species = "psi"
  ) 
dfForall <- rbind(tseASOrth,psiASOrth) %>%
  filter(!(V4 %in% c("brain", "gonad", "kidney")))
colnames(dfForall) <- c("type", "tissues", "number", "species")
df_ratio <- dfForall %>% 
  group_by(type, tissues) %>%
  pivot_wider(names_from = species, values_from = number) %>%
  mutate(ratio = tse / psi) %>%
  ungroup()
df_avg_ratio <- df_ratio %>%
  group_by(tissues) %>%
  summarise(weighted_ratio = sum(ratio * tse) / sum(tse)) %>%
  ungroup()
typecolorpanel <- c("A3SS"= "#cd0000","A5SS"="#0000cd","MXE"="#4daf4a","SE"="#ff7f00","RI"="#984ea3")
tissuecolorpanel <- c("femaleBrain"= "#cd0000","femaleKidney"="#cd0000","ovary"="#cd0000",
                      "maleBrain"= "#0000cd","maleKidney"="#0000cd","testis"="#0000cd")
#Absolute number of AS events in each tissue
ggplot(dfForall, aes(x = tissues, y = number, fill = type)) +
  geom_bar(stat = "identity") + scale_fill_manual(values = typecolorpanel) +
  labs(x = "Species x Tissue", y = "AS number") + theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  facet_wrap( ~ species)  
ggplot(df_avg_ratio, aes(x = tissues, y = weighted_ratio, fill = tissues)) +
  geom_bar(stat = "identity") + geom_hline(yintercept = 1)+
  scale_fill_manual(values = tissuecolorpanel) +
  labs(x = "Species x Tissue", y = "AS Ratio") + theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +ylim(0,2.8)
#Calculate by orthologs
tseASbyGene <- tseAS %>% filter(V1 %in% orthlogs$V2) %>%
  filter(V3 %in% c("SE", "RI")) %>% 
  group_by(V1, V4) %>%
  summarise(V2_sum = sum(as.numeric(V2)), .groups = "drop") %>% mutate(
    species = "tse"
  ) %>% mutate(
    V1 = orthlogs$V1[match(V1,orthlogs$V2)],
    factor = paste(V1, V4)
  )  %>%
  filter(!(V4 %in% c("brain", "gonad", "kidney"))) 
psiASbyGene <- psiAS %>% filter(V1 %in% orthlogs$V1) %>%
  filter(V3 %in% c("SE", "RI")) %>% 
  group_by(V1, V4) %>% 
  summarise(V2_sum = sum(as.numeric(V2)), .groups = "drop") %>% mutate(
    species = "psi",
    factor = paste(V1, V4)
  )  %>%
  filter(!(V4 %in% c("brain", "gonad", "kidney"))) 
dfASbyGene <- merge(tseASbyGene, psiASbyGene, by = "factor", all = TRUE) %>% 
  dplyr::select(1, 4, 8) %>%       
  mutate_all( ~ replace_na(., 0)) %>%
  separate(factor, into = c("GeneID", "Tissue"), sep = " ")  
colnames(dfASbyGene) <- c("GeneID", "Tissue", "ASnumTse", "ASnumPsi")
dfASbyGeneRatio <- dfASbyGene %>% filter(ASnumTse > 0 & ASnumPsi > 0)
tissuecolorpanel <- c("femaleBrain"= "#cd0000","femaleKidney"="#cd0000","ovary"="#cd0000",
                      "maleBrain"= "#0000cd","maleKidney"="#0000cd","testis"="#0000cd")
ggplot(dfASbyGeneRatio, aes(x = Tissue, y = log2((ASnumTse) / (ASnumPsi)), color = Tissue)) +
  geom_boxplot(notch = T) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_color_manual(values = tissuecolorpanel) +
  stat_compare_means(method = "wilcox.test",
                     ref.group = "testis",
                     label = "p.signif") +  
  theme_minimal()





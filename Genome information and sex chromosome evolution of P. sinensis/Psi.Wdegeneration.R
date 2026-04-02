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
#W gametolog degeneration
zBestHit <- read.table(file = "psi.z.besthit.txt")
geneloci <- read.table(file = "psi.geneID.finalinfo.txt")
zBestHit <- zBestHit %>% mutate(
  st = geneloci$V2[match(V1, geneloci$V5)],
  ed = geneloci$V3[match(V1, geneloci$V5)],
  prod = geneloci$V6[match(V1, geneloci$V5)],
  prod = gsub("_[0-9]*","",prod)
)
wgeneORF <- read.table(file = "psi.w.disrupted.txt") #disrupted orf
expGene <- read.table(file = "psi.tpmAtLeast1.txt") #exp or silence
zgeneInfo <- zBestHit %>% mutate(
  strata = case_when(ed <= 1.5e6 ~ "PAR",
                     st >= 1.5e6 & ed <= 3.2e6 ~ "S3",
                     st >= 3.2e6 & ed <= 10.5e6 ~ "S1",
                     st >= 10.5e6 & ed <= 12.5e6 ~ "S2",
                     st >= 12.5e6 ~ "S0"
  ),
  orf = ifelse(V2 %in% wgeneORF$V1, "Intact", "Disrupted"),
  tpm = ifelse(V2 %in% expGene$x, "Expressed", "Silenced"),
  status = ifelse(V3 < 50 | (V4/V5 < .5), "Deleted", "Retained")
) %>% select(2,8,11:14)
zgeneAll <- geneloci %>% filter(V1 == "chrZ") %>% mutate(
  strata = case_when(V3 <= 1.5e6 ~ "PAR",
                     V2 >= 1.5e6 & V3 <= 3.0e6 ~ "S3",
                     V2 >= 3.0e6 & V3 <= 10.5e6 ~ "S1",
                     V2 >= 10.5e6 & V3 <= 12.5e6 ~ "S2",
                     V2 >= 12.5e6 ~ "S0"
  )
)
dfForBar <- data.frame(
  Region = rep(c("S0", "S1", "S2", "S3"), each = 4),
  Status = rep(c("Fun", "Silence", "ORF disruption", "Deleted"), 4),
  Count = c(11, 13, 35-13, 219, 2, 4, 16-4, 57, 1, 1, 2, 32, 1, 0, 2, 39)
) %>%
  group_by(Region) %>%
  mutate(Proportion = Count / sum(Count))
ggplot(dfForBar, aes(x = Region, y = Proportion, fill = Status)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(y = "Proportion", x = NULL, fill = "Gene Status") +
  theme_minimal() +
  theme(
    text = element_text(size = 14),
    legend.position = "right"
  )
 
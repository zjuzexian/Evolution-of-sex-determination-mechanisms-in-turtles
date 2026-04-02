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
#Tissue specificity comparison between turtles 
tseGeneInfo <- read.table(file = "tse.wgcna.module.txt")
psiGeneInfo <- read.table(file = "psi.wgcna.module.txt.new")
orthologs <- read.table("psi-red.orth-1to1.50")
geneRename <- read.table("psi.geneID.finalinfo.txt")
tseSub <- tseGeneInfo %>% filter(gene %in% orthologs$V2) %>%
  mutate(product = 
           geneRename$V6[match(orthologs$V1[match(gene, orthologs$V2)],geneRename$V4)]
  )
psiSub <- psiGeneInfo %>% filter(gene %in% orthologs$V1) %>%
  mutate(product = geneRename$V6[match(gene, geneRename$V4)]
  )
dfmerge <- merge(tseSub, psiSub, by = "product" ,all = T) %>% dplyr::select(5,7,2,4,1) %>% 
  mutate(product = gsub("_[0-9]*", "", product))
colnames(dfmerge) <- c("psiID","psiType","tseID","tseType","Product")
dfmerge <- dfmerge %>%
  mutate(psiID = ifelse(is.na(psiID), 
                        orthologs$V1[match(tseID, orthologs$V2)], 
                        psiID)) %>%
  dplyr::select(1, 2, 4, 5) %>%
  mutate(across(everything(), ~replace(., is.na(.), "Not tissue-biased")))
df_freq <- dfmerge %>%
  filter(!is.na(psiType) & !is.na(tseType)) %>% 
  group_by(psiType, tseType) %>% 
  summarise(freq = n(), .groups = "drop") %>% # Count tissue specificity
  mutate(type = ifelse(psiType != tseType, "Shifted", "Conserved"))
summary_table <- df_freq %>% filter(psiType != "Not tissue-biased" & tseType != "Not tissue-biased") %>%
  group_by(tseType, type) %>%
  summarise(n = sum(freq), .groups = "drop") %>%
  rename(tseType = "tissue")
#Gonad had higher proportion of tissue-shifted genes
ggplot(summary_table, aes(x = tissue, y = n, fill = type)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_minimal() +
  labs(x = "Tissue", y = "Count", fill = "Type", title = "") +
  scale_fill_manual(values = c("Conserved" = "#bdbdbd", "Shifted" = "#737373"))
#Gonad - brain : X-squared = 600.52, df = 1, p-value < 2.2e-16
#Gonad - kidney : X-squared = 5.42, df = 1, p-value = 0.020
library(ggalluvial)
#Lots of genes undergo transcript reprogramming during TSD-GSD transition
ggplot(df_freq, aes(axis1 = tseType, axis2 = psiType, y = freq)) +
  geom_alluvium(aes(fill = psiType),width = 1/8, curve_type = "sine")+ 
  geom_stratum(width = 1/8) +
  geom_text(stat="stratum",aes(label=after_stat(stratum))) +
  scale_x_discrete(limits = c("psi", "tse"), expand = c(0.15, 0.15)) +
  scale_fill_manual(values = c("Gonad" = "#984ea3", "Brain" = "#4daf4a",
                               "Kidney" = "#ff7f00", "Not tissue-biased" = "#bdbdbd")) +
  theme_void()  +
  labs(title = "Sankey Plot of psiType and tseType Frequencies")
#Sex chromosome recruitment on Micro chr that are already enriched with gonad genes
geneTissueType <- dfmerge %>% mutate(
  tissueTP = case_when(
    tseType != "Gonad" & psiType == "Gonad" ~ "Psi",
    tseType == "Gonad" & psiType != "Gonad" ~ "Tse",
    tseType == "Gonad" & psiType == "Gonad" ~ "Shared",
    TRUE ~ NA
  ),
  chr = geneRename$V1[match(psiID, geneRename$V4)]) %>% na.omit()
tissueFreq <- geneTissueType %>% group_by(tissueTP, chr) %>% 
  summarise(freq = n(), .groups = "drop") 
Var1 <- geneRename$V1[!grepl("HiC_scaffold_", geneRename$V1)]
all_freq <-  data.frame(table(Var1))
tissueRatioChr <- merge(tissueFreq, all_freq, by.x = "chr", by.y = "Var1") %>%
  rename(Tar = freq, All = Freq) %>% mutate(
    chr = factor(chr, levels = c(paste0("chr", 1:32), "chrZ", "chrW")),
    type = case_when(
      chr %in% c(paste0("chr", 1:8)) ~ "Macro", 
      chr == "chrW" ~ "chrW",
      chr == "chrZ" ~ "chrZ",
      TRUE ~ "Micro"
    ), 
    chrCol = case_when(
      type == "Macro" ~ "Macro",
      TRUE ~ "Micro"
    ),
    ratio = (Tar/All) *100
  ) 
ggplot(data = tissueRatioChr, aes(x= chrCol, y = ratio)) + 
  geom_boxplot() + geom_jitter(aes(color = type), size = 2, width = 0.2) + 
  facet_wrap( ~ tissueTP, scales = "free_y") + scale_color_manual(values = c("Macro" = "#737373", "Micro" = "#bdbdbd", "chrZ" = "#0000cd", "chrW" = "#cd0000")) +
  stat_compare_means(comparisons = list(c("Macro", "Micro")),
                     method="wilcox.test", method.args = list(alternative = "greater")) + theme_bw()
#Psi gonad-gain genes had higher dN/dS compared to gonad-shared & tse-gonad genes
dndsttGga <- read.table(file = "psi.tse.dnds.GgaOP.txt")
#Gonad genes
colnames(dndsttGga) <- c("pID", "dNpsi", "dSpsi", "dNtse", "dStse")
dndst2pByTissue <-  dndsttGga %>% filter(dSpsi >= 0.01 & dSpsi <= 2) %>% 
  filter(dStse >= 0.01 & dStse <= 2) %>% mutate(
    gene = geneRename$V4[match(pID, geneRename$V5)],
    pT = dfmerge$psiType[match(gene, dfmerge$psiID)],
    tT = dfmerge$tseType[match(gene, dfmerge$psiID)]
  ) %>% na.omit() %>%
  filter(tT != "Gonad" & pT == "Gonad") %>% mutate(
    wPsi = dNpsi/dSpsi,
    wTse = dNtse/dStse
  ) %>% select(wPsi, wTse) %>%
  pivot_longer(cols = c(wPsi, wTse), names_to = "Species", values_to = "w")
ggplot(dndst2pByTissue, aes(x = Species, y = w, fill = Species)) +
  geom_violin(trim = FALSE) +
  stat_summary(fun = mean, geom = "crossbar", 
               width = 0.2, color = "black", 
               fatten = 2) +
  stat_compare_means(method = "wilcox.test", 
                     comparisons = list(c("wPsi", "wTse"),
                                        method.args = list(alternative = "greater"))) +
  theme_test() + ylim(0,2)+
  labs(x = "Gonad Type", y = "Value", title = "Boxplot of wPsi and wTse by Gonad Type")
#TSE DE genes
dndst2pByTseDE <- dndst2pByTissue <-  dndsttGga %>% filter(dSpsi >= 0.01 & dSpsi <= 2) %>% 
  filter(dStse >= 0.01 & dStse <= 2) %>% filter(dNtse > 0 & dNpsi > 0) %>%  mutate(
    gene = geneRename$V4[match(pID, geneRename$V5)],
    tseID = orthologs$V2[match(gene, orthologs$V1)],
    tseType = tseDEgene$group[match(tseID, tseDEgene$gene)]) %>% na.omit() %>% mutate(
      wPsi = dNpsi/dSpsi,
      wTse = dNtse/dStse
    ) 
dfVisualDNDS <- dndst2pByTseDE %>% select(wPsi, wTse, tseType) %>%
  pivot_longer(cols = c(wPsi, wTse), names_to = "Species", values_to = "w")
ggplot(dfVisualDNDS, aes(x = Species, y = w, fill = Species)) +
  geom_violin(trim = FALSE) + facet_wrap( ~ tseType) +
  stat_summary(fun = mean, geom = "crossbar", 
               width = 0.2, color = "black", 
               fatten = 2) +
  stat_compare_means(method = "wilcox.test", 
                     comparisons = list(c("wPsi", "wTse"),
                                        method.args = list(alternative = "greater"))) +
  theme_test() + ylim(0,2)+
  labs(x = "Gonad Type", y = "Value", title = "Boxplot of wPsi and wTse by Gonad Type")
#What are these TSD lost gonad genes
psiGainGonad <- dfmerge %>% filter(tseType != "Gonad" & psiType == "Gonad")
tseLossGonad <- dfmerge %>% filter(tseType == "Gonad" & psiType != "Gonad")
library(clusterProfiler)
library(org.Hs.eg.db)
targetGenes <- bitr(unique(tseLossGonad$Product),
                    fromType = "SYMBOL",
                    toType = "ENTREZID",
                    OrgDb = "org.Hs.eg.db")
ego <- enrichGO(targetGenes$ENTREZID,
                keyType = "ENTREZID",
                OrgDb = "org.Hs.eg.db",
                ont = "BP",
                pvalueCutoff = 0.05,
                pAdjustMethod = "fdr",
                readable = T)
#write.table(ego@result,file = "turtleTSElostGonad.GO.txt")
tseLossGonadGO <- read.table(file = "turtleTSElostGonad.GO.txt")
tseLossGonadGO <- tseLossGonadGO %>% filter(p.adjust <= 0.05)
#Overlap with TSE pathway genes
CaAndHSPgenes <- read.table("TRP_Ca_HSP.txt",sep ="\t")
filtered_rows <- tseLossGonadGO[sapply(tseLossGonadGO$geneID, function(gene_ids) {
  any(sapply(CaAndHSPgenes$V1, function(gene) grepl(gene, gene_ids)))
}), ] %>% mutate(
  geneRatio = sapply(strsplit(as.character(GeneRatio), "/"), function(x) as.numeric(x[1]) / as.numeric(x[2])),
  Description = factor(Description, levels = Description[order(-p.adjust)])
)
ggplot(filtered_rows, aes(x = geneRatio, y = Description, size = Count, fill = -log10(p.adjust))) +
  geom_point(shape = 21) + scale_fill_gradient(low = "white", high = "#cd0000") + theme_bw()
TSDdfSUB <- dfmerge %>% filter(Product %in% CaAndHSPgenes$V1) %>% mutate(
  func = CaAndHSPgenes$V2[match(Product, CaAndHSPgenes$V1)]
)
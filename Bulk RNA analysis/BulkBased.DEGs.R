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
{
dat = read.table("redcounts.txt")
rownames(dat) = dat[,1]
dat = dat[,2*(1:30)]
nam = read.table("redname.txt")
colnames(dat) = nam$V1
for (m in 1:30) {
  dat[,m] = as.numeric(dat[,m])
}
TSPcount <- dat %>% select(matches("15|16|17|18|19|21"))
library(DESeq2)
sampleInfo <- data.frame(sp = colnames(TSPcount)) %>% mutate(
  condition = substr(sp, 1, 1)
) #Input
dds <- DESeqDataSetFromMatrix(countData = TSPcount,
                              colData = sampleInfo,
                              design = ~ condition) #Turn to deseq dat.
keep <- rowSums(counts(dds)) >= 6  #Filter low count genes
dds <- dds[keep,]
dds <- DESeq(dds)
res <- results(dds, contrast=c("condition","F","M"))
resLFC <- lfcShrink(dds, coef="condition_M_vs_F", type="apeglm") # lfcShrinking optional
resOrdered <- resLFC[order(resLFC$pvalue),]
res05 <- results(dds, alpha = 0.05)
#write.table(as.data.frame(resOrdered), file="tse.gonad.DEgene.txt")
tseDEgene <- read.table(file = "tse.gonad.DEgene.txt")
#Only all tseDEgenes volcano
library(ggrastr)
tseDEgene$gene <- rownames(tseDEgene)
tseDEgene$group <- "ns"
tseDEgene$group[tseDEgene$padj <= 0.05 & tseDEgene$log2FoldChange >= 1] <- "male"
tseDEgene$group[tseDEgene$padj <= 0.05 & tseDEgene$log2FoldChange <= -1] <- "female"
# volcano plot
tseDEgene$gene <- rownames(tseDEgene)
tseDEgene$group <- "ns"
tseDEgene$group[tseDEgene$padj <= 0.05 & tseDEgene$log2FoldChange >= 1] <- "male"
tseDEgene$group[tseDEgene$padj <= 0.05 & tseDEgene$log2FoldChange <= -1] <- "female"
ns_genes <- tseDEgene[tseDEgene$group == "ns", ]
sig_genes <- tseDEgene[tseDEgene$group != "ns", ]
p <- ggplot() +
  geom_point_rast(
    data = ns_genes,
    aes(x = log2FoldChange, y = -log10(padj)),
    fill = "#dcdcdc", color= "#dcdcdc", size = 0.5
  ) +
  geom_point(
    data = sig_genes,
    aes(x = log2FoldChange, y = -log10(padj), color = group),
    size = 1.2, alpha = 0.8
  ) +
  scale_color_manual(values = c("male" = "#8cb7df", "female" = "#ed7d61")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  theme_bw(base_size = 14) +
  labs(x = "log2(Fold Change)", y = "-log10(p-value)", color = "Group")
label_genes <- tseDEgene[tseDEgene$gene %in% GenesToshowSub1$tseID, ]
p + geom_text_repel(
  data = label_genes,
  aes(x = log2FoldChange, y = -log10(pvalue), label = gene, color = group),
  size = 3,
  box.padding = 0.4,
  point.padding = 0.3,
  max.overlaps = Inf,
  show.legend = FALSE
)
#How about the expr. in somatic tissues?
tseExp <- read.table("tse.tpm.mean.allTissue.txt", sep = ",", header = T, row.names = 1)
tseExpSoma_ratio <- tseExp[rownames(tseExp) %in% subset(sig_genes, group == "male")$gene, ] %>%
  select(grep("FB|MB|FM|MM", colnames(.))) %>%
  select(grep("12|15|16|18|21", colnames(.))) %>%
  filter_all(any_vars(. >= 1))
FB_cols <- grep("FB", colnames(tseExpSoma_ratio))
MB_cols <- grep("MB", colnames(tseExpSoma_ratio))
FM_cols <- grep("FM", colnames(tseExpSoma_ratio))
MM_cols <- grep("MM", colnames(tseExpSoma_ratio))
FB_MB_ratio <- (tseExpSoma_ratio[, FB_cols] + 0.1) / (tseExpSoma_ratio[, MB_cols])
FM_MM_ratio <- (tseExpSoma_ratio[, FM_cols] + 0.1) / (tseExpSoma_ratio[, MM_cols])
colnames(FB_MB_ratio) <- paste0("FB_MB_", c(12,15,16,18,21))
colnames(FM_MM_ratio) <- paste0("FM_MM_", c(12,15,16,18,21))
tseExp_ratio_all <- cbind(Gene = rownames(tseExpSoma_ratio), FB_MB_ratio, FM_MM_ratio)
rownames(tseExp_ratio_all) <- tseExp_ratio_all$Gene
classify_bias <- function(x) {
  ifelse(x > 2, "bias",
         ifelse(x < 0.5, "reverse", "non-bias"))
}
bias_class <- as.data.frame(apply(tseExp_ratio_all[,-1], 2, classify_bias)) 
} #Call DE genes in Tse gonad merged stages
{
pmalecount <- read.table(file="psi.male.count",header = T,row.names = 1) 
pfemalecount <- read.table(file="psi.female.count",header = T,row.names = 1)
pcount <- merge(pmalecount, pfemalecount, by = 0) 
rownames(pcount) <- pcount$Row.names
pcount <- pcount[,-1]
pSDcount <- pcount %>% select(matches('FG|MG')) %>% select(matches("14|15|16|18|21|24"))
library(DESeq2)
PsampleInfo <- data.frame(sp = colnames(pSDcount)) %>% mutate(
  condition = substr(sp, 1, 1)
) #Input
dds <- DESeqDataSetFromMatrix(countData = pSDcount,
                              colData = PsampleInfo,
                              design = ~ condition) #Turn to deseq dat.
keep <- rowSums(counts(dds)) >= 6  #Filter low count genes
dds <- dds[keep,]
dds <- DESeq(dds)
res <- results(dds, contrast=c("condition","F","M"))
resultsNames(dds)
resLFC <- lfcShrink(dds, coef="condition_M_vs_F", type="apeglm") # lfcShrinking optional
resOrdered <- resLFC[order(resLFC$pvalue),]
res05 <- results(dds, alpha = 0.05)
#write.table(as.data.frame(resOrdered), file="psi.gonad.DEgene.txt")
psiDEgene <- read.table(file = "psi.gonad.DEgene.txt")
psiDEgeneOrthSub <- psiDEgene %>% mutate(
  type = case_when(log2FoldChange >= 1 & padj <= 0.05 ~ "M", 
                   log2FoldChange <= -1 & padj <= 0.05 ~ "F",
                   TRUE ~ "E")
) %>% filter(rownames(.) %in% orthologs$V1)
#DEGs transitions of two turtles pieDonut
DEtype <- orthologs %>% mutate(
  tseType = tseDEgeneOrthSub$type[match(V2, rownames(tseDEgeneOrthSub))],
  psiType = psiDEgeneOrthSub$type[match(V1, rownames(psiDEgeneOrthSub))],
  GeneProduct = psiGeneProduct$product[match(V1, psiGeneProduct$V4)]
) %>% na.omit()
dfFreq <- DEtype %>% group_by(tseType, psiType) %>%
  summarise(freq = n(), .groups = "drop") %>% filter(tseType != "E")
library(webr)
PieDonut(dfFreq, aes(tseType, psiType, count=freq),  r0 = 0.8,title = "DEGs")
#
write.table(DEtype, file = "psitse.detype.txt")
#DEtypes of TSD genes
DEtype <- read.table("psitse.detype.txt")
DEtypeTSD <- DEtype %>% filter(GeneProduct %in% CaAndHSPgenes$V1) 
DEtypeTSDFreq <- DEtypeTSD %>% group_by(tseType, psiType) %>%
  summarise(freq = n(), .groups = "drop") 
PieDonut(DEtypeTSDFreq, aes(tseType, psiType, count=freq),  r0 = 0.8,title = "DEGs")
#W of sex-biased loss gene
dndsttGga <- read.table(file = "psi.tse.dnds.GgaOP.txt")
#
colnames(dndsttGga) <- c("pID", "dNpsi", "dSpsi", "dNtse", "dStse")
dndst2pByTissue <-  dndsttGga %>% filter(dSpsi > 0.01 & dSpsi < 2) %>% 
  filter(dStse > 0.01 & dStse < 2) %>% mutate(
    gene = geneRename$V4[match(pID, geneRename$V5)],
    pT = DEtype$psiType[match(gene, DEtype$V1)],
    tT = DEtype$tseType[match(gene, DEtype$V1)]
  ) %>% na.omit() %>%
  filter(tT == "E" & pT == "M") %>% mutate(
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
  theme_test() + ylim(0,1)+
  labs(x = "Gonad Type", y = "Value", title = "Boxplot of wPsi and wTse by Gonad Type")

} #Call DE genes in Psi gonad merged stages
{
dat <- read.table("redcounts.txt")
rownames(dat) <- dat[,1]
dat <- dat[, 2*(1:30)]
nam <- read.table("redname.txt")
colnames(dat) <- nam$V1
dat[] <- lapply(dat, as.numeric)
{
# By Stages 
stages <- c("12", "15", "16", "17", "18", "19", "21")
# Form list
DEG.list <- list()
for (stage in stages) {
  message("Running DESeq2 for stage ", stage, " ...")
  stage_pattern <- paste0("s", stage, "(_|$)")
  stage.counts <- dat %>% select(matches(stage_pattern))
  sampleInfo <- data.frame(sp = colnames(stage.counts)) %>%
    mutate(
      condition = substr(sp, 1, 1),                # F or M
      stage = gsub(".*s(\\d+).*", "\\1", sp)       # Stages
    )
  dds <- DESeqDataSetFromMatrix(countData = stage.counts,
                                colData = sampleInfo,
                                design = ~ condition)
  keep <- rowSums(counts(dds)) >= 6
  dds <- dds[keep, ]
  dds <- DESeq(dds)
  res <- results(dds, contrast = c("condition", "F", "M"))
  resLFC <- lfcShrink(dds, coef = "condition_M_vs_F", type = "apeglm")
  resOrdered <- resLFC[order(resLFC$pvalue), ]
  DEG.list[[stage]] <- as.data.frame(resOrdered)
} # Calculate
all_genes <- unique(unlist(lapply(DEG.list, rownames)))
bias_matrix <- matrix("non-bias",
                      nrow = length(all_genes),
                      ncol = length(DEG.list),
                      dimnames = list(all_genes, names(DEG.list)))
for (stage in names(DEG.list)) {
  t <- DEG.list[[stage]]
  t$bias <- "non-bias"
  t$bias[t$padj <= 0.05 & t$log2FoldChange >= 1] <- "male"
  t$bias[t$padj <= 0.05 & t$log2FoldChange <= -1] <- "female"
  bias_matrix[rownames(t), stage] <- t$bias
}
bias_df <- as.data.frame(bias_matrix)
#write.table(bias_df, file="tse.gonad.DEgene.byStage.txt") 
} # Tse gonad DEGs by stage
dat = read.table("tse.all.tissue.count.txt", header = T, row.names = 1)
datBr <- dat %>% select(matches("F_B|M_B")) 
datMm <- dat %>% select(matches("F_M|M_M")) 
{
new_colnamesBr <- sapply(colnames(datBr), function(x){
  m <- str_match(x, "Tse_(\\d+)(F|M)_B(\\d)")[,2:4]
  stage <- m[1]
  sex <- m[2]
  rep <- m[3]
  paste0(sex, "B", stage, "_", rep)
})
colnames(datBr) <- new_colnamesBr
datBr <- datBr[, order(
  substr(colnames(datBr), 1, 1),                             
  as.numeric(sub("^[FM]B?(\\d+)_\\d+", "\\1", colnames(datBr))),  
  as.numeric(sub("^[FM]B?\\d+_(\\d+)", "\\1", colnames(datBr)))  
)]
} #rename br
{
new_colnamesMm <- sapply(colnames(datMm), function(x){
  m <- str_match(x, "Tse_(\\d+)(F|M)_M(\\d)")[,2:4]
  stage <- m[1]
  sex <- m[2]
  rep <- m[3]
  paste0(sex, "M", stage, "_", rep)
})
colnames(datMm) <- new_colnamesMm
datMm <- datMm[, order(
  substr(colnames(datMm), 1, 1),                             
  as.numeric(sub("^[FM]M?(\\d+)_\\d+", "\\1", colnames(datMm))),  
  as.numeric(sub("^[FM]M?\\d+_(\\d+)", "\\1", colnames(datMm)))  
)]
} #rename mm
{
stages <- c("12", "15", "16", "18", "21")
DEG.list <- list()
for (stage in stages) {
  message("Running DESeq2 for stage ", stage, " ...")
    stage_pattern <- paste0(stage, "(_|$)")
  stage.counts <- datBr %>% select(matches(stage_pattern))
  sampleInfo <- data.frame(sp = colnames(stage.counts)) %>%
    mutate(
      condition = substr(sp, 1, 1),   # F or M
      stage = stage
    )
  dds <- DESeqDataSetFromMatrix(
    countData = round(stage.counts),
    colData = sampleInfo,
    design = ~ condition
  )
  keep <- rowSums(counts(dds)) >= 6
  dds <- dds[keep, ]
  dds <- DESeq(dds)
  res <- results(dds, contrast = c("condition", "F", "M"))
  resLFC <- lfcShrink(dds, coef = "condition_M_vs_F", type = "apeglm")
  resOrdered <- resLFC[order(resLFC$pvalue), ]
  DEG.list[[stage]] <- as.data.frame(resOrdered)
}
all_genes <- unique(unlist(lapply(DEG.list, rownames)))
bias_matrix <- matrix("non-bias",
                      nrow = length(all_genes),
                      ncol = length(DEG.list),
                      dimnames = list(all_genes, names(DEG.list)))

for (stage in names(DEG.list)) {
  t <- DEG.list[[stage]]
  t$bias <- "non-bias"
  t$bias[t$padj <= 0.05 & t$log2FoldChange >= 1] <- "male"
  t$bias[t$padj <= 0.05 & t$log2FoldChange <= -1] <- "female"
  bias_matrix[rownames(t), stage] <- t$bias
}
bias_df <- as.data.frame(bias_matrix)
write.table(bias_df, file="tse.brain.DEgene.byStage.txt") 
} # Tse brain DEGs by stage
{
  stages <- c("12", "15", "16", "18", "21")
  DEG.list <- list()
  for (stage in stages) {
    message("Running DESeq2 for stage ", stage, " ...")
    stage_pattern <- paste0(stage, "(_|$)")
    stage.counts <- datMm %>% select(matches(stage_pattern))
    sampleInfo <- data.frame(sp = colnames(stage.counts)) %>%
      mutate(
        condition = substr(sp, 1, 1),   # F or M
        stage = stage
      )
    dds <- DESeqDataSetFromMatrix(
      countData = round(stage.counts),
      colData = sampleInfo,
      design = ~ condition
    )
    keep <- rowSums(counts(dds)) >= 6
    dds <- dds[keep, ]
    dds <- DESeq(dds)
    res <- results(dds, contrast = c("condition", "F", "M"))
    resLFC <- lfcShrink(dds, coef = "condition_M_vs_F", type = "apeglm")
    resOrdered <- resLFC[order(resLFC$pvalue), ]
    DEG.list[[stage]] <- as.data.frame(resOrdered)
  }
  all_genes <- unique(unlist(lapply(DEG.list, rownames)))
  bias_matrix <- matrix("non-bias",
                        nrow = length(all_genes),
                        ncol = length(DEG.list),
                        dimnames = list(all_genes, names(DEG.list)))
  
  for (stage in names(DEG.list)) {
    t <- DEG.list[[stage]]
    t$bias <- "non-bias"
    t$bias[t$padj <= 0.05 & t$log2FoldChange >= 1] <- "male"
    t$bias[t$padj <= 0.05 & t$log2FoldChange <= -1] <- "female"
    bias_matrix[rownames(t), stage] <- t$bias
  }
  bias_df <- as.data.frame(bias_matrix)
  write.table(bias_df, file="tse.meso.DEgene.byStage.txt") 
} # Tse mesonephros DEGs by stage


# Visualize
dfStGd <- read.table(file = "tse.gonad.DEgene.byStage.txt",header=TRUE, check.names=FALSE)
dfStBr <- read.table(file = "tse.brain.DEgene.byStage.txt",header=TRUE, check.names=FALSE)
dfStMm <- read.table(file = "tse.meso.DEgene.byStage.txt",header=TRUE, check.names=FALSE)

df_list <- list(dfStGd, dfStBr, dfStMm)

Reduce(merge, lapply(l, function(x) data.frame(x, rn = row.names(x))))



deg_summary <- dfStBr %>%
  tidyr::pivot_longer(cols = everything(),
                      names_to = "stage",
                      values_to = "bias") %>%
  group_by(stage, bias) %>%
  summarise(count = n(), .groups = "drop") %>%
  filter(bias %in% c("male", "female"))
deg_summary$stage <- factor(deg_summary$stage,
                            levels = c("12","15","16","17","18","19","21")) 
deg_summary <- deg_summary %>% filter(stage %in% c("12","15","16","18","21"))
# 
p2 <- ggplot(deg_summary, aes(x = stage, y = count, fill = bias)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  scale_fill_manual(values = c("female" = "#ED7D61", "male" = "#8CB7DF")) +
  labs(
    x = "Stage",
    y = "Number of sex-biased DEGs",
    fill = "Bias type",
    title = "Sex-biased DEGs across developmental stages"
  ) +
  theme_classic() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
ggarrange(p1,p2,p3, ncol = 1, align ="h",common.legend = T)
} #Call DE genes of Tse 12-21 (tissue x sex x stage)
{
pmalecount <- read.table("psi.male.count", header = TRUE, row.names = 1)
pfemalecount <- read.table("psi.female.count", header = TRUE, row.names = 1)
pcount <- merge(pmalecount, pfemalecount, by = 0)
rownames(pcount) <- pcount$Row.names
pcount <- pcount[, -1]
# Select stages and tissues
pSDcount <- pcount %>%
  select(matches("FG|MG")) %>%
  select(matches("12|13|14|15|16|18|21|24"))
stages <- c("12","13","14","15","16","18","21","24")
# Calculate DEGs by stages
DEG.list <- list()
for (stage in stages) {
  message("Running DESeq2 for stage ", stage, " ...")
  stage_pattern <- paste0("(MC)?", stage, "(_|$)")
  stage.counts <- pSDcount %>% select(matches(stage_pattern))
  # Sample list
  PsampleInfo <- data.frame(sp = colnames(stage.counts)) %>%
    mutate(
      condition = substr(sp, 1, 1) # F/M
    )
  # form DESeq object
  dds <- DESeqDataSetFromMatrix(countData = stage.counts,
                                colData = PsampleInfo,
                                design = ~ condition)
  # Low expr. filt
  keep <- rowSums(counts(dds)) >= 6
  dds <- dds[keep, ]
  dds <- DESeq(dds)
  # F vs M
  res <- results(dds, contrast = c("condition", "F", "M"))
  # optional LFC shrink
  resLFC <- lfcShrink(dds, coef = "condition_M_vs_F", type = "apeglm")
  # sort and save
  resOrdered <- resLFC[order(resLFC$pvalue), ]
  DEG.list[[stage]] <- as.data.frame(resOrdered)
}
# Define male/female/non-bias
all_genes <- unique(unlist(lapply(DEG.list, rownames)))
bias_matrix <- matrix("non-bias",
                      nrow = length(all_genes),
                      ncol = length(DEG.list),
                      dimnames = list(all_genes, names(DEG.list)))
for (stage in names(DEG.list)) {
  t <- DEG.list[[stage]]
  t$bias <- "non-bias"
  t$bias[t$padj <= 0.05 & t$log2FoldChange >= 1] <- "male"
  t$bias[t$padj <= 0.05 & t$log2FoldChange <= -1] <- "female"
  bias_matrix[rownames(t), stage] <- t$bias
}
bias_df <- as.data.frame(bias_matrix)
#write.table(bias_df, file="psi.gonad.DEgene.byStage.txt")
#
dfStGd <- read.table(file = "psi.gonad.DEgene.byStage.txt",header=TRUE, check.names=FALSE)
dfStBr <- read.table(file = "psi.brain.DEgene.byStage.txt",header=TRUE, check.names=FALSE)
deg_summary <- dfStGd %>%
  tidyr::pivot_longer(cols = everything(),
                      names_to = "stage",
                      values_to = "bias") %>%
  group_by(stage, bias) %>%
  summarise(count = n(), .groups = "drop") %>%
  filter(bias %in% c("male", "female"))
deg_summary$stage <- factor(deg_summary$stage,
                            levels = c("12","13","14","15","16","18","21","24")) 
deg_summary <- deg_summary %>% filter(stage %in% c("12","15","16","18","21"))

# 
p1 <- ggplot(deg_summary, aes(x = stage, y = count, fill = bias)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  scale_fill_manual(values = c("female" = "#ED7D61", "male" = "#8CB7DF")) +
  labs(
    x = "Stage",
    y = "Number of sex-biased DEGs",
    fill = "Bias type",
    title = "Sex-biased DEGs across developmental stages"
  ) +
  theme_classic() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
ggarrange(p1,p2, ncol = 1, align ="h",common.legend = T)

} #Call DE genes of Psi 12-21 (tissue x sex x stage)
{
psiDEGgd <- read.table(file = "psi.gonad.DEgene.byStage.txt",header=TRUE, check.names=FALSE)
tseDEGgd <- read.table(file = "tse.gonad.DEgene.byStage.txt",header=TRUE, check.names=FALSE)
orthologs <- read.table("psi-red.orth-1to1.50")
# Only overlap stages
psi_sel <- psiDEGgd[, c("12", "15", "16", "18", "21")]
tse_sel <- tseDEGgd[, c("12", "15", "16", "18", "21")]
orth_clean <- orthologs %>%
  filter(!is.na(V1), !is.na(V2)) %>%
  distinct() %>%
  filter(V1 %in% rownames(psi_sel), V2 %in% rownames(tse_sel))
# Define first DE stage
get_first_DE_stage <- function(x) {
  stages <- colnames(x)
  first_stage <- apply(x, 1, function(v) {
    de_stages <- which(v != "non-bias")
    if(length(de_stages)==0) return(NA)
    return(as.numeric(stages[min(de_stages)]))
  })
  first_bias <- apply(x, 1, function(v) {
    de_stages <- which(v != "non-bias")
    if(length(de_stages)==0) return(NA)
    return(v[min(de_stages)])
  })
  return(list(stage = first_stage, bias = first_bias))
}

psi_first <- get_first_DE_stage(psi_sel)
tse_first <- get_first_DE_stage(tse_sel)
# Merge info.
merged <- orth_clean %>%
  rename(V1 = "psi", V2 = "tse") %>%
  mutate(
    psi_first_stage = psi_first$stage[psi],
    psi_first_bias  = psi_first$bias[psi],
    tse_first_stage = tse_first$stage[tse],
    tse_first_bias  = tse_first$bias[tse]
  )
result <- merged %>%
  filter(!is.na(psi_first_stage), !is.na(tse_first_stage),
         psi_first_stage >= 16,
         tse_first_stage <= 15,
         !is.na(psi_first_bias), !is.na(tse_first_bias),
         psi_first_bias == tse_first_bias) %>%
  mutate(first_bias = psi_first_bias) %>%
  select(psi, tse,
         psi_first_stage, psi_first_bias,
         tse_first_stage, tse_first_bias,
         first_bias)
#write.table(result,file = "tse15DEbutPsiLaterDE.gene.txt")
} #Do Tse had more DEGs than Psi at early stages

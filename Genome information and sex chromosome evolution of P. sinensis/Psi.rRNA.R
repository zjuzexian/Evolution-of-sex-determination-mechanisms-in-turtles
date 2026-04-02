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
#rDNA amplification on W
#5S rDNA
{
#rRNA expression
all <- read.table(file = "fall.txt",header = T) 
fsub <- all  %>%  dplyr::select(1,7:50)
rownames(fsub) <- fsub[,1]
fsub <- fsub[,-1]
#get mean of 2 reps
for (i in seq(from=1, to=44, by=2)) {
  new <- rowMeans(fsub %>% dplyr::select(i,i+1))
  fsub[,ncol(fsub) + 1] <- new
  colnames(fsub)[ncol(fsub)] <- paste0("new", i)
}
tpmmean <- fsub[45:66]
colnames(tpmmean) <- gsub("_1.filt.bam","",colnames(fsub %>% dplyr::select(1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,
                                                                           33,35,37,39,41,43)))
library(pheatmap)
pheatmap(tpmmean,cluster_rows = T, cluster_cols = F,
         cellwidth = 10, cellheight = 10, fontsize = 10,
         color = c(colorRampPalette(c("#8cb7df", "white"))(50),colorRampPalette(c("white", "#ed7d61"))(50)))  

} 
#45S rDNA
{
library(IRanges)
all <- read.table(file = "reeves.rRNA.loci")
z <- subset(all,all$V1=="chrZ")
w <- subset(all,all$V1=="chrW")
#Get signals by length
z <- z %>% filter(V4 > .8) %>% 
  mutate(start = pmin(V2, V3), end = pmax(V2, V3))
merged <- z %>%
  group_by(V1) %>%
  group_modify(~ {
    ir <- IRanges(start = .x$start, end = .x$end)
    reduced <- reduce(ir)
    # Find overlaps
    hits <- findOverlaps(reduced, ir)
    merged_data <- tibble(
      V2 = start(reduced),
      V3 = end(reduced)
    )
    # Combind V5 and connect
    merged_data$V5 <- tapply(.x$V5[subjectHits(hits)], queryHits(hits), function(x) paste(unique(x), collapse = ","))
    merged_data
  }) %>%
  ungroup()
merged$V1 <- z$V1[1] 
ggplot(merged,aes(x=V2/1e6, y = (V3-V2)))+
  geom_bar(stat = "identity") + theme_classic() + labs(x="", y="Length (bp)")
#Count input
fcount <- read.table(file = "rdna.f.count.txt",header = T) 
mcount <- read.table(file = "rdna.m.count.txt",header = T) 
#TPM
tarsub <- mcount  %>%  dplyr::select(1,6:50)
gene_lengths <- tarsub$Length 
count_matrix <- as.matrix(tarsub[, 3:ncol(tarsub)]) 
rpk <- count_matrix / (gene_lengths / 1000)
scaling_factors <- colSums(rpk)
tpm <- t( t(rpk) / scaling_factors ) * 1e6
tpm_df <- cbind(Geneid = tarsub$Geneid, as.data.frame(tpm))
rownames(tpm_df) <- tpm_df$Geneid
tpm_df <- tpm_df[, -1]
#get mean of 2 reps
for (i in seq(from=1, to=44, by=2)) {
  new <- rowMeans(tpm_df %>% dplyr::select(i,i+1))
  tpm_df[,ncol(tpm_df) + 1] <- new
  colnames(tpm_df)[ncol(tpm_df)] <- paste0("new", i)
}
tpmmean <- tpm_df[45:66]
colnames(tpmmean) <- gsub("_1.filt.bam","",colnames(tpm_df %>% dplyr::select(1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,
                                                                           33,35,37,39,41,43)))
ftpm <- tpmmean
mtpm <- tpmmean
ftpm$Geneid <- rownames(ftpm)
mtpm$Geneid <- rownames(mtpm)
rDNAtpm <- merge(ftpm, mtpm, by = "Geneid", all = TRUE)
rownames(rDNAtpm) <- rDNAtpm$Geneid
rDNAtpm$Geneid <- NULL
#write.table(rDNAtpm, file = "psi.45SrDNA.tpm.txt")
#From sRNA mapping
srnacount <- read.table(file="psi.nor.sRNA.count.txt",header = T) 
srnacount <- read.table(file="psi.nor.sRNA.count.uniq.txt",header = T) 
tarsub <- srnacount  %>%  dplyr::select(1,6:94)
gene_lengths <- tarsub$Length 
count_matrix <- as.matrix(tarsub[, 3:ncol(tarsub)]) 
rpk <- count_matrix / (gene_lengths / 1000)
scaling_factors <- colSums(rpk)
tpm <- t( t(rpk) / scaling_factors ) * 1e6
tpm_df <- cbind(Geneid = tarsub$Geneid, as.data.frame(tpm))
rownames(tpm_df) <- tpm_df$Geneid
tpm_df <- tpm_df[, -1]
colnames(tpm_df) <- c(
  "FB12_1", "FB12_2", "FB13_1", "FB13_2", "FB14_1", "FB14_2", "FB15_1", "FB15_2",
  "FB16_1", "FB16_2", "FB18_1", "FB18_2", "FB21_1", "FB21_2", "FB24_1", "FB24_2",
  "FG12_1", "FG12_2", "FG14_1", "FG14_2", "FG15_1", "FG15_2", "FG16_1", "FG16_2",
  "FG18_1", "FG18_2", "FG21_1", "FG21_2", "FG24_1", "FG24_2", "FG13_1", "FG13_2",
  "FM14_1", "FM14_2", "FM15_1", "FM15_2", "FM16_1", "FM16_2", "FM18_1", "FM18_2",
  "FM21_1", "FM21_2", "FM24_1", "FM24_2", "MB12_1", "MB12_2", "MB13_1", "MB13_2",
  "MB14_1", "MB14_2", "MB15_1", "MB15_2", "MB16_1", "MB16_2", "MB18_1", "MB18_2",
  "MB21_1", "MB21_2", "MB24_1", "MB24_2", "MG12_1", "MG12_2", "MG14_1", "MG14_2",
  "MG15_1", "MG15_2", "MG16_1", "MG16_2", "MG18_1", "MG18_2", "MG21_1", "MG21_2",
  "MG24_1", "MG24_2", "MG13_1", "MG13_2", "MM14_1", "MM14_2", "MM15_1", "MM15_2",
  "MM16_1", "MM16_2", "MM18_1", "MM18_2", "MM21_1", "MM21_2", "MM24_1", "MM24_2"
)
tpm_df <- tpm_df %>% dplyr::select(1:16,45:60,17:18,31:32,19:30,61:62,75:76,63:74,33:44,77:88)
for (i in seq(from=1, to=88, by=2)) {
  new <- rowMeans(tpm_df %>% dplyr::select(i,i+1))
  tpm_df[,ncol(tpm_df) + 1] <- new
  colnames(tpm_df)[ncol(tpm_df)] <- paste0("new", i)
}
tpmmean <- tpm_df[89:132]
colnames(tpmmean) <- gsub("_1","",colnames(tpm_df %>% dplyr::select(1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,
                                                                     33,35,37,39,41,43,45,47,49,51,53,55,57,59,61,
                                                                     63,65,67,69,71,73,75,77,79,81,83,85,87)))
#Heatmap
rDNAtpm <- read.table(file = "psi.45SrDNA.sRNA.tpm.txt")
rDNAtpm <- read.table(file = "psi.45SrDNA.sRNA.tpm.uniq.txt")
rDNAtpmsub <- rDNAtpm %>%
  select(which(grepl("FG|MG", colnames(.)))) %>%
  filter_all(.,any_vars(.>=1)) 
annoRow <- fcount %>% select(1,2) %>% mutate(
  type = ifelse( Chr %in% c("chrW", "chrZ"), Chr, "Autosome")
) 
rownames(annoRow) <- annoRow$Geneid
annoRow <- annoRow %>% select(type) 
library(pheatmap)
pheatmap(rDNAtpm %>% filter(rownames(.) %in% paste(rep("rDNA",22),c(1:22), sep = "_")) %>%
           select(grep("F", colnames(.))) %>%
           filter_all(.,any_vars(.>=1)) ,
         cluster_rows = T, cluster_cols = F, annotation_row = annoRow, scale = "row",
         color = c(colorRampPalette(c("#8cb7df", "white"))(50),colorRampPalette(c("white", "#ed7d61"))(50)))  
#For boxplot
rDNAtpm$Geneid <- rownames(rDNAtpm)
dfForboxplot <- melt(rDNAtpm, id.vars = "Geneid", variable.name = "key", value.name = "TPM")
dfForboxplot <- dfForboxplot %>% mutate(
  sex = gsub("[0-9].*", "", key),
  stage = gsub("^[A-Z]+", "", key),
  chr = fcount$Chr[match(Geneid, fcount$Geneid)],
  st = fcount$Start[match(Geneid, fcount$Geneid)],
  ed = fcount$End[match(Geneid, fcount$Geneid)]
) %>% filter(sex %in% c("FG", "MG", "FGMC", "MGMC")) %>% 
  mutate(sex = gsub("G.*","", sex)) %>% filter(TPM >1)
ggplot(dfForboxplot, aes(x = chr, y = log2(TPM+1), fill = sex)) +
  geom_violin(trim = FALSE) + theme_bw() + 
  stat_summary(fun = mean, geom = "crossbar", 
               width = 0.2, color = "black", 
               fatten = 2) 
}
#R2 location
{
R2w <- read.table(text = "chrW	3505823	3505918
chrW	3505972	3509086
chrW	7059112	7061034
chrW	7065309	7066620
chrW	14776212	14779326
chrW	14779380	14779478
chrW	14791606	14794653
chrW	14794707	14794804", header = FALSE)
R2z <- read.table(text = "chrZ	37072359	37073787
chrZ	37073859	37073961
chrZ	37106081	37108336
chrZ	39195597	39195651
chrZ	39195720	39195823
chrZ	39204921	39205016
chrZ	39205074	39208189", header = FALSE)
p1 <- ggplot(R2w, aes(x = (V2+V3)/2000000,y=abs(V3-V2))) +geom_histogram(stat = "identity")+xlim(0,15)
p2 <- ggplot(R2z, aes(x = (V2+V3)/2000000,y=abs(V3-V2))) +geom_histogram(stat = "identity")+xlim(0,40)
p1+p2


}
#R2 expression
{
TEexp <- read.table(file="psiTE.loci.tpmMean.txt",sep = ",",
                  row.names = 1,header = T) %>% dplyr::select(1:8,23:30,15:16,9:14,37:38,31:36,17:22,39:44)
rownames(TEexp) <- gsub("linear;","linear",rownames(TEexp))
tetype <- read.table(file="TEtype.txt")
R2sub <- tetype %>% filter(grepl("LINE/R2", V2)) %>% filter(!grepl("NeSL", V2)) 
R2exp <- TEexp[rownames(TEexp) %in% R2sub$V1,] %>%
  select(which(grepl("FG|MG", colnames(.)))) %>%  
  filter_all(.,any_vars(.> 0)) 
pheatmap(R2exp,cluster_rows = T, cluster_cols = F,  scale = "row",
         color = c(colorRampPalette(c("#8cb7df", "white"))(50),colorRampPalette(c("white", "#ed7d61"))(50)))  
}
#R2 kimura
{
tarKim <- read.table(file = "psi.tar.kimura.txt")
dfKim <- tarKim %>% mutate(
  chr = ifelse(V1 %in% c("chrZ", "chrW"), V1, "Autosome")
)  %>% filter(V5=="LINE/R2") 
dfKim$kimura <- as.numeric(gsub("Kimura", "", dfKim$V6)) 
dfKim <- dfKim %>% filter(kimura > 0)
ggplot(dfKim, aes(x = chr, y = kimura, fill = chr)) +
  geom_violin(trim = FALSE) + theme_bw() + 
  stat_summary(fun = mean, geom = "crossbar", 
               width = 0.2, color = "white", 
               fatten = 1)  +theme_test() +
  stat_compare_means(
    aes(group = chr),
    method = "wilcox.test",
    ref.group = "Autosome",
    method.args = list(alternative = "less"),
    label = "p.format"  
  ) + ylim(0, 60)
}
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
orthologs <- read.table("psi-red.orth-1to1.50")
expTse <- read.table(file = "tse.tpm.mean.allTissue.txt", header = T,sep = ',',row.names = 1)
expPsi <- read.table(file = "psi.tpmmean.NewFilt.txt", header = T,sep = ',',row.names = 1)
expPsi <- expPsi %>% select(15:16, 9:14, 37:38, 31:36, 1:8, 23:30, 17:22, 39:44)
colnames(expPsi) <- gsub("GMC", "G", colnames(expPsi)) 
geneTse <- c("ZNRF3", "TBX1")
genePsi <- c("TBX1A_ENSCPBP00000014334-D1", "ZNRF3_ENSCPBP00000039495-D1", "ZNRF3_ENSCPBP00000039495-D2")
expPsiSub <- expPsi[rownames(expPsi) %in% genePsi,]
expTseSub <- expTse[rownames(expTse) %in% geneTse,]
expPsiSub_longer <- expPsiSub %>% 
  rownames_to_column(var = "gene_id") %>% 
  tidyr::gather(key = "key", value = "value", -gene_id) %>%  
  mutate(
    stage = substr(key, 3, 4),
    sex = substr(key, 1, 1),
    tissue = substr(key, 2, 2),
    species = "Psi"
  )
ggplot(data = expPsiSub_longer, aes(x = stage, y = value, group = sex, color = sex)) +
  geom_point() + geom_line() + theme_test() +
  facet_wrap(~ gene_id + tissue, scales = "free_y") +
  scale_color_manual(values = c("F" = "#ed7d61", "M" = "#8cb7df"))
expTseSub_longer <- expTseSub %>% 
  rownames_to_column(var = "gene_id") %>% 
  tidyr::gather(key = "key", value = "value", -gene_id) %>%  
  mutate(
    stage = substr(key, 3, 4),
    sex = substr(key, 1, 1),
    tissue = substr(key, 2, 2),
    species = "Tse"
  ) 
ggplot(data = expTseSub_longer, aes(x = stage, y = value, group = sex, color = sex)) +
  geom_point() + geom_line() + theme_bw() +
  facet_wrap(~ gene_id + tissue, scales = "free_y", ncol = 2) +
  scale_color_manual(values = c("F" = "#ed7d61", "M" = "#8cb7df"))
} #Fast check of Expression of TBX1/ZNRF3 in turtles
{
expMir <- read.table(file="psi.sRNA.tpmmean.final.shortstack.txt") #data used
#Z or W fitering
targetexpZ <- expMir[grep("chrZ", rownames(expMir)), ] %>%
  select(all_of(colnames(.)[grepl("FG|MG", colnames(.))])) %>% filter_all(.,any_vars(.>=1))
targetexpW <- expMir[grep("chrW",rownames(expMir)),] %>%
  select(all_of(colnames(.)[grepl("FG", colnames(.))])) %>% filter_all(.,any_vars(.>=1))
#Heat map of all expressed w-linked micro RNAs
library(pheatmap)
pheatmap(log2(targetexpW+1),cluster_rows = T,cluster_cols = F,
         color = c(colorRampPalette(c("#8cb7df", "white"))(50),
                   colorRampPalette(c("white", "#ed7d61"))(50)),
         border_color= "white")
#Who target ZNRF3-W
wtarlist <- read.table('psi.wmir.target.txt')
wtarlistZNRF3 <- wtarlist %>%
  filter(V5 == "ZNRF3") %>%
  select(V2) %>%
  mutate(chr = str_extract(V2, "^[^_]+_[^_]+")) %>%
  pull(chr) %>%
  unique()
candidateMir <- c(intersect(rownames(targetexpW), wtarlistZNRF3), "chrW_42977")
candidateMirExp <- expMir[rownames(expMir) %in% candidateMir,] %>% 
  rownames_to_column(var = "gene_id") %>% 
  tidyr::gather(key = "key", value = "value", -gene_id) %>%  
  mutate(
    stage = substr(key, 3, 4),
    sex = substr(key, 1, 1),
    tissue = substr(key, 2, 2),
    species = "Psi"
  ) %>% filter(tissue == "G") %>% filter(sex == "F")
ggplot(candidateMirExp, aes(x = stage, y = value, color = sex, group = sex)) + 
  geom_line() + geom_point() + theme_test() + facet_wrap(~ gene_id, scale = "free_y")
#Get brain tissue dat.
znrf3Tar <- c("chrW_42977", "chrW_42949","chrW_43028")
znrf3TarExp <- expMir[rownames(expMir) %in% znrf3Tar,] %>% 
  rownames_to_column(var = "gene_id") %>% 
  tidyr::gather(key = "key", value = "value", -gene_id) %>%  
  mutate(
    stage = substr(key, 3, 4),
    sex = substr(key, 1, 1),
    tissue = substr(key, 2, 2),
    species = "Psi"
  ) %>% filter(tissue == "B") %>% filter(sex == "F")
ggplot(znrf3TarExp, aes(x = stage, y = value, color = sex, group = sex)) + 
  geom_line() + geom_point() + theme_test() + facet_wrap(~ gene_id, scale = "free_y")
#Expression of their Z counterpart 
znrf3TarZ <- c("chrZ_43834")
znrf3TarZExp <- expMir[rownames(expMir) %in% znrf3TarZ,] %>% 
  rownames_to_column(var = "gene_id") %>% 
  tidyr::gather(key = "key", value = "value", -gene_id) %>%  
  mutate(
    stage = substr(key, 3, 4),
    sex = substr(key, 1, 1),
    tissue = substr(key, 2, 2),
    species = "Psi"
  ) %>% filter(tissue == "G") 
ggplot(znrf3TarZExp, aes(x = stage, y = value, color = sex, group = sex)) + 
  geom_line() + geom_point() + theme_test() + facet_wrap(~ gene_id, scale = "free_y")
#Exp of AutoMir target ZNRF3
znrf3TarA <- read.table("mir.tarZNRF3.txt")
znrf3TarAexp <- expMir[rownames(expMir) %in% znrf3TarA$V1,] %>% 
  filter_all(.,any_vars(.>=1))
pheatmap(znrf3TarAexp,cluster_rows = T,cluster_cols = F, scale = "row",
         color = c(colorRampPalette(c("#8cb7df", "white"))(50),
                   colorRampPalette(c("white", "#ed7d61"))(50)),
         border_color= "white") 
#Who target KDM6B
wtarlistKDM6B <- wtarlist %>%
  filter(V5 == "KDM6B") %>%
  select(V2) %>%
  mutate(chr = str_extract(V2, "^[^_]+_[^_]+")) %>%
  pull(chr) %>%
  unique()
candidateMir <- c(intersect(rownames(targetexpW), wtarlistKDM6B), "chrW_42977") 
candidateMirExp <- expMir[rownames(expMir) %in% candidateMir,] %>% 
  rownames_to_column(var = "gene_id") %>% 
  tidyr::gather(key = "key", value = "value", -gene_id) %>%  
  mutate(
    stage = substr(key, 3, 4),
    sex = substr(key, 1, 1),
    tissue = substr(key, 2, 2),
    species = "Psi"
  ) %>% filter(tissue == "G") %>% filter(sex == "F")
ggplot(candidateMirExp, aes(x = stage, y = value, color = sex, group = sex)) + 
  geom_line() + geom_point() + theme_test() + facet_wrap(~ gene_id, scale = "free_y")
  } #Fast check of Expression of micro RNA
{
psiExp <- read.table(file = "psi.tpmmean.NewFilt.txt", sep=',', row.names=1, header=T) 
tseExp <- read.table("tse.tpm.mean.allTissue.txt", header = T, row.names = 1)
rownames(tseExp) <- tseExp$Gene
tseExp <- tseExp[, -1]
geneidproduct <- read.table("psi.geneID.finalinfo.txt")
orthologs <- read.table("psi-red.orth-1to1.50") } #Data input for gene check
{
psiOnlyExpr <- psiExp %>%
  dplyr::select(grep("FG|MG", colnames(.))) %>% 
  dplyr::select(7:8, 1:6, 15:16, 9:14)
psiOnlyExpr <- psiExp %>%
  dplyr::select(grep("FB|MB|FG|MG|FM|MM", colnames(.)))
psiOnlyExpr$MM12 <- psiOnlyExpr$MGMC12
psiOnlyExpr$MM13 <- psiOnlyExpr$MGMC13
psiOnlyExpr$FM12 <- psiOnlyExpr$FGMC12
psiOnlyExpr$FM13 <- psiOnlyExpr$FGMC13
psiOnlyExpr_longer <- psiOnlyExpr %>% tibble::rownames_to_column("gene") %>% 
  gather(key = "variables", value = "tpm",- gene) %>% mutate(variables =  gsub("GMC","G",variables)) %>% mutate(
    sex = substr(variables, 1, 1),   
    tissue = substr(variables, 2, 2),  
    stage = as.numeric(substr(variables, 3, 4)),
    chr = geneidproduct$V1[match(gene, geneidproduct$V4)],
    factor = paste(sex, tissue)) 
tarGene <- psiOnlyExpr_longer %>% filter(gene == "DHB1_ENSCPBP00000022668-D1") # Select one gene
tarGene_ratio <- tarGene %>%
  group_by(gene, tissue, stage, chr) %>%
  summarise(
    tpm_F = mean(tpm[sex == "F"], na.rm = TRUE),
    tpm_M = mean(tpm[sex == "M"], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    tpmR = tpm_F / tpm_M
  ) 
ovlpstage <- tarGene %>% filter(stage %in% c("12","15","16","18","21"))
ggplot(ovlpstage, aes(x = as.character(stage), y = tpm, group = factor, color = sex)) +
  geom_line(size = 1)  + facet_wrap( ~ tissue) + scale_color_manual(values = c("F" = "#ed7d61", "M" = "#8cb7df")) +
  theme_classic() + labs(title = geneidproduct$V6[match(tarGene$gene, geneidproduct$V4)],
                         x="",y="TPM")
ggplot(tarGene_ratio %>% filter(stage %in% c("12","15","16","18","21")) %>% filter(tissue =="G"), aes(x = as.character(stage), y = log2(tpmR), group = gene)) +
  geom_line(size = 1)  + geom_point(size =2) +geom_hline(yintercept = 1)+
  facet_wrap( ~ gene, scales="free_y")  + scale_color_manual(values = c("F" = "#ed7d61", "M" = "#8cb7df")) +
  theme_classic() + labs(title = geneidproduct$V6[match(tarGene$gene, geneidproduct$V4)],
                         x="",y="TPM") 
subset1 <- tarGene_ratio %>% filter(stage %in% c("12","15","16","18","21")) %>% filter(tissue =="G") %>% mutate(species = "Psi") } #Expression in Psi no normalization
{
tseExpsub <- tseExp[rownames(tseExp) %in% orthologs$V2, ] %>%
  dplyr::select(grep("FG|MG|FB|MB|FM|MM", colnames(.))) %>% 
  dplyr::select(grep("12|15|16|18|21", colnames(.)))
tseExpsub$orthname <- orthologs$V1[match(row.names(tseExpsub),orthologs$V2)]
tseExpsub <- tseExpsub[tseExpsub$orthname %in% orthologs$V1,]
row.names(tseExpsub) <- tseExpsub$orthname
tseExpsub <- tseExpsub[,-31]
tseOnlyExpr <- tseExpsub %>%
  dplyr::select(grep("FG|MG", colnames(.)))
colnames(tseOnlyExpr) <- c(paste("tseFG",c("12","15","16","18","21"), sep = ""),
                           paste("tseMG",c("12","15","16","18","21"), sep = ""))
tseOnlyExpr_longer <- tseOnlyExpr %>% tibble::rownames_to_column("gene") %>% 
  gather(key = "variables", value = "tpm",- gene) %>% mutate(
    species = substr(variables, 1, 3), 
    sex = substr(variables, 4, 4),   
    tissue = substr(variables, 5, 5),  
    stage = as.numeric(substr(variables, 6, 7)),
    factor = paste(species, sex)) 
tseExpr_longer <- tseExpsub %>% tibble::rownames_to_column("gene") %>% 
  gather(key = "variables", value = "tpm",- gene) %>% mutate(
    sex = substr(variables, 1, 1),   
    tissue = substr(variables, 2, 2),  
    stage = as.numeric(substr(variables, 3, 4)),
    chr = geneidproduct$V1[match(gene, geneidproduct$V4)],
    factor = paste(sex, tissue)) 
tarGene <- tseExpr_longer %>% filter(gene == "DHB1_ENSCPBP00000022668-D1") # Select one gene
tarGene_ratio <- tarGene %>%
  group_by(gene, tissue, stage, chr) %>%
  summarise(
    tpm_F = mean(tpm[sex == "F"], na.rm = TRUE),
    tpm_M = mean(tpm[sex == "M"], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    tpmR = tpm_F / tpm_M
  ) 
ggplot(tarGene, aes(x = as.character(stage), y = tpm, group = factor, color = sex)) +
  geom_line(size = 1)  + facet_wrap( ~ tissue) + scale_color_manual(values = c("F" = "#ed7d61", "M" = "#8cb7df")) +
  theme_classic() + labs(title = tarGene$gene,
                         x="",y="TPM")
ggplot(tarGene_ratio %>% filter(stage %in% c("12","15","16","18","21")) %>% filter(tissue =="G"), aes(x = as.character(stage), y = log2(tpmR), group = gene)) +
  geom_line(size = 1)  + geom_point(size =2) +geom_hline(yintercept = 1)+
  facet_wrap( ~ gene, scales="free_y")  + scale_color_manual(values = c("F" = "#ed7d61", "M" = "#8cb7df")) +
  theme_classic() + labs(title = geneidproduct$V6[match(tarGene$gene, geneidproduct$V4)],
                         x="",y="TPM") 
subset2 <- tarGene_ratio %>% filter(stage %in% c("12","15","16","18","21")) %>% filter(tissue =="G") %>% mutate(species = "Tse") 
} #Expression in Tse no normalization
{
{
hormone_df <- data.frame(
  geneID = c(
    # Receptors
    "ESR1","ESR2","AR","PGR","NR3C1","NR3C2","GPER1",
    
    # Membrane / non-classical
    "PGRMC1","PGRMC2","PAQR5","PAQR6","PAQR7","PAQR8","PAQR9","GPRC6A",
    
    # Steroidogenesis
    "STAR","CYP11A1","CYP17A1","CYP21A2","CYP19A1",
    
    # HSD metabolism
    "HSD3B","HSD17B1","HSD17B2","HSD17B3","HSD17B6","HSD17B7",
    
    # Androgen metabolism
    "SRD5A1","SRD5A2",
    
    # Conjugation
    "SULT1E1","UGT1A1","UGT2B7","UGT2B15",
    
    # Coactivators
    "NCOA1","NCOA2","NCOA3","CREBBP","EP300","MED1",
    
    # Corepressors
    "NCOR1","NCOR2",
    
    # Signaling
    "MAPK1","MAPK3","PIK3CA","AKT1","PTEN",
    
    # Growth factors
    "IGF1","IGF1R","FGF2","FGFR1","FGFR2","FGFR3","FGFR4",
    
    # Ovary pathway
    "FOXL2","WNT4","RSPO1","BMP15","FST",
    
    # Testis pathway
    "SOX9","DMRT1","AMH",
    
    # Supporting TFs
    "GATA4","ZFPM2","NR5A1","NR0B1",
    
    # Endocrine axis
    "GNRH1","GNRHR","LHB","FSHB","LHCGR","FSHR","SHBG"
  ),
  module_function = c(
    # Receptors
    rep("Steroid receptor",7),
    
    # Membrane
    rep("Non-genomic steroid signaling",8),
    
    # Steroidogenesis
    rep("Steroid biosynthesis",5),
    
    # HSD
    rep("Steroid metabolism",6),
    
    # Androgen metabolism
    rep("Androgen metabolism",2),
    
    # Conjugation
    rep("Hormone inactivation",4),
    
    # Coactivators
    rep("Nuclear receptor coactivator",6),
    
    # Corepressors
    rep("Nuclear receptor corepressor",2),
    
    # Signaling
    rep("Downstream signaling pathway",5),
    
    # Growth
    rep("Growth factor signaling",7),
    
    # Ovary
    rep("Ovary differentiation program",5),
    
    # Testis
    rep("Testis differentiation program",3),
    
    # Supporting TF
    rep("Gonadal supporting TF",4),
    
    # Endocrine axis
    rep("Hypothalamic-pituitary-gonadal axis",7)
  ),
  stringsAsFactors = FALSE
)
} #Hormone related gene list
geneidproduct$prod <- gsub("_[0-9]*","", geneidproduct$V6)
hormone_df <- hormone_df %>% mutate(
  psiId = geneidproduct$V4[match(geneID, geneidproduct$prod)]
)
targetGenes <- hormone_df$psiId
#in Psi
id_map <- geneidproduct %>%
  mutate(prod = gsub("_[0-9]*","", V6)) %>%
  select(psiId = V4, geneSymbol = prod)
mat <- psiExp %>% select(1:8,23:30,17:22,39:44,15:16,9:14,37:38,31:36) %>%
  filter(rownames(.) %in% hormone_df$psiId)
mat <- mat %>%
  rownames_to_column("psiId") %>%
  left_join(id_map, by = "psiId") %>%
  mutate(geneSymbol = ifelse(is.na(geneSymbol), psiId, geneSymbol)) %>%
  column_to_rownames("geneSymbol") %>%
  select(-psiId)
#in Tse
id_map <- geneidproduct %>%
  mutate(prod = gsub("_[0-9]*","", V6)) %>%
  select(psiId = V4, geneSymbol = prod) %>%
  mutate(tseID = orthologs$V2[match(psiId,orthologs$V1)]) %>% na.omit()
mat <- tseExp %>% select(1:5,11:15,6:10,14:20,21:30) %>% 
  filter(rownames(.) %in% (orthologs[orthologs$V1 %in% hormone_df$psiId,] %>% pull(V2)))
mat <- mat %>%
  rownames_to_column("tseID") %>%
  left_join(id_map, by = "tseID") %>%
  mutate(geneSymbol = ifelse(is.na(geneSymbol), tseID, geneSymbol)) %>%
  column_to_rownames("geneSymbol") %>%
  select(-tseID) %>%
  select(-psiId)
##heatmap
anno_row <- hormone_df %>%
  filter(geneID %in% rownames(mat)) %>%
  left_join(id_map, by = "psiId") %>%
  mutate(geneSymbol = ifelse(is.na(geneSymbol), tseID, geneSymbol)) %>%
  distinct(geneSymbol, module_function) %>%
  column_to_rownames("geneSymbol")
mat <- mat[rownames(anno_row), ]
color_palette <- colorRampPalette(c("#8cb7df","white","#ed7d61"))(100)
pheatmap::pheatmap(
  mat, 
  treeheight_row = 0,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  scale = "row",
  color = color_palette,
  annotation_row = anno_row,
  show_rownames = TRUE,
  fontsize_row = 8,
  border_color = NA,
  main = "Hormone signaling related gene expression"
)
} #Hormone related genes expression
{
targetGenes <- c("ESR1_ENSGEVP00005018606-D1",
                 "PRGR_ENSGEVP00005011508-D1",
                 "GREB1_ENSCPBP00000019364-D1",
                 "NRIP1_ENSCPBP00000041311-D1") 
{
  psiOnlyExpr <- psiExp %>%
    select(grep("FB|MB|FG|MG|FM|MM", colnames(.)))
  
  psiOnlyExpr$MM12 <- psiOnlyExpr$MGMC12
  psiOnlyExpr$MM13 <- psiOnlyExpr$MGMC13
  psiOnlyExpr$FM12 <- psiOnlyExpr$FGMC12
  psiOnlyExpr$FM13 <- psiOnlyExpr$FGMC13
  
  psi_long <- psiOnlyExpr %>%
    tibble::rownames_to_column("gene") %>%
    gather(key="variables", value="tpm",-gene) %>%
    mutate(
      variables = gsub("GMC","G",variables),
      sex = substr(variables,1,1),
      tissue = substr(variables,2,2),
      stage = as.numeric(substr(variables,3,4)),
      species = "Ps",
      factor = paste(sex,tissue)
    ) %>%
    filter(stage %in% c(12,15,16,18,21)) %>%
    filter(gene %in% targetGenes)
} #process psi expr.
{
tseExpsub <- tseExp[rownames(tseExp) %in% orthologs$V2,] %>%
  select(grep("FG|MG|FB|MB|FM|MM",colnames(.))) %>%
  select(grep("12|15|16|18|21",colnames(.)))

tseExpsub$orthname <- orthologs$V1[match(row.names(tseExpsub),orthologs$V2)]

tseExpsub <- tseExpsub[tseExpsub$orthname %in% orthologs$V1,]
row.names(tseExpsub) <- tseExpsub$orthname
tseExpsub <- tseExpsub[,-31]

tse_long <- tseExpsub %>%
  tibble::rownames_to_column("gene") %>%
  gather(key="variables", value="tpm",-gene) %>%
  mutate(
    sex = substr(variables,1,1),
    tissue = substr(variables,2,2),
    stage = as.numeric(substr(variables,3,4)),
    species = "Tse",
    factor = paste(sex,tissue)
  ) %>%
  filter(stage %in% c(12,15,16,18,21)) %>%
  filter(gene %in% targetGenes)
} #process tse expr.
plotData <- bind_rows(psi_long, tse_long)
ggplot(plotData,
       aes(x=factor(stage),
           y=tpm,
           color=sex,
           group=factor)) +
  
  geom_line(size=1) +
  geom_point(size=1) +
  
  facet_grid(gene ~ species + tissue,
             scales="free_y") +
  
  scale_color_manual(values=c(
    "F"="#ed7d61",
    "M"="#8cb7df"
  )) +
  
  theme_classic() +
  
  labs(
    x="Stage",
    y="TPM",
    color="Sex"
  )
} #Expression in both turtles no normalization
{
psiExpsub <- psiExp[rownames(psiExp) %in% orthologs$V1, ] %>%
  select(grep("FG|MG|FB|MB", colnames(.))) %>% 
  select(grep("12|15|16|18|21", colnames(.))) %>%
  select(10, 6:9, 20, 16:19, 1:5, 11:15)
tseExpsub <- tseExp[rownames(tseExp) %in% orthologs$V2, ] %>%
  select(grep("FG|MG|FB|MB", colnames(.))) %>% 
  select(grep("12|15|16|18|21", colnames(.)))
tseExpsub$orthname <- orthologs$V1[match(row.names(tseExpsub),orthologs$V2)]
tseExpsub <- tseExpsub[tseExpsub$orthname %in% orthologs$V1,]
row.names(tseExpsub) <- tseExpsub$orthname
tseExpsub <- tseExpsub[,-21]
psiExpsub <- psiExpsub[rownames(psiExpsub) %in% orthologs$V1,]
mergedInfo <-  merge(psiExpsub, tseExpsub, by='row.names', all=TRUE)
row.names(mergedInfo) <- mergedInfo$Row.names
mergedInfo <- mergedInfo[,-1]
#UQ normalization
ei <- as.matrix(mergedInfo)
er2 <- as.data.frame(UQ_FN(ei))
psisub <- er2[1:20]
colnames(psisub) <- c(paste("psiFG",c("12","15","16","18","21"), sep = ""),
                      paste("psiMG",c("12","15","16","18","21"), sep = ""),
                      paste("psiFB",c("12","15","16","18","21"), sep = ""),
                      paste("psiMB",c("12","15","16","18","21"), sep = ""))
tsesub <- er2[21:40]
colnames(tsesub) <- c(paste("tseFB",c("12","15","16","18","21"), sep = ""),
                      paste("tseMB",c("12","15","16","18","21"), sep = ""),
                      paste("tseFG",c("12","15","16","18","21"), sep = ""),
                      paste("tseMG",c("12","15","16","18","21"), sep = ""))
#only gonad
psiGsub <- er2[1:10]
colnames(psiGsub) <- c(paste("psiFG",c("12","15","16","18","21"), sep = ""),
                      paste("psiMG",c("12","15","16","18","21"), sep = ""))
tseGsub <- er2[31:40]
colnames(tseGsub) <- c(paste("tseFG",c("12","15","16","18","21"), sep = ""),
                       paste("tseMG",c("12","15","16","18","21"), sep = ""))
#Merge
allExpr <- cbind(psiGsub, tseGsub) %>%
  filter_all(.,any_vars(.> 1))
allExpr_longer <- allExpr %>% tibble::rownames_to_column("gene") %>% 
  gather(key = "variables", value = "tpm",- gene) %>% mutate(
    species = substr(variables, 1, 3), 
    sex = substr(variables, 4, 4),   
    tissue = substr(variables, 5, 5),  
    stage = as.numeric(substr(variables, 6, 7)),
    chr = geneidproduct$V1[match(gene, geneidproduct$V4)],
    factor = paste(species, sex)) 
tarGene <- allExpr_longer %>% filter(gene == "NRIP1_ENSCPBP00000041311-D1")
tarGene <- allExpr_longer %>% filter(gene %in% targetGenes)
tarGene_ratio <- tarGene %>%
  group_by(gene, species, stage) %>%
  summarise(
    tpm_F = mean(tpm[sex == "F"], na.rm = TRUE),
    tpm_M = mean(tpm[sex == "M"], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    tpmR = tpm_F / tpm_M
  ) 
ggplot(tarGene, aes(x = as.character(stage), y = tpm, group = factor, color = sex, linetype = species)) +
  geom_line(size = 1)  + geom_point()+
  facet_wrap( ~ gene, scales = "free") + scale_color_manual(values = c("F" = "#ed7d61", "M" = "#8cb7df")) +
  theme_classic()
ggplot(tarGene_ratio, aes(x = as.character(stage), y = log2(tpmR), group = species, color = species, linetype = species)) +
  geom_hline(yintercept = 1)+
  geom_line(size = 1)  + geom_point()+
  facet_wrap( ~ gene, scales = "free") + scale_color_manual(values = c("F" = "#ed7d61", "M" = "#8cb7df")) +
  theme_classic() 
} #Expression in both turtles normalized by UQ



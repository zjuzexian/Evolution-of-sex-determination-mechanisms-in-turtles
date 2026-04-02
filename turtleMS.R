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
} # lib
{
#DNA coverage
{
depthall <- read.table(file = "psi.dna.fm.100k.noN.txt")
cptDf <- subset(depthall,depthall$V1=="chrZ") %>% select(2, 3, 4, 6) %>% mutate(
  x = (V2+V3)/2e6,
  y = V6/V4
) %>% select(5,6)
cpt <- cpt.meanvar(cptDf$y, method = "PELT") 
depthall %>% filter(V1 == "chrZ") %>% filter(row_number(.) %in% cpts(cpt))
zdepth <- depthall %>% filter(V1 == "chrZ") %>% mutate(
  strata = case_when(V3 <= 1517136 ~ "PAR",
                     V2 >= 1517136 & V3 <= 3178082 ~ "S2",
                     V2 >= 3178082 & V3 <= 10.4e6 ~ "S0.1",
                     V2 >= 10.4e6 & V3 <= 12.5 ~ "S0.1",
                     V2 >= 12.5e6 ~ "S0.2"
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
  geom_vline(xintercept = c(1.51, 3.17, 10.4, 12.5), linetype = 2) +
  guides(color = "none") +
  ylim(0, 2.5) +
  geom_hline(data = strata_means, aes(yintercept = mean_ratio, color = strata), linetype = "dashed")

}
#W-linked scaffolds
{
depthHiC <- read.table(file = "hic-nonasb-mf-cov-mpsite.ratio")
ggplot(depthHiC,aes(x=V4,y=V3,size=V2/1e6))+ coord_cartesian(xlim = c(0,1),ylim = c(0,1)) + 
  geom_point() + guides(size=guide_legend("Scaffold size / Mb")) +
  labs(y="M/F coverage",x="M/F mappable sites")+ scale_size(range = c(0.0005,8.75)) +theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_line(size=.1, color="black",linetype = "dashed"),
        panel.grid.major.y = element_line(size=.1, color="black",linetype = "dashed"),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        axis.text.x= element_text(size=15),
        axis.text.y= element_text(size=15)) 
  
wScafs <- subset(depthHiC,depthHiC$V3<0.1&depthHiC$V4<0.1&depthHiC$V2>5000)
#write.table(wScafs$V1,file="C:/Users/13719/Desktop/r_file/hic-nonasb-W-candidate.name", row.names =FALSE, col.names =FALSE, quote =FALSE,sep = "\t")

}
#SNP density
{
fsnp <- read.table("psi.snp.f.50k.txt")
msnp <- read.table("psi.snp.m.50k.txt")
fsnp <- fsnp %>% mutate(
  sex = "F"
)
msnp <- msnp %>% mutate(
  sex = "M"
)
Ratiosnp <- merge.data.frame(fsnp, msnp, by = "V4") %>% 
  select(chr = 2, st = 3, ed = 4, f = 5, m = 10) %>% 
  mutate(
    r = (f+1)/(m+1)
  )
Zsnp <- rbind(fsnp, msnp) %>% filter(V1 == "chrZ") %>% mutate(
  strata = case_when(V3 <= 1517136 ~ "PAR",
                     V2 >= 1517136 & V3 <= 3178082 ~ "S2",
                     V2 >= 3178082 & V3 <= 14051980 ~ "S1",
                     V2 >= 14051980 ~ "S0"
  )) %>% filter(
    sex == "M" | (sex == "F" & strata %in% c("PAR", "S2"))
  )
parf <- subset(Zsnp, (Zsnp$strata == "PAR" | Zsnp$strata == "S2") & Zsnp$sex=="F")
parm <- subset(Zsnp, (Zsnp$strata == "PAR" | Zsnp$strata == "S2") & Zsnp$sex=="M")
otherm <- subset(Zsnp, (Zsnp$strata != "PAR" & Zsnp$strata != "S2") & Zsnp$sex=="M")
allm <- subset(Zsnp, Zsnp$sex=="M")
ggplot(Zsnp,aes(x=(V2+V3)/2e6,y=V5,group=sex,color=sex)) +xlim(0,40) +theme_bw() + 
  scale_color_manual(values=c("#ed7d61","#8cb7df")) + 
  geom_point(alpha = .5, size = 1) + 
  geom_smooth(data = parf, method = "loess", span = .5) + 
  geom_smooth(data = allm, method = "loess", span = .1) + 
  geom_vline(aes(xintercept=1.51),linetype=2) +
  geom_vline(aes(xintercept=3.17),linetype=2) +
  geom_vline(aes(xintercept=6),linetype=2) +
  geom_vline(aes(xintercept=10.4),linetype=2) +
  geom_vline(aes(xintercept=13.1),linetype=2) +
  guides(color="none") + theme_test()
#SNP density by chromosome Boxplot
suball <- fsnp
color_panel <- c("grey90","grey40","grey90","grey40","grey90","grey40","grey90","grey40",
                 "grey90","grey40","grey90","grey40","grey90","grey40","grey90","grey40",
                 "grey90","grey40","grey90","grey40","grey90","grey40","grey90","grey40",
                 "grey90","grey40","grey90","grey40","grey90","grey40","grey90","grey40",
                 "#0000cd","#cd0000")
suball$V1 <- fct_relevel(suball$V1,"chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8",
                         "chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16",
                         "chr17","chr18","chr19","chr20","chr21","chr22","chr23","chr24",
                         "chr25","chr26","chr27","chr28","chr29","chr30","chr31","chr32",
                         "chrZ","chrW"); levels(suball$V1)
ggplot(suball,aes(x=V1,y=V5,fill=V1)) +
  geom_boxplot(outlier.colour = NULL,
               width = .4, outlier.shape = NA,
               show.legend  = FALSE) +
  scale_fill_manual(values = color_panel)+ 
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        axis.text.x= element_text(size=10),
        axis.text.y= element_text(size=15))+
  labs(x="",y="female-snp")
}
#DNA polymorphism (pai)
{
pai <- read.table(file = "pai500k-bialle.windowed.pi",header=T)
cptDf <- subset(pai,pai$CHROM=="chrZ"&pai$N_VARIANTS>200) %>% select(2, 3, 5) %>% mutate(
  x = (BIN_START+BIN_END)/2e6,
  y = PI
) %>% select(4,5)
cpt <- cpt.meanvar(cptDf$y, method = "PELT") 
pai %>% filter(CHROM == "chrZ") %>% filter(row_number(.) %in% cpts(cpt))
subpai <- subset(pai,pai$CHROM=="chrZ"&pai$N_VARIANTS>0) %>% mutate(
  strata = case_when(BIN_END <= 1.5e6 ~ "PAR",
                     BIN_START >= 1.5e6 & BIN_END <= 3.0e6 ~ "S2",
                     BIN_START >= 3.0e6 & BIN_END <= 10.5e6 ~ "S0.2",
                     BIN_START >= 10.5e6 & BIN_END <= 12.5e6 ~ "S1",
                     BIN_START >= 12.5e6 ~ "S0.1"
  )) 
ggplot(subpai,aes(x=(BIN_START)/1e6,y=PI*1e4))+ 
  geom_point(aes(color = strata))  + guides(size="none") +
  scale_color_manual(values=c("#a7d49e","#ed7d61","#fede8e","#8cb7df")) +
  geom_smooth(method = "loess", span = .4)+ theme_bw()+xlim(0,40)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        axis.text.x= element_text(size=15),
        axis.text.y= element_text(size=15)) + 
  geom_vline(aes(xintercept=1.51),linetype=2) +
  geom_vline(aes(xintercept=3.17),linetype=2) +
  geom_vline(aes(xintercept=8),linetype=2) +
  geom_vline(aes(xintercept=12.5),linetype=2) +
  geom_hline(aes(yintercept = mean(subset(subpai, strata == "S1")$PI*1e4)))+
  geom_hline(aes(yintercept = mean(subset(subpai, strata == "S0")$PI*1e4)))+
  guides(color="none") + theme_test() +
  labs(x="",y="") + ylim(0, 4)
ggplot(subpai %>% filter(strata %in% c("S0.1", "S0.2", "S1")),
       aes(x = strata, y = PI * 1e4)) +
  geom_violin(aes(fill = strata), trim = FALSE) +
  stat_summary(fun = mean, geom = "crossbar", 
               width = 0.2, color = "black", 
               fatten = 2) +
  stat_compare_means(method = "wilcox.test", 
                     comparisons = list(c("S0.2", "S0.1"),
                                        c("S1", "S0.1"),
                                        c("S1", "S0.2")), 
                     method.args = list(alternative = "greater")) +
  theme_test()

}
#ZW gametologs pair dS 
{
synZW <- read.table(file = "zwgame.identity.alignPerc") # with gametolog pairs
dndsZW <- read.table('psi.zw.game.dnds.info',header = F,row.names = 1)
geneloci <- read.table(file="psi.geneID.finalinfo.txt")
geneloci <- geneloci %>% filter(V1 == "chrZ" | V1 == "chrW")
synZW <- synZW %>% mutate(
  zID = geneloci$V5[match(synZW$V5,geneloci$V2)],
  dN = dndsZW$V2[match(zID,rownames(dndsZW))],
  dS = dndsZW$V3[match(zID,rownames(dndsZW))],
  w = dndsZW$V4[match(zID,rownames(dndsZW))]
) %>% mutate(
  strata = case_when(V6 <= 1517136 ~ "PAR",
                     V5 >= 1517136 & V6 <= 3178082 ~ "S2",
                     V5 >= 3178082 & V6 <= 10e6 ~ "S0.2",
                     V5 >= 10e6 & V6 <= 12.4e6 ~ "S1",
                     V5 >= 12.4e6 ~ "S0.1"
  )) 
ggplot(synZW, aes(x = V5/1e6, y = dS))+
  geom_point(aes(color = strata),size =3) + theme_test()+ 
  scale_color_manual(values=c("#a7d49e","#ed7d61","#fede8e","#8cb7df")) +
  geom_vline(aes(xintercept=1.51),linetype=2)+ guides(color = F)+
  geom_vline(aes(xintercept=10),linetype=2) +
  geom_vline(aes(xintercept=3.17),linetype=2) + geom_vline(aes(xintercept=12.4),linetype=2) + ylim(0,1)

ggplot(synZW,
       aes(x = strata, y = dS)) +
  geom_violin(aes(fill = strata), trim = FALSE) + geom_point() +
  stat_summary(fun = mean, geom = "crossbar", 
               width = 0.2, color = "black", 
               fatten = 2) +
  stat_compare_means(method = "wilcox.test", 
                     comparisons = list(c("S0.1", "S0.2"),
                                        c("S0.2", "S1"),
                                        c("S0.1", "S1")), 
                     method.args = list(alternative = "greater")) + 
  theme_test() + ylim(0,0.5)





}
#ZW similarity
{
zw.simi <- read.table(file = "zw.simi.lastz.final100k")
colnames(zw.simi)<-c("Chr","Pos","Mis","All","Similarity")
zw.simi <- zw.simi %>%
    filter(All>5000) %>% mutate(strata = 
                                  case_when(
                                    Pos <= 1.51e6 ~ "PAR",
                                    Pos > 1.51e6 & Pos < 3.17e6 ~ "S2",
                                    Pos > 3.17e6 & Pos < 10.6e6 ~ "S0.1",
                                    Pos > 10.6e6 & Pos < 13e6 ~ "S1",
                                    Pos > 13e6 ~ "S0.2"
                                  )) %>% mutate(
                                    color_group = ifelse(Similarity > 89.8, "high",
                                                         ifelse(Similarity >= 87 & Similarity <= 89.8, "mid", "low"))
                                  )
ggplot(zw.simi,aes(x=Pos/1e6,y=Similarity)) + 
    geom_point(aes(size = All,color= strata)) + 
    geom_smooth(method = "loess", span = .15, se =F, color = "#dcdcdc") +
    labs(x="Position(Mbp)",y="Z-W Similarity")+ theme_bw() + theme(legend.position = "none") +
  scale_color_manual(values = c("PAR"="#add69d","S2"="#dc6a67","S1"="#f3d68a",
                                "S0.1"="#8ac4eb","S0.2"="#8ac4eb")) + theme_test()+ ylim(78,94)+
  geom_vline(xintercept = c(1.51,3.17,10.4,13.1),linetype =2)+
  scale_size(range = c(.5, 4))
ggplot(zw.simi,aes(x=strata,y=Similarity)) + 
  geom_boxplot(aes(color= strata)) +
  stat_summary(fun = mean, geom = "crossbar", 
               width = 0.2, color = "black", 
               fatten = 2) +
  stat_compare_means(method = "wilcox.test", 
                     comparisons = list(c("S0.1", "S0.2"),
                                        c("S0.2", "S1"),
                                        c("S0.1", "S1")), 
                     method.args = list(alternative = "less"))

cptDf <- subset(zw.simi %>% filter(Similarity >= 78 & Similarity <= 94)) %>% select(x = 2, y = 5) 
cpt <- cpt.mean(cptDf$y, method = "PELT") 
cpts(cpt)
zw.simi %>% filter(Similarity >= 78 & Similarity <= 94) %>% filter(row_number(.) %in% cpts(cpt))
}
#RC rate
{
rc <- read.table(file = "chrZ-rcrate-200k",header=T)
rc1 <- subset(rc,rc$CIL>=0) %>%
  mutate(strata = case_when(
    End/1e6 <= 1.51 ~ "par",
    End/1e6 > 1.51 & End/1e6 <= 3.17 ~ "s2",
    End/1e6 > 3.17 & End/1e6 <= 8 ~ "s0.1",
    End/1e6 > 8 & End/1e6 <= 12.5 ~ "s1",
    End/1e6 > 12.5 ~ "s0.2"
  ))
rc1 <- rc1 %>%
  mutate(pos = (Start + End) / 2 / 1e6)  
ggplot(rc1, aes(x = Start/1e6)) +
  geom_ribbon(aes(ymin = CIL, ymax = CIR), fill = "#d0d9e8", alpha = 0.6) +
  geom_step(aes(y = Rho), color = "#2a6fba", size = 1) +
  labs(
    x = "Position on Chromosome (Mb)",
    y = "Recombination rate (ρ)"
  ) +
  theme_test() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  )+ geom_vline(aes(xintercept=1.51),linetype=2)+ geom_vline(aes(xintercept=10.4),linetype=2)+
  geom_vline(aes(xintercept=3.17),linetype=2) + geom_vline(aes(xintercept=12.5),linetype=2) 
ggplot(data = rc1, aes(x = strata, y = Rho, color = strata),shape =21) +
  geom_boxplot()  +
  theme_bw() +
  stat_compare_means(
    comparisons = list(c("s0.1", "s0.2"),c("s2","s1")),
    method = "wilcox.test",
    method.args = list(alternative = "greater")
  ) + 
  scale_color_manual(
    values = c("par" = "#a7d49e",
               "s2" = "#ed7d61",
               "s1" = "#fede8e",
               "s0.1" = "#8cb7df",
               "s0.2" = "#8cb7df"))
cptDf <- rc1 %>% select(x=7, y=3)
cpt <- cpt.mean(cptDf$y, method = "PELT") 
rc1 %>% filter(row_number(.) %in% cpts(cpt))
}
#Chromosome synteny   
{
psi2GgaSC <- read.table(file = "psisex-chick.delta.filt.coords", skip = 4)  
psi2GgaZ <- psi2GgaSC %>% filter(V13 == "chrZ" &V5 >= 500)
p1 <- ggplot(data = psi2GgaZ) +
  geom_segment(aes(x=V3/1000000,xend=V1/1000000,y=1,yend=2)) + theme_test() + xlim(0,40) + 
  geom_vline(aes(xintercept=1.51),linetype=2)+
  geom_vline(aes(xintercept=3.17),linetype=2)+
  geom_vline(aes(xintercept=6),linetype=2)+
  geom_vline(aes(xintercept=10.4),linetype=2)+
  geom_vline(aes(xintercept=13.05),linetype=2) +
  guides(color="none") 
psizw <- read.table("psizw.nucmer.coords", skip = 4)
psizw <- psizw %>% filter(
  V7 > 60 & V7 < 90 & V5 >= 1000
) %>% mutate(
  strata = case_when(V2 <= 1517136 ~ "PAR",
                     V1 >= 1517136 & V2 <= 3178082 ~ "S2",
                     V1 >= 3178082 & V2 <= 8e6 ~ "S1.1",
                     V1 >= 8e6 & V2 <= 13e6 ~ "S1.2",
                     V1 >= 13e6 ~ "S0"
  ))
p2 <- ggplot(data = psizw) +
  geom_segment(aes(x=V4/1000000,xend=V2/1000000,y=1,yend=2,color=strata)) + theme_test() + xlim(0,40)+ 
  geom_vline(aes(xintercept=1.51),linetype=2)+
  geom_vline(aes(xintercept=3.17),linetype=2)+
  geom_vline(aes(xintercept=6),linetype=2)+
  geom_vline(aes(xintercept=10.4),linetype=2)+
  geom_vline(aes(xintercept=13.05),linetype=2) +
  guides(color="none") 
ggplot(data = psizw %>% filter(V1>3e6&V2<10e6)) +
  geom_segment(aes(x=V4/1000000,xend=V2/1000000,y=1,yend=2,color=strata)) + theme_test() + xlim(0,15)+ 
  geom_vline(aes(xintercept=1.51),linetype=2)+
  geom_vline(aes(xintercept=3.17),linetype=2)+
  geom_vline(aes(xintercept=6),linetype=2)+
  geom_vline(aes(xintercept=10.4),linetype=2)+
  geom_vline(aes(xintercept=13.05),linetype=2) +
  guides(color="none") 

ggarrange(p1,p2, ncol =1) 
}
#Gene density
{
geneDens <- read.table(file = "psi.geneDens.100k.txt")
geneDensZ <- geneDens %>% filter(V1 == "chrZ") %>% mutate(
  strata = case_when(V3 <= 1517136 ~ "PAR",
                     V2 >= 1517136 & V3 <= 3178082 ~ "S2",
                     V2 >= 3178082 & V3 <= 10.4e6 ~ "S0.1",
                     V2 >= 10.4e6 & V3 <= 13e6 ~ "S1",
                     V2 >= 13e6 ~ "S0.2"
  ))
geneDensW <- geneDens %>% filter(V1 == "chrW") %>% mutate(
  strata = case_when(V3 <= 1517136 ~ "PAR",
                     V2 >= 1517136 & V3 <= 3178082 ~ "S2",
                     V2 >= 3178082 & V3 <= 10.4e6 ~ "S0.1",
                     V2 >= 10.4e6 & V3 <= 13e6 ~ "S1",
                     V2 >= 13e6 ~ "S0.2"
  ))
ggplot(geneDensW, aes(x = (V2+V3)/2e6, y = V5)) +
  geom_bar(stat = "identity", aes(fill = strata)) + theme_test()
}
#Dosage compensation
{
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
      chr == "chrZ" & st > 1517136 & ed <= 3178082 ~ "S2",
      chr == "chrZ" & st >= 3178082 & ed <= 10.4e6 ~ "S0.1",
      chr == "chrZ" & st >= 10.4e6 & ed <= 12.5e6 ~ "S1",  # 12.5 corrected to 12.5e6
      chr == "chrZ" & st >= 12.5e6 ~ "S0.2",
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






}
#W gametolog degeneration
{
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
 
}
#Faster Z
{
dndsAll <- read.table("psi2tse.dnds.pair.txt")
geneloci <- read.table("trans.loci")
zsub <- dndsAll %>% filter(V1 %in% subset(geneloci, V1 == "chrZ")$V4) %>% 
  filter(V4 >= 0.01 & V4 <= 2) %>% mutate(
    chr = case_when(
      V7 <= 1.5e6 ~ "Z-PAR",
      V6 > 1.5e6 & V7 < 3.17e6 ~ "Z-S3",
      V6 > 3.17e6 & V7 < 10.4e6 ~ "Z-S1",
      V6 > 10.4e6 & V7 < 12.5e6 ~ "Z-S2",
      V6 > 12.5e6 ~ "Z-S0"
    ))
wsub <- dndsAll %>% filter(V1 %in% subset(geneloci, V1 == "chrW")$V4) %>% 
  filter(V4 >= 0.01 & V4 <= 2) %>%
  mutate(chr = "chrW") 
asub <- dndsAll %>% 
  filter(V1 %in% subset(geneloci, V1 != "chrZ" & V1 != "chrW" & !grepl("Hic", V1))$V4) %>%
  filter(V4 >= 0.01 & V4 <= 2) %>% 
  mutate(chr = gsub("chr","", V5),
         chr = ifelse(chr < 9, "micro", "macro")) 
dfsub <- rbind(asub,zsub,wsub) %>% mutate(
  V2 = ifelse(V2 < 0 ,0, V2)
)
ggplot(dfsub,aes(x = chr, y=V2,fill=chr))+
  geom_violin(trim = FALSE) + theme_bw() + 
  stat_summary(fun = mean, geom = "crossbar", 
               width = 0.2, color = "white", 
               fatten = 2)  +
  stat_compare_means(method = "wilcox.test", 
                     comparisons = list(c("Z-S0", "macro"),
                                        c("Z-S1", "macro"),
                                        c("chrW", "macro")), 
                     method.args = list(alternative = "greater"))+
  theme_test() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        axis.text.x= element_text(size=10),
        axis.text.y= element_text(size=15))+
  labs(x="",y="dN/dS") + ylim(-0.5,1)



}
#Masculinization of Z-linked genes   
{
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
testisZsub <- psiDEgeneSub %>% filter(chr == "chrZ") %>% select(6:9) %>% 
  mutate(
    strata = case_when(ed <= 1517136 ~ "PAR",
                       st >= 1517136 & ed <= 3178082 ~ "S2",
                       st >= 3178082 & ed <= 14051980 ~ "S1",
                       st >= 14051980 ~ "S0"
    )) 
testisStrataProp <- testisZsub %>% group_by(type, strata) %>%
  summarise(freq = n(), .groups = "drop") %>%
  group_by(strata) %>%
  mutate(perc = freq / sum(freq) * 100)
ggplot(data = testisStrataProp, aes(x = strata, y = perc, fill = type)) +
  geom_bar(stat = "identity") + labs(x = "", y = "% Genes") +
  theme_test() + scale_fill_manual(values=c("#a7d49e","#ed7d61","#8cb7df"))
{
  testisZsub <- psiDEgeneSub %>% filter(chr == "chrZ") %>% select(6:9)
  #100k win
  window_size <- 1000000
  
  testisZsub <- testisZsub %>%
    mutate(window = floor(st / window_size) * window_size)
  
  window_summary <- testisZsub %>%
    group_by(window, type) %>%
    summarise(count = n(), .groups = "drop")
  
  window_plot <- window_summary %>%
    group_by(window) %>%
    mutate(total = sum(count),
           prop = count / total)
  
  ggplot(window_summary, aes(x = as.numeric(window / 1e6), y = count, fill = type)) +
    geom_bar(stat = "identity") +
    labs(x = "Window start (Mb)", y = "Proportion", fill = "Type") +
    theme_test() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12)
    ) + scale_fill_manual(values = c("white", "#8cb7df"))
  } # position on Z
}
#rDNA amplification on W
{
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
#write.table(tpmmean, file = "psi.45SrDNA.sRNA.tpm.uniq.txt")
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
}


  





synZW <-  synZW %>% mutate(
  strata = case_when(V6 <= 1517136 ~ "PAR",
                     V5 >= 1517136 & V6 <= 3178082 ~ "S2",
                     V5 >= 3178082 & V6 <= 14051980 ~ "S1",
                     V5 >= 14051980 ~ "S0"
  )) %>% filter(V7 > 70)
ggplot(synZW,aes(x=V2/1000000,xend=V5/1000000,y=1,yend=2))+
  geom_segment(stat = "identity", aes(color = strata),size =1) + theme_test()+
  scale_color_manual(values=c("#a7d49e","#ed7d61","#fede8e","#8cb7df")) +
  geom_vline(aes(xintercept=1.51),linetype=2)+ guides(color = F)+
  geom_vline(aes(xintercept=3.17),linetype=2) + geom_vline(aes(xintercept=14.05),linetype=2)+
  geom_vline(aes(xintercept=6),linetype=2)+
  geom_vline(aes(xintercept=10.4),linetype=2)



nucmerZW <- read.table(file = "z2w.delta.filt.coords", skip = 4)

nucmerZW <- nucmerZW %>% filter(V7 >=60 & V7 <=96) 
ggplot(nucmerZW %>% filter(V5 > 2000),aes(x=V1/1000000,xend=V3/1000000,y=1,yend=2))+
  geom_segment(stat = "identity",size =1) + theme_test()+
  scale_color_manual(values=c("#a7d49e","#ed7d61","#fede8e","#8cb7df")) +
  geom_vline(aes(xintercept=1.51),linetype=2)+ guides(color = F)+
  geom_vline(aes(xintercept=3.17),linetype=2) + geom_vline(aes(xintercept=14.05),linetype=2)



synZW <- read.table(file = "psi.zwlast.txt") # with nucleotide alignments 
synZW <-  synZW %>% filter(V4 == "chrW" & V7 > 1500) %>% mutate(
  strata = case_when(V3 <= 1517136 ~ "PAR",
                     V2 >= 1517136 & V3 <= 3178082 ~ "S2",
                     V2 >= 3178082 & V3 <= 14051980 ~ "S1",
                     V2 >= 14051980 ~ "S0" ))
ggplot(synZW,aes(x=V5/1000000,xend=V2/1000000,y=1,yend=2))+
  geom_segment(stat = "identity", aes(color = strata),size =1) + theme_minimal()+
  geom_vline(aes(xintercept=1.51),linetype=2)+
  geom_vline(aes(xintercept=3.17),linetype=2) + geom_vline(aes(xintercept=14.05),linetype=2)

ggplot(synZW,aes(x=V2/1E6, y=V8))+
  geom_point(aes(color = strata),size =1) + theme_minimal()+
  geom_vline(aes(xintercept=1.51),linetype=2)+
  geom_vline(aes(xintercept=3.17),linetype=2) + geom_vline(aes(xintercept=14.05),linetype=2)

#ZW similarity test ver
simiZW <- read.table(file = "psi.zwlast.txt")
simiZWsub <-  simiZW %>% mutate(
  strata = case_when(V6 <= 1517136 ~ "PAR",
                     V5 >= 1517136 & V6 <= 3178082 ~ "S2",
                     V5 >= 3178082 & V6 <= 14051980 ~ "S1",
                     V5 >= 14051980 ~ "S0"
  )) %>% filter(V7 > 2000)
ggplot(simiZWsub,aes(x=V2/1e6,xend=V3/1e6,y=1,yend=1))+
  geom_segment(stat = "identity", aes(color = V8),size = 10)+
  scale_color_gradient2(low = "#8cb7df", mid = "#fede8e", high = "#ed7d61", midpoint = 90) + theme_test()

ggplot(simiZWsub, aes(x = V2/1e6, y = V8)) +
  geom_point(aes(size = V7, color = V8)) + theme_test() + 
  geom_vline(aes(xintercept=1.51),linetype=2)+
  geom_smooth(span = .3) +
  geom_vline(aes(xintercept=3.17),linetype=2) + geom_vline(aes(xintercept=13.05),linetype=2) 




##this now
simiZW <- read.table(file = "psi.test.last.txt")
simiZWsub <-  simiZW %>% mutate(
  strata = case_when(V3 <= 1517136 ~ "PAR",
                     V2 >= 1517136 & V3 <= 3178082 ~ "S2",
                     V2 >= 3178082 & V3 <= 13e6 ~ "S1",
                     V2 >= 13e6 ~ "S0"
  ))  %>% filter(V5 > 10000) %>% mutate(
    color_group = ifelse(V6 > 92, "high",
                                        ifelse(V6 >= 89 & V6 <= 92, "mid", "low"))
  )
ggplot(simiZWsub, aes(x = V2/1e6, y = V6)) + scale_size(range = c(.4,8))+
  geom_point(aes(size = V5, color = color_group)) + theme_test() + 
  geom_vline(aes(xintercept=1.51),linetype=2)+
  geom_smooth(span = .1) +
  geom_vline(aes(xintercept=3.17),linetype=2) + geom_vline(aes(xintercept=13),linetype=2) 

ggplot(simiZWsub,aes(x=V2/1e6,xend=V3/1e6,y=1,yend=1))+
  geom_segment(stat = "identity", aes(color = color_group),size = 10) +
  scale_color_manual(values = c(
    "low" = "#8cb7df",
    "mid" = "#fede8e",
    "high" = "#ed7d61"
  )) + theme_test()

simiZWsub <- simiZW %>%
  filter(V5 > 2000) %>%
  mutate(
    group = case_when(
      V6 < 88 ~ "low",
      V6 >= 88 & V6 <= 90 ~ "mid",
      V6 > 90 ~ "high"
    )
  ) %>%
  group_by(group) %>%
  mutate(
    bin = case_when(
      group == "low"  ~ ntile(V6, 3),
      group == "mid"  ~ ntile(V6, 3) + 3,
      group == "high" ~ ntile(V6, 4) + 6
    ),
    color_group = paste0("g", bin)
  ) %>%
  ungroup()

color10 <- colorRampPalette(c("#8cb7df", "#fede8e", "#ed7d61"))(10)
ggplot(simiZWsub, aes(x = V2/1e6, xend = V3/1e6, y = 1, yend = 1)) +
  geom_segment(aes(color = color_group), size = 10) +
  scale_color_manual(values = setNames(color10, paste0("g", 1:10))) +
  theme_test() +
  labs(color = "Similarity Level")


cptDf <- simiZWsub %>% select(2, 3, 6) %>% mutate(
  x = (V2+V3)/2e6,
  y = V6
) %>% select(4, 5) %>% 
  arrange(x)
cpt <- cpt.mean(cptDf$y, method = "PELT") 

cpts(cpt)

cpWin <- simiZWsub %>% filter(row_number(.) %in% cpts(cpt)) %>% 
  arrange(V2)

cpWin


#ZW similarity
all <- read.table(file = "zw.simi.lastz.final.200k")
suball <- subset(all,all$V4>5000)
par <- subset(suball,suball$V2<=1517136) %>% mutate(strata = "par")
s2 <- subset(suball,suball$V2>1517136&suball$V2<=3178082) %>% mutate(strata = "s2")
s1 <- subset(suball,suball$V2>3178082&suball$V2<=14051980) %>% mutate(strata = "s1")
s0 <- subset(suball,suball$V2>14051980) %>% mutate(strata = "s0")
allsub <- rbind(s0,s1,s2,par)
ggplot(allsub,aes(x=strata,y=V5,color=strata)) + 
  geom_boxplot() + 
  stat_compare_means(comparisons = list(c("s0", "s1")), label="p.signif",
                     method="wilcox.test", method.args = list(alternative = "less"))  + theme_bw()

suball$color_group <- with(suball, ifelse(V5 > 88, "high",
                                          ifelse(V5 >= 85.5 & V5 <= 88, "mid", "low")))
ggplot() + 
  geom_point(data = suball, aes(x = V2 / 1e6, y = V5, size = V4, color = color_group)) + 
  xlim(0, 40) + ylim(80, 95) +
  stat_smooth(data = suball, aes(x = V2 / 1e6, y = V5), method = "loess", span = 0.3) +
  scale_color_manual(values = c(
    "low" = "#8cb7df",
    "mid" = "#fede8e",
    "high" = "#ed7d61"
  )) +
  labs(x = "", y = "") + 
  guides(size = "none") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 15)
  ) +
  geom_vline(xintercept = c(1.51, 3.17, 14.05), linetype = 2)







#GC level
colfunc <- colorRampPalette(c("#ffffff", "#8cb7df"))
gc <- read.table(file = "gc.100k")
gcsub <- subset(gc,gc$V1=="chrZ")
ggplot(gcsub, aes(x=V2/1000000, xend=V3/1000000, y=V4, yend=V4)) +
  geom_step()+
  theme_bw() + xlab("Psi chrZ") +ylab("GC Prop.") 

mid = median(gcsub$V4)
ggplot()+geom_segment(data=gcsub,aes(x=V2/1000000,y=1,xend=V3/1000000,yend=1,color=V4),size = 10)+
  labs(x="",y="") + theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x= element_blank(),
        axis.text.y= element_blank())+ geom_vline(aes(xintercept=3.5),linetype=2)+
  geom_vline(aes(xintercept=12),linetype=2) + geom_vline(aes(xintercept=26),linetype=2)+ xlim(0,40)+
  scale_colour_gradientn(colours = colfunc(50))
#repeat level
colfunc <- colorRampPalette(c("#ffffff", "#ed7d61"))
rep <- read.table(file = "repeat.proc.100k")
repall <- read.table(file = "repeat.info.all.denovo", sep='\t', header=T)
repsub <- subset(rep,rep$V1=="chrZ")
ggplot(repsub, aes(x=V2/1000000, xend=V3/1000000, y=V4, yend=V4)) +
  geom_step()+ 
  theme_bw() + xlab("Psi chrZ") +ylab("Repeat Prop.")+
  geom_vline(aes(xintercept=1.51),linetype=2)+
  geom_vline(aes(xintercept=3.17),linetype=2) + geom_vline(aes(xintercept=14.05),linetype=2)
ggplot()+geom_point(data=repsub,aes(x=V2/1000000,y=V4,color=V4))+
  labs(x="",y="") + xlim(0,40)+ theme_bw()+ 
  scale_color_gradient2(low = "white",mid="#ed7d61",midpoint = 0.8,high="red3")
ggplot()+geom_segment(data=repsub,aes(x=V2/1000000,y=1,xend=V3/1000000,yend=1,color=V4),size = 10)+
  labs(x="",y="") + theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x= element_blank(),
        axis.text.y= element_blank())+ geom_vline(aes(xintercept=1.5),linetype=2)+
  geom_vline(aes(xintercept=3.7),linetype=2) + geom_vline(aes(xintercept=14.5),linetype=2)+ xlim(0,40)+
  scale_color_gradient2(low = "white",mid="#ed7d61",midpoint = 0.7,high="#cd0000")
#Dosage compensation with Z-excluded par gene
exp <- read.table(file = "psi.tpmmean.NewFilt.txt",header = T,sep = ',',row.names = 1)
zexpar <- read.table(file = "psi.zexpar.genelist")
expz <- exp[rownames(exp) %in% zexpar$V1,] %>% dplyr::select(1:8,23:30) #brain
expz <- exp[rownames(exp) %in% zexpar$V1,] %>% dplyr::select(15:22,37:44) #kidney
expz <- exp[rownames(exp) %in% zexpar$V1,] %>% dplyr::select(15:16,9:14,37:38,31:36) #gonad
expz <- expz %>% filter_all(.,any_vars(.>=1))
test <- tidyr::gather(expz)
test$stage <- gsub("FB|MB","",test$key) #brain
test$sex <- gsub("B.*","",test$key)
test$stage <- gsub("FGMC|MGMC|FM|MM","",test$key) #kidney
test$sex <- gsub("GMC[0-9].*|M[0-9].*","",test$key)
test$stage <- gsub("FGMC|MGMC|FG|MG","",test$key) #gonad
test$sex <- gsub("GMC[0-9].*|G[0-9].*","",test$key)
expauto <- exp[!rownames(exp) %in% zexpar$V1,] %>% dplyr::select(1:8,23:30) #brain
expauto <- exp[!rownames(exp) %in% zexpar$V1,] %>% dplyr::select(15:22,37:44) #kidney
expauto <- exp[!rownames(exp) %in% zexpar$V1,] %>% dplyr::select(15:16,9:14,37:38,31:36) #gonad
expauto <- expauto %>% filter_all(.,any_vars(.>=1))
genchrinfo <- read.table(file = "gene.chrinfo") 
chr15 <- subset(genchrinfo,genchrinfo$V1=="chr15") #pick one micro Chr with similar gene density
expautosub <- expauto[rownames(expauto) %in% chr15$V2,]
test2 <- tidyr::gather(expautosub)
test2$stage <- gsub("FB|MB","",test2$key) #brain
test2$sex <- gsub("B.*","A",test2$key)
test2$stage <- gsub("FGMC|MGMC|FM|MM","",test2$key) #kidney
test2$sex <- gsub("GMC.*|M[0-9].*","A",test2$key)
test2$stage <- gsub("FGMC|MGMC|FG|MG","",test2$key) #gonad
test2$sex <- gsub("GMC.*|G[0-9].*","A",test2$key)
testallpsi <- rbind(test,test2)
testallpsi$Species <- "Psi"
ggplot(data = testallpsi,aes(x=stage,y=value,fill=factor(sex,levels = c("F","FA","M","MA")),
                             color=factor(sex,levels = c("F","FA","M","MA"))))+
  geom_boxplot(outlier.colour = "white", position=position_dodge(0.8),
               width = .6, outlier.shape = NA,
               show.legend  = T) + coord_cartesian(ylim = c(-2,160))+
  scale_fill_manual(values = c("#cd0000","white","#0000cd","white"),name="")+ 
  scale_color_manual(values = c("black","#cd0000","black","#0000cd"),name="")+
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        axis.text.x= element_text(size=10),
        axis.text.y= element_text(size=15))+
  labs(x="",y="Z-linked genes")
#W gene degeneration
wgenehomo <- read.table(file = "wgene.homoType") #zw or zwa or za
wgeneORF <- read.table(file = "psi.w.disrupted.txt") #disrupted orf
expGene <- read.table(file = "psi.tpmAtLeast1.txt") #exp or silence
wgeneBestGame <- read.table(file = "psi.zwORzwa.bestZ.txt") #zw or zwa gametologs(best one)
geneloci <- read.table(file = "trans.loci")
wgeneBestGame <- wgeneBestGame %>% mutate(
  strata = case_when(
    V4 <= 1.5e6 ~ "PAR",
    V3 > 1.5e6 & V4 < 3.7e6 ~ "S2",
    V3 > 3.7e6 & V4 < 14.5e6 ~ "S1",
    V3 > 14.5e6 ~ "S0"
  )
)
wgeneStatus <- wgenehomo %>% mutate(
  chr = geneloci$V1[match(V1, geneloci$V4)],
  st = geneloci$V2[match(V1, geneloci$V4)],
  ed = geneloci$V3[match(V1, geneloci$V4)]) %>% filter(chr == "chrW") %>% mutate(
    orf = ifelse(V1 %in% wgeneORF$V1, "Disrupted", "Intact"),
    exp = ifelse(V1 %in% expGene$x, "Expressed", "Silenced"),
    game = wgeneBestGame$strata[match(V1, wgeneBestGame$V1)]
  ) 
#What are those genes with autosome paralogs
geneProduct <- read.table("psi.geneID.finalinfo.txt")
targetAuto <- geneProduct %>% filter(V5 %in% subset(wgeneStatus, V2 %in% c("WA", "ZWA"))$V1) %>%
  mutate(prod = gsub("_[0-9]*", "" , V6))
library(clusterProfiler)
library(org.Hs.eg.db)
library(cols4all)
targetGenes <- bitr(unique(targetAuto$prod),
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
#Go
GOresult <- ego@result %>% filter(p.adjust <= 0.05 | qvalue <= 0.05) %>% head(10) %>%
  mutate(
    GeneRatio_numeric = as.numeric(sapply(strsplit(as.character(GeneRatio), "/"), function(x) as.numeric(x[1])/as.numeric(x[2]))),
  ) 
ggplot(data = GOresult, aes(x = GeneRatio_numeric, y = factor(Description, levels = rev(GOresult$Description)), fill = -log10(p.adjust))) +
  geom_bar(stat = 'identity', width = 0.8, alpha = 0.9) +  
  scale_fill_gradient(low = "#fff6f1",high = "#ed7d61")+
  labs(x = 'Gene Ratio', y = '') +
  geom_text(aes(x = 0.005, 
                label= Description),
            hjust= 0) +
  theme_classic() + theme(axis.text.y = element_blank()) 
#Kegg
ekegg <- enrichKEGG(gene = targetGenes$ENTREZID, 
                    organism = 'hsa',
                    pvalueCutoff = 0.05,
                    pAdjustMethod = "fdr")
keggresult <- ekegg@result
keggresult$FoldEnrichment <- apply(keggresult,1,function(x){
  GeneRatio=eval(parse(text=x["GeneRatio"]))
  BgRatio=eval(parse(text=x["BgRatio"]))
  enrichment_fold=round(GeneRatio/BgRatio,2)
  enrichment_fold
})
keggresult <- mutate(keggresult, RichFactor= Count / as.numeric(sub("/\\d+", "", BgRatio)))
geneID_list <- strsplit(keggresult$geneID, "/")
gene_symbols <- lapply(geneID_list, function(ids) {
  mapped_genes <- bitr(ids, fromType = "ENTREZID", toType = "SYMBOL", OrgDb = org.Hs.eg.db)
  return(mapped_genes$SYMBOL) 
})
keggresult$geneSymbols <- sapply(gene_symbols, function(symbols) paste(symbols, collapse = "/"))
keggresult$Description <- factor(keggresult$Description,levels= rev(keggresult$Description))
kegg1 <- keggresult[order(keggresult$Count, keggresult$pvalue, decreasing = T),] %>% 
  filter(p.adjust <= 0.05 | qvalue <= 0.05)
ggplot(data = kegg1[1:10,], aes(x = Count, y = Description, fill = -log10(p.adjust))) +
  geom_bar(stat = 'identity', width = 0.8, alpha = 0.9) +  
  scale_fill_gradient(low = "#f3f9fd",high = "#0000cd")+
  labs(x = 'Number of Genes', y = '') +
  geom_text(aes(x = .5, 
                label= Description),
            hjust= 0) +
  theme_classic() + theme(axis.text.y = element_blank()) 
#Define by strata
wgeneStatusOnlyZWorZWA <- wgeneStatus %>% filter(grepl("ZW",V2)) %>% 
  select(game, orf, exp) %>% group_by(game, orf, exp) %>% 
  summarise(freq = n(), .groups = "drop") 
tardf <- wgeneStatusOnlyZWorZWA
game_data <- tardf %>%
  group_by(game) %>%
  summarise(width = sum(freq), .groups = "drop") %>%
  arrange(game) %>%
  mutate(ymid = 1,
         ymax = ymid + 0.5, ymin = ymid - 0.5,
         xmin = c(0, head(cumsum(width), -1)),
         xmax = cumsum(width),
         xmid = (xmax + xmin) / 2)
#cal orf width based on game
orf_data <- tardf %>%
  group_by(game, orf) %>%
  summarise(width = sum(freq), .groups = "drop") %>%
  left_join(game_data %>% select(game, xmin, xmax), by = "game") %>%
  arrange(game, orf) %>%
  group_by(game) %>%
  mutate(ymid = 2,
         ymax = ymid + 0.5, ymin = ymid - 0.5,
         xmin = xmin + c(0, head(cumsum(width), -1)),
         xmax = xmin + width,
         xmid = (xmax + xmin) / 2)
# cal exp width based on game + orf
exp_data <- tardf %>%
  group_by(game, orf, exp) %>%
  summarise(width = sum(freq), .groups = "drop") %>%
  left_join(orf_data %>% select(game, orf, xmin, xmax), by = c("game", "orf")) %>%
  arrange(game, orf, exp) %>%
  group_by(game, orf) %>%
  mutate(ymid = 3,
         ymax = ymid + 0.5, ymin = ymid - 0.5,
         xmin = xmin + c(0, head(cumsum(width), -1)),
         xmax = xmin + width,
         xmid = (xmax + xmin) / 2)
df_long <- bind_rows(game_data, orf_data, exp_data) #merge layers
ggplot(df_long, aes(xmid, ymid, fill = as.factor(game))) +
  geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, 
                alpha = as.factor(ymid), color = as.factor(game))) +
  geomtextpath::geom_textpath(aes(y = ymid + 0.25, label = paste(game, orf, exp, sep = " - "), 
                                  group = paste(game, orf, exp, sep = " - ")), size = 3)   +
  scale_y_continuous(limits = c(-0.5, 3.6)) +
  coord_polar() +
  theme_void() +
  theme(legend.position = "none")
#Define by Homologs
tardf <- wgeneStatus %>%
  select(V2, orf, exp) %>% rename(c(homo = V2)) %>% group_by(homo, orf, exp) %>% 
  summarise(freq = n(), .groups = "drop") 
homo_data <- tardf %>%
  group_by(homo) %>%
  summarise(width = sum(freq), .groups = "drop") %>%
  arrange(homo) %>%
  mutate(ymid = 1,
         ymax = ymid + 0.5, ymin = ymid - 0.5,
         xmin = c(0, head(cumsum(width), -1)),
         xmax = cumsum(width),
         xmid = (xmax + xmin) / 2)
#cal orf width based on homo
orf_data <- tardf %>%
  group_by(homo, orf) %>%
  summarise(width = sum(freq), .groups = "drop") %>%
  left_join(homo_data %>% select(homo, xmin, xmax), by = "homo") %>%
  arrange(homo, orf) %>%
  group_by(homo) %>%
  mutate(ymid = 2,
         ymax = ymid + 0.5, ymin = ymid - 0.5,
         xmin = xmin + c(0, head(cumsum(width), -1)),
         xmax = xmin + width,
         xmid = (xmax + xmin) / 2)
# cal exp width based on game + orf
exp_data <- tardf %>%
  group_by(homo, orf, exp) %>%
  summarise(width = sum(freq), .groups = "drop") %>%
  left_join(orf_data %>% select(homo, orf, xmin, xmax), by = c("homo", "orf")) %>%
  arrange(homo, orf, exp) %>%
  group_by(homo, orf) %>%
  mutate(ymid = 3,
         ymax = ymid + 0.5, ymin = ymid - 0.5,
         xmin = xmin + c(0, head(cumsum(width), -1)),
         xmax = xmin + width,
         xmid = (xmax + xmin) / 2)
df_long <- bind_rows(homo_data, orf_data, exp_data) #merge layers
ggplot(df_long, aes(xmid, ymid, fill = as.factor(homo))) +
  geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, 
                fill = as.factor(homo)), color = "white") +
  geomtextpath::geom_textpath(aes(y = ymid + 0.25, label = paste(homo, orf, exp, sep = " - "), 
                                  group = paste(homo, orf, exp, sep = " - ")), size = 3)   +
  scale_y_continuous(limits = c(-0.5, 3.6)) +
  coord_polar() +
  theme_void() +
  theme(legend.position = "none")
#Autosome distribution of WA and ZWA genes
wgeneBestPara <- read.table(file = "psi.waORzwa.bestA.txt") #wa or zwa autosomal paralogs(best one)
wgeneBestPara <- wgeneBestPara %>%  mutate(
  chr = geneloci$V1[match(V1, geneloci$V4)],
  st = geneloci$V2[match(V1, geneloci$V4)],
  ed = geneloci$V3[match(V1, geneloci$V4)]) %>% filter(chr == "chrW") %>% mutate(
    orf = ifelse(V1 %in% wgeneORF$V1, "Disrupted", "Intact"),
    exp = ifelse(V1 %in% expGene$x, "Expressed", "Silenced")
  ) 
dfFreq <- wgeneBestPara %>% mutate(factor = paste(orf, exp, sep = "_")) %>% 
  group_by(V2, factor) %>% 
  summarise(freq = n(), .groups = "drop") 
chrorder <- c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8",
              "chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16",
              "chr17","chr18","chr19","chr20","chr21","chr22","chr23","chr24",
              "chr25","chr26","chr27","chr28","chr29","chr30","chr31","chr32",
              "chrZ","chrW")

ggplot(dfFreq, aes(x = factor(V2, levels = chrorder), y = freq, fill = factor)) +
  geom_bar(stat = "identity", position = "stack") + theme_bw() +labs(x="")
#Fastz with tse RBH
dndsAll <- read.table("psi2tse.dnds.pair.txt")
geneloci <- read.table("trans.loci")
zsub <- dndsAll %>% filter(V1 %in% subset(geneloci, V1 == "chrZ")$V4) %>% 
  filter(V4 >= 0.01 & V4 <= 1) %>% mutate(
    chr = case_when(
      V7 <= 1.5e6 ~ "Z-PAR",
      V6 > 1.5e6 & V7 < 3.7e6 ~ "Z-S2",
      V6 > 3.7e6 & V7 < 14.5e6 ~ "Z-S1",
      V6 > 14.5e6 ~ "Z-S0"
    ))
wsub <- dndsAll %>% filter(V1 %in% subset(geneloci, V1 == "chrW")$V4) %>% 
  filter(V4 >= 0.01 & V4 <= 1) %>%
  mutate(chr = "W chromosome") 
asub <- dndsAll %>% 
  filter(V1 %in% subset(geneloci, V1 != "chrZ" & V1 != "chrW" & !grepl("Hic", V1))$V4) %>%
  filter(V4 >= 0.01 & V4 <= 1) %>% 
  mutate(chr = gsub("chr","", V5),
         chr = ifelse(chr < 9, "micro", "macro")) 
dfsub <- rbind(asub,zsub,wsub)
ggplot(dfsub,aes(x = chr, y=V2,fill=chr))+
  geom_boxplot(outlier.colour = NULL,
               width = .4) + theme_bw()+ 
  coord_cartesian(ylim = c(0,2)) + 
  stat_compare_means(ref.group = "Z-S2", method = "wilcox.test",
                     method.args = list(alternative = "less"))  +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        axis.text.x= element_text(size=10),
        axis.text.y= element_text(size=15))+
  labs(x="",y="dN/dS")
#rDNA clusters with reeves rRNA seq
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
  geom_step() + xlim(0,40) + theme_classic() + labs(x="", y="Length (bp)")
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

} # Genome information and sex chromosome evolution of Chinese soft shell turtle
{
#PCA of chinese softshell turtles
library(plot3D)
tt <- read.csv("rpkm_gonad_add1213.csv", row.names = 1, header = TRUE) 
tt_order <- tt[,order(colnames(tt))]
annotation_col<- read.csv("sample_list.csv",row.names = 5,header = TRUE)
annotation_col <- annotation_col[,-1]
data.ttB <- tt_order[,grep("FG_*|MG_*|FM12*|FM13*|MM12*|MM13*",colnames(tt_order))]
sort.ttB <- data.ttB[,sort(colnames(data.ttB))]
filt.ttB <- sort.ttB[sapply(1:nrow(sort.ttB), function(x) 
  length(which(sort.ttB[x,c(1:32)] >=1))>=1)==T,]
log.ttB <- log2(filt.ttB+0.01)
annotation_colB <- annotation_col[grep("FG_*|MG_*|FM12_*|FM13_*|MM12_*|MM13_*",rownames(annotation_col)),]
pca_dataB=prcomp(t(log.ttB))
summary(pca_dataB)
plot(pca_dataB, type="lines")
pca_dataB_perc=round(100*pca_dataB$sdev^2/sum(pca_dataB$sdev^2),1)
df_pca_B = data.frame(PC1 = pca_dataB$x[,1], PC2 = pca_dataB$x[,2],PC3=pca_dataB$x[,3],
                      sample = colnames(log.ttB), 
                      Sex = annotation_colB$Sex,
                      Stage = annotation_colB$Stage,
                      tissue=annotation_colB$Tissue,
                      labell=rownames(annotation_colB))
df_pca_B$gender <- rep(c(18,17),c(16,16))
scatter3D(df_pca_B$PC1,df_pca_B$PC2,df_pca_B$PC3,bty="u",col.axis="grey",axes=FALSE,pch=16,colvar = as.integer(df_pca_B$gender), 
          col=c('#cd0000','#0000cd'), main = "tse gonad",cex=df_pca_B$Stage/5, 
          colkey = FALSE, theta =-50,phi=30, xlab = paste0("PC1 (",pca_dataB_perc[1],")"),
          ylab =paste0("PC2 (",pca_dataB_perc[2],")"), zlab = paste0("PC3 (",pca_dataB_perc[3],")"))
text3D(df_pca_B$PC1,df_pca_B$PC2,df_pca_B$PC3,labels = df_pca_B$Stage,add = TRUE,cex=2)
#PCA of red-eared slider turtles gonad (8 stages)
tt <- read.table(file = "red.rpkm.txt", sep=",", row.names=1, header=T,
                 quote="", comment="", check.names=F)

tt_order <- tt[,order(colnames(tt))]
annotation_colB<- data.frame(row.names =  colnames(tt_order),Sex=rep(c("female","male"),c(15,15)),Stage=as.numeric(gsub("_1|_2","",gsub("Fs|Ms", "", colnames(tt_order)))))
annotation_colB$Tissue <- "Gonad"
data.ttB <- tt_order[,grep("Ms_*",colnames(tt_order))]
sort.ttB <- data.ttB[,sort(colnames(data.ttB))]
filt.ttB <- sort.ttB[sapply(1:nrow(sort.ttB), function(x) 
  length(which(sort.ttB[x,c(1:15)] >=1))>=1)==T,]
log.ttB <- log2(filt.ttB+0.01)
annotation_colB <- annotation_colB[grep("Ms_*",rownames(annotation_colB)),]
pca_dataB=prcomp(t(log.ttB))
summary(pca_dataB)
plot(pca_dataB, type="lines")
pca_dataB_perc=round(100*pca_dataB$sdev^2/sum(pca_dataB$sdev^2),1)
df_pca_B = data.frame(PC1 = pca_dataB$x[,1], PC2 = pca_dataB$x[,2],PC3=pca_dataB$x[,3],
                      sample = colnames(log.ttB), 
                      Sex = annotation_colB$Sex,
                      Stage = annotation_colB$Stage,
                      tissue=annotation_colB$Tissue,
                      labell=rownames(annotation_colB))
df_pca_B$gender <- rep(c(7,8),c(7,8))
scatter3D(df_pca_B$PC1,df_pca_B$PC2,df_pca_B$PC3,bty="u",col.axis="grey",axes=FALSE,pch=16,colvar = as.integer(df_pca_B$gender), 
          col=c('#cd0000','#0000cd'), main = "tse gonad",cex=as.numeric(df_pca_B$Stage)/5, 
          colkey = FALSE, theta =-90,phi=90, xlab = paste0("PC1 (",pca_dataB_perc[1],")"),
          ylab =paste0("PC2 (",pca_dataB_perc[2],")"), zlab = paste0("PC3 (",pca_dataB_perc[3],")"))
text3D(df_pca_B$PC1,df_pca_B$PC2,df_pca_B$PC3,labels = df_pca_B$Stage,add = TRUE,cex=2)
#PCA of tse simatic tissues(As controls with only 5 stages)
exprMat <- "tse.new.tpm.txt"
dataExpr <- read.table(file = exprMat, sep=' ', row.names=1, header=T)
colnames(dataExpr) <- colnames(dataExpr) %>%
  str_split("_", simplify = TRUE) %>%
  apply(1, function(x) paste0(x[1], x[3], x[2], x[4]))

data2 <- dataExpr %>% select(matches("FM|MM")) %>%  select(!matches("FM15rep1|FM12rep3")) %>% filter_all(.,any_vars(.>1)) 
dataExprsub <- log2(data2+0.1)
dataExprsub <- as.data.frame(t(dataExprsub))
nGenes = ncol(dataExprsub)
nSamples = nrow(dataExprsub)
dim(dataExprsub)
sampleTree = hclust(dist(dataExprsub), method = "average")
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="")
scaled_data <- scale(dataExprsub) 


pca_result <- prcomp(scaled_data, center = TRUE, scale. = TRUE)
pca_scores <- as.data.frame(pca_result$x) %>%
    mutate(
      sex = substr(rownames(.), 1, 1),   
      tissue = substr(rownames(.), 2, 2),  
      stage = substr(rownames(.), 3, 4),
      gender = case_when( sex == "F" ~ "1",
                          sex == "M" ~ "2")
    )
ggplot(pca_scores, aes(x = PC1, y = PC3, color = sex, shape = tissue, size = stage)) +
    geom_point() + geom_text(label = gsub("FM|MM","",rownames(pca_scores)),size = 5)+
    labs(x = "PC1", y = "PC2", title = "PCA of Gene Expression Matrix") +
    theme_minimal() + scale_color_manual(values = c("F"="#cd0000","M"="#0000cd"))

scatter3D(pca_scores$PC1,pca_scores$PC2,pca_scores$PC3,bty="u",col.axis="grey",
          pch=16,colvar = as.integer(pca_scores$gender),col=c('#cd0000','#0000cd'),
          cex=as.numeric(pca_scores$stage)/5,
          theta = 70,phi= 45)
text3D(pca_scores$PC1,pca_scores$PC2,pca_scores$PC3,labels = pca_scores$stage,add = TRUE,cex=2)

} # Data quality control of two turtles
{
{
psiExp <- read.table(file = "psi.tpmmean.NewFilt.txt", sep=',', row.names=1, header=T) 
tseExp <- read.table(file = "tse.tpm.mean.allTissue.txt", header = T, row.names = 1)
rownames(tseExp) <- tseExp$Gene
tseExp <- tseExp[, -1]
geneidproduct <- read.table("psi.geneID.finalinfo.txt")
orthologs <- read.table("psi-red.orth-1to1.50")
#Exp. normalization between species by UQ
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
ei <- as.matrix(mergedInfo)
er2 <- as.data.frame(UQ_FN(ei))
psisub <- er2[1:20]
colnames(psisub) <- c(paste("FG",c("12","15","16","18","21"), sep = ""),
                      paste("MG",c("12","15","16","18","21"), sep = ""),
                      paste("FB",c("12","15","16","18","21"), sep = ""),
                      paste("MB",c("12","15","16","18","21"), sep = ""))
tsesub <- er2[21:40]
colnames(tsesub) <- colnames(psisub)
#Psi Stage 15 Exp. increase fitting curves
colnames(psiExp) <- gsub("GMC", "G", colnames(psiExp))
colnames_data <- colnames(psiExp)
psiExpGd <- psiExp %>% select(grep("FG|MG", colnames(.))) %>%
  select(grep("12|15|16|18|21", colnames(.))) %>%
  filter_all(.,any_vars(.> 1))
targetdata <- psiExpGd %>% 
  filter(FG15/FG12 >= 2 & FG15/FG16 >= 2) 
targetdata_scaled <- data.frame(t(scale(t(targetdata)))) %>% tibble::rownames_to_column("gene") %>% 
  gather(key = "variables", value = "tpm",- gene) %>% mutate(
    sex = substr(variables, 1, 1),   
    stage = as.numeric(substr(variables, 3, 4)),
    chr = geneidproduct$V1[match(gene, geneidproduct$V4)],
    type = ifelse(chr %in% c("chrZ", "chrW"), "Sex Chr", "Other"),
    spe = "Psi"
  )
tseExp <- tseExp %>% select(grep("FG|MG", colnames(.))) %>%
  select(grep("12|15|16|18|21", colnames(.))) %>%
  filter_all(.,any_vars(.> 1)) %>% filter(rownames(.) %in% orthologs$V2)
rownames(tseExp) <- orthologs$V1[match(row.names(tseExp),orthologs$V2)]
targetdata2 <- tseExp %>% filter(rownames(.) %in% rownames(targetdata))
targetdata2_scaled <- data.frame(t(scale(t(targetdata2)))) %>% tibble::rownames_to_column("gene") %>% 
  gather(key = "variables", value = "tpm",- gene) %>% mutate(
    sex = substr(variables, 1, 1),   
    stage = as.numeric(substr(variables, 3, 4)),
    chr = geneidproduct$V1[match(gene, geneidproduct$V4)],
    type = ifelse(chr %in% c("chrZ", "chrW"), "Sex Chr", "Other Chr"),
    spe = "Tse"
  )
tt_trajectory <- rbind(targetdata_scaled, targetdata2_scaled) %>% mutate(
  factor = paste(spe, type, sep = "_"))
ggplot(tt_trajectory, aes(x = as.character(stage), y = tpm, group = sex, color = sex)) +
  geom_smooth(method = "loess", span = .5) + facet_wrap(~ factor) +
  labs(x = "Stage", y = "TPM after Z-scaled", title = "Relative expression trajectory") +
  theme(legend.position = "none") + scale_color_manual(values = c("F" = "#ed7d61", "M" = "#8cb7df")) +
  theme_classic()
ggplot(tt_trajectory, aes(x = as.character(stage), y = tpm, group = sex, color = sex)) +
  geom_smooth(method = "loess", span = .5) + facet_wrap(~ spe) +
  labs(x = "Stage", y = "TPM after Z-scaled", title = "Relative expression trajectory") +
  theme(legend.position = "none") + scale_color_manual(values = c("F" = "#ed7d61", "M" = "#8cb7df")) +
  theme_classic()
#Controls in brains and kidney
psiExpCt <- psiExp %>% select(grep("FB|MB|FM|MM|FG12|MG12", colnames(.))) %>%
  select(grep("12|15|16|18|21", colnames(.))) %>%
  filter_all(.,any_vars(.> 1))
targetdataCt <- psiExpCt %>% 
  filter(rownames(.) %in% rownames(targetdata)) 
targetdataCt_scaled <- data.frame(t(scale(t(targetdataCt)))) %>% tibble::rownames_to_column("gene") %>% 
  gather(key = "variables", value = "tpm",- gene) %>% mutate(
    sex = substr(variables, 1, 1),   
    stage = as.numeric(substr(variables, 3, 4)),
    tissue = substr(variables, 2, 2),
    tissue = gsub("G", "M", tissue),
    chr = geneidproduct$V1[match(gene, geneidproduct$V4)],
    type = ifelse(chr %in% c("chrZ", "chrW"), "Sex Chr", "Other")
  )
ggplot(targetdataCt_scaled, aes(x = as.character(stage), y = tpm, group = sex, color = sex)) +
  geom_smooth(method = "loess", span = .5) + facet_wrap( ~ tissue) +
  labs(x = "Stage", y = "TPM after Z-scaled", title = "Relative expression trajectory") +
  theme(legend.position = "none") + scale_color_manual(values = c("F" = "#ed7d61", "M" = "#8cb7df")) +
  theme_classic()
#How many male genes and female genes in SD list
sdlist <- read.csv("newSDRelaetedGenes.csv", row.names = 1)
geneForGo <- unique(gsub("_[0-9]+$", "", geneidproduct %>% filter(V4 %in% targetdata_scaled$gene) %>% pull(V6)))
sdlistSub <- sdlist %>% filter(Gene %in% geneForGo)

#Gene location by chromosomes
allcounts <- data.frame(table(subset(geneidproduct, V4 %in% rownames(psiExpGd))$V1)) %>%
  filter(grepl("chr", Var1)) # Only expressed genes!
s15counts <- data.frame(table(geneidproduct$V1[geneidproduct$V4 %in% rownames(targetdata)])) %>%
  filter(grepl("chr", Var1))
tarprop <- merge(allcounts, s15counts, by = "Var1", all = TRUE) %>%
  mutate(
    type = case_when(
      Var1 %in% paste("chr", rep(1:8), sep = "") ~ "macro",  
      Var1 %in% paste("chr", rep(9:32), sep = "") ~ "micro", 
      TRUE ~ "Sex Chr"  
    )
  ) %>%
  group_by(type) %>%
  summarise(
    all = sum(Freq.x, na.rm = TRUE), 
    s15 = sum(Freq.y, na.rm = TRUE)  
  )
ggplot(tarprop, aes(x=type,y=s15/all)) + geom_bar(stat = "identity")
#Fast chi-sq test
data <- matrix(c(44, 357-44, 558, 8002-558), nrow = 2, byrow = TRUE)
rownames(data) <- c("Group A", "Group B")
colnames(data) <- c("Increase", "Other")
chisq.test(data)
#Merged UQ expression curve of two turtles
psidata_long <- psisub %>% tibble::rownames_to_column("gene") %>% 
  filter(gene %in% rownames(targetdata)) %>% 
  gather(key = "variables", value = "tpm",- gene) %>% mutate(
    sex = substr(variables, 1, 1),   
    stage = as.numeric(substr(variables, 3, 4)),
    spe = "Psi"
  )
#Tse Dat.
tsedata_long <- tsesub %>% tibble::rownames_to_column("gene") %>%
  filter(gene %in% rownames(targetdata)) %>% 
  gather(key = "variables", value = "tpm",- gene) %>% mutate(
    sex = substr(variables, 1, 1),   
    stage = as.numeric(substr(variables, 3, 4)),
    spe = "Tse"
  )
#Plot
ttdata_long <- rbind(psidata_long,tsedata_long) %>% mutate(
  log2tpm = log(tpm+1),
  factor =  paste(spe, sex,sep = "_")
)
ggplot(ttdata_long, aes(x = factor, y = log2tpm)) +
  geom_boxplot(
    notch = TRUE, width = .7,
    aes(
      fill = ifelse(spe == "Psi" & sex == "F", "#cd0000", 
                    ifelse(spe == "Psi" & sex == "M", "#0000cd", "white")),
      color = ifelse(spe == "Psi" & sex == "F", "#cd0000", 
                     ifelse(spe == "Psi" & sex == "M", "#0000cd", 
                            ifelse(spe == "Tse" & sex == "F", "#cd0000", "#0000cd")))
    ),
    outlier.shape = 16, outlier.size = 1,
    alpha = 0.6
  ) + 
  facet_wrap(~ stage, nrow = 1) +
  theme_minimal(base_size = 14) +
  labs(x = "", y = "UQ (TPM)", title = "") + 
  scale_color_identity() +  # Ensure that color is treated as actual color
  scale_fill_identity() +  # Ensure that fill is treated as actual color
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    strip.text = element_text(size = 16),
    plot.title = element_text(hjust = 0.5, size = 18),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  ) +
  stat_compare_means(
    method = "wilcox.test", 
    label = "p.signif", 
    comparisons = list(c("Psi_F", "Psi_M"),
                       c("Tse_F", "Tse_M")), 
    method.args = list(alternative = "greater"),
    size = 5,
    label.x.npc = "center",
    label.y.npc = "top"
  ) 
#Specific location of 1187 genes
fg15prop <- geneidproduct %>% 
  mutate(
    type = ifelse(V4 %in% unique(targetdata_scaled %>% pull(gene)), "Up", "Equal"),
    strata = case_when(
      V1 == "chrZ" & V3 > 12.5e6 ~ "ZS0",
      V1 == "chrZ" & V3 <= 12.5e6 ~ "ZOther",
      V1 %in% paste0("chr", 1:8) ~ "macro",
      V1 %in% paste0("chr", 9:32) ~ "micro",
      TRUE ~ "Other"
    )
  ) %>% filter(strata != "Other") %>% filter(V4 %in% rownames(psiExpGd)) %>%
  group_by(type, strata) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(strata) %>%
  mutate(ratio = count / sum(count)) %>% filter(type == "Up")
ggplot(fg15prop, aes(x = strata, y = ratio, fill = type)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("Up" = "#ed7d61", "Equal" = "#8cb7df")) +
  labs(x = "Strata", y = "Proportion", fill = "Type") +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
fg15table <- geneidproduct %>%
  mutate(
    type = ifelse(V4 %in% unique(targetdata_scaled %>% pull(gene)), "Up", "Equal"),
    strata_group = case_when(
      V1 == "chrZ" & V3 > 12.5e6 ~ "ZS0",
      TRUE ~ "Other"
    )
  )
contingency_table <- table(fg15table$strata_group, fg15table$type)
chisq.test(contingency_table)
#Specific location of sex chromosomes
scGenes <- unique(targetdata_scaled %>% filter(type == "Sex Chr") %>% select(gene))
geneZsub <- geneidproduct %>% filter(V1 == "chrZ") %>%
  mutate(
    type = ifelse(V4 %in% scGenes$gene, "Up", "Equal"),
    strata = ifelse(V3 > 12.5e6,"S0", "Other")
  )
freq_data <- geneZsub  %>%
  group_by(type, strata) %>%
  summarise(count = n(), .groups = "drop")
PieDonut(freq_data, aes(strata,type, count=count), title = "",
         r0 = 0.6, r1 = 0.95)
#Fast chi-sq test
data <- matrix(c(30, 180, 14, 247), nrow = 2, byrow = TRUE)
rownames(data) <- c("Group A", "Group B")
colnames(data) <- c("Increase", "Other")
chisq.test(data)
#What are those genes GO enrichments
library(clusterProfiler)
library(org.Hs.eg.db)
library(cols4all)
sdlist <- read.csv("newSDRelaetedGenes.csv", row.names = 1)
geneForGo <- unique(gsub("_[0-9]+$", "", geneidproduct %>% filter(V4 %in% targetdata_scaled$gene) %>% pull(V6)))
notsdlistsub <- setdiff(geneForGo, unique(sdlist$Gene))
sdlistsub <- intersect(geneForGo, unique(sdlist$Gene))
targetGenes <- bitr(sdlistsub,
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
#Go
GOresult <- ego@result %>% head(10) %>%
  mutate(
    GeneRatio_numeric = as.numeric(sapply(strsplit(as.character(GeneRatio), "/"), function(x) as.numeric(x[1])/as.numeric(x[2]))),
  )
ggplot(data = GOresult, aes(x = GeneRatio_numeric, y = factor(Description, levels = rev(GOresult$Description)), fill = -log10(p.adjust))) +
  geom_bar(stat = 'identity', width = 0.8, alpha = 0.9) +  
  scale_fill_gradient(low = "#fff6f1",high = "#cd0000")+
  labs(x = 'Gene Ratio', y = '') +
  geom_text(aes(x = 0.005, 
                label= Description),
            hjust= 0) +
  theme_classic() + theme(axis.text.y = element_blank()) 
#Kegg
ekegg <- enrichKEGG(gene = targetGenes$ENTREZID, 
                    organism = 'hsa',
                    pvalueCutoff = 0.05,
                    pAdjustMethod = "fdr")
keggresult <- ekegg@result
keggresult$FoldEnrichment <- apply(keggresult,1,function(x){
  GeneRatio=eval(parse(text=x["GeneRatio"]))
  BgRatio=eval(parse(text=x["BgRatio"]))
  enrichment_fold=round(GeneRatio/BgRatio,2)
  enrichment_fold
})
keggresult <- mutate(keggresult, RichFactor= Count / as.numeric(sub("/\\d+", "", BgRatio)))
geneID_list <- strsplit(keggresult$geneID, "/")
gene_symbols <- lapply(geneID_list, function(ids) {
  mapped_genes <- bitr(ids, fromType = "ENTREZID", toType = "SYMBOL", OrgDb = org.Hs.eg.db)
  return(mapped_genes$SYMBOL) 
})
keggresult$geneSymbols <- sapply(gene_symbols, function(symbols) paste(symbols, collapse = "/"))
keggresult$Description <- factor(keggresult$Description,levels= rev(keggresult$Description))
kegg1 <- keggresult[order(keggresult$Count, keggresult$pvalue, decreasing = T),]
ggplot(data = kegg1[1:10,], aes(x = Count, y = Description, fill = -log10(p.adjust))) +
  geom_bar(stat = 'identity', width = 0.8, alpha = 0.9) +  
  scale_fill_gradient(low = "#f3f9fd",high = "#0000cd")+
  labs(x = 'Number of Genes', y = '') +
  geom_text(aes(x = .5, 
                label= Description),
            hjust= 0) +
  theme_classic() + theme(axis.text.y = element_blank()) 
#check FG15 down micros
expmicro <- read.table(file="psi.sRNA.tpmmean.final.shortstack.txt")
expmicroGd <- expmicro %>% select(grep("FG|MG|FB|MB",colnames(.))) %>% 
  select(grep("12|15|16|18|21", colnames(.))) %>% filter_all(.,any_vars(.>=1))
expmicrosub <- subset(expmicroGd, expmicroGd$FG15/expmicroGd$FG12 <= .5 
                      & expmicroGd$FG15/expmicroGd$FG16 <= .5) # Expression fold change <= 0.5 against nearby stages
expmicrosub_scaled <- data.frame(t(scale(t(expmicrosub)))) %>% tibble::rownames_to_column("gene") %>% 
  gather(key = "variables", value = "tpm",- gene) %>% mutate(
    sex = substr(variables, 1, 1),   
    stage = as.numeric(substr(variables, 3, 4)),
    tissue = substr(variables, 2, 2),
  )
ggplot(expmicrosub_scaled, aes(x = as.character(stage), y = tpm, group = sex, color = sex)) +
  geom_smooth(method = "loess", span = .5) + facet_wrap( ~ tissue) +
  labs(x = "Stage", y = "TPM after Z-scaled", title = "Relative expression trajectory") +
  theme(legend.position = "none") + scale_color_manual(values = c("F" = "#ed7d61", "M" = "#8cb7df")) +
  theme_classic()
{
##check % 8mer/multisites of those genes of micros
targenelist <- read.table(file="psi.automir.tar.txt")
targenelist$mir <- paste(data.frame(do.call('rbind', strsplit(as.character(targenelist$V2),'_',fixed=TRUE)))$X1,
                         data.frame(do.call('rbind', strsplit(as.character(targenelist$V2),'_',fixed=TRUE)))$X2,
                         sep="_")
targenesub <- subset(targenelist,targenelist$mir %in% rownames(expmicrosub))
targenesub$factor <- interaction(targenesub$V2,targenesub$V4,targenesub$V5,sep=":") 

targenesub$type <- ifelse(!(targenesub$V1 %in% rownames(expsub2)), "Irrelevant","Related")
##split gene sets 
nomalsub <- subset(targenesub,targenesub$type=="Irrelevant")
sdsub <- subset(targenesub,targenesub$type=="Related")
##irrelevant genes
df <- data.frame(factor=unique(nomalsub$factor)) # all bind sites sum up
df$mir <- data.frame(do.call('rbind', strsplit(as.character(df$factor),':',fixed=TRUE)))$X1
dfnew <- data.frame(table(df$mir)) 
sub8mer <- subset(nomalsub,grepl("8mer", V3)) # 8mer sites
df1.1 <- data.frame(factor=unique(sub8mer$factor))
submulti <- nomalsub #multi-sites check
duplicates <- duplicated(submulti$factor) | duplicated(submulti$factor, fromLast = TRUE)
duplicate_rows <- submulti[duplicates, ]
df1.2 <- data.frame(factor=unique(duplicate_rows$factor))
df1 <- unique(rbind(df1.1,df1.2)) # intersect 8mer/multi-sites target
df1$mir <- data.frame(do.call('rbind', strsplit(as.character(df1$factor),':',fixed=TRUE)))$X1
df1new <- data.frame(table(df1$mir)) # 8mer/multi-sites sum up
merged_df <- merge(dfnew, df1new, by = "Var1", all = TRUE)
colnames(merged_df) <- c("mir","AllsiteTar","MultiSitesOr8MerTar")
merged_df$type <- "Irrelevant"
##Important SD genes
df <- data.frame(factor=unique(sdsub$factor))
df$mir <- data.frame(do.call('rbind', strsplit(as.character(df$factor),':',fixed=TRUE)))$X1
dfnew <- data.frame(table(df$mir))
sub8mer <- subset(sdsub,grepl("8mer", V3)) # 8mer sites
df1.1 <- data.frame(factor=unique(sub8mer$factor))
submulti <- sdsub #multi-sites check
duplicates <- duplicated(submulti$factor) | duplicated(submulti$factor, fromLast = TRUE)
duplicate_rows <- submulti[duplicates, ]
df1.2 <- data.frame(factor=unique(duplicate_rows$factor))
df1 <- unique(rbind(df1.1,df1.2)) # intersect 8mer/multi-sites target
df1$mir <- data.frame(do.call('rbind', strsplit(as.character(df1$factor),':',fixed=TRUE)))$X1
df1new <- data.frame(table(df1$mir)) # 8mer/multi-sites sum up
merged_df2 <- merge(dfnew, df1new, by = "Var1", all = TRUE)
colnames(merged_df2) <- c("mir","AllsiteTar","MultiSitesOr8MerTar")
merged_df2$type <- "SD related"
##Visualize
merged_dfall <- rbind(merged_df,merged_df2)
ggplot(merged_dfall, aes(x = mir, y = (MultiSitesOr8MerTar/AllsiteTar)*100,fill=type)) +
  geom_bar(stat = "identity", position = position_dodge(.7), width = 0.5) + 
  labs(title = "Grouped Barplot",
       x = "mir",
       y = "value",
       fill = "chr") +theme_classic()+scale_fill_manual(values=c("#4daf4a","#0000cd"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
##Boxplot of all
ggplot(merged_dfall, aes(x = type, y = (MultiSitesOr8MerTar/AllsiteTar)*100)) +
  geom_boxplot()
##Significance check
merged_dfall <- merged_dfall[order(merged_dfall$mir), ]
merged_dfall[is.na(merged_dfall)] <- 0
results_list <- list()
for (i in seq(1, nrow(merged_dfall), by = 2)) {
  values <- as.matrix(merged_dfall[i:(i+1), c("AllsiteTar", "MultiSitesOr8MerTar")])
  matrix <- matrix(c((values[1,1]-values[1,2]),values[1,2],
                     (values[2,1]-values[2,2]),values[2,2]),nrow=2,ncol=2)
  result <- fisher.test(matrix) ##use fisher here
  results_list[[length(results_list) + 1]] <- list(matrix_values = values, chisq_result = result)
}
for (res in results_list) {
  print(res$chisq_result)
}

} # May due to Micro targeting
} # Psi S15 ovary up-regulation
{
psiExp <- read.table(file = "psi.tpmmean.NewFilt.txt", sep=',', row.names=1, header=T) 
tseExp <- read.table(file = "tse.tpm.mean.allTissue.txt", header = T, row.names = 1)
rownames(tseExp) <- tseExp$Gene
tseExp <- tseExp[, -1]
geneidproduct <- read.table("psi.geneID.finalinfo.txt")
orthologs <- read.table("psi-red.orth-1to1.50")
#Exp. normalization between species by UQ
psiExpsub <- psiExp[rownames(psiExp) %in% orthologs$V1, ] %>%
  select(grep("FG|MG", colnames(.))) %>% 
  select(grep("12|15|16|18|21", colnames(.))) %>% 
  select(5, 1:4, 10, 6:9)
tseExpsub <- tseExp[rownames(tseExp) %in% orthologs$V2, ] %>%
  select(grep("FG|MG", colnames(.))) %>% 
  select(grep("12|15|16|18|21", colnames(.)))
tseExpsub$orthname <- orthologs$V1[match(row.names(tseExpsub),orthologs$V2)]
tseExpsub <- tseExpsub[tseExpsub$orthname %in% orthologs$V1,]
row.names(tseExpsub) <- tseExpsub$orthname
tseExpsub <- tseExpsub[,-11]
psiExpsub <- psiExpsub[rownames(psiExpsub) %in% orthologs$V1,]
mergedInfo <-  merge(psiExpsub, tseExpsub, by='row.names', all=TRUE)
row.names(mergedInfo) <- mergedInfo$Row.names
mergedInfo <- mergedInfo[,-1]
ei <- as.matrix(mergedInfo)
er2 <- as.data.frame(UQ_FN(ei))
psisub <- er2[1:10]
colnames(psisub) <- c(paste("FG",c("12","15","16","18","21"), sep = ""),paste("MG",c("12","15","16","18","21"), sep = ""))
tsesub <- er2[11:20]
colnames(tsesub) <- c(paste("FG",c("12","15","16","18","21"), sep = ""),paste("MG",c("12","15","16","18","21"), sep = ""))
#Tse Stage 15 Exp. decrease fitting curves
tseExpGd <- tseExp %>% select(grep("FG|MG", colnames(.))) %>%
  select(grep("12|15|16|18|21", colnames(.))) %>%
  filter_all(.,any_vars(.> 1))
targetdata <- tseExpGd %>% 
  filter(FG15/FG12 <= .5 & FG15/FG16 <= .5) 
targetdata_scaled <- data.frame(t(scale(t(targetdata)))) %>% tibble::rownames_to_column("gene") %>% 
  gather(key = "variables", value = "tpm",- gene) %>% mutate(
    sex = substr(variables, 1, 1),   
    stage = as.numeric(substr(variables, 3, 4)),
    spe = "Tse"
  )
colnames(psiExp) <- gsub("GMC", "G", colnames(psiExp))
colnames_data <- colnames(psiExp)
psiExp <- psiExp %>% select(grep("FG|MG", colnames(.))) %>%
  select(grep("12|15|16|18|21", colnames(.))) %>%
  filter_all(.,any_vars(.> 1)) %>% filter(rownames(.) %in% orthologs$V1)
rownames(psiExp) <- orthologs$V2[match(row.names(psiExp),orthologs$V1)]
targetdata2 <- psiExp %>% filter(rownames(.) %in% rownames(targetdata))
targetdata2_scaled <- data.frame(t(scale(t(targetdata2)))) %>% tibble::rownames_to_column("gene") %>% 
  gather(key = "variables", value = "tpm",- gene) %>% mutate(
    sex = substr(variables, 1, 1),   
    stage = as.numeric(substr(variables, 3, 4)),
    spe = "Psi"
  )
tt_trajectory <- rbind(targetdata_scaled, targetdata2_scaled) 
ggplot(tt_trajectory, aes(x = as.character(stage), y = tpm, color = sex, group = sex)) +
  geom_smooth(method = "loess", span = .3) + facet_wrap(~ spe) + theme_bw()+
  labs(x = "Stage", y = "TPM after Z-scaled", title = "Relative expression trajectory") +
  theme(legend.position = "none") + scale_color_manual(values = c("F" = "#cd0000", "M" = "#0000cd")) 
#Controls in brains and kidney
tseExpCt <- tseExp %>% select(grep("FB|MB|FM|MM", colnames(.))) %>%
  select(grep("12|15|16|18|21", colnames(.))) %>%
  filter_all(.,any_vars(.> 1))
targetdataCt <- tseExpCt %>% 
  filter(rownames(.) %in% rownames(targetdata)) 
targetdataCt_scaled <- data.frame(t(scale(t(targetdataCt)))) %>% tibble::rownames_to_column("gene") %>% 
  gather(key = "variables", value = "tpm",- gene) %>% mutate(
    sex = substr(variables, 1, 1),   
    stage = as.numeric(substr(variables, 3, 4)),
    tissue = substr(variables, 2, 2),
    tissue = gsub("G", "M", tissue),
    chr = geneidproduct$V1[match(gene, geneidproduct$V4)],
    type = ifelse(chr %in% c("chrZ", "chrW"), "Sex Chr", "Other")
  )
ggplot(targetdataCt_scaled, aes(x = as.character(stage), y = tpm, group = sex, color = sex)) +
  geom_smooth(method = "loess", span = .3) + facet_wrap( ~ tissue) + ylim(-1.5,0.5)+
  labs(x = "Stage", y = "TPM after Z-scaled", title = "Relative expression trajectory") +
  theme(legend.position = "none") + scale_color_manual(values = c("F" = "#ed7d61", "M" = "#8cb7df")) +
  theme_classic()

#How many male genes and female genes in SD list
sdlist <- read.csv("newSDRelaetedGenes.csv", row.names = 1)
geneForGo <- unique(targetdata_scaled$gene)
sdlistSub <- sdlist %>% filter(Gene %in% geneForGo)









#Location of genes
tsegeneloci <- read.table(file = "tse.gene.loci") 
tsechrlen <- read.table(file = "tse.chrlen") 
tsegeneloci$type <- ifelse(tsegeneloci$V4 %in% rownames(targetdata), "Suppression", "Other")
freq_data <- tsegeneloci  %>%
  group_by(V1, type) %>%
  summarise(count = n(), .groups = "drop") 
freq_dataMerge <- freq_data %>% mutate(
  chrType = ifelse(V1 %in% paste(rep("chr",13),c(1:13), sep = ""),
                   "Macro", "Micro")
) %>%
  group_by(chrType, type) %>%
  summarise(count = n(), .groups = "drop")
freq_ratio <- freq_data %>%
  spread(key = type, value = count, fill = 0) %>%  
  mutate(ratio = Suppression / Other)
ggplot(freq_data, aes(x= factor(V1, levels = c(paste("chr",rep(1:25),sep = ""))), y=count,fill=type)) +
  geom_bar(stat = "identity") + theme_bw() + labs(x="",y="Gene number") + 
  scale_fill_manual(values = c("Suppression" = "black", "Other" = "#BDBDBD"))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggplot(freq_dataMerge, aes(x= chrType, y=count,fill=type)) +
  geom_bar(stat = "identity") + theme_bw() + labs(x="",y="Gene number") + 
  scale_fill_manual(values = c("Suppression" = "black", "Other" = "#BDBDBD"))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggplot(freq_ratio, aes(x= factor(V1, levels = c(paste("chr",rep(1:25),sep = ""))), y=ratio*100)) +
  geom_bar(stat = "identity") + theme_bw() + labs(x="",y="% Supression genes") 
#Location on chromosomes
suppression_genes <- tsegeneloci %>%
  filter(type == "Suppression")
ggplot(tsegeneloci) +
  geom_segment(aes(x = 0, xend = V2, y = V1, yend = V1)) +
  geom_segment(data = suppression_genes, aes(x = V2, xend = V3, y = V1), size = 2, color = "red") + 
  theme_minimal() +
  labs(x = "Position", y = "Chromosome", title = "Gene Locations for Suppression with Chromosome Lengths") 
#What are those genes GO enrichments
library(clusterProfiler)
library(org.Hs.eg.db)
library(cols4all)
sdlist <- read.csv("newSDRelaetedGenes.csv", row.names = 1)
geneForGo <- unique(targetdata_scaled$gene)
notsdlistsub <- setdiff(geneForGo, unique(sdlist$Gene))
sdlistsub <- intersect(geneForGo, unique(sdlist$Gene))
targetGenes <- bitr(geneForGo,
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
#Go
GOresult <- ego@result %>% head(10) %>%
  mutate(
    GeneRatio_numeric = as.numeric(sapply(strsplit(as.character(GeneRatio), "/"), function(x) as.numeric(x[1])/as.numeric(x[2]))),
  )
ggplot(data = GOresult, aes(x = GeneRatio_numeric, y = factor(Description, levels = rev(GOresult$Description)), fill = -log10(p.adjust))) +
  geom_bar(stat = 'identity', width = 0.8, alpha = 0.9) +  
  scale_fill_gradient(low = "#fff6f1",high = "#cd0000")+
  labs(x = 'Gene Ratio', y = '') +
  geom_text(aes(x = 0.005, 
                label= Description),
            hjust= 0) +
  theme_classic() + theme(axis.text.y = element_blank()) 
#Kegg
ekegg <- enrichKEGG(gene = targetGenes$ENTREZID, 
                    organism = 'hsa', keyType = 'kegg' ,
                    pvalueCutoff = 0.05,
                    pAdjustMethod = "fdr")
keggresult <- ekegg@result
keggresult$FoldEnrichment <- apply(keggresult,1,function(x){
  GeneRatio=eval(parse(text=x["GeneRatio"]))
  BgRatio=eval(parse(text=x["BgRatio"]))
  enrichment_fold=round(GeneRatio/BgRatio,2)
  enrichment_fold
})
keggresult <- mutate(keggresult, RichFactor= Count / as.numeric(sub("/\\d+", "", BgRatio)))
geneID_list <- strsplit(keggresult$geneID, "/")
gene_symbols <- lapply(geneID_list, function(ids) {
  mapped_genes <- bitr(ids, fromType = "ENTREZID", toType = "SYMBOL", OrgDb = org.Hs.eg.db)
  return(mapped_genes$SYMBOL) 
})
keggresult$geneSymbols <- sapply(gene_symbols, function(symbols) paste(symbols, collapse = "/"))
keggresult$Description <- factor(keggresult$Description,levels= rev(keggresult$Description))
kegg1 <- keggresult[order(keggresult$Count, keggresult$pvalue, decreasing = T),]
ggplot(data = kegg1[1:10,], aes(x = Count, y = Description, fill = -log10(pvalue))) +
  scale_fill_continuous_c4a_seq('ag_sunset') +
  geom_bar(stat = 'identity', width = 0.8) + 
  labs(x = 'Number of Genes', y = '') + 
  theme_bw() + mytheme
ggplot(data = kegg1[1:10,], aes(x = Count, y = Description, fill = -log10(p.adjust))) +
  geom_bar(stat = 'identity', width = 0.8, alpha = 0.9) +  
  scale_fill_gradient(low = "#f3f9fd",high = "#0000cd")+
  labs(x = 'Number of Genes', y = '') +
  geom_text(aes(x = .5, 
                label= Description),
            hjust= 0) +
  theme_classic() + theme(axis.text.y = element_blank()) 
#Do they overlap with TSD-related genes
targetgenes <- read.table("TRP_Ca_HSP.txt",sep ="\t")
tarsub1 <- orthologs[orthologs$V2 %in% targetgenes$V1,]##pick
calsiumgenepsi <- read.table("calsium.genes.txt")
tarsub2 <- orthologs[orthologs$V1 %in% calsiumgenepsi$V1,]##pick
tarsub <- unique(rbind(tarsub1,tarsub2)) %>% mutate(
  type = ifelse(V2 %in% geneForGo, "In", "Out")
)

} # Tse S15 ovary down-regulation
} # Stage 15 things
{
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
#dndsTTs <- read.table(file = "turtle.dnds.free.txt", header = T)
#dndst2p <- read.table(file = "psi2tse.dnds.pair.txt") # pair
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
#Call DE genes in tse gonad
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


##All
orthologs <- read.table(file = "psi-red.orth-1to1.50")
tseDEgeneOrthSub <- tseDEgene %>% mutate(
  type = case_when(log2FoldChange >= 1 & padj <= 0.05 ~ "M", 
                   log2FoldChange <= -1 & padj <= 0.05 ~ "F",
                   TRUE ~ "E")
) %>% filter(rownames(.) %in% orthologs$V2)
#What about those TSE gonad - Psi brain genes?
tseDEgeneG2Bsub <- tseDEgeneOrthSub %>% filter(
  rownames(.) %in% c("CIRBP","CL4K","DNAJA4",
                     "KDM1A","PCGF6","ATP2B1",
                     "TXNDC11","STAT4")
)
#Genes to show labels
top20Genes <- head(tseDEgeneOrthSub, 20)
CaAndHSPgenes <- read.table("TRP_Ca_HSP.txt",sep ="\t")
psiGeneProduct <- read.table(file = "psi.geneID.finalinfo.txt")
psiGeneProduct <- psiGeneProduct %>% mutate(
  product = gsub("_[0-9]*", "", V6)
) %>% filter(V4 %in% orthologs$V1) %>%
  mutate(
    tseID = orthologs$V2[match(V4, orthologs$V1)]
  )
CaAndHSPgenes <- CaAndHSPgenes %>% filter(V1 %in% psiGeneProduct$product) # Calsium genes and HSP genes
GenesToshowSub1 <- psiGeneProduct[match(CaAndHSPgenes$V1,psiGeneProduct$product),] %>% filter(tseID %in% rownames(tseDEgeneOrthSub)) %>% select(tseID)
GenesToshowSub2 <- rownames(top20Genes %>% filter(rownames(.) %in% orthologs$V2))

tseDEgeneToshow <- tseDEgeneG2Bsub %>% filter(rownames(.) %in% c(GenesToshowSub1$tseID,GenesToshowSub2)) %>% 
  filter(type != "E") %>% mutate(
    ID = rownames(.),
    ID = case_when(ID == "LOC117884452" ~ "CYP19A1",
                   ID == "LOC117867547" ~ "TRPV2",
                   TRUE ~ ID),
    color = ifelse(ID %in% CaAndHSPgenes$V1, "TSD", "Other")
  )
library(ggrepel)
ggplot(tseDEgeneOrthSub, aes(x = log2FoldChange, y = -log10(padj), color = type)) + coord_cartesian(xlim = c(-6,6)) + 
  scale_color_manual(values = c("F" = "#CD0000", "M" = "#0000cd", "E" = "#bdbdbd", "TSD" = "#4daf4a"))+
  geom_point(alpha = .3) + geom_hline(yintercept = c(1.30103),linetype = 2) + theme_bw() +
  geom_label_repel(data = tseDEgeneToshow, aes(label = ID, colour = color)) 
#Call DE genes in psi gonad
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

#Expression profiles by Time
expTse <- read.table(file = "tse.tpm.mean.allTissue.txt", header = T,sep = ',',row.names = 1)
expPsi <- read.table(file = "psi.tpmmean.NewFilt.txt", header = T,sep = ',',row.names = 1)
expTseSub <- expTse %>% select(matches("FG|MG")) %>% filter(rownames(.) %in% orthologs$V2)
expPsiSub <- expPsi %>% select(matches("FG|MG")) %>% select(7, 2:5, 15, 10:13) %>% filter(rownames(.) %in% orthologs$V1)  # overlapped stages and orthologs
colnames(expTseSub) <- paste(colnames(expTseSub), "Tse",sep = "_")
colnames(expPsiSub) <- gsub("GMC", "G", colnames(expPsiSub)) 
colnames(expPsiSub) <- paste(colnames(expPsiSub), "Psi",sep = "_") # rename col. names
expTseSub$pid <- orthologs$V1[match(rownames(expTseSub), orthologs$V2)]
expTseSub$pid <- orthologs$V1[match(rownames(expTseSub), orthologs$V2)]
row.names(expTseSub) <- expTseSub$pid
expTseSub <- expTseSub[,-11]
mergedExp <-  merge(expTseSub, expPsiSub, by= 0, all=TRUE)
row.names(mergedExp) <- mergedExp$Row.names
mergedExp <- mergedExp[,-1] # merge exp by psi ID
library(scone)
ei <- as.matrix(mergedExp)
er2 <- as.data.frame(UQ_FN(ei))
#write.table(er2, file = "TseAndPsiexpUQ.txt")
ptUQexp <- read.table(file = "TseAndPsiexpUQ.txt")
geneToshow <- c("DNAJA4", "HSPB6", "TRPV2", "KDM6B", "JARID2") 
ptUQsub <- ptUQexp[rownames(ptUQexp) %in% psiGeneProduct$V4[psiGeneProduct$product %in% geneToshow],]
ptUQsub_longer <- ptUQsub %>% 
  rownames_to_column(var = "gene_id") %>% 
  tidyr::gather(key = "key", value = "value", -gene_id) %>%  
  mutate(
    stage = substr(key, 3, 4),
    sex = substr(key, 1, 1),
    species = substr(key, 6, 8),
    product = psiGeneProduct$product[match(gene_id, psiGeneProduct$V4)]
  ) 
ptUQsub_longerRatio <- ptUQsub_longer %>%
  group_by(product, stage, species) %>%
  summarise(
    F_value = value[sex == "F"],
    M_value = value[sex == "M"],
    Ratio = F_value / M_value,
    .groups = "drop"
  ) %>%
  select(product, stage, Ratio, species)
p1 <- ggplot(data = ptUQsub_longer, aes(x = stage, y = log2(value), group = factor(paste(sex, species)), color = sex, shape = species)) +
  geom_line(aes(linetype = species)) + geom_point() + facet_wrap( ~ product, scales = "free_y", ncol =1) + 
  scale_shape_manual(values = c("Psi"= 16, "Tse" = 21)) + 
  scale_color_manual(values = c("F" = "#cd0000", "M"= "#0000cd")) + theme_bw() + labs(x="", y="Log2 (expr.)")
p2 <- ggplot(data = ptUQsub_longerRatio, aes(x = stage, y = log2(Ratio), group = species, shape = species)) +
  geom_line(aes(linetype = species)) + geom_point() + facet_wrap( ~ product, scales = "free_y", ncol =1) + 
  scale_shape_manual(values = c("Psi"= 16, "Tse" = 21)) + theme_bw() + ylim(-3,3) + 
  geom_hline(yintercept = c(-1, 1)) + labs(x="", y="Log2 (F/M expr.)")
ggarrange(p1, p2)
#Expr. M/F ratio of all TSD genes
ptUQTSDsub <- ptUQexp[rownames(ptUQexp) %in% psiGeneProduct$V4[psiGeneProduct$product %in% targetgenes$V1],] %>%
  rownames_to_column(var = "gene_id") %>% 
  tidyr::gather(key = "key", value = "value", -gene_id) %>%  
  mutate(
    stage = substr(key, 3, 4),
    sex = substr(key, 1, 1),
    species = substr(key, 6, 8),
    product = psiGeneProduct$product[match(gene_id, psiGeneProduct$V4)]
  ) %>% group_by(product, stage, species) %>%
  summarise(
    F_value = value[sex == "F"]+1,
    M_value = value[sex == "M"]+1,
    Ratio = M_value / F_value,
    .groups = "drop"
  ) %>%
  select(product, stage, Ratio, species)
ggplot(data = ptUQTSDsub, aes(x = stage, y = log2(Ratio), group = species, shape = species)) +
  geom_line(aes(linetype = species)) + geom_point() + facet_wrap( ~ product, scales = "free_y") + 
  scale_shape_manual(values = c("Psi"= 16, "Tse" = 21)) + theme_bw() + ylim(-3,3) + 
  geom_hline(yintercept = c(-1, 1)) + labs(x="", y="Log2 (F/M expr.)")



#Do gonad
{
{
psiExp <- read.table("psi.tpmmean.NewFilt.txt", sep = ",", header = T, row.names = 1)
psiGonadGenes <- read.table("psi-wgcna-gonad-gene-after-rd2.txt") 
psiExpGonad <- psiExp[rownames(psiExp) %in% psiGonadGenes$ID, ] %>%
  select(grep("FG|MG", colnames(.))) %>% 
  select(grep("12|15|16|18|21", colnames(.))) %>% 
  select(5, 1:4, 10, 6:9)
psitpmF <- psiExpGonad %>%
  select(grep("FG", colnames(.))) 
psitpmM <- psiExpGonad %>%
  select(grep("MG", colnames(.))) 
} # Psi data
{
tseExp <- read.table("tse.tpm.mean.allTissue.txt", sep = ",", header = T, row.names = 1)
hubinfos <- read.table(file = "tse.wgcna.module.txt")
tseGonadGenes <- hubinfos %>% filter(tissue == "Gonad") %>% dplyr::select(gene)
tseExpGonad <- tseExp[rownames(tseExp) %in% tseGonadGenes$gene, ] %>%
  select(grep("FG|MG", colnames(.))) %>% 
  select(grep("12|15|16|18|21", colnames(.)))
#
tsetpmF <- tseExpGonad %>%
  select(grep("FG", colnames(.))) 
tsetpmM <- tseExpGonad %>%
  select(grep("MG", colnames(.))) 
} # Tse data
#
pca=PCA(t(tsetpmM),scale.unit=T,ncp=5,graph=F)
res=pca$ind$coord
res=res[,1:2]
dis=dist(res)
dis=as.matrix(dis)
#利用相邻样本点的欧式距离计算5个样本点原始的时间线
point1=0
point2=point1+dis[1,2]
point3=point2+dis[2,3]
point4=point3+dis[3,4]
point5=point4+dis[4,5]
raw.timeline=c(point1,point2,point3,point4,point5)
#通过scales包中的rescale函数将原始的时间线经过线性变换为0-10之间的新的时间线
new.timeline=rescale(raw.timeline,to=c(0,10))

##
#对tpm值进行zscale标准化
tpm.z=as.data.frame(t(apply(tsetpmM,1,scale))) %>% drop_na()
names(tpm.z)= names(tsetpmM)
#将5个时间点的zcale之后的表达值变成500个点的表达值
grid=data.frame(time=seq(0,10,length.out=500))
pseudotime.model.fun=function(value){time=new.timeline;data = tibble(value=value,time =time);model=loess(value ~ time,data);predict=grid %>% add_predictions(model);return(predict)}
res=apply(tpm.z,1,pseudotime.model.fun)
results=res %>% reduce(bind_cols)
exp500=as.data.frame(t(results[,seq(0,ncol(results),2)]))
names(exp500)=results$time...1
row.names(exp500)=row.names(tpm.z)

#对基因进行PCA分析，得到每个基因的Dim1和Dim2
gene.pca=PCA(exp500,scale.unit=T,ncp=5,graph=F)
res=gene.pca$ind$coord
res=res[,1:2]
exp500=cbind(exp500,res)
#用基因的Dim1和Dim2，进行排列组合，计算四种atan2值
exp500$atan2.1=atan2(exp500$Dim.1,exp500$Dim.2)
exp500$atan2.2=atan2(exp500$Dim.1,-exp500$Dim.2)
exp500$atan2.3=atan2(-exp500$Dim.1,exp500$Dim.2)
exp500$atan2.4=atan2(-exp500$Dim.1,-exp500$Dim.2)
#分别根据四种atan2值从低到高进行排序，然后用ComplexHeatmap可视化结果选定排序方式
order1=arrange(exp500,exp500$atan2.1)
order2=arrange(exp500,exp500$atan2.2)
order3=arrange(exp500,exp500$atan2.3)
order4=arrange(exp500,exp500$atan2.4)
col_fun <- colorRamp2(c(-2, 0, 2), c("#0000cd", "white", "#cd0000"))
p1=Heatmap(order1[,1:500],raster_quality = 1, cluster_rows = F,cluster_columns = F,show_row_names = F,show_column_names = F,column_title = "Order1",heatmap_legend_param = list(title="Order1",legend_height=unit(2,"cm")),col=col_fun)
p2=Heatmap(order2[,1:500],raster_quality = 1, cluster_rows = F,cluster_columns = F,show_row_names = F,show_column_names = F,column_title = "Order2",heatmap_legend_param = list(title="Order2",legend_height=unit(2,"cm")),col=col_fun)
p3=Heatmap(order3[,1:500],raster_quality = 1, cluster_rows = F,cluster_columns = F,show_row_names = F,show_column_names = F,column_title = "Order3",heatmap_legend_param = list(title="Order3",legend_height=unit(2,"cm")),col=col_fun)
p4=Heatmap(order4[,1:500],raster_quality = 1, cluster_rows = F,cluster_columns = F,show_row_names = F,show_column_names = F,column_title = "Order4",heatmap_legend_param = list(title="Order4",legend_height=unit(2,"cm")),col=col_fun)
p1+p2+p3+p4
##
datadf <- order1[,1:500]
htTse <- Heatmap(datadf,
          cluster_rows = F,
          cluster_columns = F,
          show_row_names = F,
          show_column_names = F,
          border = TRUE,
          use_raster = TRUE, raster_quality = 5,
          show_row_dend = TRUE,
          column_title = "Psi Ovary Exprssion",
          heatmap_legend_param = list(title="Z score"),col=col_fun)
ordergenes <- rownames(order1[,1:500])[row_order(htTse)]
#
breaks <- new.timeline
line_data <- data.frame(x = breaks, y = 0)
line_plot <- ggplot(line_data, aes(x = x, y = y)) +
  geom_line(color = "black", size = 1) + 
  geom_point(size = 3, color = "red") + 
  scale_x_continuous(breaks = breaks, labels = round(breaks, 2)) + 
  theme_minimal() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  labs(x = "Value", y = "")
# 将 ComplexHeatmap 转换为 ggplot 图形
heatmap_plot <- grid.grabExpr(draw(htTse))
combined_plot <- ggarrange(
  ggplotify::as.ggplot(heatmap_plot),
  line_plot,
  ncol = 1,
  heights = c(4, 1)
)
print(combined_plot)
#Original exp
tpmF <- tpmF[match(ordergenes, rownames(tpmF)), ]
zccale_df <- t(scale(t(tpmF)))
Heatmap(zccale_df,
        cluster_rows = F,
        cluster_columns = F,
        show_row_names = F,
        show_column_names = F,
        border = TRUE,
        use_raster = TRUE, raster_quality = 5,
        show_row_dend = TRUE,
        column_title = "Psi Ovary Exprssion",
        heatmap_legend_param = list(title="Z score"),col=col_fun)
#Gene Orders
psiFGorder <- ordergenes
psiMGorder <- ordergenes
tseFGorder <- ordergenes
tseMGorder <- ordergenes
turtletimecourseDF <- data.frame(
  Gene = c(psiFGorder, psiMGorder, tseFGorder, tseMGorder),
  type = c(
    rep("psiFG", length(psiFGorder)),
    rep("psiMG", length(psiMGorder)),
    rep("tseFG", length(tseFGorder)),
    rep("tseMG", length(tseMGorder))
  ),
  rank = c(
    seq_along(psiFGorder),
    seq_along(psiMGorder),
    seq_along(tseFGorder),
    seq_along(tseMGorder)
  )
)
#write.table(turtletimecourseDF, file = "turtletimecourse.txt")
} # Get all gonad genes expression order
genehub <- subset(dfmerge, psiType == "Gonad" & tseType == "Gonad") # select both gonad genes
timedf <- read.table("turtletimecourse.txt")
psiDEgene <- read.table(file = "psi.gonad.DEgene.txt")
psiDEgeneOrthSub <- psiDEgene %>% mutate(
  type = case_when(log2FoldChange >= 1 & padj <= 0.05 ~ "M", 
                   log2FoldChange <= -1 & padj <= 0.05 ~ "F",
                   TRUE ~ "E")
) %>% filter(rownames(.) %in% orthologs$V1)
tseDEgene <- read.table(file = "tse.gonad.DEgene.txt")
tseDEgeneOrthSub <- tseDEgene %>% mutate(
  type = case_when(log2FoldChange >= 1 & padj <= 0.05 ~ "M", 
                   log2FoldChange <= -1 & padj <= 0.05 ~ "F",
                   TRUE ~ "E")
) %>% filter(rownames(.) %in% orthologs$V2) %>%  mutate(psiID = orthologs$V1[match(rownames(.), orthologs$V2)])  
genehub <- genehub %>% mutate(
  DEpsi = psiDEgeneOrthSub$type[match(psiID, rownames(psiDEgeneOrthSub))],
  DEtse = tseDEgeneOrthSub$type[match(psiID, tseDEgeneOrthSub$psiID)]
) # Assign DEGs
gonadDEfreq <- genehub %>% group_by(DEpsi, DEtse) %>% 
  summarise(freq = n(), .groups = "drop") 
gonadBothSex <- timedf %>% mutate(
  species = gsub("FG|MG", "", type),
  tissue = gsub("psi|tse", "", type),
  psiID = case_when(
    species == "psi" ~ Gene,
    species == "tse" ~ orthologs$V1[match(Gene,orthologs$V2)]
  )) %>% filter(psiID %in% orthologs$V1) %>% mutate(
    product = gsub("_[0-9]*", "", geneidproduct$V6[match(psiID, geneidproduct$V4)])
  ) %>% filter(psiID %in% subset(genehub, DEpsi == "E" & DEtse == "E")$psiID)
timedfsubG <- gonadBothSex %>%
  filter(tissue %in% c("FG", "MG")) %>%              
  select(psiID, species, rank, tissue) %>%         
  pivot_wider(names_from = species,        
              values_from = rank) %>%     
  rename(psi_rank = psi, tse_rank = tse) %>% 
  na.omit() %>%
  mutate(
    psi_rank = psi_rank / max(psi_rank) * 100,
    tse_rank = tse_rank / max(tse_rank) * 100,
    tissue = ifelse(tissue == "FG", "Ovary", "Testis")
  ) %>%
  pivot_longer(cols = c(psi_rank, tse_rank), 
               names_to = "species", 
               values_to = "rank")
ggplot(timedfsubG, aes(x = species, y = rank)) +
  geom_boxplot(notch = T) + facet_wrap(~ tissue) +
  stat_compare_means(method = "wilcox.test", paired = TRUE,
                     method.args = list(alternative = "greater")) +
  theme_minimal()
wilcox.test(subset(timedfsubG, tissue == "Ovary")$rank, subset(timedfsubG, tissue == "Testis")$rank,
            paired = T, alternative = "greater")
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






} # Shut down of TSD pathways

{
# 1. Define gene modules --------------------------------------------------
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
# 2. Define DEGs ----------------------------------------------------------
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
dat = read.table("tse.all.tissue.count.txt", header = T, row.names = 1)
datBr <- dat %>% select(matches("F_B|M_B")) 
datMm <- dat %>% select(matches("F_M|M_M")) 
new_colnamesBr <- sapply(colnames(datBr), function(x){
  m <- str_match(x, "Tse_(\\d+)(F|M)_B(\\d)")[,2:4]
  stage <- m[1]
  sex <- m[2]
  rep <- m[3]
  paste0(sex, "B", stage, "_", rep)
})
new_colnamesMm <- sapply(colnames(datMm), function(x){
  m <- str_match(x, "Tse_(\\d+)(F|M)_M(\\d)")[,2:4]
  stage <- m[1]
  sex <- m[2]
  rep <- m[3]
  paste0(sex, "M", stage, "_", rep)
})
colnames(datBr) <- new_colnamesBr
colnames(datMm) <- new_colnamesMm
TSPcountBr <- datBr %>% select(matches("15|16|17|18|19|21"))
TSPcountMm <- datMm %>% select(matches("15|16|17|18|19|21"))
TSPcountBr <- TSPcountBr[, order(
  substr(colnames(TSPcountBr), 1, 1),                             
  as.numeric(sub("^[FM]B?(\\d+)_\\d+", "\\1", colnames(TSPcountBr))),  
  as.numeric(sub("^[FM]B?\\d+_(\\d+)", "\\1", colnames(TSPcountBr)))  
)] 
TSPcountMm <- TSPcountMm[, order(
  substr(colnames(TSPcountMm), 1, 1),                             
  as.numeric(sub("^[FM]M?(\\d+)_\\d+", "\\1", colnames(TSPcountMm))),  
  as.numeric(sub("^[FM]M?\\d+_(\\d+)", "\\1", colnames(TSPcountMm)))  
)]
library(DESeq2)
TSPcount <- TSPcountBr
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
#write.table(as.data.frame(resOrdered), file="tse.brain.DEgene.txt")
} # DEseq2 of tse 15-21
{
pmalecount <- read.table(file="psi.male.count",header = T,row.names = 1) 
pfemalecount <- read.table(file="psi.female.count",header = T,row.names = 1)
pcount <- merge(pmalecount, pfemalecount, by = 0) 
rownames(pcount) <- pcount$Row.names
pcount <- pcount[,-1]
pSDcount <- pcount %>% select(matches('FM|MM')) %>% select(matches("15|16|18|21|24"))
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
#write.table(as.data.frame(resOrdered), file="psi.meso.DEgene.txt")
} # DEseq2 of psi 15-24
# 3. Tse gonad DEGs volcano -----------------------------------------------
library(ggrastr)
tseDEgene <- read.table("tse.gonad.DEgene.txt", header=TRUE, row.names=1)
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
label_genes <- tseDEgene[tseDEgene$gene %in% c("KDM6B","JARID2"), ]
p + geom_text_repel(
  data = label_genes,
  aes(x = log2FoldChange, y = -log10(pvalue), label = gene, color = group),
  size = 3,
  box.padding = 0.4,
  point.padding = 0.3,
  max.overlaps = Inf,
  show.legend = FALSE
)
# 4. DEGs shift in turtles ------------------------------------------------
process_DEgene <- function(df){
  df$gene <- rownames(df)
  df$group <- "ns"
  df$group[df$padj <= 0.05 & df$log2FoldChange >= 1] <- "male"
  df$group[df$padj <= 0.05 & df$log2FoldChange <= -1] <- "female"
  return(df)
}
tseDEgeneGd <- read.table("tse.gonad.DEgene.txt", header=TRUE, row.names=1)
tseDEgeneBr <- read.table("tse.brain.DEgene.txt", header=TRUE, row.names=1)
tseDEgeneMm <- read.table("tse.meso.DEgene.txt", header=TRUE, row.names=1)
##Tse DEGs Venn
tseDEgeneGd <- process_DEgene(tseDEgeneGd)
tseDEgeneBr <- process_DEgene(tseDEgeneBr)
tseDEgeneMm <- process_DEgene(tseDEgeneMm)
Mm_list <- tseDEgeneMm %>% filter(group != "ns") %>% pull(gene)
Br_list <- tseDEgeneBr %>% filter(group != "ns") %>% pull(gene)
Gd_list <- tseDEgeneGd %>% filter(group != "ns") %>% pull(gene)
library(ggVennDiagram)
gene_sets <- list(
  Mm = Mm_list,
  Br = Br_list,
  Gd = Gd_list
)
ggVennDiagram(gene_sets, label_alpha = 0) +
  scale_fill_gradient(low = "white", high = "skyblue") +
  ggtitle("Non-ns Gene Overlap Across Tissues")
##Tse KDM6B no gonad specific sex-bias, no sertoli specificity 
tseExp <- read.table("tse.tpm.mean.allTissue.txt", header = T)
candigenelist <- c("KDM6B","DNAJA4","JARID2","HSPB6","RBM20","DDX4")
df_long <- tseExp %>% filter(Gene == "LRRFIP1") %>% 
  pivot_longer(
    cols = -Gene,
    names_to = "Sample",
    values_to = "TPM"
  ) %>%
  mutate(
    sex = substr(Sample, 1, 1),       # sex
    tissue = substr(Sample, 2, 2),    # tissue
    stage = as.numeric(str_extract(Sample, "\\d+")),  # stage
  )
ggplot(df_long, aes(x = as.character(stage), y = TPM, color = sex, group = interaction(sex, Gene))) +
  geom_point(alpha = 0.7) +
  geom_line(alpha = 0.7) +
  facet_wrap(~ tissue, scales = "free_y") +
  theme_bw(base_size = 14) +
  labs(
    x = "Stage",
    y = "TPM",
    color = "Sex"
  ) +
  theme(
    strip.background = element_rect(fill = "gray90", color = NA),
    panel.grid.minor = element_blank()
  )
##Psi DEGs
psiDEgeneGd <- read.table("psi.gonad.DEgene.txt", header=TRUE, row.names=1)
psiDEgeneBr <- read.table("psi.brain.DEgene.txt", header=TRUE, row.names=1)
psiDEgeneMm <- read.table("psi.meso.DEgene.txt", header=TRUE, row.names=1)
orthologs <- read.table("psi-red.orth-1to1")
tseDEgeneGd <- tseDEgeneGd %>% mutate(ID = orthologs$V1[match(gene, orthologs$V2)], tissue = "G") %>% select(ID, group, tissue)
tseDEgeneBr <- tseDEgeneBr %>% mutate(ID = orthologs$V1[match(gene, orthologs$V2)], tissue = "B") %>% select(ID, group, tissue)
tseDEgeneMm <- tseDEgeneMm %>% mutate(ID = orthologs$V1[match(gene, orthologs$V2)], tissue = "M") %>% select(ID, group, tissue)
psiDEgeneGd <- process_DEgene(psiDEgeneGd) %>% filter(gene %in% orthologs$V1) %>% mutate(ID = gene, tissue = "G") %>% select(ID, group, tissue)
psiDEgeneBr <- process_DEgene(psiDEgeneBr) %>% filter(gene %in% orthologs$V1) %>% mutate(ID = gene, tissue = "B") %>% select(ID, group, tissue)
psiDEgeneMm <- process_DEgene(psiDEgeneMm) %>% filter(gene %in% orthologs$V1) %>% mutate(ID = gene, tissue = "M") %>% select(ID, group, tissue)
tse_all <- bind_rows(
  tseDEgeneBr,
  tseDEgeneGd,
  tseDEgeneMm
) %>%
  rename(group = "groupT")
psi_all <- bind_rows(
  psiDEgeneBr,
  psiDEgeneGd,
  psiDEgeneMm
) %>%
  rename(group = "groupP")
merged_df <- tse_all %>%
  inner_join(psi_all, by = c("ID", "tissue")) %>%
  select(ID, groupT, groupP, tissue) %>% mutate(type = ifelse(groupT == groupP, "Conserve", "Shift"))

group_change_freq <- merged_df %>% filter(groupT != "ns") %>% 
  group_by(type, tissue) %>%
  summarise(freq = n()) %>%
  arrange(desc(freq))
group_change_freq <- group_change_freq %>%
  mutate(tissue_name = case_when(
    tissue == "G" ~ "Gonad",
    tissue == "B" ~ "Brain",
    tissue == "M" ~ "Muscle",
    TRUE ~ tissue
  ))
ggplot(group_change_freq, aes(x = change, y = freq, fill = tissue_name)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_bw() + facet_wrap(~ tissue)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Group Change", y = "Frequency", fill = "Tissue",
       title = "Group Change Frequency Across Tissues")


df_plot <- group_change_freq %>%
  group_by(tissue) %>%
  mutate(percentage = freq / sum(freq) * 100,
         ypos = cumsum(percentage) - 0.5*percentage)  # label位置

# 绘制donut图
ggplot(df_plot, aes(x = 2, y = percentage, fill = change)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar(theta = "y") +
  facet_wrap(~tissue) +  # 每个tissue一个donut
  xlim(0.5, 2.5) +       # 留空心形成donut
  geom_text(aes(y = ypos, label = paste0(round(percentage,1), "%")),
            color = "white", size = 3) +
  theme_void() +
  theme(legend.position = "right") +
  ggtitle("Group Change Proportion by Tissue")







#2. check gonad DEGs expr in somtatic tissues
tseExp <- read.table("tse.tpm.mean.allTissue.txt", sep = ",", header = T, row.names = 1)
candigenelist <- c("KDM6B","DNAJA4","JARID2","HSPB6","RBM20")
tseExpsub <- tseExp %>% filter(rownames(.) %in% candigenelist)



tarExp <- tseExpsub %>% select(grep("FM|MM", colnames(.)))
pheatmap(log2(tarExp),scale = "row",
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         order_color = NA)

{
# 5. 循环每个基因画热图
for(g in candigenelist){
  df <- tseExpsub[g, , drop = FALSE]   # 取单基因
  # 转置后更直观：列=组织-阶段
  df_t <- t(df)
  df_t <- as.data.frame(df_t)
  df_t$sample <- rownames(df_t)
  
  # 提取 tissue 和 stage 信息
  df_t$tissue <- substr(df_t$sample, 2, 2)   # G/B/M
  df_t$stage <- gsub("^[FM][GBM]", "", df_t$sample)
  
  # 按 tissue 和 stage 排序
  df_t <- df_t[order(df_t$tissue, as.numeric(df_t$stage)), ]
  
  # 重新整理成矩阵（行为 tissue-stage，列为表达）
  mat <- t(as.matrix(df_t[,1,drop=FALSE]))
  colnames(mat) <- df_t$sample
  
  # 绘图
  pheatmap(mat,
           scale = "row",
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           color = col_fun,
           main = g,
           border_color = NA)
}
}

tissue_list <- c("B","G","M")   # B:brain, G:gonad, M:meso
stages <- c("12","15","16","18","21")
ratio_list <- list()
for(tiss in tissue_list){
  for(stage in stages){
    F_col <- paste0("F", tiss, stage)
    M_col <- paste0("M", tiss, stage)
    
    if(F_col %in% colnames(tseExp) & M_col %in% colnames(tseExp)){
      ratio_list[[paste0(tiss, stage)]] <- data.frame(
        gene = rownames(tseExp),
        stage = stage,
        tissue = tiss,
        F_expr = tseExp[,F_col],
        M_expr = tseExp[,M_col],
        ratio = log2((tseExp[,F_col]+0.1)/(tseExp[,M_col]+0.1))   # log2(F/M)
      )
    }
  }
}
ratio_df <- bind_rows(ratio_list)
ratio_df <- ratio_df %>% filter(gene %in% candigenelist) %>%
  mutate(
    bias = case_when(
      ratio >= 1 ~ "female",
      ratio <= -1 ~ "male",
      TRUE ~ "non-bias"
    )
  ) 
#
ggplot(ratio_df, aes(x=stage, y=ratio, color = tissue)) +
  geom_point(alpha=0.7) + geom_line(aes(group = tissue))+ facet_wrap(~ gene) +
  scale_color_manual(values=c("G"="#ad9ad5","B"="#a7d49e","M"="#f7ab60")) + theme_classic()

ratio_long <- ratio_df %>% filter(tissue == "B")%>%
  select(gene, stage, tissue, F_expr, M_expr) %>%
  pivot_longer(cols = c(F_expr, M_expr),
               names_to = "sex", values_to = "TPM") %>%
  mutate(sex = ifelse(sex == "F_expr", "Female", "Male"),
         stage = as.numeric(stage))  # stage改为数值方便连线排序

ggplot(ratio_long, aes(x = as.character(stage), y = TPM, color = sex, group = sex)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  facet_wrap(~ gene + tissue, scales = "free_y") +   # 每个gene一个facet，也保留tissue信息
  scale_color_manual(values = c("Female" = "#ed7d61", "Male" = "#8cb7df")) +
  theme_bw(base_size = 14) +
  labs(x = "Stage", y = "TPM", color = "Sex") +
  theme(
    strip.text = element_text(face = "bold", size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank()
  )



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
p + geom_label_repel(
  data = label_genes,
  aes(x = log2FoldChange, y = -log10(padj), label = gene, color = group),
  size = 3,
  box.padding = 0.4,
  point.padding = 0.3,
  max.overlaps = Inf,
  show.legend = FALSE
) + geom_hline(yintercept = -log10(0.001399718)) + geom_vline(xintercept = 1.238515)
  
  
  
# 5. Gonad genes shift in turtles ----------------------------------------
# Lot of tse-gonad genes lost tissue specificity in psi
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
df_freq <- dfmerge %>%
  filter(!is.na(psiType) & !is.na(tseType)) %>%  filter(tseType == "Gonad" | tseType == "Brain") %>% 
  group_by(psiType, tseType) %>%  
  summarise(freq = n(), .groups = "drop") %>% # Count tissue specificity
  mutate(type = ifelse(psiType != tseType, "Shifted", "Conserved"))
summary_table <- df_freq %>% filter(psiType != "Not tissue-biased" & tseType != "Not tissue-biased") %>%
  group_by(tseType, type) %>%
  summarise(n = sum(freq), .groups = "drop") %>%
  rename(tseType = "tissue") %>%
  group_by(tissue) %>%
  mutate(percent = n / sum(n) * 100)
ggplot(summary_table, aes(x = tissue, y = percent, fill = type)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_minimal() +
  labs(
    x = "Tissue",
    y = "Percentage (%)",
    fill = "Type",
    title = ""
  ) +
  scale_fill_manual(values = c("Conserved" = "#bdbdbd", "Shifted" = "#737373")) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 100)) #Gonad - brain : X-squared = 600.52, df = 1, p-value < 2.2e-16
#TSD gonad lost genes show RNA splicing related patterns
dfmerge <- dfmerge %>% mutate(tseID = orthologs$V2[match(psiID, orthologs$V1)])
tseLossGonadGO <- read.table(file = "turtleTSElostGonad.GO.txt")
tseLossGonadGO <- tseLossGonadGO %>% filter(p.adjust <= 0.05)
tsdSplice <- tseLossGonadGO %>% 
  filter(grepl("splicing", Description, ignore.case = TRUE))
gene_list <- tsdSplice$geneID %>%
  strsplit(split = "/") %>%
  unlist() %>%
  trimws() %>%
  unique()
spliceDEGs <- tseDEgeneGd %>% filter(rownames(.) %in% gene_list)

tseGonadGenes <- tseGeneInfo %>% filter(tissue == "Gonad") %>% mutate(
  DE = tseDEgeneGd$group[match(gene, tseDEgeneGd$gene)]
) 


#Overlap with TSE pathway genes
CaAndHSPgenes <- read.table("TRP_Ca_HSP.txt",sep ="\t")
filtered_rows <- tseLossGonadGO[sapply(tseLossGonadGO$geneID, function(gene_ids) {
  any(sapply(CaAndHSPgenes$V1, function(gene) grepl(gene, gene_ids)))
}), ] %>% mutate(
  geneRatio = sapply(strsplit(as.character(GeneRatio), "/"), function(x) as.numeric(x[1]) / as.numeric(x[2])),
  Description = factor(Description, levels = Description[order(-p.adjust)])
)
ggplot(filtered_rows, aes(x = geneRatio, y = Description, size = Count, fill = -log10(p.adjust))) +
  geom_point(shape = 21) + scale_fill_gradient(low = "white", high = "#ED7D61") + theme_bw()
TSDdfSUB <- dfmerge %>% filter(Product %in% CaAndHSPgenes$V1) %>% mutate(
  func = CaAndHSPgenes$V2[match(Product, CaAndHSPgenes$V1)]
)
} # TSD shut off


{
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
} # DEseq2 of tse 12-21
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

} # DEseq2 of psi 12-21
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
} # Do Tse DEGs earlier than Psi
} # Call DEGs by stages




{
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
} #Expression of USDGs in turtles

  
{
#data used
expMir <- read.table(file="psi.sRNA.tpmmean.final.shortstack.txt")
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
candidateMir <- c(intersect(rownames(targetexpW), wtarlistZNRF3), "chrW_42977") #Note gga-miR-1687-3p is chrW_42977
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
#WGCNA check if FG12-14 significant micros?



#Who target KDM6B
wtarlistKDM6B <- wtarlist %>%
  filter(V5 == "KDM6B") %>%
  select(V2) %>%
  mutate(chr = str_extract(V2, "^[^_]+_[^_]+")) %>%
  pull(chr) %>%
  unique()
candidateMir <- c(intersect(rownames(targetexpW), wtarlistKDM6B), "chrW_42977") #Note gga-miR-1687-3p is chrW_42977
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
#Do chrW_43028 target most ZW-gametolog Z copy?
chrW43028 <- wtarlist %>%
  filter(V2 == "chrW_43028_gga-miR-130b-3p") %>%
  select(V3, V4, V5) %>%
  unique() 



#find zw game - > test 8mer% multi %



  } # micro RNA part


} # USDG in GSD turtles


{
psiAR <- read.table(file = "psi.ar.filt.txt")  
psiACf <- read.table(file = "psi.27ac.f.Peak.txt") 
psiACm <- read.table(file = "psi.27ac.m.Peak.txt")  
psiMe3f <- read.table(file = "psi.27me3.f.Peak.txt")  
psiMe3m <- read.table(file = "psi.27me3.m.Peak.txt")  
geneloci <- read.table(file = "psi.geneID.finalinfo.txt")
exonloci <- read.table(file = "psi.exon.region")
#Turn genes into GR ranges
library(GenomicRanges)
grAR <- GRanges(seqnames = gsub("Psi.","",psiAR$V1),
               ranges = IRanges(start = psiAR$V2, end = psiAR$V3))
grGene <- GRanges(seqnames = geneloci$V1,
               ranges = IRanges(start = geneloci$V2, end = geneloci$V3))
grExon <- GRanges(seqnames = exonloci$V1,
                  ranges = IRanges(start = exonloci$V2, end = exonloci$V3))
grACf <- GRanges(seqnames = psiACf$V1,
                 ranges = IRanges(start = psiACf$V2, end = psiACf$V3))
grACm <- GRanges(seqnames = psiACm$V1,
                 ranges = IRanges(start = psiACm$V2, end = psiACm$V3))
#More than 70% ARs falls into intergenenic regions(IGR) / introns
ARexon_hits <- findOverlapPairs(grAR, grExon)
intersect_ranges <- pintersect(ARexon_hits)
overlap_lengths <- width(intersect_ranges)
dfARexon <- data.frame(AR = ARexon_hits@first, Exon = ARexon_hits@second, ovlp = overlap_lengths) 
dfARexonFilt.8 <- unique(dfARexon %>% filter((ovlp/AR.width) >= .8) %>% select(1:3))
psiAR <- psiAR %>% mutate(
  type = ifelse(V2 %in% dfARexonFilt.8$AR.start, "Exon", "InterGenetic/Intron")
)
#Density of AR lengths
ggplot(psiAR, aes(x = log2(V3-V2),fill = type)) +
  geom_density(alpha = .8) + theme_bw() +labs(x= "Log2(Length)", y = "") + 
  ylim(0, .6) + scale_fill_manual(values = c("Exon" = "#4daf4a", "InterGenetic/Intron" = "#bdbdbd"))
type_counts <- psiAR %>%
  count(type) %>%
  mutate(percentage = n / sum(n) * 100,  
         label = paste0(type, " (", round(percentage, 1), "%)")) 
ggplot(type_counts, aes(x = "", y = n, fill = type)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") + 
  labs(title = "", x = NULL, y = NULL) +
  geom_text(aes(label = label), position = position_stack(vjust = 0.5)) +  
  theme_void() + scale_fill_manual(values = c("Exon" = "#4daf4a", "InterGenetic/Intron" = "#bdbdbd"))
#Find Active ARs
overlaps1 <- findOverlaps(grAR, grACf) 
fActiveARs <- unique(psiAR[overlaps1@from,])
overlaps2 <- findOverlaps(grAR, grACm) 
mActiveARs <- unique(psiAR[overlaps2@from,])
funiqARs <- psiAR %>% filter(V4 %in% setdiff(fActiveARs$V4, mActiveARs$V4))
muniqARs <- psiAR %>% filter(V4 %in% setdiff(mActiveARs$V4, fActiveARs$V4))
#write.table(funiqARs, file = "psi.f.ar.txt")
#write.table(muniqARs, file = "psi.m.ar.txt")
#Find active AR overlap genes
grARf <- GRanges(seqnames = gsub("Psi.","",funiqARs$V1),
                 ranges = IRanges(start = funiqARs$V2, end = funiqARs$V3))
grARm <- GRanges(seqnames = gsub("Psi.","",muniqARs$V1),
                 ranges = IRanges(start = muniqARs$V2, end = muniqARs$V3))
genesWithARf <- findOverlaps(grARf, grGene) 
genesWithARf <- unique(geneloci[genesWithARf@to,])
genesWithARm <- findOverlaps(grARm, grGene) 
genesWithARm <- unique(geneloci[genesWithARm@to,])
#TFBS enrichment
tfbs_all <- read_tsv("psi.ar.all.ame.tsv")
tfbs_f <- read_tsv("psi.f.ar.ame.tsv")
tfbs_m <- read_tsv("psi.m.ar.ame.tsv")
tfbs_filt <- tfbs_all %>%
  filter(`adj_p-value` <= 5e-2) %>%
  mutate(
    geneID = if_else(
      is.na(motif_alt_ID),
      str_extract(motif_ID, "^[^_]+"),
      str_extract(motif_alt_ID, "^[^_]+")
    ),
    geneID = toupper(geneID)
  )
targetTF <- tfbs_filt %>%
  mutate(geneID_split = str_split(geneID, "::")) %>%
  unnest(geneID_split) %>%
  pull(geneID_split) %>%
  unique()
psiTFsub <- geneloci %>% mutate(
  geneID = gsub("_[0-9]*", "", V6)
) %>% filter(
  geneID %in% targetTF
) 
#Genes in SD list
sdlist <- read.csv("newSDRelaetedGenes.csv", row.names = 1)
psiTFsd <- psiTFsub %>% filter(
  geneID %in% sdlist$Gene
) 
#GO enrichment
library(clusterProfiler)
library(org.Hs.eg.db)
library(cols4all)
geneForGo <- psiTFsub$geneID
targetGenes <- bitr(geneForGo,
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
#Go
GOresult <- ego@result %>% head(20) %>%
  mutate(
    GeneRatio_numeric = as.numeric(sapply(strsplit(as.character(GeneRatio), "/"), function(x) as.numeric(x[1])/as.numeric(x[2]))),
  )
ggplot(data = GOresult, aes(x = GeneRatio_numeric, y = factor(Description, levels = rev(GOresult$Description)), fill = -log10(p.adjust))) +
  geom_bar(stat = 'identity', width = 0.8, alpha = 0.9) +  
  scale_fill_gradient(low = "#fff6f1",high = "#cd0000")+
  labs(x = 'Gene Ratio', y = '') +
  geom_text(aes(x = 0.005, 
                label= Description),
            hjust= 0) +
  theme_classic() + theme(axis.text.y = element_blank()) 
#Kegg
ekegg <- enrichKEGG(gene = targetGenes$ENTREZID, 
                    organism = 'hsa',
                    pvalueCutoff = 0.05,
                    pAdjustMethod = "fdr")
keggresult <- ekegg@result
keggresult$FoldEnrichment <- apply(keggresult,1,function(x){
  GeneRatio=eval(parse(text=x["GeneRatio"]))
  BgRatio=eval(parse(text=x["BgRatio"]))
  enrichment_fold=round(GeneRatio/BgRatio,2)
  enrichment_fold
})
keggresult <- mutate(keggresult, RichFactor= Count / as.numeric(sub("/\\d+", "", BgRatio)))
geneID_list <- strsplit(keggresult$geneID, "/")
gene_symbols <- lapply(geneID_list, function(ids) {
  mapped_genes <- bitr(ids, fromType = "ENTREZID", toType = "SYMBOL", OrgDb = org.Hs.eg.db)
  return(mapped_genes$SYMBOL) 
})
keggresult$geneSymbols <- sapply(gene_symbols, function(symbols) paste(symbols, collapse = "/"))
keggresult$Description <- factor(keggresult$Description,levels= rev(keggresult$Description))
kegg1 <- keggresult[order(keggresult$Count, keggresult$pvalue, decreasing = T),]
ggplot(data = kegg1[1:20,], aes(x = Count, y = Description, fill = -log10(p.adjust))) +
  geom_bar(stat = 'identity', width = 0.8, alpha = 0.9) +  
  scale_fill_gradient(low = "#f3f9fd",high = "#cd0000")+
  labs(x = 'Number of Genes', y = '') +
  geom_text(aes(x = .5, 
                label= Description),
            hjust= 0) +
  theme_classic() + theme(axis.text.y = element_blank()) 


{
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
  } # Data of gene tissue specificity


} # AR in Psi



{
psiExp <- read.table(file = "psi.tpmmean.NewFilt.txt", sep=',', row.names=1, header=T) 
tseExp <- read.table("tse.tpm.mean.allTissue.txt", header = T, row.names = 1)
rownames(tseExp) <- tseExp$Gene
tseExp <- tseExp[, -1]
geneidproduct <- read.table("psi.geneID.finalinfo.txt")
orthologs <- read.table("psi-red.orth-1to1.50")
#Single species no normalization
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
tarGene <- psiOnlyExpr_longer %>% filter(gene == "DHB1_ENSCPBP00000022668-D1")
tarGene <- psiOnlyExpr_longer %>% filter(gene %in% c(
  "DHB1_ENSCPBP00000022668-D1",
  "ANDR_ENSCPBP00000009353-D1",
  "ESR1_ENSGEVP00005018606-D1",
  "PRGR_ENSGEVP00005011508-D1",
  "GCR_ENSGEVP00005016354-D2",
  "FOXO1_ENSGEVP00005004416-D1"))
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
subset1 <- tarGene_ratio %>% filter(stage %in% c("12","15","16","18","21")) %>% filter(tissue =="G") %>% mutate(species = "Psi")
#Single species no normalization Tse
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
tarGene <- tseExpr_longer %>% filter(gene == "DHB1_ENSCPBP00000022668-D1")
tarGene <- psiOnlyExpr_longer %>% filter(gene %in% c(
  "DHB1_ENSCPBP00000022668-D1",
  "ANDR_ENSCPBP00000009353-D1",
  "ESR1_ENSGEVP00005018606-D1",
  "PRGR_ENSGEVP00005011508-D1",
  "GCR_ENSGEVP00005016354-D2",
  "FOXO1_ENSGEVP00005004416-D1"))
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
#Plot candidate hormone genes
ggplot(rbind(subset1,subset2), aes(x = as.character(stage), y = log2(tpmR), group = factor(paste(gene,species)), color = species)) +
  geom_line(size = 1)  + geom_point(size =2) +geom_hline(yintercept = 1)+
  facet_wrap( ~ gene, scales="free_y")  + scale_color_manual(values = c("Psi" = "#ed7d61", "Tse" = "#8cb7df")) +
  theme_classic() + labs(title = "",
                         x="",y="Log2(F/M TPM)") 

#Hub genes expr. merge
targetGenes <- c(
  "KDM6B_ENSPMRP00000034777-D1",
  "JARD2_ENSCPBP00000031845-D1",
  "KCNK3_ENSCPBP00000024946-D3",
  "GCR_ENSGEVP00005016354-D2",
  "DMRT1_ENSCPBP00000036668-D1",
  "FOXL2_ENSCPBP00000027202-D1",
  "CRYAB_ENSCPBP00000011979-D1",
  "WT1A_ENSGALP00000019761-D1",
  "CTNB1_ENSGALP00000045406-D1"
)
targetGenes <- c(
  "DHB1_ENSCPBP00000022668-D1",
  "ANDR_ENSCPBP00000009353-D1",
  "ESR1_ENSGEVP00005018606-D1",
  "PRGR_ENSGEVP00005011508-D1",
  "GCR_ENSGEVP00005016354-D2")
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
} # PGR-AR-ESR genes
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
targetGenes <- c("ESR1_ENSGEVP00005018606-D1",
                 "PRGR_ENSGEVP00005011508-D1",
                 "GREB1_ENSCPBP00000019364-D1",
                 "NRIP1_ENSCPBP00000041311-D1") # esr1 target genes
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
#Exp. normalization between species by UQ
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
} # Gene expr. check


{
library(pheatmap)
library(dplyr)
exp <- read.table(file="psi.sRNA.tpmmean.final.shortstack.txt")
#exp <- read.table(file="C:/Users/13719/Desktop/r_file/psi.sRNA.tpmmean.final.uniq.shortstack.txt")
targetexp <- exp[grep("chrZ_43837|chrZ_43932",rownames(exp)),]
targetexp <- targetexp[,grep("G",colnames(exp))]
targetexp <- targetexp %>% filter_all(.,any_vars(.>0))
library(pheatmap)
pheatmap(targetexp,cluster_rows = T,cluster_cols = F,scale="row",
         color = c(colorRampPalette(c("#0000cd", "white"))(50),
                   colorRampPalette(c("white", "#cd0000"))(50)),
         border_color=NA)
} # Micro expr. check

{

{
library(data.table)
library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)  
library(enrichplot)
library(ggplot2)

psiF15 <- fread("psi_female.gene_spectra_score.k_15.dt_0_5.txt")
tseF19 <- fread("tse_female.gene_spectra_score.k_19.dt_0_5.txt")
psiM13 <- fread("psi_male.gene_spectra_score.k_13.dt_0_5.txt")
tseM16 <- fread("tse_male.gene_spectra_score.k_16.dt_0_5.txt")

psiMF18 <- fread("02.psi.gene_spectra_score.k_18.dt_0_5.txt")
tseMF16 <- fread("02.tse.gene_spectra_score.k_16.dt_0_5.txt")
  
orthologs <- read.table("psi-red.orth-1to1.50")
orthologs$V1 <- gsub("_","-",orthologs$V1)

module_overlap_pipeline <- function(psi_df, tse_df, orthologs, top_n=100){
  
  convert_matrix <- function(dt){
    
    modules <- dt[[1]]
    mat <- as.matrix(dt[,-1])
    
    mat <- t(mat)
    
    colnames(mat) <- modules
    rownames(mat) <- colnames(dt)[-1]
    
    mat
  }
  
  psi_mat <- convert_matrix(psi_df)
  tse_mat <- convert_matrix(tse_df)
  
  # 限制到1:1 ortholog
  orth_df <- orthologs[
    orthologs$V1 %in% rownames(psi_mat) &
      orthologs$V2 %in% rownames(tse_mat),
  ]
  
  psi_mat <- psi_mat[orth_df$V1,]
  tse_mat <- tse_mat[orth_df$V2,]
  
  # 统一 gene ID（TSE）
  rownames(psi_mat) <- orth_df$V2
  rownames(tse_mat) <- orth_df$V2
  
  get_top <- function(mat,n){
    
    lapply(colnames(mat),function(m){
      
      rownames(mat)[order(mat[,m],decreasing=TRUE)][1:n]
      
    }) |> setNames(colnames(mat))
    
  }
  
  psi_top <- get_top(psi_mat,top_n)
  tse_top <- get_top(tse_mat,top_n)
  
  background <- nrow(psi_mat)
  
  psi_modules <- names(psi_top)
  tse_modules <- names(tse_top)
  
  overlap <- matrix(0,length(psi_modules),length(tse_modules))
  pval <- matrix(1,length(psi_modules),length(tse_modules))
  jaccard <- matrix(0,length(psi_modules),length(tse_modules))
  
  rownames(overlap) <- psi_modules
  colnames(overlap) <- tse_modules
  
  rownames(pval) <- psi_modules
  colnames(pval) <- tse_modules
  
  rownames(jaccard) <- psi_modules
  colnames(jaccard) <- tse_modules
  
  for(i in psi_modules){
    
    g1 <- psi_top[[i]]
    
    for(j in tse_modules){
      
      g2 <- tse_top[[j]]
      
      inter <- intersect(g1,g2)
      union <- union(g1,g2)
      
      k <- length(inter)
      
      overlap[i,j] <- k
      jaccard[i,j] <- k/length(union)
      
      pval[i,j] <- phyper(
        k-1,
        length(g2),
        background-length(g2),
        length(g1),
        lower.tail=FALSE
      )
      
    }
    
  }
  
  padj <- matrix(p.adjust(pval,"BH"),nrow=nrow(pval))
  
  rownames(padj) <- rownames(pval)
  colnames(padj) <- colnames(pval)
  
  list(
    overlap=overlap,
    pvalue=pval,
    padj=padj,
    jaccard=jaccard,
    psi_top=psi_top,
    tse_top=tse_top,
    psi_mat=psi_mat,
    tse_mat=tse_mat
  )
  
}
res <- module_overlap_pipeline(
  psiMF18,
  tseMF16,
  orthologs,
  top_n=50
)

stars <- ifelse(res$padj < 0.05,"*","")
pheatmap(
  res$jaccard,
  fontsize_number=10,
  cluster_rows=TRUE,
  cluster_cols=TRUE,
  color=colorRampPalette(c("white","#ED7D61"))(100),
  display_numbers=stars,
  number_color="black",
  border_color="grey"
)
#gene cell specificity
library(data.table)
psiCS <- fread("0005.sc-all-cell.psi")
tseSC <- fread("0005.sc-all-cell.tse")
psiCS <- psiCS[, .(
  psi_gene = V1,
  psi_cell = V4,
  psi_sex  = V6
)]
tseSC <- tseSC[, .(
  tse_gene = V1,
  tse_cell = V4,
  tse_sex  = V6
)]
setnames(orthologs, c("psi_gene","tse_gene"))
orthoCS <- merge(
  orthologs,
  psiCS,
  by = "psi_gene",
  all.x = TRUE
)
orthoCS <- merge(
  orthoCS,
  tseSC,
  by = "tse_gene",
  all.x = TRUE
)
#sup GEP module top genes
library(UpSetR)
psi_list <- list(
  psi_m2  = res$psi_top[["2"]],
  psi_m11 = res$psi_top[["11"]],
  psi_m8  = res$psi_top[["8"]],
  psi_m5  = res$psi_top[["5"]]
) # psi modules
tse_list <- list(
  tse_m10 = res$tse_top[["10"]],
  tse_m2  = res$tse_top[["2"]],
  tse_m6  = res$tse_top[["6"]]
) # tse modules
psi_input <- fromList(psi_list)
upset(
  psi_input,
  nsets = 4,
  order.by = "freq"
)
tse_input <- fromList(tse_list)
upset(
  tse_input,
  nsets = 3,
  order.by = "freq"
)
psi_genes <- unique(unlist(psi_list))
psi_df <- data.frame(
  gene = psi_genes,
  psi_m2  = psi_genes %in% psi_list$psi_m2,
  psi_m11 = psi_genes %in% psi_list$psi_m11,
  psi_m8  = psi_genes %in% psi_list$psi_m8,
  psi_m5  = psi_genes %in% psi_list$psi_m5
)
rownames(psi_df) <- psi_df$gene
psi_df$gene <- NULL
tse_genes <- unique(unlist(tse_list))
tse_df <- data.frame(
  gene = tse_genes,
  tse_m10 = tse_genes %in% tse_list$tse_m10,
  tse_m2  = tse_genes %in% tse_list$tse_m2,
  tse_m6  = tse_genes %in% tse_list$tse_m6
)
rownames(tse_df) <- tse_df$gene
tse_df$gene <- NULL
#The interaction of tse-3 to psi-8
tse3 <- res$tse_top[["3"]]
psi8 <- res$psi_top[["8"]]
intersect(tse3,psi8)

#GO enrich for 3 sets
geneIdproduct <- fread("psi.geneID.finalinfo.txt",header = F)
geneIdproduct$V4 <- gsub("_","-",geneIdproduct$V4)
psiFsupGEP <-  gsub("_[0-9]*","", geneIdproduct$V6[match(orthologs %>% filter(tse_gene %in% psi_list$psi_m2) %>% pull(psi_gene), geneIdproduct$V4)])
psiMsupGEP <-  gsub("_[0-9]*","", geneIdproduct$V6[match(orthologs %>% filter(tse_gene %in% c(psi_list$psi_m11,psi_list$psi_m8)) %>% pull(psi_gene), geneIdproduct$V4)])
psiFMsupGEP <- gsub("_[0-9]*","", geneIdproduct$V6[match(orthologs %>% filter(tse_gene %in% psi_list$psi_m5) %>% pull(psi_gene), geneIdproduct$V4)])
psiMearlyGEP <- gsub("_[0-9]*","", geneIdproduct$V6[match(orthologs %>% filter(tse_gene %in% psi_list$psi_m8) %>% pull(psi_gene), geneIdproduct$V4)]) 
tseFMsupLateGEP <- gsub("_[0-9]*","", geneIdproduct$V6[match(orthologs %>% filter(tse_gene %in% tse_list$tse_m10) %>% pull(psi_gene), geneIdproduct$V4)])
tseFMsupEarlyGEP <- gsub("_[0-9]*","", geneIdproduct$V6[match(orthologs %>% filter(tse_gene %in% tse_list$tse_m2) %>% pull(psi_gene), geneIdproduct$V4)])
tseFsupGEP <- gsub("_[0-9]*","", geneIdproduct$V6[match(orthologs %>% filter(tse_gene %in% tse_list$tse_m6) %>% pull(psi_gene), geneIdproduct$V4)])
tsePreSteroGEP <- gsub("_[0-9]*","", geneIdproduct$V6[match(orthologs %>% filter(tse_gene %in% res$tse_top[["3"]]) %>% pull(psi_gene), geneIdproduct$V4)])

runGO <- function(genes, label){
  genes <- unique(genes)
  genes <- genes[!is.na(genes)]
  eg <- bitr(genes,
             fromType = "SYMBOL",
             toType = "ENTREZID",
             OrgDb = org.Hs.eg.db)
  ego <- enrichGO(eg$ENTREZID,
                  OrgDb = org.Hs.eg.db,
                  keyType = "ENTREZID",
                  ont = "ALL",
                  pvalueCutoff = 0.05,
                  pAdjustMethod = "fdr",
                  readable = TRUE)
  ego <- simplify(ego, cutoff=0.9, by="p.adjust", select_fun=min)
  res <- ego@result %>%
    filter(qvalue < 0.05) %>%
    mutate(
      GeneRatio_numeric = as.numeric(sapply(strsplit(GeneRatio, "/"),
                                            function(x) as.numeric(x[1]) / as.numeric(x[2]))),
      set = label
    )
  return(res)
} # GO function
runKEGG <- function(genes, label){
  genes <- unique(genes)
  genes <- genes[!is.na(genes)]
  eg <- bitr(genes,
             fromType = "SYMBOL",
             toType = "ENTREZID",
             OrgDb = org.Hs.eg.db)
  ekegg <- enrichKEGG(
    gene = eg$ENTREZID,
    organism = "hsa",
    pvalueCutoff = 0.05,
    pAdjustMethod = "fdr"
  )
  res <- ekegg@result
  if(nrow(res) == 0) return(NULL)
  res <- res %>%
    mutate(
      GeneRatio_num = as.numeric(sapply(strsplit(GeneRatio, "/"),
                                        function(x) as.numeric(x[1]) / as.numeric(x[2]))),
      BgRatio_num = as.numeric(sapply(strsplit(BgRatio, "/"),
                                      function(x) as.numeric(x[1]) / as.numeric(x[2]))),
      RichFactor = Count / as.numeric(sub("/\\d+", "", BgRatio)),
      FoldEnrichment = GeneRatio_num / BgRatio_num,
      set = label
    ) 
  return(res)
} # kegg function
res_list <- list(
  psiF = runGO(psiFsupGEP, "psi_F"),
  psiM = runGO(psiMsupGEP, "psi_M"),
  psiFM = runGO(psiFMsupGEP, "psi_FM"),
  tseLate = runGO(tseFMsupLateGEP, "tse_late"),
  tseEarly = runGO(tseFMsupEarlyGEP, "tse_early")
)
shift_list <- list(
  psiEM <- runGO(psiMearlyGEP, "psi_EM"),
  tseSt <- runGO(tsePreSteroGEP, "tse_St")
)
GO_all <- bind_rows(shift_list)
top_terms <- GO_all  %>%
  group_by(Description) %>%
  summarise(min_q = min(qvalue)) %>%
  arrange(min_q) %>%
  slice(1:20) %>%
  pull(Description)
GO_plot <- GO_all %>%
  filter(Description %in% top_terms)
ggplot(GO_plot,
       aes(x = set,
           y = Description)) +
  
  geom_point(aes(size = GeneRatio_numeric,
                 color = -log10(qvalue))) +
  
  facet_grid(ONTOLOGY ~ ., scales = "free_y", space = "free") +
  
  scale_color_gradient(low = "#8CB7DF", high = "#ED7D61") +
  
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 8),
    strip.text = element_text(face = "bold")
  ) +
  
  labs(
    x = NULL,
    y = NULL,
    color = "-log10(qvalue)",
    size = "GeneRatio"
  )
kegg_list <- list(
  psiF = runKEGG(psiFsupGEP, "psi_F"),
  psiM = runKEGG(psiMsupGEP, "psi_M"),
  psiFM = runKEGG(psiFMsupGEP, "psi_FM"),
  tseLate = runKEGG(tseFMsupLateGEP, "tse_late"),
  tseEarly = runKEGG(tseFMsupEarlyGEP, "tse_early")
)
KEGG_all <- bind_rows(kegg_list)
top_kegg <- KEGG_all %>%
  group_by(Description) %>%
  summarise(min_q = min(p.adjust)) %>%
  arrange(min_q) %>%
  slice(1:20) %>%
  pull(Description)
KEGG_plot <- KEGG_all %>%
  filter(Description %in% top_kegg)
ggplot(KEGG_plot,
       aes(x = set,
           y = Description)) +
  
  geom_point(aes(size = RichFactor,
                 color = -log10(p.adjust))) +
  
  scale_color_gradient(low = "#8CB7DF", high = "#ED7D61") +
  
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 8)
  ) +
  
  labs(
    x = NULL,
    y = NULL,
    color = "-log10(FDR)",
    size = "RichFactor"
  )
#Binding by which TFs infer by ChEA3
psiftf <- fread("psiFTF_meanRank.tsv")
psimtf <- fread("psiMTF_meanRank.tsv")
psifmtf <- fread("psiFMTF_meanRank.tsv")
tseetf <- fread("tseEarly_meanRank.tsv")
tseltf <- fread("tseLate_meanRank.tsv")

geneIdproduct <- fread("psi.name2product.txt",header = F)
geneIdproduct <- geneIdproduct %>% mutate(
  psiID = gsub("_","-",geneIdproduct$V1),
  prod = gsub("_[0-9]*","", geneIdproduct$V2)
)

processTF <- function(tf_df, label, geneIdproduct){
  tf_df <- tf_df %>%
    mutate(
      psiID = geneIdproduct$psiID[match(TF, geneIdproduct$prod)],
      set = label
    ) %>%
    slice(1:20)
  return(tf_df)
} 

psiftf2  <- processTF(psiftf,  "psi_F", geneIdproduct)
psimtf2  <- processTF(psimtf,  "psi_M", geneIdproduct)
psifmtf2 <- processTF(psifmtf, "psi_FM", geneIdproduct)
tseetf2  <- processTF(tseetf,  "tse_early", geneIdproduct)
tseltf2  <- processTF(tseltf,  "tse_late", geneIdproduct)

tf_all <- bind_rows(psiftf2, psimtf2, psifmtf2, tseetf2, tseltf2)

orthoCSsup <- as.data.frame(orthoCS)
orthoCSsup <- orthoCSsup %>%
  mutate(
    psi_cell = ifelse(str_detect(psi_cell, "^[0-9.]+$"), "None", psi_cell),
    tse_cell = ifelse(str_detect(tse_cell, "^[0-9.]+$"), "None", tse_cell)
  ) %>%
  mutate(
    psi_cell = ifelse(str_detect(psi_sex, "None"), "None", psi_cell),
    tse_cell = ifelse(str_detect(tse_sex, "None"), "None", tse_cell)
  )

plot_df <- tf_all %>%
  left_join(orthoCSsup, by = c("psiID" = "psi_gene"))

plot_long <- plot_df %>%
  dplyr::select(TF, set, Score, psi_cell, tse_cell) %>%
  tidyr::pivot_longer(
    cols = c(psi_cell, tse_cell),
    names_to = "species",
    values_to = "cell_type"
  ) %>%
  mutate(
    species = recode(species,
                     psi_cell = "PSI",
                     tse_cell = "TSE"),
    score_plot = -log10(Score + 1e-10)
  )

ggplot(plot_long,
       aes(x = set,
           y = TF)) +
  
  geom_point(aes(
    size = score_plot,
    fill = cell_type
  ),
  shape = 21,
  color = "black") +
  
  facet_wrap(~species, nrow = 1)  +
  
  theme_bw() +
  
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 12, face = "bold")
  ) +
  
  labs(
    x = NULL,
    y = "TF",
    fill = "Cell specificity",
    size = "-log10(score)"
  )

# links
links <- orthoCSsup %>% filter(tse_gene %in% psi_list$psi_m11) %>%
  group_by(psi_cell, tse_cell) %>%
  summarise(value = n(), .groups="drop") %>%
  mutate(
    flow_type = ifelse(psi_cell == tse_cell, "conserved", "changed")
  )
library(ggalluvial)
ggplot(
  links,
  aes(axis1 = psi_cell,
      axis2 = tse_cell,
      y = value)
) +
  geom_alluvium(
    aes(fill = flow_type),
    width = 0.2,
    alpha = 0.8
  ) +
  geom_stratum(
    aes(fill = after_stat(stratum)),
    width = 0.25,
    color = "black"
  ) +
  
  geom_text(
    stat = "stratum",
    aes(label = after_stat(stratum)),
    size = 4
  ) +
  scale_x_discrete(
    limits = c("PSI", "TSE"),
    expand = c(.1, .1)
  ) +
  scale_fill_manual(
    values = c(
      "conserved" = "grey70",
      "changed" = "#ED7D61"
    )
  ) +
  labs(x = NULL, y = "Gene count") +
  theme_classic() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks = element_blank()
  )







orthoCSsup <- orthoCS %>% filter(psi_gene %in% psimtf$psiID)








  

psim2Genes <- orthoCSsup %>% dplyr::select(psi_gene, psi_cell, tse_cell) %>%
  mutate(prod = gsub("_[0-9]*","", geneIdproduct$V2[match(psi_gene, geneIdproduct$V1)]))



#Sup-GEPs cell specificity?
orthoCSsup <- orthoCS %>% filter(tse_gene %in% psi_list$psi_m8)
orthoCSsup <- as.data.frame(orthoCSsup)
orthoCSsup <- orthoCSsup %>%
  mutate(
    psi_cell = ifelse(str_detect(psi_cell, "^[0-9.]+$"), "None", psi_cell),
    tse_cell = ifelse(str_detect(tse_cell, "^[0-9.]+$"), "None", tse_cell)
  )



psim2Genes <- unique(gsub("_[0-9]*","",geneIdproduct %>% filter(V1 %in% orthoCSsup$psi_gene) %>% pull(V2)))



orthoCSsup$psi_gene





orthoCSsup <- orthoCSsup %>%
  filter(
    str_detect(psi_cell, "Sertoli|Granulosa|Supporting|OCLN") |
      str_detect(tse_cell, "Sertoli|Granulosa|Supporting|OCLN")
  )
links <- orthoCSsup %>%
  group_by(psi_cell, tse_cell) %>%
  summarise(value = n(), .groups="drop") %>%
  mutate(
    flow_type = ifelse(psi_cell == tse_cell, "conserved", "changed")
  )
library(ggalluvial)
ggplot(
  links,
  aes(axis1 = psi_cell,
      axis2 = tse_cell,
      y = value)
) +
  geom_alluvium(
    aes(fill = flow_type),
    width = 0.2,
    alpha = 0.8
  ) +
  geom_stratum(
    aes(fill = after_stat(stratum)),
    width = 0.25,
    color = "black"
  ) +
  
  geom_text(
    stat = "stratum",
    aes(label = after_stat(stratum)),
    size = 4
  ) +
  scale_x_discrete(
    limits = c("PSI", "TSE"),
    expand = c(.1, .1)
  ) +
  scale_fill_manual(
    values = c(
      "conserved" = "grey70",
      "changed" = "#ED7D61"
    )
  ) +
  labs(x = NULL, y = "Gene count") +
  theme_classic() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks = element_blank()
  )





psi6 <- unique(c(res$psi_top[["8"]]))
tse8 <- res$tse_top[["3"]]
conserved_genes <- intersect(psi6,tse8)

psi_specific <- setdiff(psi6,tse8)

tse_specific <- setdiff(tse8,psi6)



} # psi x tse alone heatmap
  
 


  
  
  
} # GEP
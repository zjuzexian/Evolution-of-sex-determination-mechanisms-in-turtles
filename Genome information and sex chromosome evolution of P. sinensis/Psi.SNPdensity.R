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
#SNP density
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
                     V2 >= 1517136 & V3 <= 3178082 ~ "S3",
                     V2 >= 3178082 & V3 <= 10.4e6 ~ "S1",
                     V2 >= 10.4e6 & V3 <= 12.5 ~ "S2",
                     V2 >= 12.5e6 ~ "S0"),
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
  guides(color="none") + theme_test()
#SNP density by chromosome 
suball <- fsnp
color_panel <- c("grey90","grey40","grey90","grey40","grey90","grey40","grey90","grey40",
                 "grey90","grey40","grey90","grey40","grey90","grey40","grey90","grey40",
                 "grey90","grey40","grey90","grey40","grey90","grey40","grey90","grey40",
                 "grey90","grey40","grey90","grey40","grey90","grey40","grey90","grey40",
                 "#8cb7df","#ed7d61")
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

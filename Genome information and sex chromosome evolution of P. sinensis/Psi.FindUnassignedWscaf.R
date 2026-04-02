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
wScafs <- subset(depthHiC,depthHiC$V3<0.1&depthHiC$V4<0.1&depthHiC$V2>5000) # both Rdepth < 0.1 & LenScaf > 5k

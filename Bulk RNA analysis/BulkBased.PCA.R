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

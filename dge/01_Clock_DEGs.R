# Distance analysis (not in the server since you need some github library impossible to install in the server)
# Library to load for this step
suppressPackageStartupMessages({
library(tidyverse)
library(ggrepel)
library(BiocParallel)
library(ggpubr)
library(magrittr)
library(broom)
library(data.table)
library(cowplot)
library(bsseq)
library(methylSig)
library(dmrseq)
library(methylKit)
library(annotatr)
library(factoextra)
library(cluster)
library(ComplexHeatmap)
library(circlize)
library(ade4)
library(DESeq)
})


RPKM=read.table("RPKM_ALL.txt",sep="\t",header=T)

filter=apply(RPKM, 1, function(x) (all(x[1:3] >= 0.5) | all(x[4:6] >= 0.5) | all(x[7:9] >= 0.5) | all(x[10:12] >= 0.5) | all(x[13:15] >= 0.5) | all(x[16:18] >= 0.5)| all(x[19:21] >= 0.5)| all(x[22:24] >= 0.5)| all(x[25:27] >= 0.5)| all(x[28:30] >= 0.5) | all(x[31:33] >= 0.5) | all(x[34:36] >= 0.5) | all(x[37:39] >= 0.5) | all(x[40:42] >= 0.5) | all(x[43:45] >= 0.5) | all(x[46:48] >= 0.5)| all(x[49:51] >= 0.5)| all(x[52:54] >= 0.5)| all(x[55:57] >= 0.5)| all(x[58:60] >= 0.5) | all(x[61:63] >= 0.5) | all(x[64:66] >= 0.5) | all(x[67:69] >= 0.5) | all(x[70:72] >= 0.5)))

cnt=RPKM[,grep("ctrl",names(RPKM))]
kd=RPKM[,grep("kd",names(RPKM))]
RPKM=cbind(kd,cnt)


mat=t(log2(RPKM[filter,]+1))
mat_scaled = t(apply(mat, 1, scale))

# Distance analysis
res.dist <- dist(mat_scaled, method = "euclidean")
res.hc <- hclust(res.dist, method = "ward.D2")
pdf("Miles_RPKM_hc.pdf")
plot(res.hc, cex = 0.5)
dev.off()

res.pca <- dudi.pca(mat_scaled, scannf = FALSE, nf = 5)
Class=c(rep("KD",42),rep("CTL",42))
pdf("PCA_Clock.pdf")
fviz_pca_ind(res.pca, label="none", habillage=Class,
             addEllipses=TRUE, ellipse.level=0.95)+theme_minimal()
dev.off()

### Differential analysis
# Do for cont vs kd.
# Filter RPKM and Count by cutoff and treatments
# Change the parameter accordingly
RPKM=read.table("RPKM_ALL.txt",sep="\t",header=T)
COUNT=read.table("COUNT_ALL.txt",sep="\t",header=T)

filter=apply(RPKM, 1, function(x) (all(x[1:3] >= 0.5) | all(x[4:6] >= 0.5) | all(x[7:9] >= 0.5) | all(x[10:12] >= 0.5) | all(x[13:15] >= 0.5) | all(x[16:18] >= 0.5)| all(x[19:21] >= 0.5)| all(x[22:24] >= 0.5)| all(x[25:27] >= 0.5)| all(x[28:30] >= 0.5) | all(x[31:33] >= 0.5) | all(x[34:36] >= 0.5) | all(x[37:39] >= 0.5) | all(x[40:42] >= 0.5) | all(x[43:45] >= 0.5) | all(x[46:48] >= 0.5)| all(x[49:51] >= 0.5)| all(x[52:54] >= 0.5)| all(x[55:57] >= 0.5)| all(x[58:60] >= 0.5) | all(x[61:63] >= 0.5) | all(x[64:66] >= 0.5) | all(x[67:69] >= 0.5) | all(x[70:72] >= 0.5)))
COUNT0.5 <- COUNT[filter,]
cnt=COUNT0.5[,grep("ctrl",names(COUNT0.5))]
kd=COUNT0.5[,grep("kd",names(COUNT0.5))]
COUNT0.5=cbind(cnt,kd)
Design <- data.frame(row.names=colnames(COUNT0.5), condition = c(rep("CNT",42), rep("KD",42)))
cds <- newCountDataSet(COUNT0.5, Design$condition)
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds, method="per-condition",sharingMode = "fit-only")
deg = nbinomTest( cds, "CNT", "KD" )
deg$abs=abs(deg$log2FoldChange)
sign=deg[deg$padj < 0.05 & deg$abs > 0.3,]
write.table(deg,"CTRvsKD_DESeq_ALL.txt",sep="\t",quote=F)
write.table(sign,"CTRvsKD_DESeq_FC_FDR.txt",sep="\t",quote=F)

# Diff Expression with DESeq2 for TimeSeries
COUNT=read.table("COUNT_ALL.txt",sep="\t",header=T)
RPKM=read.table("RPKM_ALL.txt",sep="\t",header=T)
filter=apply(RPKM, 1, function(x) (all(x[1:3] >= 0.5) | all(x[4:6] >= 0.5) | all(x[7:9] >= 0.5) | all(x[10:12] >= 0.5) | all(x[13:15] >= 0.5) | all(x[16:18] >= 0.5)| all(x[19:21] >= 0.5)| all(x[22:24] >= 0.5)| all(x[25:27] >= 0.5)| all(x[28:30] >= 0.5) | all(x[31:33] >= 0.5) | all(x[34:36] >= 0.5) | all(x[37:39] >= 0.5) | all(x[40:42] >= 0.5) | all(x[43:45] >= 0.5) | all(x[46:48] >= 0.5)| all(x[49:51] >= 0.5)| all(x[52:54] >= 0.5)| all(x[55:57] >= 0.5)| all(x[58:60] >= 0.5) | all(x[61:63] >= 0.5) | all(x[64:66] >= 0.5) | all(x[67:69] >= 0.5) | all(x[70:72] >= 0.5)))
tab=COUNT[filter,]

DESIGN <- data.frame(row.names=colnames(tab), time=gsub("ZT|.CLKkd\\.[1-3]|.ctrl\\.[1-3]", "", colnames(tab)), treat = c(rep(c("KD","KD","KD","WT","WT","WT"),14)))
dds <- DESeqDataSetFromMatrix(countData = tab,colData = DESIGN,design = ~ time + treat + time:treat)
dds <- DESeq(dds, test="LRT", reduced = ~ time + treat)
DEGtime <- results(dds)
write.table(DEGtime,"DESeq_CLOCK_TimeCourse.txt",sep="\t",quote=F)
resultsNames(dds)

pdf("DESeq_TimeCourse.pdf")
data <- plotCounts(dds, which.min(DEGtime$padj < 0.05), 
                   intgroup=c("time","treat"), returnData=TRUE)
ggplot(data, aes(x=time, y=count, color=treat, group=treat)) + 
  geom_point() + stat_smooth(se=FALSE,method="loess") +  scale_y_log10()+theme_minimal()
dev.off()



### Differential analysis 24-52
RPKM=read.table("RPKM_ALL.txt",sep="\t",header=T)
COUNT=read.table("COUNT_ALL.txt",sep="\t",header=T)

RPKM=RPKM[,grep("24|28|32|36|40|44|48|52",names(RPKM))]
COUNT=COUNT[,grep("24|28|32|36|40|44|48|52",names(COUNT))]

filter=apply(RPKM, 1, function(x) (all(x[1:3] >= 0.5) | all(x[4:6] >= 0.5) | all(x[7:9] >= 0.5) | all(x[10:12] >= 0.5) | all(x[13:15] >= 0.5) | all(x[16:18] >= 0.5)| all(x[19:21] >= 0.5)| all(x[22:24] >= 0.5)| all(x[25:27] >= 0.5)| all(x[28:30] >= 0.5) | all(x[31:33] >= 0.5) | all(x[34:36] >= 0.5) | all(x[37:39] >= 0.5) | all(x[40:42] >= 0.5) | all(x[43:45] >= 0.5) | all(x[46:48] >= 0.5)))
COUNT0.5 <- COUNT[filter,]
cnt=COUNT0.5[,grep("ctrl",names(COUNT0.5))]
kd=COUNT0.5[,grep("kd",names(COUNT0.5))]
COUNT0.5=cbind(cnt,kd)
Design <- data.frame(row.names=colnames(COUNT0.5), condition = c(rep("CNT",24), rep("KD",24)))
cds <- newCountDataSet(COUNT0.5, Design$condition)
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds, method="per-condition",sharingMode = "fit-only")
deg = nbinomTest( cds, "CNT", "KD" )
deg$abs=abs(deg$log2FoldChange)
sign=deg[deg$padj < 0.05 & deg$abs > 0.3,]
write.table(deg,"CTRvsKD_DESeq_ALL_24_52.txt",sep="\t",quote=F)
write.table(sign,"CTRvsKD_DESeq_FC_FDR_24_52.txt",sep="\t",quote=F)

# Diff Expression with DESeq2 for TimeSeries
COUNT=read.table("COUNT_ALL.txt",sep="\t",header=T)
RPKM=read.table("RPKM_ALL.txt",sep="\t",header=T)


RPKM=RPKM[,grep("24|28|32|36|40|44|48|52",names(RPKM))]
COUNT=COUNT[,grep("24|28|32|36|40|44|48|52",names(COUNT))]

filter=apply(RPKM, 1, function(x) (all(x[1:3] >= 0.5) | all(x[4:6] >= 0.5) | all(x[7:9] >= 0.5) | all(x[10:12] >= 0.5) | all(x[13:15] >= 0.5) | all(x[16:18] >= 0.5)| all(x[19:21] >= 0.5)| all(x[22:24] >= 0.5)| all(x[25:27] >= 0.5)| all(x[28:30] >= 0.5) | all(x[31:33] >= 0.5) | all(x[34:36] >= 0.5) | all(x[37:39] >= 0.5) | all(x[40:42] >= 0.5) | all(x[43:45] >= 0.5) | all(x[46:48] >= 0.5)))
tab=COUNT[filter,]

DESIGN <- data.frame(row.names=colnames(tab), time=gsub("ZT|.CLKkd\\.[1-3]|.ctrl\\.[1-3]", "", colnames(tab)), treat = c(rep(c("KD","KD","KD","WT","WT","WT"),8)))
dds <- DESeqDataSetFromMatrix(countData = tab,colData = DESIGN,design = ~ time + treat + time:treat)
dds <- DESeq(dds, test="LRT", reduced = ~ time + treat)
DEGtime <- results(dds)
write.table(DEGtime,"DESeq_CLOCK_TimeCourse_24_52.txt",sep="\t",quote=F)
resultsNames(dds)

pdf("DESeq_TimeCourse_24_52.pdf")
data <- plotCounts(dds, which.min(DEGtime$padj < 0.05), 
                   intgroup=c("time","treat"), returnData=TRUE)
ggplot(data, aes(x=time, y=count, color=treat, group=treat)) + 
  geom_point() + stat_smooth(se=FALSE,method="loess") +  scale_y_log10()+theme_minimal()
dev.off()

















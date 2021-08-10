# Distance analysis (not in the server since you need some github library impossible to install in the server)
# Library to load for this step
library(DESeq2)

# Deg By Time
# CT17 vs 30 min
rm(list=ls())
COUNT=read.table("CLOCK_ATAC_HP_PEAKCOUNTS.txt")
RPKM=read.table("CLOCK_ATAC_HP_PEAKRPKM.txt")
first=apply(RPKM, 1, function(x) (all(x[1:3] > 0.5) | all(x[4:6] > 0.5)))
COUNT=COUNT[first,]

DESIGN=data.frame(row.names = colnames(COUNT), Treatment=c(rep("ctrl",3),rep("hp",3)))
dds <- DESeqDataSetFromMatrix(countData = COUNT,colData = DESIGN,design = ~ Treatment)
dds$Treatment <- factor(dds$Treatment, levels=c("ctrl", "hp"))
dds <- estimateSizeFactors(dds)
dds <- DESeq(dds, full=design(dds))
res <- as.data.frame(results(dds,contrast=c("Treatment","ctrl","hp")))
write.table(res,"DESeq_CLOCK_ATAC_HP.txt",sep="\t",quote=F)






















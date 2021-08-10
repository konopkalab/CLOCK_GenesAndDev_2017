suppressPackageStartupMessages({
library(Cairo)
library(ChIPseeker)
library(clusterProfiler)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
})

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
files <- list.files(".", pattern = ".bed", full.names = TRUE)
promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
tagMatrix <- lapply(files, getTagMatrix, windows=promoter)
tagMatrixList <- setNames(tagMatrix,c(files))
peakAnno <- lapply(files, annotatePeak, TxDb=txdb, tssRegion=c(-3000, 3000), verbose=FALSE)
peakAnnoList <- setNames(peakAnno,c(files))
genes = lapply(peakAnnoList, function(i) as.data.frame(i)$geneId)
names(genes) = sub("_", "\n", names(genes))

pdf("signal.compilation.pdf")
par(mfrow=c(10,10), mar=c(1,1,1,1))
plotAvgProf(tagMatrixList, xlim=c(-3000, 3000),facet="row")
dev.off()

pdf("BA10_UpSetHist.pdf",width=7,height=6)
upsetplot(peakAnno[[1]], vennpie=F)
dev.off()

pdf("BA40_UpSetHist.pdf",width=7,height=6)
upsetplot(peakAnno[[2]], vennpie=F)
dev.off()

pdf("signal.compilation_all.pdf")
par(mfrow=c(10,10), mar=c(1,1,1,1))
plotAvgProf(tagMatrixList, xlim=c(-3000, 3000))
dev.off()

pdf("tag.heatmap.pdf")
par(mfrow=c(10,10), mar=c(1,1,1,1))
tagHeatmap(tagMatrixList, xlim=c(-3000, 3000), color=NULL)
dev.off()

pdf("feature.distribution.pdf")
plotAnnoBar(peakAnnoList)
dev.off()

pdf("PIEs.pdf")
plotAnnoPie(peakAnno[[1]])
plotAnnoPie(peakAnno[[2]])
dev.off()

pdf("feature.distribution_acrossTss.pdf")
plotDistToTSS(peakAnnoList)
dev.off()

pdf("relative2tss.pdf")
plotAnnoBar(peakAnnoList)
dev.off()

# Overlap by consensus
ba10=genes[1:5]
ba10=unique(unlist(ba10))
ba40=genes[6:10]
ba40=unique(unlist(ba40))


list=list(ba10,ba40)
names(list)=c("BA10","BA40")

pdf("venn.pdf")
vennplot(list)
dev.off()

pdf("venn_oe.pdf")
vennplot(list1)
dev.off()


one=hp[!(hp%in%ctr)]
two=oe[!(oe%in%mch)]

library(org.Hs.eg.db)
library(annotate)

HP=as.data.frame(getSYMBOL(one, data='org.Hs.eg'))
OE=as.data.frame(getSYMBOL(two, data='org.Hs.eg'))


write.table(HP,"HP_SPECIFIC_GENES.txt",sep="\t",quote=F)
write.table(OE,"OE_SPECIFIC_GENES.txt",sep="\t",quote=F)















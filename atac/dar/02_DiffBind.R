# Call the differential peaks by diffbind
library(DiffBind)

tx=dba(sampleSheet="diffbind_ATAC_HP.csv",peakFormat='bed')
tx=dba.count(tx, minOverlap = 5,score=DBA_SCORE_TMM_READS_FULL,insertLength=200,bParallel=TRUE,bRemoveDuplicates=TRUE)
tx = dba.contrast(tx, categories=DBA_CONDITION)
tx=dba.analyze(tx, method=DBA_DESEQ, bTagwise = TRUE,bParallel=TRUE,bFullLibrarySize=TRUE)
tx.DB = dba.report(tx,method=DBA_DESEQ)
df=as.data.frame(tx.DB)
save(tx,file="CLOCK_ATAC_DiffBind_HP.RData")
write.table(df, "CLOCK_ATAC_DiffBind_HP.txt",sep="\t",quote=F)

peaksets=paste(tx$peaks[[1]]$Chr,tx$peaks[[1]]$Start,tx$peaks[[1]]$End,sep="_")
data=data.frame(row.names=peaksets,ctrhp_1=tx$peaks[[1]]$Reads,ctrhp_2=tx$peaks[[2]]$Reads,ctrhp_3=tx$peaks[[3]]$Reads,CLKhp_1=tx$peaks[[4]]$Reads,CLKhp_2=tx$peaks[[5]]$Reads,CLKhp_3=tx$peaks[[6]]$Reads)
write.table(data,"CLOCK_ATAC_HP_PEAKCOUNTS.txt",sep="\t",quote=F)

data2=data.frame(row.names=peaksets,ctrhp_1=tx$peaks[[1]]$RPKM,ctrhp_2=tx$peaks[[2]]$RPKM,ctrhp_3=tx$peaks[[3]]$RPKM,CLKhp_1=tx$peaks[[4]]$RPKM,CLKhp_2=tx$peaks[[5]]$RPKM,CLKhp_3=tx$peaks[[6]]$RPKM)
write.table(data2,"CLOCK_ATAC_HP_PEAKRPKM.txt",sep="\t",quote=F)

# OE
rm(list=ls())
tx=dba(sampleSheet="diffbind_ATAC_OE.csv",peakFormat='bed')
tx=dba.count(tx, minOverlap = 5,score=DBA_SCORE_TMM_READS_FULL,insertLength=200,bParallel=TRUE,bRemoveDuplicates=TRUE)
tx = dba.contrast(tx, categories=DBA_CONDITION)
tx=dba.analyze(tx, method=DBA_DESEQ, bTagwise = TRUE,bParallel=TRUE,bFullLibrarySize=TRUE)
tx.DB = dba.report(tx,method=DBA_DESEQ)
df=as.data.frame(tx.DB)
save(tx,file="CLOCK_ATAC_DiffBind_OE.RData")
write.table(df, "CLOCK_ATAC_DiffBind_OE.txt",sep="\t",quote=F)

peaksets=paste(tx$peaks[[1]]$Chr,tx$peaks[[1]]$Start,tx$peaks[[1]]$End,sep="_")
data=data.frame(row.names=peaksets,ctrhp_1=tx$peaks[[1]]$Reads,ctrhp_2=tx$peaks[[2]]$Reads,ctrhp_3=tx$peaks[[3]]$Reads,CLKhp_1=tx$peaks[[4]]$Reads,CLKhp_2=tx$peaks[[5]]$Reads,CLKhp_3=tx$peaks[[6]]$Reads)
write.table(data,"CLOCK_ATAC_OE_PEAKCOUNTS.txt",sep="\t",quote=F)

data2=data.frame(row.names=peaksets,ctrhp_1=tx$peaks[[1]]$RPKM,ctrhp_2=tx$peaks[[2]]$RPKM,ctrhp_3=tx$peaks[[3]]$RPKM,CLKhp_1=tx$peaks[[4]]$RPKM,CLKhp_2=tx$peaks[[5]]$RPKM,CLKhp_3=tx$peaks[[6]]$RPKM)
write.table(data2,"CLOCK_ATAC_OE_PEAKRPKM.txt",sep="\t",quote=F)

suppressPackageStartupMessages({
library(ComplexHeatmap)
library(circlize)
})

RPKM=read.table("RPKM_ALL.txt",sep="\t",header=T)

filter=apply(RPKM, 1, function(x) (all(x[1:3] >= 0.5) | all(x[4:6] >= 0.5) | all(x[7:9] >= 0.5) | all(x[10:12] >= 0.5) | all(x[13:15] >= 0.5) | all(x[16:18] >= 0.5)| all(x[19:21] >= 0.5)| all(x[22:24] >= 0.5)| all(x[25:27] >= 0.5)| all(x[28:30] >= 0.5) | all(x[31:33] >= 0.5) | all(x[34:36] >= 0.5) | all(x[37:39] >= 0.5) | all(x[40:42] >= 0.5) | all(x[43:45] >= 0.5) | all(x[46:48] >= 0.5)| all(x[49:51] >= 0.5)| all(x[52:54] >= 0.5)| all(x[55:57] >= 0.5)| all(x[58:60] >= 0.5) | all(x[61:63] >= 0.5) | all(x[64:66] >= 0.5) | all(x[67:69] >= 0.5) | all(x[70:72] >= 0.5)))

cnt=RPKM[,grep("ctrl",names(RPKM))]
kd=RPKM[,grep("kd",names(RPKM))]
RPKM=cbind(kd,cnt)

mat=log2(RPKM[filter,]+1)
gene=read.table("CTRvsKD_DESeq_FC_FDR.txt")
mat=mat[rownames(mat) %in% gene$id,]
mat_scaled = t(apply(mat, 1, scale))

ha = HeatmapAnnotation(df = DESIGN,col=list(Time = c("20" = "olivedrab4","24"="palevioletred3","28" = "cadetblue4","32"="lightcyan","36" = "orchid","40" = "cornsilk4","44"="dodgerblue1","48" = "darkorchid4","52"="thistle1","56"="darkorange4","60"="maroon3","64"="navajowhite3","68" = "lavenderblush1","72"="darkorchid"),Treat=c("KD" = "red","CTL" = "black")))

pdf("Complex_Heatmap_Clok.pdf")
Heatmap(mat_scaled, name = "expression", km = 2, col = colorRamp2(c(-2, 0, 2), c("cornflowerblue", "white", "gold")),
     top_annotation = ha, top_annotation_height = unit(4, "mm"), 
     show_row_names = FALSE, show_column_names = FALSE,cluster_columns = FALSE)
dev.off()

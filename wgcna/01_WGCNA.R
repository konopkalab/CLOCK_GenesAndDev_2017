library(WGCNA);
library(cluster);
library(ggplot2)
library(reshape2)
library(RColorBrewer)
options(stringsAsFactors = FALSE);
enableWGCNAThreads()

RPKM=read.table("RPKM_ALL.txt",sep="\t",header=T)
cnt=RPKM[,grep("ctrl",names(RPKM))]
kd=RPKM[,grep("kd",names(RPKM))]
RPKM=cbind(kd,cnt)

filter=apply(RPKM, 1, function(x) (all(x[1:42] >= 0.5) | all(x[43:84] >= 0.5)))
tab=RPKM[filter,]
tab=log2(tab+1)
datExpr <- as.data.frame(t(tab));
names(datExpr) <- rownames(tab);
rownames(datExpr) <- names(tab);

## Powers analysis

powers = c(seq(2,30,2))
sft=pickSoftThreshold(datExpr,powerVector=powers,verbose = 5, blockSize=12000, networkType = "signed") 
pdf("SoftThresholdingPower_signed.pdf")
par(mfrow = c(1,2), mar=c(5.1,5.1,4.1,2.1));
cex1 = 0.9;
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit, signed R^2",type="n",main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers,cex=cex1,col="red");
abline(h=0.5,col="red"); abline(h=0.8,col="blue")
plot(sft$fitIndices[,1], sft$fitIndices[,5],xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

######################################################################################################################
############################################ MODULE construction #####################################################
############################################    signed network   #####################################################
PWR=22
net = blockwiseModules(datExpr,corType="pearson", maxBlockSize = 12000, networkType="signed",power=PWR, minModuleSize=50,nThreads=15,
TOMType = "signed",TOMDenom = "mean",deepSplit=2,verbose=5,mergeCutHeight=0.15,detectCutHeight = 0.999,reassignThreshold = 1e-10,numericLabels=TRUE,saveTOMs=TRUE,pamStage=TRUE, pamRespectsDendro=TRUE,saveTOMFileBase="TOM_SIGNED")
moduleLabelsAutomatic=net$colors
moduleColorsAutomatic = labels2colors(moduleLabelsAutomatic)
MEsAutomatic=net$MEs
unique(moduleColorsAutomatic)
save(net,file="CLOCK_015_pw16_ModSize50.RData")
write.table(moduleColorsAutomatic, "CLOCK_colors.txt",sep="\t",quote=F)

#KMEs
KMEs<-signedKME(datExpr, net$MEs,corFnc = "cor", corOptions = "use = 'p'")
kme=data.frame(rownames(tab), moduleColorsAutomatic, KMEs)
kme$rownames.tab.=NULL
write.table(kme,"KMEs_CLOCK.txt",sep="\t",quote=F)


## Network Dendogram
# Load the tables we need
# Single column tables containing the gene sets needed

load("CLOCK_Degs_GeneSets.RData")
ZT20=degs[[1]]
ZT24=degs[[2]]
ZT28=degs[[3]]
ZT32=degs[[4]]
ZT36=degs[[5]]
ZT40=degs[[6]]
ZT44=degs[[7]]
ZT48=degs[[8]]
ZT52=degs[[9]]
ZT56=degs[[10]]
ZT60=degs[[11]]
ZT64=degs[[12]]
ZT68=degs[[13]]
ZT72=degs[[14]]

pdf("NetworkDendrogram.pdf")
par(mfrow=c(2,1))
ZT20colors=rep(rgb(0,0,1,0.01), nrow(tab)); 
ZT24colors=rep(rgb(0,0,1,0.01), nrow(tab)); 
ZT28colors=rep(rgb(0,0,1,0.01), nrow(tab)); 
ZT32colors=rep(rgb(0,0,1,0.01), nrow(tab)); 
ZT36colors=rep(rgb(0,0,1,0.01), nrow(tab));
ZT40colors=rep(rgb(0,0,1,0.01), nrow(tab)); 
ZT44colors=rep(rgb(0,0,1,0.01), nrow(tab)); 
ZT48colors=rep(rgb(0,0,1,0.01), nrow(tab)); 
ZT52colors=rep(rgb(0,0,1,0.01), nrow(tab)); 
ZT56colors=rep(rgb(0,0,1,0.01), nrow(tab)); 
ZT60colors=rep(rgb(0,0,1,0.01), nrow(tab)); 
ZT64colors=rep(rgb(0,0,1,0.01), nrow(tab)); 
ZT68colors=rep(rgb(0,0,1,0.01), nrow(tab)); 
ZT72colors=rep(rgb(0,0,1,0.01), nrow(tab)); 

ZT20colors[colnames(datExpr)%in%ZT20$Gene]="red"
ZT24colors[colnames(datExpr)%in%ZT24$Gene]="red"
ZT28colors[colnames(datExpr)%in%ZT28$Gene]="red"
ZT32colors[colnames(datExpr)%in%ZT32$Gene]="red"
ZT36colors[colnames(datExpr)%in%ZT36$Gene]="red"
ZT40colors[colnames(datExpr)%in%ZT40$Gene]="red"
ZT44colors[colnames(datExpr)%in%ZT44$Gene]="red"
ZT48colors[colnames(datExpr)%in%ZT48$Gene]="red"
ZT52colors[colnames(datExpr)%in%ZT52$Gene]="red"
ZT56colors[colnames(datExpr)%in%ZT56$Gene]="red"
ZT60colors[colnames(datExpr)%in%ZT60$Gene]="red"
ZT64colors[colnames(datExpr)%in%ZT64$Gene]="red"
ZT68colors[colnames(datExpr)%in%ZT68$Gene]="red"
ZT72colors[colnames(datExpr)%in%ZT72$Gene]="red"

plotColors=cbind(moduleColorsAutomatic, ZT20colors, ZT24colors, ZT28colors, ZT32colors, ZT36colors, ZT40colors, ZT44colors, ZT48colors, ZT52colors,ZT56colors, ZT60colors, ZT64colors, ZT68colors, ZT72colors); 
colnames(plotColors)=c("module", "ZT20","ZT24","ZT28","ZT32","ZT36","ZT40","ZT44","ZT48","ZT52","ZT56","ZT60","ZT64","ZT68","ZT72")
plotDendroAndColors(net$dendrograms[[1]],plotColors, dendroLabels=FALSE )
dev.off()

# Reclusterd without the grey
restGenes=(moduleColorsAutomatic != "grey")
diss1=1-TOMsimilarityFromExpr(datExpr[,restGenes], power = PWR,corType = "pearson",networkType="signed",TOMType="signed",TOMDenom = "mean",nThreads = 15,verbose = 5, indent = 0)
colnames(diss1) = rownames(diss1) = rownames(tab)[restGenes]
hier1=flashClust(as.dist(diss1), method="average" )
plotColors =as.data.frame(plotColors)
nogrey= plotColors[plotColors$module != "grey",]
pdf("NetworkDendrogram_nogrey.pdf")
par(mfrow=c(2,1))
plotDendroAndColors(hier1, nogrey, dendroLabels=FALSE )
dev.off()

# Make Ggplot2 input for WGCNA eigengene values
nGenes = ncol(datExpr) 
nSamples = nrow(datExpr)
MEs0 = moduleEigengenes(datExpr,moduleColorsAutomatic,softPower = PWR,impute = TRUE)$eigengenes
MEs0$Rows=names(tab)
MEs0$Time=gsub(".CLKkd\\.[1-3]|.ctrl\\.[1-3]", "", colnames(tab))
MEs0$Class = c(rep(c("KD","KD","KD","WT","WT","WT"),14))
write.table(MEs0, "Matrix_module_correlation.txt",sep="\t",quote=F)

# Define the association between genes and detected modules.
Adj = adjacency(datExpr, power = PWR,type="signed",corFnc = "cor", corOptions = "use = 'p'")
moduleOutput <- data.frame(rownames(tab))
moduleOutput[,2]<- moduleColorsAutomatic
intraCon <- intramodularConnectivity(Adj, moduleColorsAutomatic)
moduleOutput[,3]<-intraCon$kWithin
colnames(moduleOutput) <- c("Gene", "ModuleColor", "kWithin")
head(moduleOutput)
moduleOutput=moduleOutput[which(moduleOutput$ModuleColor != "grey"),]
write.table(moduleOutput, "ModuleOutput_CLOCK.txt", sep="\t", quote=F)

TOM = TOMsimilarityFromExpr(datExpr, power= PWR,corType = "pearson",networkType="signed",TOMType="signed",TOMDenom = "mean",nThreads = 15,verbose = 5, indent = 0)
colnames(TOM)=rownames(TOM)=colnames(datExpr)
save(TOM,file="TOM_CLOCK_SIGNED.RData")

## Create a database for Table S6
load("CLOCK_GeneSets.RData")
GeneSets$Mod=moduleOutput
LIST=append(GeneSets,degs) 
DB <- Reduce(function(x, y) {
    merge(x, y, all=TRUE, by="Gene")
}, LIST)
DB=DB[which(DB$kWithin != "NA"),]
save(LIST,file="GeneSets_WGCNAmodule_CLOCK.RData")
write.table(DB,"WGCNA_DataBase_CLOCK.txt",sep="\t",quote=F)

# CytoScape output for detected modules
dir.create("Cyto")
setwd("Cyto")
for(module in unique(moduleColorsAutomatic)){
inModule <- is.finite(match(moduleColorsAutomatic, module))
modTOM <- TOM[inModule, inModule]
cyt = exportNetworkToCytoscape(modTOM, edgeFile=paste("CytoEdge",paste(module,collapse="-"),".txt",sep=""), nodeFile=paste("CytoNode",paste(module,collapse="-"),".txt",sep=""), weighted = TRUE, threshold = 0, nodeAttr = moduleColorsAutomatic[inModule], nodeNames = names(datExpr)[inModule])
}

# WGCNA eigengenes output
library(ggplot2)
library(reshape2)
library(RColorBrewer)
df=melt(MEs0)
dir.create("Plot")
setwd("Plot")
df$Rows=factor(df$Rows, levels = colnames(tab))
df=df[!df$variable == "MEgrey", ]
# Function to make the barplot and saving as pdf according to the module name
cols <- c("red","blue")
doPlot = function(sel_name) 
{
    df = subset(df, variable == sel_name)
    PLOT=ggplot(data=df, aes(x=Rows, y=value)) +
 geom_bar(aes(fill=Class),stat="identity",position = "dodge")+
 scale_y_continuous(limits=c(-1,+1))+
 theme_bw()+
 theme(strip.text.x = element_text(size=12, face="bold"),
 strip.background = element_rect(colour="black", fill="#CCCCFF"))+
 scale_fill_manual(values = cols)+
 theme(axis.title.x = element_blank(),
            axis.text.x  = element_text(face="bold", size=6,angle = 45, hjust = 1))+
 theme(axis.title.y = element_blank(),
            axis.text.y  = element_text(face="bold", size=6))
    print(PLOT)
    ggsave(sprintf("%s.pdf", sel_name))
 }

lapply(unique(df$variable), doPlot)












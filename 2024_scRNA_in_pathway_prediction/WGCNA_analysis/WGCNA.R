library('preprocessCore')
library('WGCNA')
library('cluster')
library('dplyr')

Exp <- read.csv("example_data.csv", row.names=1)
Exp <- Exp[apply(Exp!=0, 1, any),]
dataExpr <- as.matrix(Exp)
dataExprVar <- Exp
type = "signed"
corType = "pearson"
corFnc = ifelse(corType=="pearson", cor, bicor)
dataExpr <- as.data.frame(t(dataExprVar))
sampleTree = hclust(dist(dataExpr), method = "average")
pdf('Sample_clustering_to_detect_outliers.pdf', 22, 22)
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="")
dev.off()

powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(dataExpr, powerVector=powers, networkType=type, verbose=5)
pdf('WGCNA_power.pdf', 10, 10)
par(mfrow = c(2,1))
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
abline(h=0.85, col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1, col="red")
dev.off()

best_beta = sft$powerEstimate
enableWGCNAThreads()
softPower = best_beta
adjacency = adjacency(dataExpr, power = softPower, type = "signed")
TOM = TOMsimilarity(adjacency, TOMType="signed")
dissTOM = 1-TOM

library(flashClust)
library(infotheo)
geneTree = flashClust(as.dist(dissTOM), method="average")
pdf('Gene_Clustering_on_TOM-based_dissimilarity.pdf')
plot(geneTree, xlab="", sub="", main="Gene Clustering on TOM-based dissimilarity", labels=FALSE, hang=0.04)
dev.off()

minModuleSize = 30
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM, deepSplit = 4, pamRespectsDendro = FALSE, minClusterSize = minModuleSize)
# table(dynamicMods)
dynamicColors = labels2colors(dynamicMods)
# table(dynamicColors)
MEList = moduleEigengenes(dataExpr, colors = dynamicColors, softPower = softPower)
MEs = MEList$eigengenes
MEDiss = 1-cor(MEs, use = "pairwise.complete.obs")
MEDiss[is.na(MEDiss)] <- 0
METree = flashClust(as.dist(MEDiss), method= "average")
pdf('Clustering_of_module_eigengenes.pdf')
plot(METree, main = "Clustering of module eigengenes", xlab= "", sub= "")
dev.off()

MEDissThres = 0.1
# abline(h=MEDissThres, col = "red")
merge = mergeCloseModules(dataExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
mergedColors = merge$colors
mergedMEs = merge$newMEs
geneordered <-colnames(dataExpr)[geneTree$order]
colorlabel <- cbind(colnames(dataExpr),mergedColors)
rownames(colorlabel) <- colorlabel[,1]
color_order <- colorlabel[geneordered,2]
write.csv(color_order, 'exp_WGCNA_MEDissThres_colors_0.1.csv', row.names=T, quote=F)

pdf(file="exp_WGCNA_MEDissThres_cluster_0.1.pdf", 8, 6)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), c("Dynamic Tree Cut", "Merged dynamic"), dendroLabels=FALSE, hang=0.03, addGuide=TRUE, guideHang=0.05)
moduleColors = mergedColors
colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder)-1
MEs = mergedMEs
dev.off()

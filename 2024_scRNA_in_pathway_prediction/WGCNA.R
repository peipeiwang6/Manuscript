library(preprocessCore)
library('WGCNA', warn.conflicts = FALSE)
library(cluster)
filepath <- getwd()
setwd(filepath)

args = commandArgs(TRUE)
input_csv <- args[1]
output_name <- args[2]

dataExpr <- as.data.frame(read.csv(input_csv,head=T,row.names=1))
NotExp <- dataExpr[apply(dataExpr==0, 1, all),]
dataExpr <- dataExpr[apply(dataExpr!=0, 1, any),]
log_data <- log10(dataExpr+1)
dataExpr <- t(log_data)
type = "signed"
corType = "pearson"
corFnc = ifelse(corType=="pearson", cor, bicor)
powers = c(c(1:10), seq(from = 12, to=70, by=2))
sft = pickSoftThreshold(dataExpr, powerVector=powers, networkType=type, verbose=5)
par(mfrow = c(2,1))
cex1 = 0.9

plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
abline(h=0.85,col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

softPower = sft$powerEstimate

TOM = TOMsimilarityFromExpr(dataExpr, power = softPower)

module = output_name
cyt = exportNetworkToCytoscape(TOM,
                               edgeFile = paste("CytoscapeInput_edges_", paste(module, collapse="_"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0,
                               nodeNames = colnames(dataExpr))
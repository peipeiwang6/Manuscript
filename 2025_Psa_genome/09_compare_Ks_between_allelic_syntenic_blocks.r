###################################
# kaks for reciprocal best match conserved genes
setwd('D:\\鲲鹏院工作\\课题\\03_Phalaenopsis\\25_allelic_genes')
library(gplots)
library(ggplot2)
library(gridExtra)
df <- read.table('KaKs_and_expression_dist/Psa_RBM_KaKs.txt',head=T,sep='\t')
allelic <- read.table('KaKs_and_expression_dist/Psa_winner_allelic_gene_pairs_KaKs.txt',head=T,sep='\t')
other <- read.table('KaKs_and_expression_dist/Psa_other_syntenic_gene_pairs_KaKs.txt',head=T,sep='\t')
allelic$Align <- 'syntenic'
colnames(allelic)[1] <- 'RBM_type'
other$Align <- 'syntenic'
colnames(other)[1] <- 'RBM_type'
df <- rbind(df, allelic)
df <- rbind(df, other)
pdf('Psa_RBM_dS.pdf')
ggplot(df, aes(x=dS, color=RBM_type)) + 
		geom_density(alpha=0.3, size=1)+ 
		#scale_x_log10() +
		xlim(0,3) +
		labs(x= "dS", subtitle="dS")
dev.off()

pdf('Psa_RBM_dNdS.pdf')
ggplot(df, aes(x=Omega, color=RBM_type)) + 
		geom_density(alpha=0.3, size=1)+ 
		#scale_x_log10() +
		xlim(0,2) +
		labs(x= "Omega", subtitle="Omega")
dev.off()

pdf('Psa_RBM_dN.pdf')
ggplot(df, aes(x=dN, color=RBM_type)) + 
		geom_density(alpha=0.3, size=1)+ 
		#scale_x_log10() +
		xlim(0,2) +
		labs(x= "dN", subtitle="dN")
dev.off()
















#########################################################
# dS and exp dist between inversed genes and syntenic genes
setwd('D:\\鲲鹏院工作\\课题\\03_Phalaenopsis\\25_allelic_genes')
Inversed <- read.table('Psa_inversed_genes.txt',head=F,sep='\t')
df <- read.table('KaKs_and_expression_dist/Psa_winner_allelic_gene_pairs_KaKs.txt',head=T,sep='\t')
tem <- df[df$Gene1 == Inversed[i, 1] | df$Gene2 == Inversed[i, 1], ]

















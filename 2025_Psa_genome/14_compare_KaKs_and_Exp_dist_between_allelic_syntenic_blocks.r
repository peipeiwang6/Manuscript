setwd('D:\\鲲鹏院工作\\课题\\03_Phalaenopsis\\25_allelic_genes\\KaKs_and_expression_dist')
library(gplots)
library(ggplot2)
library(gridExtra)
df <- read.table('Psa_winner_allelic_gene_pairs_KaKs.txt',head=T,sep='\t')
df <- cbind(df,'Chr' = substr(df[,2], 1, 2))
pdf('Psa_winner_allelic_gene_pairs_Ks_boxplot.pdf', onefile = TRUE, height = 5, width=10)
for(chr in unique(df$Chr)){
	chr1 <- df[df$Chr==chr,]
	mycol <- c("#BD6263","#8EA325","#84CAC0","#F5AE6B","#BCB8D3","#4387B5")
	chr1$dS <- log10(chr1$dS+0.0001)
	p <- ggplot(chr1, aes(x = Chr_pair, y=dS, color=Chr_pair)) + 
		geom_boxplot(alpha=0.3, size=1,outlier.shape = NA)+ 
		geom_point(position = position_jitter(0.2),color = "gray") + 
		labs(x= "Ks", subtitle=paste("Ks of allelic gene pairs on Chr",chr,sep=''))
	print(p)
	}
dev.off()

pdf('Psa_winner_allelic_gene_pairs_Ks_violin.pdf', onefile = TRUE, height = 5, width=10)
for(chr in unique(df$Chr)){
	chr1 <- df[df$Chr==chr,]
	mycol <- c("#BD6263","#8EA325","#84CAC0","#F5AE6B","#BCB8D3","#4387B5")
	chr1$dS <- log10(chr1$dS+0.0001)
	p <- ggplot(chr1, aes(x = Chr_pair, y=dS, color=Chr_pair)) + 
		geom_point(position = position_jitter(0.2),color = "gray") + 
		geom_violin(alpha=0.3, size=1)+ 
		stat_summary(fun='median', geom='crossbar', color='red') +
		labs(x= "Ks", subtitle=paste("Ks of allelic gene pairs on Chr",chr,sep=''))
	print(p)
	}
dev.off()


pdf('Psa_winner_allelic_gene_pairs_Ks_density_plot.pdf', onefile = TRUE, height = 5, width=10)
for(chr in unique(df$Chr)){
	chr1 <- df[df$Chr==chr,]
	mycol <- c("#BD6263","#8EA325","#84CAC0","#F5AE6B","#BCB8D3","#4387B5")
	chr1$dS <- log10(chr1$dS)
	p <- ggplot(chr1, aes(x=dS, color=Chr_pair)) + 
		geom_density(alpha=0.3, size=1)+ 
		labs(x= "Ks", subtitle=paste("Ks of allelic gene pairs on Chr",chr,sep=''))
	print(p)
	}
dev.off()

# dS
setwd('D:\\鲲鹏院工作\\课题\\03_Phalaenopsis\\25_allelic_genes\\')
df <- read.table('KaKs_and_expression_dist/Psa_winner_allelic_gene_pairs_KaKs.txt',head=T,sep='\t')
df2 <- read.table('KaKs_and_expression_dist/Psa_other_syntenic_gene_pairs_KaKs.txt',head=T,sep='\t')
df2 <- unique(df2)
df <- cbind(df, 'Chr1'=substr(df$Chr_pair,1,3), 'Chr2'=substr(df$Chr_pair,5,7), 'Type'='allelic_chr')
df2 <- cbind(df2, 'Chr1'=substr(df2$Chr_pair,1,3), 'Chr2'=substr(df2$Chr_pair,5,7), 'Type'='other')
#df2 <- cbind(df2, 'Chr1'=substr(df2$Chr_pair,1,3), 'Chr2'=substr(df2$Chr_pair,5,7), 'Type'='same_chr')
#df2[df2$Chr1 != df2$Chr2,]$Type = 'diff_chr'
res <- rbind(df, df2)
pdf('dS_of_syntenic_genes_on_allelic_and_others.pdf')
ggplot(res, aes(x=dS, color=Type)) + 
		geom_density(alpha=0.3, size=1)+ 
		scale_x_log10() +
		labs(x= "dS", subtitle="Median dS of syntenic gene pairs")
dev.off()

pdf('Omega_of_syntenic_genes_on_allelic_and_others.pdf')
ggplot(res, aes(x=Omega, color=Type)) + 
		geom_density(alpha=0.3, size=1)+ 
		scale_x_log10() +
		labs(x= "Omega", subtitle="Median Omega of syntenic gene pairs")
dev.off()


# median dS of block
df <- read.table('Psa_winner_allelic_gene_pairs_KaKs.txt',head=T,sep='\t')
df_median <- data.frame()
for(align in unique(df$Align)){	
	tem <- c('Align'=align, 'Median_dS' = summary(df[df$Align == align,]$dS)[3])
	df_median <- rbind(df_median, tem)
	}
df_median <- cbind(df_median, 'Type' = 'allelic')
colnames(df_median) <- c('Align', 'Median_dS','Type')

df2 <- read.table('Psa_other_syntenic_gene_pairs_KaKs.txt',head=T,sep='\t')
df2 <- unique(df2)
df_median2 <- data.frame()
for(align in unique(df2$Align)){	
	tem <- c('Align'=align, 'Median_dS' = summary(df2[df2$Align == align,]$dS)[3])
	df_median2 <- rbind(df_median2, tem)
	}
df_median2 <- cbind(df_median2, 'Type' = 'other')
colnames(df_median2) <- c('Align', 'Median_dS','Type')

res <- rbind(df_median, df_median2)

pdf('Median_dS_between_allelic_and_other_syntenic_blocks.pdf')
ggplot(res, aes(x=Median_dS, color=Type)) + 
		geom_density(alpha=0.3, size=1)+ 
		scale_x_log10() +
		labs(x= "Ks", subtitle="Median Ks of syntenic gene pairs")
dev.off()


# omega
df <- read.table('Psa_winner_allelic_gene_pairs_KaKs.txt',head=T,sep='\t')
df_median <- data.frame()
for(align in unique(df$Align)){	
	tem <- c('Align'=align, 'Median_Omega' = summary(df[df$Align == align,]$Omega)[3])
	df_median <- rbind(df_median, tem)
	}
df_median <- cbind(df_median, 'Type' = 'allelic')
colnames(df_median) <- c('Align', 'Median_Omega','Type')

df2 <- read.table('Psa_other_syntenic_gene_pairs_KaKs.txt',head=T,sep='\t')
df2 <- unique(df2)
df_median2 <- data.frame()
for(align in unique(df2$Align)){	
	tem <- c('Align'=align, 'Median_Omega' = summary(df2[df2$Align == align,]$Omega)[3])
	df_median2 <- rbind(df_median2, tem)
	}
df_median2 <- cbind(df_median2, 'Type' = 'other')
colnames(df_median2) <- c('Align', 'Median_Omega','Type')

res <- rbind(df_median, df_median2)

pdf('Median_Omega_between_allelic_and_other_syntenic_blocks.pdf')
ggplot(res, aes(x=Median_Omega, color=Type)) + 
		geom_density(alpha=0.3, size=1)+ 
		scale_x_log10() +
		labs(x= "Omega", subtitle="Median Omega of syntenic gene pairs")
dev.off()


setwd('D:\\鲲鹏院工作\\课题\\03_Phalaenopsis\\25_allelic_genes\\KaKs_and_expression_dist')
library(gplots)
library(ggplot2)
library(gridExtra)
##################################
# expression dist
df <- read.table('Psa_winner_allelic_gene_pair_expression_euclid_dist.txt',head=T,sep='\t')
df <- cbind(df,'Chr' = substr(df[,2], 1, 2))
pdf('Psa_winner_allelic_gene_pairs_expression_dist_violin.pdf', onefile = TRUE, height = 5, width=10)
for(chr in unique(df$Chr)){
	chr1 <- df[df$Chr==chr,]
	mycol <- c("#BD6263","#8EA325","#84CAC0","#F5AE6B","#BCB8D3","#4387B5")
	chr1$Euclid_dist <- log10(chr1$Euclid_dist+0.0001)
	p <- ggplot(chr1, aes(x = Chr_pair, y=Euclid_dist, color=Chr_pair)) + 
		geom_point(position = position_jitter(0.2),color = "gray") + 
		geom_violin(alpha=0.3, size=1,outlier.shape = NA)+ 
		stat_summary(fun='median', geom='crossbar', color='red') +
		labs(x= "Euclid_dist", subtitle=paste("Euclid_dist of allelic gene pairs on Chr",chr,sep=''))
	print(p)
	}
dev.off()


# median expression dist across syntenic blocks, to remove genes that are totally not expressed. 20250312
df <- read.table('Psa_winner_allelic_gene_pair_expression_euclid_dist.txt',head=T,sep='\t')
df_median <- data.frame()
for(align in unique(df$Align)){	
	tem <- c('Align'=align, 'Median_Euclid_dist' = summary(df[df$Align == align,]$Euclid_dist)[3])
	df_median <- rbind(df_median, tem)
	}
df_median <- cbind(df_median, 'Type' = 'allelic')
colnames(df_median) <- c('Align', 'Median_Euclid_dist','Type')

df2 <- read.table('Psa_other_syntenic_gene_pair_expression_euclid_dist.txt',head=T,sep='\t')
df2 <- unique(df2)
df_median2 <- data.frame()
for(align in unique(df2$Align)){	
	tem <- c('Align'=align, 'Median_Euclid_dist' = summary(df2[df2$Align == align,]$Euclid_dist)[3])
	df_median2 <- rbind(df_median2, tem)
	}
df_median2 <- cbind(df_median2, 'Type' = 'other')
colnames(df_median2) <- c('Align', 'Median_Euclid_dist','Type')

res <- rbind(df_median, df_median2)

pdf('Median_expression_Euclid_dist_between_allelic_and_other_syntenic_blocks.pdf')
# ggplot(res, aes(x=log10(Median_Euclid_dist+0.0001), color=Type)) + 
		# geom_density(alpha=0.3, size=1)+ 
		# scale_x_log10() +
		# labs(x= "Euclid_dist", subtitle="Median Euclid_dist of syntenic gene pairs")
# dev.off()
ggplot(res, aes(x=Median_Euclid_dist, color=Type)) + 
		geom_density(alpha=0.3, size=1)+ 
		scale_x_log10() +
		labs(x= "Euclid_dist", subtitle="Median Euclid_dist of syntenic gene pairs")
dev.off()


df <- read.table('Psa_winner_allelic_gene_pair_expression_euclid_dist.txt',head=T,sep='\t')
df2 <- read.table('Psa_other_syntenic_gene_pair_expression_euclid_dist.txt',head=T,sep='\t')
df2 <- unique(df2)
df <- cbind(df, 'Type'='allelic')
df2 <- cbind(df2, 'Type'='other')
res <- rbind(df, df2)

pdf('Expression_dist_between_allelic_and_other_syntenic_blocks.pdf')
ggplot(res, aes(x=Euclid_dist, color=Type)) + 
		geom_density(alpha=0.3, size=1)+ 
		scale_x_log10() +
		labs(x= "Euclid_dist", subtitle="Median Euclid_dist of syntenic gene pairs")
dev.off()

# expression dist, three types
df <- read.table('Psa_winner_allelic_gene_pair_expression_euclid_dist.txt',head=T,sep='\t')
df2 <- read.table('Psa_other_syntenic_gene_pair_expression_euclid_dist.txt',head=T,sep='\t')
df2 <- unique(df2)
df <- cbind(df, 'Chr1'=substr(df$Chr_pair,1,3), 'Chr2'=substr(df$Chr_pair,5,7), 'Type'='allelic_chr')
df2 <- cbind(df2, 'Chr1'=substr(df2$Chr_pair,1,3), 'Chr2'=substr(df2$Chr_pair,5,7), 'Type'='others')
#df2 <- cbind(df2, 'Chr1'=substr(df2$Chr_pair,1,3), 'Chr2'=substr(df2$Chr_pair,5,7), 'Type'='same_chr')
#df2[df2$Chr1 != df2$Chr2,]$Type = 'diff_chr'
res <- rbind(df, df2)
pdf('Expression_dist_of_syntenic_genes_on_allelic_and_others.pdf')
ggplot(res, aes(x=Euclid_dist, color=Type)) + 
		geom_density(alpha=0.3, size=1)+ 
		scale_x_log10() +
		labs(x= "Euclid_dist", subtitle="Median Euclid_dist of syntenic gene pairs")
dev.off()


# correlation betweed dS and expression dist
df <- read.table('Psa_winner_allelic_gene_pairs_KaKs.txt',head=T,sep='\t')
df2 <- read.table('Psa_other_syntenic_gene_pairs_KaKs.txt',head=T,sep='\t')
df2 <- unique(df2)
df <- cbind(df, 'Chr1'=substr(df$Chr_pair,1,3), 'Chr2'=substr(df$Chr_pair,5,7), 'Type'='allelic_chr')
df2 <- cbind(df2, 'Chr1'=substr(df2$Chr_pair,1,3), 'Chr2'=substr(df2$Chr_pair,5,7), 'Type'='other')
res_dS <- rbind(df, df2)
res_dS <- cbind(res_dS, 'Gene_pair' = paste(res_dS$Gene1, '__', res_dS$Gene2,sep=''))

df <- read.table('Psa_winner_allelic_gene_pair_expression_euclid_dist.txt',head=T,sep='\t')
df2 <- read.table('Psa_other_syntenic_gene_pair_expression_euclid_dist.txt',head=T,sep='\t')
df2 <- unique(df2)
df <- cbind(df, 'Type'='allelic')
df2 <- cbind(df2, 'Type'='other')
res_exp_dist <- rbind(df, df2)
res_exp_dist <- cbind(res_exp_dist, 'Gene_pair' = paste(res_exp_dist$Gene1, '__', res_exp_dist$Gene2,sep=''))

res <- merge(res_dS, res_exp_dist, by.x='Gene_pair', by.y = 'Gene_pair')
library(pvclust)
library(gplots)
pdf('Correlation_between_dS_and_exp_dist.pdf')
smoothScatter(log10(res$dS), log10(res$Euclid_dist), xlab="dS", ylab="exp_Euclid_dist",xlim=c(-3,2),ylim=c(-3,4))
at.y <- log10(outer(1:9, 10^-(4:-1)))
lab.y <- 10^ifelse(at.y %% 1 == 0, at.y, NA)
at.x <- log10(outer(1:9, 10^-(4:-1)))
lab.x <- 10^ifelse(at.x %% 1 == 0, at.x, NA)
axis(2, at=at.y, labels=lab.y, las=1) 
axis(1, at=at.x, labels=lab.x, las=1) 
dev.off()


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

















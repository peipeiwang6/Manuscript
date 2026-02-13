library(ggplot2)
setwd('D:\\鲲鹏院工作\\课题\\03_Phalaenopsis\\18_P.sa_transcriptome')
dat <- read.csv('Psa_normalized_mean_TPM_clean_20260120-rename.csv',head=T,row.names=1,sep=',')
allelic <- read.table('../25_allelic_genes/Psa_syri_allelic_fillout_using_McScanX_and_inversion_20250328_combine_tandem_duplicates.txt',head=T,sep='\t')

Sum_exp <- function(dat, allelic, gene, type, pair){
	genes <- allelic[allelic[,1] == gene | allelic[,2] == gene | allelic[,3] == gene | allelic[,4] == gene,]
	subdat <- c()
	for(i in 1:4){
		if(genes[1,i]!=''){
			subdat <- rbind(subdat, dat[genes[1,i],])
			}
		}
	summed_exp <- t(as.data.frame(colSums(subdat)))
	rownames(summed_exp)[1] <- paste(c(rownames(subdat)),collapse = ';')
	summed_exp <- cbind(summed_exp, 'Type' = type)
	summed_exp <- cbind(summed_exp, 'Pair' = pair)
	return(summed_exp)
	}

# C1
Pairs <- read.table('Expression_dosage_20260109/Reciprocal_best_hist_Type6_set3_and_set4_genes_20260109.txt',head=F,sep='\t')
Res <- c()
for(i in 1:nrow(Pairs)){
	Res <- rbind(Res, Sum_exp(dat, allelic, Pairs[i,1], '3/3', i))
	Res <- rbind(Res, Sum_exp(dat, allelic, Pairs[i,2], '4/4', i))
	}
for(i in 1:38){
	subRes <- data.frame('Exp'=log10(as.numeric(Res[,i])+1),'Type' = Res[,'Type'], 'Pair' = Res[,'Pair'])
	ggplot(subRes, aes(x = Type, y = Exp)) +
	  geom_line(aes(group = Pair), color = "gray", alpha = 0.5) +  # 连线
	  geom_point(aes(color = Type), size = 1) +                  # 散点
	  geom_boxplot(alpha = 0.5, color = c('blue', 'red'), fill = c('blue', 'red'), width=0.2, outlier.shape = NA)+
	  scale_color_manual(values = c("3/3" = "blue", "4/4" = "red")) +
	  annotate('text', x = 1.5, y = 3, label = paste('p-value = ',wilcox.test(subRes[subRes$Type == '3/3',]$Exp, subRes[subRes$Type == '4/4',]$Exp, paired = TRUE)[3], sep='')) +
	  labs(x = "Type", y = paste('log10(TPM_',colnames(Res)[i],' + 1)',sep=''), color = "Type") +
	  theme_minimal()
	ggsave(paste('Expression_dosage_20260109/C1_3of3_vs_4of4/Expression_dosage_C1_3of3_vs_4of4_', colnames(Res)[i],'.pdf',sep=''))
	}

# C1, average
Pairs <- read.table('Expression_dosage_20260109/Reciprocal_best_hist_Type6_set3_and_set4_genes_20260109.txt',head=F,sep='\t')
Res <- c()
for(i in 1:nrow(Pairs)){
	Res <- rbind(Res, Sum_exp(dat, allelic, Pairs[i,1], '3/3', i))
	Res <- rbind(Res, Sum_exp(dat, allelic, Pairs[i,2], '4/4', i))
	}
Res[Res[,39]=='3/3',1:38] <- as.numeric(Res[Res[,39]=='3/3',1:38])/3
Res[Res[,39]=='4/4',1:38] <- as.numeric(Res[Res[,39]=='4/4',1:38])/4
for(i in 1:38){
	subRes <- data.frame('Exp'=log10(as.numeric(Res[,i])+1),'Type' = Res[,'Type'], 'Pair' = Res[,'Pair'])
	ggplot(subRes, aes(x = Type, y = Exp)) +
	  geom_line(aes(group = Pair), color = "gray", alpha = 0.5) +  # 连线
	  geom_point(aes(color = Type), size = 1) +                  # 散点
	  geom_boxplot(alpha = 0.5, color = c('blue', 'red'), fill = c('blue', 'red'), width=0.2, outlier.shape = NA)+
	  scale_color_manual(values = c("3/3" = "blue", "4/4" = "red")) +
	  annotate('text', x = 1.5, y = 3, label = paste('p-value = ',wilcox.test(subRes[subRes$Type == '3/3',]$Exp, subRes[subRes$Type == '4/4',]$Exp, paired = TRUE)[3], sep='')) +
	  labs(x = "Type", y = paste('log10(TPM_',colnames(Res)[i],' + 1)',sep=''), color = "Type") +
	  theme_minimal()
	ggsave(paste('Expression_dosage_20260109/C1_3of3_vs_4of4_average/Expression_dosage_C1_3of3_vs_4of4_average_', colnames(Res)[i],'.pdf',sep=''))
	}


# C2
Pairs <- read.table('Expression_dosage_20260109/Reciprocal_best_hist_2_of_3_Type2_vs_3n_Type6_20260109.txt',head=F,sep='\t')
Res <- c()
for(i in 1:nrow(Pairs)){
	Res <- rbind(Res, Sum_exp(dat, allelic, Pairs[i,1], '2/3', i))
	Res <- rbind(Res, Sum_exp(dat, allelic, Pairs[i,2], '3/3', i))
	}
for(i in 1:38){
	subRes <- data.frame('Exp'=log10(as.numeric(Res[,i])+1),'Type' = Res[,'Type'], 'Pair' = Res[,'Pair'])
	ggplot(subRes, aes(x = Type, y = Exp)) +
	  geom_line(aes(group = Pair), color = "gray", alpha = 0.5) +  # 连线
	  geom_point(aes(color = Type), size = 1) +                  # 散点
	  geom_boxplot(alpha = 0.5, color = c('blue', 'red'), fill = c('blue', 'red'), width=0.2, outlier.shape = NA)+
	  scale_color_manual(values = c("2/3" = "blue", "3/3" = "red")) +
	  annotate('text', x = 1.5, y = 3, label = paste('p-value = ',wilcox.test(subRes[subRes$Type == '2/3',]$Exp, subRes[subRes$Type == '3/3',]$Exp, paired = TRUE)[3], sep='')) +
	  labs(x = "Type", y = paste('log10(TPM_',colnames(Res)[i],' + 1)',sep=''), color = "Type") +
	  theme_minimal()
	ggsave(paste('Expression_dosage_20260109/C2_2of3_vs_3of3/Expression_dosage_C2_2of3_vs_3of3_', colnames(Res)[i],'.pdf',sep=''))
	}


# C3
Pairs <- read.table('Expression_dosage_20260109/Reciprocal_best_hist_2_of_4_Type2_vs_3n_Type6_20260109.txt',head=F,sep='\t')
Res <- c()
for(i in 1:nrow(Pairs)){
	Res <- rbind(Res, Sum_exp(dat, allelic, Pairs[i,1], '2/4', i))
	Res <- rbind(Res, Sum_exp(dat, allelic, Pairs[i,2], '3/3', i))
	}
for(i in 1:38){
	subRes <- data.frame('Exp'=log10(as.numeric(Res[,i])+1),'Type' = Res[,'Type'], 'Pair' = Res[,'Pair'])
	ggplot(subRes, aes(x = Type, y = Exp)) +
	  geom_line(aes(group = Pair), color = "gray", alpha = 0.5) +  # 连线
	  geom_point(aes(color = Type), size = 1) +                  # 散点
	  geom_boxplot(alpha = 0.5, color = c('blue', 'red'), fill = c('blue', 'red'), width=0.2, outlier.shape = NA)+
	  scale_color_manual(values = c("2/4" = "blue", "3/3" = "red")) +
	  annotate('text', x = 1.5, y = 3, label = paste('p-value = ',wilcox.test(subRes[subRes$Type == '2/4',]$Exp, subRes[subRes$Type == '3/3',]$Exp, paired = TRUE)[3], sep='')) +
	  labs(x = "Type", y = paste('log10(TPM_',colnames(Res)[i],' + 1)',sep=''), color = "Type") +
	  theme_minimal()
	ggsave(paste('Expression_dosage_20260109/C3_2of4_vs_3of3/Expression_dosage_C3_2of4_vs_3of3_', colnames(Res)[i],'.pdf',sep=''))
	}

# C4
Pairs <- read.table('Expression_dosage_20260109/Reciprocal_best_hist_2_of_3_Type2_vs_3_of_4_Type2_20260109.txt',head=F,sep='\t')
Res <- c()
for(i in 1:nrow(Pairs)){
	Res <- rbind(Res, Sum_exp(dat, allelic, Pairs[i,1], '2/3', i))
	Res <- rbind(Res, Sum_exp(dat, allelic, Pairs[i,2], '3/4', i))
	}
for(i in 1:38){
	subRes <- data.frame('Exp'=log10(as.numeric(Res[,i])+1),'Type' = Res[,'Type'], 'Pair' = Res[,'Pair'])
	ggplot(subRes, aes(x = Type, y = Exp)) +
	  geom_line(aes(group = Pair), color = "gray", alpha = 0.5) +  # 连线
	  geom_point(aes(color = Type), size = 1) +                  # 散点
	  geom_boxplot(alpha = 0.5, color = c('blue', 'red'), fill = c('blue', 'red'), width=0.2, outlier.shape = NA)+
	  scale_color_manual(values = c("2/3" = "blue", "3/4" = "red")) +
	  annotate('text', x = 1.5, y = 3, label = paste('p-value = ',wilcox.test(subRes[subRes$Type == '2/3',]$Exp, subRes[subRes$Type == '3/4',]$Exp, paired = TRUE)[3], sep='')) +
	  labs(x = "Type", y = paste('log10(TPM_',colnames(Res)[i],' + 1)',sep=''), color = "Type") +
	  theme_minimal()
	ggsave(paste('Expression_dosage_20260109/C4_2of3_vs_3of4/Expression_dosage_C4_2of3_vs_3of4_', colnames(Res)[i],'.pdf',sep=''))
	}

# C5
Pairs <- read.table('Expression_dosage_20260109/Reciprocal_best_hist_2_of_4_Type2_vs_3_of_4_Type2_20260109.txt',head=F,sep='\t')
Res <- c()
for(i in 1:nrow(Pairs)){
	Res <- rbind(Res, Sum_exp(dat, allelic, Pairs[i,1], '2/4', i))
	Res <- rbind(Res, Sum_exp(dat, allelic, Pairs[i,2], '3/4', i))
	}
for(i in 1:38){
	subRes <- data.frame('Exp'=log10(as.numeric(Res[,i])+1),'Type' = Res[,'Type'], 'Pair' = Res[,'Pair'])
	ggplot(subRes, aes(x = Type, y = Exp)) +
	  geom_line(aes(group = Pair), color = "gray", alpha = 0.5) +  # 连线
	  geom_point(aes(color = Type), size = 1) +                  # 散点
	  geom_boxplot(alpha = 0.5, color = c('blue', 'red'), fill = c('blue', 'red'), width=0.2, outlier.shape = NA)+
	  scale_color_manual(values = c("2/4" = "blue", "3/4" = "red")) +
	  annotate('text', x = 1.5, y = 3, label = paste('p-value = ',wilcox.test(subRes[subRes$Type == '2/4',]$Exp, subRes[subRes$Type == '3/4',]$Exp, paired = TRUE)[3], sep='')) +
	  labs(x = "Type", y = paste('log10(TPM_',colnames(Res)[i],' + 1)',sep=''), color = "Type") +
	  theme_minimal()
	ggsave(paste('Expression_dosage_20260109/C5_2of4_vs_3of4/Expression_dosage_C5_2of4_vs_3of4_', colnames(Res)[i],'.pdf',sep=''))
	}

# C6
Pairs <- read.table('Expression_dosage_20260109/Reciprocal_best_hist_2_of_3_Type2_vs_4n_Type6_20260109.txt',head=F,sep='\t')
Res <- c()
for(i in 1:nrow(Pairs)){
	Res <- rbind(Res, Sum_exp(dat, allelic, Pairs[i,1], '2/3', i))
	Res <- rbind(Res, Sum_exp(dat, allelic, Pairs[i,2], '4/4', i))
	}
for(i in 1:38){
	subRes <- data.frame('Exp'=log10(as.numeric(Res[,i])+1),'Type' = Res[,'Type'], 'Pair' = Res[,'Pair'])
	ggplot(subRes, aes(x = Type, y = Exp)) +
	  geom_line(aes(group = Pair), color = "gray", alpha = 0.5) +  # 连线
	  geom_point(aes(color = Type), size = 1) +                  # 散点
	  geom_boxplot(alpha = 0.5, color = c('blue', 'red'), fill = c('blue', 'red'), width=0.2, outlier.shape = NA)+
	  scale_color_manual(values = c("2/3" = "blue", "4/4" = "red")) +
	  annotate('text', x = 1.5, y = 3, label = paste('p-value = ',wilcox.test(subRes[subRes$Type == '2/3',]$Exp, subRes[subRes$Type == '4/4',]$Exp, paired = TRUE)[3], sep='')) +
	  labs(x = "Type", y = paste('log10(TPM_',colnames(Res)[i],' + 1)',sep=''), color = "Type") +
	  theme_minimal()
	ggsave(paste('Expression_dosage_20260109/C6_2of3_vs_4of4/Expression_dosage_C6_2of3_vs_4of4_', colnames(Res)[i],'.pdf',sep=''))
	}

# C7
Pairs <- read.table('Expression_dosage_20260109/Reciprocal_best_hist_2_of_4_Type2_vs_4n_Type6_20260109.txt',head=F,sep='\t')
Res <- c()
for(i in 1:nrow(Pairs)){
	Res <- rbind(Res, Sum_exp(dat, allelic, Pairs[i,1], '2/4', i))
	Res <- rbind(Res, Sum_exp(dat, allelic, Pairs[i,2], '4/4', i))
	}
for(i in 1:38){
	subRes <- data.frame('Exp'=log10(as.numeric(Res[,i])+1),'Type' = Res[,'Type'], 'Pair' = Res[,'Pair'])
	ggplot(subRes, aes(x = Type, y = Exp)) +
	  geom_line(aes(group = Pair), color = "gray", alpha = 0.5) +  # 连线
	  geom_point(aes(color = Type), size = 1) +                  # 散点
	  geom_boxplot(alpha = 0.5, color = c('blue', 'red'), fill = c('blue', 'red'), width=0.2, outlier.shape = NA)+
	  scale_color_manual(values = c("2/4" = "blue", "4/4" = "red")) +
	  annotate('text', x = 1.5, y = 3, label = paste('p-value = ',wilcox.test(subRes[subRes$Type == '2/4',]$Exp, subRes[subRes$Type == '4/4',]$Exp, paired = TRUE)[3], sep='')) +
	  labs(x = "Type", y = paste('log10(TPM_',colnames(Res)[i],' + 1)',sep=''), color = "Type") +
	  theme_minimal()
	ggsave(paste('Expression_dosage_20260109/C7_2of4_vs_4of4/Expression_dosage_C7_2of4_vs_4of4_', colnames(Res)[i],'.pdf',sep=''))
	}

# C8
Pairs <- read.table('Expression_dosage_20260109/Reciprocal_best_hist_3_of_4_Type2_vs_4n_Type6_20260109.txt',head=F,sep='\t')
Res <- c()
for(i in 1:nrow(Pairs)){
	Res <- rbind(Res, Sum_exp(dat, allelic, Pairs[i,1], '3/4', i))
	Res <- rbind(Res, Sum_exp(dat, allelic, Pairs[i,2], '4/4', i))
	}
for(i in 1:38){
	subRes <- data.frame('Exp'=log10(as.numeric(Res[,i])+1),'Type' = Res[,'Type'], 'Pair' = Res[,'Pair'])
	ggplot(subRes, aes(x = Type, y = Exp)) +
	  geom_line(aes(group = Pair), color = "gray", alpha = 0.5) +  # 连线
	  geom_point(aes(color = Type), size = 1) +                  # 散点
	  geom_boxplot(alpha = 0.5, color = c('blue', 'red'), fill = c('blue', 'red'), width=0.2, outlier.shape = NA)+
	  scale_color_manual(values = c("3/4" = "blue", "4/4" = "red")) +
	  annotate('text', x = 1.5, y = 3, label = paste('p-value = ',wilcox.test(subRes[subRes$Type == '3/4',]$Exp, subRes[subRes$Type == '4/4',]$Exp, paired = TRUE)[3], sep='')) +
	  labs(x = "Type", y = paste('log10(TPM_',colnames(Res)[i],' + 1)',sep=''), color = "Type") +
	  theme_minimal()
	ggsave(paste('Expression_dosage_20260109/C8_3of4_vs_4of4/Expression_dosage_C8_3of4_vs_4of4_', colnames(Res)[i],'.pdf',sep=''))
	}

# C18, average
Pairs <- read.table('Expression_dosage_20260109/Reciprocal_best_hist_3_of_4_Type2_vs_4n_Type6_20260109.txt',head=F,sep='\t')
Res <- c()
for(i in 1:nrow(Pairs)){
	Res <- rbind(Res, Sum_exp(dat, allelic, Pairs[i,1], '3/4', i))
	Res <- rbind(Res, Sum_exp(dat, allelic, Pairs[i,2], '4/4', i))
	}
Res[Res[,39]=='3/4',1:38] <- as.numeric(Res[Res[,39]=='3/4',1:38])/3
Res[Res[,39]=='4/4',1:38] <- as.numeric(Res[Res[,39]=='4/4',1:38])/4
for(i in 1:38){
	subRes <- data.frame('Exp'=log10(as.numeric(Res[,i])+1),'Type' = Res[,'Type'], 'Pair' = Res[,'Pair'])
	ggplot(subRes, aes(x = Type, y = Exp)) +
	  geom_line(aes(group = Pair), color = "gray", alpha = 0.5) +  # 连线
	  geom_point(aes(color = Type), size = 1) +                  # 散点
	  geom_boxplot(alpha = 0.5, color = c('blue', 'red'), fill = c('blue', 'red'), width=0.2, outlier.shape = NA)+
	  scale_color_manual(values = c("3/4" = "blue", "4/4" = "red")) +
	  annotate('text', x = 1.5, y = 3, label = paste('p-value = ',wilcox.test(subRes[subRes$Type == '3/4',]$Exp, subRes[subRes$Type == '4/4',]$Exp, paired = TRUE)[3], sep='')) +
	  labs(x = "Type", y = paste('log10(TPM_',colnames(Res)[i],' + 1)',sep=''), color = "Type") +
	  theme_minimal()
	ggsave(paste('Expression_dosage_20260109/C8_3of4_vs_4of4_average/Expression_dosage_C8_3of4_vs_4of4_average_', colnames(Res)[i],'.pdf',sep=''))
	}






setwd('E:\\Manuscripts\\2026_RNA-seq_mapping\\Figures\\Simulation')
# TPM calculation
genome <- read.table("SRR3664372_t_adjustranscriptome_Tophat2_HTSeqCount_fraction.out", header=TRUE, sep="\t", stringsAsFactors=FALSE)
colnames(genome)[1] <- "geneID"
colnames(genome)[2] <- "Counts_genome"
genome <- genome[1:(nrow(genome)-5), ]
gene_length <- read.csv("Araport11_longest_transcript_length.tsv", stringsAsFactors=FALSE,sep='\t')
merged <- merge(genome, gene_length[,1:2], by.x="geneID", by.y="gene_id")
merged$RPK <- merged$Counts_genome / (merged$length / 1000)
total_RPK <- sum(merged$RPK)
merged$TPM <- merged$RPK / total_RPK * 1e6
write.csv(merged, "SRR3664372_t_adjustranscriptome_Tophat2_HTSeqCount_fraction_TPM.csv", row.names=FALSE, quote = FALSE)


setwd("your work directory")

index <- "SRR3664372"

# read tandem gene duplicates list
tandem_file1 <- "at.tandem.pairs"
tandem_file2 <- "at.proximal.pairs"
file1 <- read.table(tandem_file1, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
file2 <- read.table(tandem_file2, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
tandem_genes <- unique(c(unique(c(file1[[1]], file1[[3]])),
                         unique(c(file2[[1]], file2[[3]]))))

# read G-guided TPM
tpm_G <- read.table(paste0(index, "_genome_Tophat2_TPM.csv"),
                      header = TRUE, sep = "\t", stringsAsFactors = FALSE)
tpm_G <- tpm_G[, c("geneID", "TPM")]
colnames(tpm_G) <- c("Gene_ID", "TPM_G")

# read T-guided TPM
tpm_T <- read.table(paste0(index, "_transcriptome_Tophat2_TPM.csv"),
                         header = TRUE, sep = "\t", stringsAsFactors = FALSE)
tpm_T <- tpm_T[, c("geneID", "TPM")]
colnames(tpm_T) <- c("Gene_ID", "TPM_T")

# combine
expr <- merge(tpm_G, tpm_T, by = "Gene_ID")
expr <- expr[-c(1:5), ]

# log10(TPM + 1)
expr$log_G    <- log10(as.numeric(expr$tpm_G) + 1)
expr$log_T <- log10(as.numeric(expr$tpm_T) + 1)

# Z-score outlier
expr$distance <- abs(expr$log_G - expr$log_T) / sqrt(2)
expr$z_score <- (expr$distance - mean(expr$distance)) / sd(expr$distance)
expr$is_outlier <- expr$z_score > 3

# mark tandem
expr$is_tandem <- sub("\\..*$", "", expr$Gene_ID) %in% tandem_genes

# ====== 4. Fisher's exact test ======
table_data <- table(
  outlier = factor(expr$is_outlier, levels = c(FALSE, TRUE)),
  tandem  = factor(expr$is_tandem,  levels = c(FALSE, TRUE))
)
cat("\n===== 2x2 matrix (outlier vs Tandem) =====\n")
print(table_data)
print(fisher.test(table_data))

expr$direction <- ifelse(expr$log_G < expr$log_T, "up", "down")

expr$up_outlier <- expr$is_outlier & expr$direction == "up"
table_up <- table(
  up_outlier = factor(expr$up_outlier, levels = c(FALSE, TRUE)),
  tandem     = factor(expr$is_tandem,  levels = c(FALSE, TRUE))
)
cat("\n===== upper outlier (T > G) vs Tandem =====\n")
print(table_up)
print(fisher.test(table_up))

expr$down_outlier <- expr$is_outlier & expr$direction == "down"
table_down <- table(
  down_outlier = factor(expr$down_outlier, levels = c(FALSE, TRUE)),
  tandem       = factor(expr$is_tandem,    levels = c(FALSE, TRUE))
)
cat("\n===== lower outlier (T < G) vs Tandem =====\n")
print(table_down)
print(fisher.test(table_down))

outliers <- expr[expr$is_outlier,
                 c("Gene_ID", "TPM_G", "TPM_T", "direction", "is_tandem")]
write.csv(outliers,
          file = paste0(index, "_Tophat2_Transcriptome_G_vs_T_outliers.csv"),
          row.names = FALSE, fileEncoding = "UTF-8", quote = FALSE)

cat("\n===== Done =====\n")

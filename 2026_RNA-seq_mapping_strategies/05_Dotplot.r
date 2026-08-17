# ============================================================================
# SRR3664372 转录组 TPM：nonunique=all vs nonunique=random 点图
# 高亮 AT1G21770 和 AT1G07520
# ============================================================================

# ====== 1. 计算 random 策略 transcriptome TPM ======
setwd("F:\\文献\\RNAseq_for_Tophat2_mapping\\Arabidopsis\\PRJNA323955")

gene_length <- read.csv("G:\\鉴定拟南芥串联重复基因\\1_拟南芥基因组文件\\Araport11_gene_length.csv",
                        stringsAsFactors = FALSE)
colnames(gene_length)[1] <- "geneID"

index <- "SRR3664372"

count_file <- paste0(index, "_t_adjustranscriptome_Tophat2_HTSeqCount_Random.out")
counts <- read.table(count_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
colnames(counts)[1] <- "geneID"
colnames(counts)[2] <- "Counts_transcriptome"
counts <- counts[1:(nrow(counts) - 5), ]
merged <- merge(counts, gene_length, by = "geneID")
merged$RPK <- merged$Counts_transcriptome / (merged$transcript_length / 1000)
merged$TPM <- merged$RPK / sum(merged$RPK) * 1e6
write.csv(merged[, c("geneID", "Counts_transcriptome", "transcript_length", "TPM")],
          paste0(index, "_transcriptome_Tophat2_TPM_Random.csv"),
          row.names = FALSE, quote = FALSE)

# ====== 2. 读取 all 和 random 转录组 TPM ======

# 读取 tandem 基因列表
tandem_file1 <- "G:\\鉴定拟南芥串联重复基因\\4_Differential_AND_Enrichment_Analysis\\MCScanX\\at.tandem.pairs"
tandem_file2 <- "G:\\鉴定拟南芥串联重复基因\\4_Differential_AND_Enrichment_Analysis\\MCScanX\\at.proximal.pairs"
file1 <- read.table(tandem_file1, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
file2 <- read.table(tandem_file2, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
tandem_genes <- unique(c(unique(c(file1[[1]], file1[[3]])),
                         unique(c(file2[[1]], file2[[3]]))))

# 读取转录组 TPM — all（在 PRJNA323955_ReadCounts 中）
tpm_all <- read.table(paste0(index, "_transcriptome_Tophat2_TPM.csv"),
                      header = TRUE, sep = ",", stringsAsFactors = FALSE)
tpm_all <- tpm_all[, c("geneID", "TPM")]
colnames(tpm_all) <- c("Gene_ID", "TPM_all")

# 读取转录组 TPM — random
tpm_random <- read.table(paste0(index, "_transcriptome_Tophat2_TPM_Random.csv"),
                         header = TRUE, sep = ",", stringsAsFactors = FALSE)
tpm_random <- tpm_random[, c("geneID", "TPM")]
colnames(tpm_random) <- c("Gene_ID", "TPM_random")

# 合并
expr <- merge(tpm_all, tpm_random, by = "Gene_ID")
expr <- expr[-c(1:5), ]

# log10(TPM + 1)
expr$log_all    <- log10(as.numeric(expr$TPM_all) + 1)
expr$log_random <- log10(as.numeric(expr$TPM_random) + 1)

# Z-score outlier
expr$distance <- abs(expr$log_all - expr$log_random) / sqrt(2)
expr$z_score <- (expr$distance - mean(expr$distance)) / sd(expr$distance)
expr$is_outlier <- expr$z_score > 3

# 标记 tandem
expr$is_tandem <- sub("\\..*$", "", expr$Gene_ID) %in% tandem_genes

# 标记高亮基因（去掉 .数字 后缀匹配）
highlight_genes <- c("AT1G21770", "AT1G07520")
expr$highlight <- sub("\\..*$", "", expr$Gene_ID) %in% highlight_genes

# ---- 统计百分比 ----
total_genes <- nrow(expr)
n_outlier  <- sum(expr$is_outlier)
n_up       <- sum(expr$is_outlier & expr$log_random > expr$log_all)   # random > all
n_down     <- sum(expr$is_outlier & expr$log_random < expr$log_all)   # random < all

cat(sprintf("\n===== 总基因数: %d =====\n", total_genes))
cat(sprintf("Outlier 总数: %d (%.2f%%)\n", n_outlier, 100 * n_outlier / total_genes))
cat(sprintf("上方 outlier (random > all): %d (占所有基因 %.2f%%)\n",
            n_up, 100 * n_up / total_genes))
cat(sprintf("下方 outlier (random < all): %d (占所有基因 %.2f%%)\n",
            n_down, 100 * n_down / total_genes))

# 打印高亮基因的 TPM 值
cat("\n===== 高亮基因 TPM =====\n")
print(expr[expr$highlight, c("Gene_ID", "TPM_all", "TPM_random")])

# ====== 3. 绘图 ======
library(ggplot2)

# 定义高亮颜色
col_AT1G21770 <- "#E69F00"  # 橙色
col_AT1G07520 <- "#009E73"  # 绿色

p <- ggplot(expr, aes(x = log_all, y = log_random)) +
  # 正常点（黑色）
  geom_point(data = expr[!expr$is_outlier & !expr$highlight, ],
             alpha = 0.5, color = "black", size = 0.8) +
  # outlier（蓝色）
  geom_point(data = expr[expr$is_outlier & !expr$highlight, ],
             alpha = 0.8, color = "#293890", size = 1.2) +
  # 高亮基因 AT1G21770
  geom_point(data = expr[expr$highlight & sub("\\..*$", "", expr$Gene_ID) == "AT1G21770", ],
             color = col_AT1G21770, size = 3, shape = 17) +
  # 高亮基因 AT1G07520
  geom_point(data = expr[expr$highlight & sub("\\..*$", "", expr$Gene_ID) == "AT1G07520", ],
             color = col_AT1G07520, size = 3, shape = 17) +
  # Y = X 参考线
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "#293890", linewidth = 0.5) +
  labs(
    x = "Transcriptome TPM log(TPM+1)  —  nonunique=all",
    y = "Transcriptome TPM log(TPM+1)  —  nonunique=random",
    title = "Tophat2 Transcriptome: nonunique=all vs random",
    subtitle = paste("Arabidopsis thaliana | Sample:", index)
  ) +
  scale_x_continuous(limits = c(0, 5), expand = c(0, 0),
                     breaks = seq(0, 5, 1),
                     labels = function(x) ifelse(x == 0, "0", x)) +
  scale_y_continuous(limits = c(0, 5), expand = c(0, 0),
                     breaks = seq(0, 5, 1),
                     labels = function(x) ifelse(x == 0, "", x)) +
  # 手动添加图例
  annotate("point", x = 3.8, y = 0.4, color = col_AT1G21770, size = 3, shape = 17) +
  annotate("text",  x = 3.9, y = 0.4, label = "AT1G21770", hjust = 0, size = 3.5) +
  annotate("point", x = 3.8, y = 0.15, color = col_AT1G07520, size = 3, shape = 17) +
  annotate("text",  x = 3.9, y = 0.15, label = "AT1G07520", hjust = 0, size = 3.5) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 13, hjust = 0.5),
    axis.title = element_text(size = 18, face = "bold"),
    axis.text = element_text(size = 12),
    axis.line = element_blank(),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16)
  ) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 1.2) +
  geom_vline(xintercept = 0, linetype = "solid", color = "black", linewidth = 1.2) +
  coord_fixed()

print(p)

ggsave(paste0(index, "_Tophat2_Transcriptome_all_vs_random_highlight.pdf"),
       p, width = 8, height = 8, dpi = 300, bg = "white")

# ====== 4. Fisher 检验 ======
table_data <- table(
  outlier = factor(expr$is_outlier, levels = c(FALSE, TRUE)),
  tandem  = factor(expr$is_tandem,  levels = c(FALSE, TRUE))
)
cat("\n===== 2x2 列联表 (outlier vs Tandem) =====\n")
print(table_data)
print(fisher.test(table_data))

expr$direction <- ifelse(expr$log_all < expr$log_random, "up", "down")

expr$up_outlier <- expr$is_outlier & expr$direction == "up"
table_up <- table(
  up_outlier = factor(expr$up_outlier, levels = c(FALSE, TRUE)),
  tandem     = factor(expr$is_tandem,  levels = c(FALSE, TRUE))
)
cat("\n===== 上方 outlier (random > all) vs Tandem =====\n")
print(table_up)
print(fisher.test(table_up))

expr$down_outlier <- expr$is_outlier & expr$direction == "down"
table_down <- table(
  down_outlier = factor(expr$down_outlier, levels = c(FALSE, TRUE)),
  tandem       = factor(expr$is_tandem,    levels = c(FALSE, TRUE))
)
cat("\n===== 下方 outlier (random < all) vs Tandem =====\n")
print(table_down)
print(fisher.test(table_down))

outliers <- expr[expr$is_outlier,
                 c("Gene_ID", "TPM_all", "TPM_random", "direction", "is_tandem")]
write.csv(outliers,
          file = paste0(index, "_Tophat2_Transcriptome_all_vs_random_outliers.csv"),
          row.names = FALSE, fileEncoding = "UTF-8", quote = FALSE)

cat("\n===== 全部完成 =====\n")

# ====SRR3664372 T-guided Tophat2, --nonunique random vs --nonunique all====
setwd("your work directory")
index <- "SRR3664372"
# read TPM files of T-guided Tophat2, under --nonunique all
tpm_all <- read.table(paste0(index, "_transcriptome_Tophat2_TPM.csv"),
                      header = TRUE, sep = "\t", stringsAsFactors = FALSE)
tpm_all <- tpm_all[, c("geneID", "TPM")]
colnames(tpm_all) <- c("Gene_ID", "TPM_all")

# read TPM files of T-guided Tophat2, under --nonunique random
tpm_random <- read.table(paste0(index, "_transcriptome_Tophat2_TPM_Random.csv"),
                         header = TRUE, sep = "\t", stringsAsFactors = FALSE)
tpm_random <- tpm_random[, c("geneID", "TPM")]
colnames(tpm_random) <- c("Gene_ID", "TPM_random")

expr <- merge(tpm_all, tpm_random, by = "Gene_ID")
expr <- expr[-c(1:5), ]

# log10(TPM + 1)
expr$log_all    <- log10(as.numeric(expr$TPM_all) + 1)
expr$log_random <- log10(as.numeric(expr$TPM_random) + 1)

# Z-score outlier
expr$distance <- abs(expr$log_all - expr$log_random) / sqrt(2)
expr$z_score <- (expr$distance - mean(expr$distance)) / sd(expr$distance)
expr$is_outlier <- expr$z_score > 3

# highlight AT1G21770 and AT1G07520
highlight_genes <- c("AT1G21770", "AT1G07520")
expr$highlight <- sub("\\..*$", "", expr$Gene_ID) %in% highlight_genes

# ---- count outlier genes and calculate percentage ----
total_genes <- nrow(expr)
n_outlier  <- sum(expr$is_outlier)
n_up       <- sum(expr$is_outlier & expr$log_random > expr$log_all)   # random > all
n_down     <- sum(expr$is_outlier & expr$log_random < expr$log_all)   # random < all

cat(sprintf("\n===== Total gene numbers: %d =====\n", total_genes))
cat(sprintf("Outlier numbers: %d (%.2f%%)\n", n_outlier, 100 * n_outlier / total_genes))
cat(sprintf("upper outlier (random > all): %d (%.2f%%)\n",
            n_up, 100 * n_up / total_genes))
cat(sprintf("lower outlier (random < all): %d (%.2f%%)\n",
            n_down, 100 * n_down / total_genes))

cat("\n===== highlighted gene TPM =====\n")
print(expr[expr$highlight, c("Gene_ID", "TPM_all", "TPM_random")])

# ====== 3. plot ======
library(ggplot2)

# define colors for AT1G21770 and AT1G07520
col_AT1G21770 <- "#E69F00"  # orange
col_AT1G07520 <- "#009E73"  # green

p <- ggplot(expr, aes(x = log_all, y = log_random)) +
  # normal dots（black）
  geom_point(data = expr[!expr$is_outlier & !expr$highlight, ],
             alpha = 0.5, color = "black", size = 0.8) +
  # outlier（blue）
  geom_point(data = expr[expr$is_outlier & !expr$highlight, ],
             alpha = 0.8, color = "#293890", size = 1.2) +
  #  AT1G21770
  geom_point(data = expr[expr$highlight & sub("\\..*$", "", expr$Gene_ID) == "AT1G21770", ],
             color = col_AT1G21770, size = 3, shape = 17) +
  #  AT1G07520
  geom_point(data = expr[expr$highlight & sub("\\..*$", "", expr$Gene_ID) == "AT1G07520", ],
             color = col_AT1G07520, size = 3, shape = 17) +
  # Y = X diagonal
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
  # add legend manually
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


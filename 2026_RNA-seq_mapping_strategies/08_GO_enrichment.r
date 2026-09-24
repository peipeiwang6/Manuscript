library(clusterProfiler)
library(org.Athaliana.eg.db)
library(ggplot2)


## ---------------------------------------------------------------------------
## 1. Settings
## ---------------------------------------------------------------------------

dataset     <- "E-MTAB-4151"            # dataset name, used in the plot title
prefix      <- "Inoculation_VS_WT"      # file name prefix of the DE analysis
direction   <- "Downregulated"          # "Upregulated" or "Downregulated"
orgdb       <- org.Athaliana.eg.db      # annotation package of the species
padj_cutoff <- 0.05                     # significance threshold
top_shared  <- 5                        # shared terms drawn above the "omitted" row


## ---------------------------------------------------------------------------
## 2. GO enrichment for both strategies
## ---------------------------------------------------------------------------

# gene IDs of the DEGs: the row names of the DE result table. Some annotation
# versions append a transcript suffix to the gene ID so that it no longer matches
# the OrgDb keys; that suffix is removed here.
clean_ids <- function(x) {
  x <- sub("^gene:", "", x)
  sub("\\.ITAG5\\.0$", "", x)
}

genes_genome <- clean_ids(rownames(read.csv(paste0(prefix, "_genome_", direction, "_genes2.csv"),
                                           row.names = 1)))
genes_trans  <- clean_ids(rownames(read.csv(paste0(prefix, "_transcriptome_", direction, "_genes2.csv"),
                                           row.names = 1)))

cat("DEGs:", length(genes_genome), "(genome) /", length(genes_trans), "(transcriptome)\n")

# pvalueCutoff = 1: keep all tested terms, significance is applied later
go_genome <- enrichGO(gene          = genes_genome,
                      OrgDb         = orgdb,
                      keyType       = "GID",
                      ont           = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 1,
                      qvalueCutoff  = 1,
                      minGSSize     = 1)

go_trans  <- enrichGO(gene          = genes_trans,
                      OrgDb         = orgdb,
                      keyType       = "GID",
                      ont           = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 1,
                      qvalueCutoff  = 1,
                      minGSSize     = 1)

go_genome <- as.data.frame(go_genome)
go_trans  <- as.data.frame(go_trans)

# odds ratio of the enrichment test: (a/b) / (c/d) with
# a = DEGs in the term, b = DEGs not in the term,
# c = background genes in the term, d = background genes not in the term
odds_ratio <- function(gene_ratio, bg_ratio) {
  a <- as.numeric(strsplit(gene_ratio, "/")[[1]])
  c <- as.numeric(strsplit(bg_ratio,   "/")[[1]])
  (a[1] / (a[2] - a[1])) / (c[1] / (c[2] - c[1]))
}

go_genome$odds_ratio <- mapply(odds_ratio, go_genome$GeneRatio, go_genome$BgRatio)
go_trans$odds_ratio  <- mapply(odds_ratio, go_trans$GeneRatio,  go_trans$BgRatio)


## ---------------------------------------------------------------------------
## 3. Shared and strategy-specific terms
## ---------------------------------------------------------------------------

sig_genome <- go_genome$ID[go_genome$p.adjust < padj_cutoff]
sig_trans  <- go_trans$ID[go_trans$p.adjust < padj_cutoff]

shared       <- intersect(sig_genome, sig_trans)
genome_only  <- setdiff(sig_genome, sig_trans)
trans_only   <- setdiff(sig_trans, sig_genome)
union_terms  <- union(sig_genome, sig_trans)

cat("significant GO terms - genome:", length(sig_genome),
    "transcriptome:", length(sig_trans), "\n")
cat("shared:", length(shared),
    "genome only:", length(genome_only),
    "transcriptome only:", length(trans_only), "\n")


## ---------------------------------------------------------------------------
## 4. Table of all union terms with the values of both strategies
## ---------------------------------------------------------------------------

# one helper: take one column of a result table, matched by GO ID
from_table <- function(table, ids, column) table[[column]][match(ids, table$ID)]

union_table <- data.frame(
  ID          = union_terms,
  Description = from_table(go_genome, union_terms, "Description"),
  group       = ifelse(union_terms %in% genome_only, "genome only",
                ifelse(union_terms %in% trans_only,  "transcriptome only", "shared")),
  padj_genome = from_table(go_genome, union_terms, "p.adjust"),
  padj_trans  = from_table(go_trans,  union_terms, "p.adjust"),
  OR_genome   = from_table(go_genome, union_terms, "odds_ratio"),
  OR_trans    = from_table(go_trans,  union_terms, "odds_ratio"),
  stringsAsFactors = FALSE)

# the more significant of the two strategies decides the position in the figure
union_table <- union_table[order(pmin(union_table$padj_genome, union_table$padj_trans,
                                      na.rm = TRUE)), ]
write.csv(union_table, "GO_overlap_down.csv", row.names = FALSE)


## ---------------------------------------------------------------------------
## 5. Rows of the figure
##    top shared terms  ->  "... n GO terms omitted ..."  ->  all specific terms
## ---------------------------------------------------------------------------

shared_drawn <- head(union_table[union_table$group == "shared", ], top_shared)
specific     <- union_table[union_table$group != "shared", ]
n_omitted    <- sum(union_table$group == "shared") - nrow(shared_drawn)

shown <- rbind(shared_drawn, specific)

# labels of the y axis, the omitted row has no data point
wrap <- function(x) vapply(x, function(s) paste(strwrap(s, width = 40), collapse = "\n"),
                           character(1))
rows <- c(wrap(shared_drawn$Description),
          if (n_omitted > 0) paste0("... ", n_omitted, " GO terms omitted ...") else NULL,
          wrap(specific$Description))
shown$label <- factor(wrap(shown$Description), levels = rev(rows))

cat("figure: ", nrow(shown), "terms drawn,", n_omitted, "shared terms omitted\n")


## ---------------------------------------------------------------------------
## 6. Figure
## ---------------------------------------------------------------------------

# both panels use the same axis, colour and size range, otherwise the two panels
# could not be compared with each other
shown$log10_genome <- -log10(shown$padj_genome)
shown$log10_trans  <- -log10(shown$padj_trans)
x_max        <- max(shown$log10_genome, shown$log10_trans, na.rm = TRUE) * 1.3
colour_range <- range(c(shown$log10_genome, shown$log10_trans), na.rm = TRUE)
size_range   <- range(c(shown$OR_genome, shown$OR_trans), na.rm = TRUE)

plot_go <- function(dat, log10_col, or_col, title, file) {
  dat$log10padj   <- dat[[log10_col]]
  dat$odds_ratio  <- dat[[or_col]]
  dat$significant <- dat$log10padj > -log10(padj_cutoff)   # p.adjust < 0.05
  nudge           <- max(dat$log10padj, na.rm = TRUE) * 0.05

  p <- ggplot(dat, aes(x = log10padj, y = label)) +
    # threshold line: p.adjust = 0.05
    geom_vline(xintercept = -log10(0.05), linetype = "dashed",
               colour = "grey40", linewidth = 1.5) +
    # not significant in this strategy: grey, but the odds ratio is still shown
    geom_point(data = subset(dat, !significant),
               aes(size = odds_ratio), colour = "grey70", alpha = 0.85) +
    # significant in this strategy: colour by -log10(p.adjust), size by odds ratio
    geom_point(data = subset(dat, significant),
               aes(size = odds_ratio, colour = log10padj), alpha = 0.85) +
    geom_text(aes(label = sprintf("%.2f", odds_ratio)),
              nudge_x = nudge, hjust = 0, size = 3.2, colour = "black") +
    scale_x_continuous(name = "-log10(p.adjust)", limits = c(0, x_max),
                       expand = expansion(mult = c(0, 0.1))) +
    # limits/drop are set so that the "omitted" row stays on the axis
    scale_y_discrete(limits = rev(rows), drop = FALSE) +
    scale_colour_gradient(low = "#FFB6C1", high = "#8B0000",
                          name = "-log10(p.adjust)", limits = colour_range) +
    scale_size_continuous(range = c(3, 8), name = "Odds Ratio", limits = size_range) +
    labs(y = NULL, title = title) +
    theme_bw(base_size = 12) +
    theme(panel.grid = element_blank(),
          plot.title = element_text(hjust = 0.5, face = "bold"))

  ggsave(file, p, width = 12, height = 12, dpi = 300)
}

plot_go(shown, "log10_genome", "OR_genome",
        paste(dataset, "genome", direction, "(p.adjust < 0.05, no simplify)"),
        "GO_overlap_down_genome.pdf")

plot_go(shown, "log10_trans", "OR_trans",
        paste(dataset, "transcriptome", direction, "(p.adjust < 0.05, no simplify)"),
        "GO_overlap_down_transcriptome.pdf")

# BiocManager::install(c("AnnotationForge","org.At.eg.db"))
library(tidyverse)
library(DESeq2)
library(dplyr)
library(clusterProfiler)
library(AnnotationForge)
library(org.At.eg.db)
setwd("E:\\Manuscripts\\2026_RNA-seq_mapping\\Figures\\Test_for_GO_enrichment_scripts")

#==== create the OrgDb GO annotation package for the target species ==== 
# input file contains two columns: gene ID and GO annotation (multiple GO annotation were deliminated with space)
geneID2GO <- read.table(file = 'SL5.0_GO_list_processed_merged_single.map', header = FALSE)
colnames(geneID2GO) <- c("GID", "GO")

# prepare gene matrix
fSym <- transform(geneID2GO, SYMBOL = GID, GENENAME = GID)[, c("GID", "SYMBOL", "GENENAME")]
fSym <- fSym[fSym[, 2] != "-", ]
fSym <- fSym[fSym[, 3] != "-", ]
fSym <- fSym[!duplicated(fSym), ]
dim(fSym)

# prepare GO matrix
fGO <- transform(geneID2GO, EVIDENCE = rep("", length(GID)))
fGO <- fGO[!duplicated(fGO), ]
dim(fGO)

# build the OrgDb GO annotation package
makeOrgPackage(
  gene_info   = fSym,
  go          = fGO,
  version     = "0.2",
  maintainer  = "jyx <873553743@qq.com>",
  author      = "jyx <873553743@qq.com>",
  outputDir   = ".",
  tax_id      = "4081",
  genus       = "Solanum",
  species     = "lycopersicum",
  goTable     = "go"
)

# install the GO annotation package
install.packages("./org.Slycopersicum.eg.db", type = "source", repos = NULL)


#==== GO enrichment for DEGs from T-guided analysis, for G-guided analysis, please modified the script correspondingly==== 
#merge read count files
filelist <- list.files(pattern="*_t_adjustranscriptome_Tophat2_HTSeqCount.out")
data <- read.table(filelist[1],header = T,sep="\t")
for (i in 2:length(filelist)) {
  df <- read.table(filelist[i],header = T,sep="\t")
  data <- merge(data,df,by="gene_id",all=T)
}
sample_name <- gsub("_t_adjustranscriptome_Tophat2_HTSeqCount.out","",filelist)
col_name <- c("gene_id",sample_name)
names(data) <- col_name
write.csv(data,file = "PRJDB8743_transcriptome.csv",quote=F,row.names = F)

transcriptome_counts <- read.csv("PRJDB8743_transcriptome.csv", header = T)
rownames(transcriptome_counts) <- transcriptome_counts[, 1]
transcriptome_counts <- transcriptome_counts[, -1]

transcriptome_control <- transcriptome_counts[, 1:3]  # control
transcriptome_S1 <- transcriptome_counts[, 4:6]      # treatment

transcriptome_combined <- cbind(transcriptome_control, transcriptome_S1)

sample <- colnames(transcriptome_combined)
condition <- factor(c(rep("control", 3), rep("S1", 3)), levels = c("control", "S1"))
colData <- data.frame(sample, condition)

transcriptome_combined_int <- round(transcriptome_combined)

# create DESeqDataSet data
dds <- DESeqDataSetFromMatrix(countData = transcriptome_combined_int, 
                              colData = colData, 
                              design = ~ condition)
dds <- DESeq(dds)

res <- results(dds)
res <- res[order(res$pvalue), ]
head(res)
summary(res)

write.csv(res, file = "PRJDB8743_inoculation_VS_mock_transcriptome_DEG.csv", quote = FALSE)

# filter significantly differentially expression genes (DEGs)
diff_deseq2 <- subset(res, padj < 0.05 & abs(log2FoldChange) > 1)
diff_deseq2 <- data.frame(diff_deseq2)
diff_deseq2 <- diff_deseq2 %>% mutate(type = if_else(log2FoldChange > 0, 'Up', 'Down'))

upregulated_genes <- diff_deseq2 %>% filter(type == 'Up')
downregulated_genes <- diff_deseq2 %>% filter(type == 'Down')

write.csv(upregulated_genes, file = "PRJDB8743_inoculation_VS_mock_transcriptome_Upregulated_genes.csv", quote = FALSE, row.names = TRUE)
write.csv(downregulated_genes, file = "PRJDB8743_inoculation_VS_mock_transcriptome_Downregulated_genes.csv", quote = FALSE, row.names = TRUE)

#GO plot
genes_up_transcriptome <- read.csv("PRJDB8743_inoculation_VS_mock_transcriptome_Upregulated_genes.csv", header = TRUE)
genes_down_transcriptome <- read.csv("PRJDB8743_inoculation_VS_mock_transcriptome_Downregulated_genes.csv", header = TRUE)
geneList_up_transcriptome <- genes_up_transcriptome[, 1]
geneList_down_transcriptome <- genes_down_transcriptome[, 1]

ego_up_transcriptome <- enrichGO(
  gene          = geneList_up_transcriptome,
  OrgDb         = org.Athaliana.eg.db,
  keyType       = "GID",            
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 1,
  minGSSize     = 1
)

ego_up_transcriptome <- clusterProfiler::simplify(
  ego_up_transcriptome,
  cutoff = 0.7,
  by = "p.adjust",
  select_fun = min,
  measure = "Wang"
)

ego_down_transcriptome <- enrichGO(
  gene          = geneList_down_transcriptome,
  OrgDb         = org.Athaliana.eg.db,
  keyType       = "GID",            
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 1,
  minGSSize     = 1
)

ego_down_transcriptome <- clusterProfiler::simplify(
  ego_down_transcriptome,
  cutoff = 0.7,
  by = "p.adjust",
  select_fun = min,
  measure = "Wang"
)

write.table(ego_up_transcriptome, "PRJDB8743_inoculation_VS_mock_transcriptome_upregulated_genes_GO_enrichment.txt", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(ego_down_transcriptome, "PRJDB8743_inoculation_VS_mock_transcriptome_downregulated_genes_GO_enrichment.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# seperate BP, CC, MF GO terms, for upregulated genes, for downregulated genes, please moditify the script correspondingly
ego_up_transcriptome_bp <- ego_up_transcriptome

ego_up_transcriptome_bp <- ego_up_transcriptome_bp %>%
  mutate(
    GeneRatio_num = as.numeric(sub("/\\d+", "", GeneRatio)),
    GeneRatio_den = as.numeric(sub("\\d+/", "", GeneRatio)),
    BgRatio_num = as.numeric(sub("/\\d+", "", BgRatio)),
    BgRatio_den = as.numeric(sub("\\d+/", "", BgRatio)),
    odds_ratio = (GeneRatio_num/(GeneRatio_den - GeneRatio_num)) / 
      (BgRatio_num/(BgRatio_den - BgRatio_num)),
    odds_ratio = ifelse(is.infinite(odds_ratio), NA, odds_ratio)
  ) %>%
  # sorted in descending order of odds ratio
  arrange(desc(odds_ratio)) %>%
  mutate(
    Description = str_wrap(Description, width = 40),
    Description = factor(Description, levels = rev(unique(Description)))
  )

ego_up_transcriptome_final <- ego_up_transcriptome_bp

write.table(ego_up_transcriptome_final, "PRJDB8743_inoculation_VS_mock_transcriptome_Upregulated_genes_GO_enrichment_ordered.txt", sep = "\t", row.names = FALSE, quote = FALSE)

#######################
# compare GO enrichment between G-guided and T-guided analysis, Biology Progress (BP)
threshold_line <- -log10(0.05)

ego_significant <- ego_up_genome_bp %>% filter(significance == "Significant")
ego_not_significant <- ego_up_genome_bp %>% filter(significance == "Not Significant")

# get the union ranges between G-guided and T-guided results
all_neg_log10_padj <- c(ego_up_genome_bp$neg_log10_padj, ego_up_transcriptome_bp$neg_log10_padj)
all_odds_ratio <- c(ego_up_genome_bp$odds_ratio, ego_up_transcriptome_bp$odds_ratio)

color_range <- range(all_neg_log10_padj, na.rm = TRUE)
size_range <- range(all_odds_ratio, na.rm = TRUE)

# first plot
p_up_genome <- ggplot() +
  
  geom_vline(xintercept = threshold_line, linetype = "dashed", color = "gray40", linewidth = 1.5) +
  
  geom_point(
    data = ego_not_significant,
    aes(x = neg_log10_padj, y = Description, size = odds_ratio),
    color = "gray70", alpha = 0.8
  ) +
  
  geom_point(
    data = ego_significant,
    aes(x = neg_log10_padj, y = Description, size = odds_ratio, color = neg_log10_padj),
    alpha = 0.8
  ) +
  
  geom_text(
    data = ego_up_genome_bp,
    aes(x = neg_log10_padj, y = Description, label = sprintf("%.2f", odds_ratio)),
    nudge_x = max(ego_up_genome_bp$neg_log10_padj, na.rm = TRUE) * 0.05,
    hjust = 0, size = 16, color = "black", fontface = "bold", vjust = 0.5
  ) +
  
  scale_x_continuous(
    name = "-log10(p.adjust)",
    limits = c(0, max(ego_up_genome_bp$neg_log10_padj, na.rm = TRUE) * 1.3),
    expand = expansion(mult = c(0, 0.1))
  ) +
  
  scale_color_gradient(
    low = "#FFB6C1",  
    high = "#8B0000", 
    name = "-log10(p.adjust)",
    limits = color_range, 
    guide = guide_colorbar(
      title.position = "top",
      title.hjust = 0.5,
      barwidth = unit(2, "cm"),
      barheight = unit(15, "cm")
    )
  ) +
  
  scale_size_continuous(
    range = c(6, 12), 
    name = "Odds Ratio",
    limits = size_range 
  ) +
  labs(y = NULL, title = "Inoculation VS Control WT genome Downregulated") +
  theme_bw(base_size = 14) +
  theme(
    panel.background = element_rect(fill = "white"),
    legend.key = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_text(size = 42),
    axis.text.x = element_text(size = 42),
    axis.text.y = element_text(
      size = 50,
      color = "black",
      lineheight = 0.9
    ),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 52),
    panel.border = element_rect(color = "black", fill = NA),
    legend.title = element_text(size = 40, face = "bold"),
    legend.text = element_text(size = 36),
    legend.key.size = unit(2, "cm"),
    legend.position = "right"
  )

# second plot - use the same scale range
p_up_transcriptome <- ggplot() +
  
  geom_vline(xintercept = threshold_line, linetype = "dashed", color = "gray40", linewidth = 1.5) +
  
  geom_point(
    data = ego_transcriptome_not_significant,
    aes(x = neg_log10_padj, y = Description, size = odds_ratio),
    color = "gray70", alpha = 0.8
  ) +
  
  geom_point(
    data = ego_transcriptome_significant,
    aes(x = neg_log10_padj, y = Description, size = odds_ratio, color = neg_log10_padj),
    alpha = 0.8
  ) +
 
  geom_text(
    data = ego_up_transcriptome_bp,
    aes(x = neg_log10_padj, y = Description, label = sprintf("%.2f", odds_ratio)),
    nudge_x = max(ego_up_transcriptome_bp$neg_log10_padj, na.rm = TRUE) * 0.05,
    hjust = 0, size = 16, color = "black", fontface = "bold", vjust = 0.5
  ) +
 
  scale_x_continuous(
    name = "-log10(p.adjust)",
    limits = c(0, max(ego_up_transcriptome_bp$neg_log10_padj, na.rm = TRUE) * 1.3),
    expand = expansion(mult = c(0, 0.1))
  ) +
  
  scale_color_gradient(
    low = "#FFB6C1",  
    high = "#8B0000", 
    name = "-log10(p.adjust)",
    limits = color_range, 
    guide = guide_colorbar(
      title.position = "top",
      title.hjust = 0.5,
      barwidth = unit(2, "cm"),
      barheight = unit(15, "cm")
    )
  ) +
 
  scale_size_continuous(
    range = c(6, 12), 
    name = "Odds Ratio",
    limits = size_range 
  ) +
  labs(y = NULL, title = "Inoculation VS Control WT transcriptome Downregulated") +
  theme_bw(base_size = 14) +
  theme(
    panel.background = element_rect(fill = "white"),
    legend.key = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_text(size = 42),
    axis.text.x = element_text(size = 42),
    axis.text.y = element_text(
      size = 50,
      color = "black",
      lineheight = 0.9
    ),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 52),
    panel.border = element_rect(color = "black", fill = NA),
    legend.title = element_text(size = 40, face = "bold"),
    legend.text = element_text(size = 36),
    legend.key.size = unit(2, "cm"),
    legend.position = "right"
  )

ggsave("Inoculation_VS_WT_genome_Downregulated_ALL_GO_terms2_Improved_adjust.pdf",
       p_up_genome,
       width = 30, 
       height = 30,
       dpi = 300)

ggsave("Inoculation_VS_WT_transcriptome_Downregulated_ALL_GO_terms2_Improved_adjust.pdf",
       p_up_transcriptome,
       width = 30, 
       height = 30,
       dpi = 300)

######################################################################
#====VennPlot-Upregulated====
library(VennDiagram)
library(grid)

genome_deg <- read.csv("Inoculation_VS_WT_genome_Upregulated_genes2.csv", row.names = 1)
transcriptome_deg <- read.csv("Inoculation_VS_WT_transcriptome_Upregulated_genes2.csv", row.names = 1)

# extract gene ID of DEGs（padj < 0.05 and |log2FoldChange| > 1）
genome_sig_genes <- rownames(genome_deg)[!is.na(genome_deg$padj) & 
                                           genome_deg$padj < 0.05 & 
                                           abs(genome_deg$log2FoldChange) > 1]

transcriptome_sig_genes <- rownames(transcriptome_deg)[!is.na(transcriptome_deg$padj) & 
                                                         transcriptome_deg$padj < 0.05 & 
                                                         abs(transcriptome_deg$log2FoldChange) > 1]

gene_lists <- list(
  Genome = genome_sig_genes,
  Transcriptome = transcriptome_sig_genes
)

# create venn plot
venn.plot <- venn.diagram(
  x = gene_lists,
  filename = NULL, 
  height = 10,  
  width = 10,  
  resolution = 300,
  col = "transparent",
  fill = c("#E0367A", "#029149"),  
  alpha = 0.50,
  label.col = c("black", "black", "black"),
  cex = 2.5, 
  fontfamily = "serif",
  fontface = "bold",
  cat.col = c("#E0367A", "#029149"),  
  cat.cex = 2.5,  
  cat.pos = c(-30, 30), 
  cat.dist = c(0.05, 0.05), 
  cat.fontfamily = "serif",
  cat.fontface = "bold",
  rotation.degree = 0,
  margin = 0.2,
  main = "Venn Diagram of Upregulated DEGs between Genome and Transcriptome",  
  main.cex = 2.0,  
  main.col = "black", 
  main.fontfamily = "serif",  
  main.fontface = "bold", 
  main.pos = c(0.5, 0.95) 
)

grid.newpage()
grid.draw(venn.plot)

pdf(file = "CRA004377_Inoculation_VS_WT_Upregulated_Genome_Transcriptome_DEGs_Venn.pdf", width = 10, height = 10)
grid.draw(venn.plot)
dev.off()

png(file = "Genome_Transcriptome_DEGs_Venn.png", width = 10, height = 10, 
    units = "in", res = 300)
grid.draw(venn.plot)
dev.off()








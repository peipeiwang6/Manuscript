# Scripts for our manuscript "Genome-guided vs transcriptome-first mapping: systematic evaluation of RNA-seq quantification bias for tandemly duplicated genes in model and crop plants"

 parameter setting for RNA-seq read mapping (Tophat2, HISAT2, STAR) and read count calculation (HTSeq-count)
## 01_RNA-seq_mapping_and_read_counts.sh

 TPM calling
## 02_TPM_calling.r

 simulate paired-end 150 bp reads from genome sequence and gff files, according to read counts calculated from real RNA-seq data
## python 03_simulate_rnaseq_from_counts.py -g genome.fa -t annotation.gff -c counts.tsv -o outdir

 assess mapping accuracy of mapping aligner (Tophat2, HISAT2, STAR) and strategy combination (G-guided or T-guided)
## python 04_calculate_benchmark_accuracy.py --real real_counts.tsv --sim simulated_counts.tsv

 draw dotplot between TPMs obtained using different aligner and parameter settings
## 05_Dotplot_between_two_TPMs.r

 Fisher's exact test for genes with different TPMs between two mapping strategies in different duplicate categories
## 06_enrichment_for_duplicate_mechanisms.r

 count overlapping genes in a genome
## python 07_count_overlap_genes.py annotation.gff -o Overlapping_genes.txt

 GO enrichment using clusterProfiler
## 08_GO_enrichment.r
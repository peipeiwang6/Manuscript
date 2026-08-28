### make genome index
	# Tophat2
	bowtie2-build genome.fas Tophat2_genome_index 

	# HISAT2
	hisat2-build genome.fas HISAT2_genome_index

	# STAR
	STAR --runThreadN 8 --runMode genomeGenerate \
		--genomeDir STAR_genome_index \
		--genomeSAindexNbases 12 \
		--genomeFastaFiles genome.fas 


# make transcriptome index
	# Tophat2
	tophat2 -G annotation.gff \
		--transcriptome-index=Tophat2_transcriptome_index \
		Tophat2_genome_index

	# HISAT2
	extract_exons.py annotation.gff > splice.txt
	extract_splice_sites.py annotation.gff > exon.txt
	hisat2-build --ss splice.txt --exon exon.txt \
		genome.fas \
		HISAT2_transcriptome_index

	# STAR
	STAR --runThreadN 8 --runMode genomeGenerate \
		--genomeDir STAR_transcriptome_index \
		--genomeFastaFiles genome.fas \
		--sjdbGTFfile annotation.gff \
		--genomeSAindexNbases 12 \
		--sjdbOverhang 99

### mapping RNA-seq reads

	# Tophat2 G-guided mapping 
	tophat2 --library-type fr-unstranded -p 20 -o Tophat2_g --max-multihits 1 -I 100000 \
		Tophat2_genome_index \
		sample_01_1.fastq.gz sample_01_2.fastq.gz

	# Tophat2 T-guided mapping 
	tophat2 --library-type fr-unstranded -p 20 -o Tophat2_t --max-multihits 1 -I 100000 \
		--transcriptome-index=Tophat2_transcriptome_index \
		Tophat2_genome_index \
		sample_01_1.fastq.gz sample_01_2.fastq.gz

	# HISAT2 G-guided mapping 
	hisat2 -p 20 --dta \
		-k 1 \ 
		-x HISAT2_genome_index \
		-1 sample_01_1.fastq.gz \
		-2 sample_01_2.fastq.gz \
		-S ./HISAT2_g/simulated_HISAT2_g.sam

	# HISAT2 T-guided mapping 
	hisat2 -p 20 --dta \
		-k 1 \ 
		-x HISAT2_transcriptome_index \
		-1 sample_01_1.fastq.gz \
		-2 sample_01_2.fastq.gz \
		-S ./HISAT2_t/simulated_HISAT2_t.sam

	# STAR G-guided mapping 
	STAR --runThreadN 20 \
		--genomeDir STAR_genome_index \
		--readFilesIn sample_01_1.fastq.gz sample_01_2.fastq.gz \
		--readFilesCommand zcat \
		--outFileNamePrefix ./STAR_g/STAR_g_ \
		--outSAMtype BAM SortedByCoordinate \
		--outFilterMultimapNmax 1

	# STAR T-guided mapping 
	STAR --runThreadN 20 \
		--genomeDir STAR_transcriptome_index \
		--readFilesIn sample_01_1.fastq.gz sample_01_2.fastq.gz \
		--readFilesCommand zcat \
		--outFileNamePrefix ./STAR_t/STAR_t_ \
		--outSAMtype BAM SortedByCoordinate \
		--outFilterMultimapNmax 1


### sort bam files for Tophat2 and HISAT2 results
	samtools sort -@ 20 ./HISAT2_t/simulated_HISAT2_t.sam -o ./HISAT2_t/simulated_HISAT2_t.sorted.bam

### call read counts using HTSeq-count
	htseq-count --format=bam -m union -s no -t gene -i ID -r pos --nonunique=fraction \
		./HISAT2_t/simulated_HISAT2_t.sorted.bam \
		annotation.gff > \
		./HISAT2_t/simulated_HISAT2_t_HTSeqCount.out

	sed -i '1s;^;gene_id \tcounts\n;' ./HISAT2_t/simulated_HISAT2_t_HTSeqCount.out
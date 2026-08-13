#!/usr/bin/env python3
"""
01_simulate_rnaseq_from_counts.py - Simulate paired-end RNA-seq reads based on real read counts.

Usage:
    python 01_simulate_rnaseq_from_counts.py -g genome.fa -t annotation.gff -c counts.tsv -o outdir

Input count file format (tab-separated, no header):
    gene_id    count

The script extracts the longest transcript for each gene from the GFF,
then generates exactly the specified number of read pairs for that gene.
"""

import argparse
import gzip
import random
import sys
from pathlib import Path
from collections import defaultdict
import numpy as np
from Bio.Seq import Seq
from pyfaidx import Fasta


def reverse_complement(seq):
    """Return the reverse complement of a DNA sequence."""
    return str(Seq(seq).reverse_complement())


def write_fastq_pair(fq1, fq2, read_id, seq1, seq2, qual):
    """Write a pair of reads to two open FASTQ files."""
    fq1.write(f"@{read_id}/1\n{seq1}\n+\n{qual}\n")
    fq2.write(f"@{read_id}/2\n{seq2}\n+\n{qual}\n")


def parse_gff(gff_file, genome, gene_type="gene", transcript_type="mRNA", feature_type="exon"):
    """
    Parse a GFF file and build gene -> transcript -> parts (exons by default).
    Selects the longest transcript per gene based on total length of parts.
    Returns dict: gene_id -> {seq, length, tx_id}
    """
    gene_children = defaultdict(list)       # gene_id -> list of transcript ids
    transcript_children = defaultdict(list) # transcript_id -> list of (chrom, start, end, strand)
    features = {}                           # feature_id -> attributes dict

    with open(gff_file, 'r', encoding='latin-1') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue
            chrom, source, ftype, start, end, score, strand, phase, attrs = parts
            start, end = int(start), int(end)
            # Parse attributes
            attr_dict = {}
            for attr in attrs.split(';'):
                if '=' in attr:
                    key, val = attr.split('=', 1)
                    attr_dict[key] = val
            if ftype == gene_type:
                gene_id = attr_dict.get('ID')
                if gene_id:
                    features[gene_id] = {'type': 'gene', 'chrom': chrom, 'start': start, 'end': end, 'strand': strand}
            elif ftype == transcript_type:
                tx_id = attr_dict.get('ID')
                parent = attr_dict.get('Parent')
                if tx_id and parent:
                    features[tx_id] = {'type': 'transcript', 'chrom': chrom, 'start': start, 'end': end,
                                       'strand': strand, 'parent': parent}
                    gene_children[parent].append(tx_id)
            elif ftype == feature_type:
                parent = attr_dict.get('Parent')
                if parent:
                    transcript_children[parent].append((chrom, start, end, strand))

    # Build gene information
    gene_info = {}
    for gene_id, tx_ids in gene_children.items():
        if gene_id not in features:
            continue
        best_tx = None
        best_len = 0
        for tx_id in tx_ids:
            if tx_id not in transcript_children:
                continue
            parts = transcript_children[tx_id]
            tx_len = sum(end - start + 1 for (_, start, end, _) in parts)
            if tx_len > best_len:
                best_len = tx_len
                best_tx = tx_id
        if best_tx is None:
            continue
        # Retrieve all parts for this transcript and sort by start coordinate
        parts = transcript_children[best_tx]
        parts.sort(key=lambda x: x[1])  # sort by start position
        seq_parts = []
        for chrom, start, end, strand in parts:
            seq = genome[chrom][start-1:end].seq
            seq_parts.append(seq)
        full_seq = "".join(seq_parts)
        if features[gene_id]['strand'] == '-':
            full_seq = reverse_complement(full_seq)
        gene_info[gene_id] = {
            'seq': full_seq,
            'length': len(full_seq),
            'tx_id': best_tx,
        }
    return gene_info


def read_counts_file(counts_file):
    """
    Read a two-column TSV file: gene_id, count.
    Returns a dict {gene_id: count}.
    """
    counts = {}
    with open(counts_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split('\t')
            if len(parts) < 2:
                continue
            gene_id = parts[0]
            try:
                count = int(parts[1])
            except ValueError:
                continue
            if count > 0:
                counts[gene_id] = count
    return counts


def main():
    parser = argparse.ArgumentParser(description="Simulate RNA-seq from real read counts")
    parser.add_argument("-g", "--genome", required=True, help="Genome FASTA (indexed)")
    parser.add_argument("-t", "--annotation", required=True, help="GFF annotation file")
    parser.add_argument("-c", "--counts", required=True, help="Tab-separated file: gene_id count")
    parser.add_argument("-o", "--outdir", default="simulated", help="Output directory")
    parser.add_argument("--readlen", type=int, default=150, help="Read length")
    parser.add_argument("--frag_mean", type=int, default=300, help="Mean fragment length")
    parser.add_argument("--frag_sd", type=int, default=50, help="Std dev fragment length")
    parser.add_argument("--seed", type=int, default=42, help="Random seed")
    parser.add_argument("--gene_type", default="gene", help="Feature type for gene")
    parser.add_argument("--transcript_type", default="mRNA", help="Feature type for transcript")
    parser.add_argument("--feature_type", default="exon",
                        help="Feature type for transcript parts (default: exon; use 'CDS' for coding sequences only)")
    args = parser.parse_args()

    random.seed(args.seed)
    np.random.seed(args.seed)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------
    # 1. Load genome index
    # ------------------------------------------------------------
    print("[1] Loading genome index...", flush=True)
    genome = Fasta(args.genome)

    # ------------------------------------------------------------
    # 2. Read counts file
    # ------------------------------------------------------------
    print("[2] Reading counts file...", flush=True)
    counts_dict = read_counts_file(args.counts)
    print(f"   Total genes with positive counts: {len(counts_dict)}", flush=True)

    # ------------------------------------------------------------
    # 3. Parse GFF and extract transcript sequences
    # ------------------------------------------------------------
    print("[3] Parsing GFF annotation file...", flush=True)
    gene_info = parse_gff(
        args.annotation,
        genome,
        gene_type=args.gene_type,
        transcript_type=args.transcript_type,
        feature_type=args.feature_type
    )
    if not gene_info:
        print("Error: No genes extracted. Please check feature type parameters.", flush=True)
        sys.exit(1)
    print(f"   Extracted {len(gene_info)} genes from GFF.", flush=True)

    # ------------------------------------------------------------
    # 4. Intersect counts with GFF genes and filter length
    # ------------------------------------------------------------
    # Keep only genes that are in counts and have length >= readlen
    valid_genes = {}
    for gid, count in counts_dict.items():
        if gid in gene_info and gene_info[gid]["length"] >= args.readlen:
            valid_genes[gid] = count
        elif gid in gene_info and gene_info[gid]["length"] < args.readlen:
            print(f"   Warning: Gene {gid} has transcript length {gene_info[gid]['length']} < readlen {args.readlen}. Skipping.", flush=True)

    print(f"   Genes with positive counts and valid length: {len(valid_genes)}", flush=True)
    if len(valid_genes) == 0:
        print("Error: No genes available for simulation.", flush=True)
        sys.exit(1)

    # Prepare lists for simulation
    gene_ids = list(valid_genes.keys())
    counts_list = np.array([valid_genes[gid] for gid in gene_ids])
    seqs = [gene_info[gid]["seq"] for gid in gene_ids]
    lengths = np.array([gene_info[gid]["length"] for gid in gene_ids])

    total_reads = counts_list.sum()
    print(f"   Total read pairs to simulate: {total_reads:,}", flush=True)

    # Write out gene length and count info
    with open(outdir / "gene_info.tsv", "w") as f:
        f.write("gene_id\tlength\tassigned_count\n")
        for gid, cnt in valid_genes.items():
            f.write(f"{gid}\t{gene_info[gid]['length']}\t{cnt}\n")
    print(f"   Gene info written to {outdir / 'gene_info.tsv'}", flush=True)

    # ------------------------------------------------------------
    # 5. Simulate reads (paired-end)
    # ------------------------------------------------------------
    print("[4] Simulating reads (this may take a while)...", flush=True)
    fq1 = gzip.open(outdir / "sample_01_1.fastq.gz", "wt")
    fq2 = gzip.open(outdir / "sample_01_2.fastq.gz", "wt")
    qual = "I" * args.readlen

    read_id = 0
    total_bp_written = 0
    skipped = 0  # should be zero as we filtered

    for i, gid in enumerate(gene_ids):
        count = counts_list[i]
        if count <= 0:
            continue
        seq = seqs[i]
        seq_len = lengths[i]

        # Random start positions and fragment lengths
        pos = np.random.randint(0, max(1, seq_len - args.frag_mean), size=count)
        frag_len = np.random.normal(args.frag_mean, args.frag_sd, size=count).astype(int)

        for j in range(count):
            # Adjust fragment length to stay within bounds
            max_len = min(seq_len - pos[j], 600)
            if frag_len[j] < args.readlen:
                frag_len[j] = args.readlen
            if frag_len[j] > max_len:
                frag_len[j] = max_len if max_len > args.readlen else args.readlen

            # Read1: from pos[j] for readlen bases
            start1 = pos[j]
            end1 = start1 + args.readlen
            if end1 > seq_len:
                start1 = seq_len - args.readlen
                end1 = seq_len
            seq1 = seq[start1:end1]
            if len(seq1) != args.readlen:
                seq1 = seq1.ljust(args.readlen, 'N')[:args.readlen]

            # Read2: from pos[j] + frag_len - readlen, reverse complemented
            start2 = pos[j] + frag_len[j] - args.readlen
            if start2 < 0:
                start2 = 0
            end2 = start2 + args.readlen
            if end2 > seq_len:
                end2 = seq_len
                start2 = end2 - args.readlen
            seq2 = seq[start2:end2]
            if len(seq2) != args.readlen:
                seq2 = seq2.ljust(args.readlen, 'N')[:args.readlen]
            seq2_rc = reverse_complement(seq2)

            write_fastq_pair(fq1, fq2, f"read_{read_id}", seq1, seq2_rc, qual)
            read_id += 1
            total_bp_written += 2 * args.readlen

            if read_id % 100000 == 0:
                print(f"   Generated {read_id:,} read pairs", flush=True)

    fq1.close()
    fq2.close()
    print(f"   Done! Generated {read_id:,} read pairs", flush=True)
    print(f"   Total bases written: {total_bp_written:,} bp ({total_bp_written/1e9:.2f} Gb)", flush=True)

    # ------------------------------------------------------------
    # 6. Compute true TPM and output table
    # ------------------------------------------------------------
    print("[5] Computing true TPM and outputting table...", flush=True)
    rpk = np.array([counts_list[i] / lengths[i] if lengths[i] > 0 else 0 for i in range(len(gene_ids))])
    tpm_real = rpk / np.sum(rpk) * 1e6 if np.sum(rpk) > 0 else np.zeros(len(gene_ids))
    tpm_dict = {gid: tpm_real[i] for i, gid in enumerate(gene_ids)}

    with open(outdir / "simulated_counts_tpm.tsv", "w") as f:
        f.write("gene_id\tlength\tcounts\tTPM\n")
        for gid in gene_ids:
            f.write(f"{gid}\t{gene_info[gid]['length']}\t{valid_genes[gid]}\t{tpm_dict[gid]:.2f}\n")

    print("All tasks completed. Output files are located in:", outdir, flush=True)


if __name__ == "__main__":
    main()
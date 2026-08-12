#!/usr/bin/env python3
"""
01_simulate_rnaseq.py - Simulate RNA-seq with expression distribution matching real data.

Expression distribution:
    The log10(TPM+1) values of genes follow a probability density proportional to:
        f(x) = x^A1 + C1 * exp( -(x - mean1)^2 / (2 * sigma1^2) )
    where x = log10(TPM+1), with parameters:
        A1 = -0.2315, C1 = 1.1590, mean1 = 1.2269, sigma1 = 0.7880.
    This ensures that the log10(frequency) of the histogram (bin width 0.01)
    matches the given mixed power-law plus Gaussian shape.

Usage:
    python 01_simulate_rnaseq.py -g genome.fa -t annotation.gff -o outdir
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
    Parse a GFF file and build a hierarchical structure: gene -> transcript -> exons/CDS.
    Uses the specified feature_type (default: exon) to define transcript parts.
    Returns a dict: gene_id -> {seq, length, tx_id}
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
                # If feature_type is 'exon', this collects all exons (including UTRs).
                # If feature_type is 'CDS', it collects only coding sequences.
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


def generate_expression_mixed_power_gaussian(n, A1, C1, mean1, sigma1, x_min=1e-6, x_max=10.0, seed=None):
    """
    Generate n log10(TPM+1) values from a mixed distribution (power-law decay + Gaussian peak).
    Uses rejection sampling.
    """
    if seed is not None:
        np.random.seed(seed)

    def f(x):
        return x ** A1 + C1 * np.exp(-(x - mean1) ** 2 / (2 * sigma1 ** 2))

    # Determine effective upper bound
    x_vals = np.linspace(x_min, x_max, 1000)
    f_vals = f(x_vals)
    max_f = np.max(f_vals)
    cutoff = None
    for xv, fv in zip(x_vals, f_vals):
        if fv < max_f * 1e-6:
            cutoff = xv
            break
    if cutoff is not None and cutoff < x_max:
        x_max = cutoff
        x_vals = np.linspace(x_min, x_max, 2000)
        f_vals = f(x_vals)
        max_f = np.max(f_vals)

    samples = []
    n_accepted = 0
    batch_size = n * 10
    while n_accepted < n:
        x_cand = np.random.uniform(x_min, x_max, batch_size)
        y_cand = np.random.uniform(0, max_f, batch_size)
        f_cand = f(x_cand)
        accepted = y_cand < f_cand
        samples.extend(x_cand[accepted])
        n_accepted += np.sum(accepted)
        if n_accepted >= n:
            break
    return np.array(samples[:n])


def main():
    parser = argparse.ArgumentParser(description="Simulate RNA-seq with mixed power+gaussian distribution")
    parser.add_argument("-g", "--genome", required=True, help="Genome FASTA (indexed)")
    parser.add_argument("-t", "--annotation", required=True, help="GFF annotation file")
    parser.add_argument("-o", "--outdir", default="simulated", help="Output directory")
    parser.add_argument("--total_bp", type=float, default=60e9, help="Total bases to sequence (bp), default 60Gb")
    parser.add_argument("--readlen", type=int, default=150, help="Read length")
    parser.add_argument("--frag_mean", type=int, default=300, help="Mean fragment length")
    parser.add_argument("--frag_sd", type=int, default=50, help="Std dev fragment length")
    parser.add_argument("--seed", type=int, default=42, help="Random seed")
    parser.add_argument("--num_genes", type=int, default=None, help="Number of genes to use (random subset)")
    parser.add_argument("--gene_type", default="gene", help="Feature type for gene")
    parser.add_argument("--transcript_type", default="mRNA", help="Feature type for transcript")
    parser.add_argument("--feature_type", default="exon",
                        help="Feature type for transcript parts (default: exon; use 'CDS' for coding sequences only)")
    # Distribution parameters (fitted from real data)
    parser.add_argument("--A1", type=float, default=-0.2315, help="Power exponent")
    parser.add_argument("--C1", type=float, default=1.1590, help="Gaussian amplitude")
    parser.add_argument("--mean1", type=float, default=1.2269, help="Gaussian mean")
    parser.add_argument("--sigma1", type=float, default=0.7880, help="Gaussian sigma")
    args = parser.parse_args()

    random.seed(args.seed)
    np.random.seed(args.seed)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------
    # 1. Load genome index
    # ------------------------------------------------------------
    print("[1] Loading genome index...")
    genome = Fasta(args.genome)

    # ------------------------------------------------------------
    # 2. Parse GFF (pure Python, no extra output)
    # ------------------------------------------------------------
    print("[2] Parsing GFF annotation file...")
    gene_info = parse_gff(
        args.annotation,
        genome,
        gene_type=args.gene_type,
        transcript_type=args.transcript_type,
        feature_type=args.feature_type
    )
    if not gene_info:
        print("Error: No genes extracted. Please check feature type parameters.")
        sys.exit(1)

    # Write out lengths of all extracted genes
    length_file_all = outdir / "all_gene_lengths.tsv"
    with open(length_file_all, "w") as f:
        f.write("gene_id\ttranscript_id\tlength\n")
        for gid, info in gene_info.items():
            f.write(f"{gid}\t{info['tx_id']}\t{info['length']}\n")
    print(f"   All gene length file: {length_file_all} (total {len(gene_info)} genes)")

    # Filter out genes shorter than read length
    original_count = len(gene_info)
    gene_info = {gid: info for gid, info in gene_info.items() if info["length"] >= args.readlen}
    print(f"   Original gene count: {original_count}, after filtering (length >= {args.readlen}): {len(gene_info)}")

    # Optionally select a random subset of genes
    if args.num_genes is not None and args.num_genes > 0:
        all_ids = list(gene_info.keys())
        if args.num_genes < len(all_ids):
            chosen = random.sample(all_ids, args.num_genes)
            gene_info = {gid: gene_info[gid] for gid in chosen}
        elif args.num_genes > len(all_ids):
            print(f"Warning: --num_genes ({args.num_genes}) exceeds available genes ({len(all_ids)}). Using all.")

    n_genes = len(gene_info)
    if n_genes == 0:
        print("Error: No genes available for simulation.")
        sys.exit(1)

    gene_ids = list(gene_info.keys())
    seqs = [gene_info[gid]["seq"] for gid in gene_ids]
    lengths = np.array([gene_info[gid]["length"] for gid in gene_ids])
    print(f"   Final number of genes: {n_genes}")
    print(f"   Transcript lengths: min {lengths.min()}, max {lengths.max()}, median {np.median(lengths):.0f} bp")

    # Write out lengths of genes used in simulation
    length_file_used = outdir / "used_gene_lengths.tsv"
    with open(length_file_used, "w") as f:
        f.write("gene_id\ttranscript_id\tlength\n")
        for gid, info in gene_info.items():
            f.write(f"{gid}\t{info['tx_id']}\t{info['length']}\n")
    print(f"   Used gene length file: {length_file_used}")

    # ------------------------------------------------------------
    # 3. Generate expression: log10(TPM+1) ~ mixed distribution (power + gaussian)
    # ------------------------------------------------------------
    print("[3] Generating gene expression (log10(TPM+1) from mixed distribution)...")
    log10_tpm1 = generate_expression_mixed_power_gaussian(
        n=n_genes,
        A1=args.A1,
        C1=args.C1,
        mean1=args.mean1,
        sigma1=args.sigma1,
        seed=args.seed
    )
    tpm_raw = 10 ** log10_tpm1 - 1
    tpm_raw = np.maximum(tpm_raw, 0.0)

    # Weight = TPM * transcript length
    weights = tpm_raw * lengths
    if np.sum(weights) == 0:
        weights += 1e-9

    total_reads = int(args.total_bp / (2 * args.readlen))
    print(f"   Target read pairs: {total_reads:,}")

    expected_counts = weights / np.sum(weights) * total_reads
    counts = np.round(expected_counts).astype(int)
    diff = total_reads - counts.sum()
    if diff != 0 and n_genes > 0:
        idx = np.random.choice(n_genes, size=abs(diff), replace=True)
        counts[idx] += np.sign(diff) * 1
    counts = np.maximum(counts, 0)
    gene_counts = {gid: int(counts[i]) for i, gid in enumerate(gene_ids)}
    print(f"   Actual assigned reads: {sum(gene_counts.values()):,}")

    # ------------------------------------------------------------
    # 4. Simulate reads (paired-end)
    # ------------------------------------------------------------
    print("[4] Simulating reads (this may take a while)...")
    fq1 = gzip.open(outdir / "sample_01_1.fastq.gz", "wt")
    fq2 = gzip.open(outdir / "sample_01_2.fastq.gz", "wt")
    qual = "I" * args.readlen

    read_id = 0
    total_bp_written = 0
    skipped = 0

    for i, gid in enumerate(gene_ids):
        count = gene_counts[gid]
        if count <= 0:
            continue
        seq = seqs[i]
        seq_len = lengths[i]
        if seq_len < args.readlen:
            skipped += count
            continue

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
                print(f"   Generated {read_id:,} read pairs")

    fq1.close()
    fq2.close()
    print(f"   Done! Generated {read_id:,} read pairs")
    print(f"   Skipped reads due to short transcripts: {skipped:,}")
    print(f"   Total bases written: {total_bp_written:,} bp ({total_bp_written/1e9:.2f} Gb)")

    # ------------------------------------------------------------
    # 5. Compute true TPM and output table
    # ------------------------------------------------------------
    print("[5] Computing true TPM and outputting table...")
    rpk = np.array([gene_counts[gid] / gene_info[gid]["length"] if gene_info[gid]["length"] > 0 else 0 for gid in gene_ids])
    tpm_real = rpk / np.sum(rpk) * 1e6 if np.sum(rpk) > 0 else np.zeros(n_genes)
    tpm_dict = {gid: tpm_real[i] for i, gid in enumerate(gene_ids)}

    with open(outdir / "simulated_counts_tpm.tsv", "w") as f:
        f.write("gene_id\tlength\tcounts\tTPM\n")
        for gid in gene_ids:
            f.write(f"{gid}\t{gene_info[gid]['length']}\t{gene_counts[gid]}\t{tpm_dict[gid]:.2f}\n")

    print("All tasks completed. Output files are located in:", outdir)


if __name__ == "__main__":
    main()
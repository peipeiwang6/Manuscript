#!/usr/bin/env python3
"""
03_simulate_rnaseq_from_counts.py - Simulate paired-end RNA-seq based on real read counts.
For transcripts shorter than read length, we use the full transcript sequence padded with 'N'
instead of cyclic repetition, as requested.

Usage:
    python 03_simulate_rnaseq_from_counts.py -g genome.fa -t annotation.gff -c counts.tsv -o outdir
    [--id_type gene|transcript] [--target_reads 200000000] [--scale 100]
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
    return str(Seq(seq).reverse_complement())


def write_fastq_pair(fq1, fq2, read_id, seq1, seq2, qual):
    fq1.write(f"@{read_id}/1\n{seq1}\n+\n{qual}\n")
    fq2.write(f"@{read_id}/2\n{seq2}\n+\n{qual}\n")


def parse_gff(gff_file, genome, gene_type="gene", transcript_type="mRNA", feature_type="exon"):
    """
    Parse GFF and return dictionaries for genes and transcripts.
    """
    gene_data = {}
    transcript_data = {}
    gene_children = defaultdict(list)
    transcript_children = defaultdict(list)
    feature_attrs = {}

    with open(gff_file, 'r', encoding='latin-1') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue
            chrom, source, ftype, start, end, score, strand, phase, attrs = parts
            start, end = int(start), int(end)
            attr_dict = {}
            for attr in attrs.split(';'):
                if '=' in attr:
                    key, val = attr.split('=', 1)
                    attr_dict[key] = val
            if ftype == gene_type:
                gene_id = attr_dict.get('ID')
                if gene_id:
                    feature_attrs[gene_id] = {'chrom': chrom, 'start': start, 'end': end, 'strand': strand,
                                              'name': attr_dict.get('Name', gene_id)}
            elif ftype == transcript_type:
                tx_id = attr_dict.get('ID')
                parent = attr_dict.get('Parent')
                if tx_id and parent:
                    feature_attrs[tx_id] = {'chrom': chrom, 'start': start, 'end': end, 'strand': strand,
                                            'parent': parent, 'name': attr_dict.get('Name', tx_id)}
                    gene_children[parent].append(tx_id)
            elif ftype == feature_type:
                parent = attr_dict.get('Parent')
                if parent:
                    transcript_children[parent].append((chrom, start, end, strand))

    # Build gene sequences
    for gene_id, tx_ids in gene_children.items():
        if gene_id not in feature_attrs:
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
        parts = transcript_children[best_tx]
        parts.sort(key=lambda x: x[1])
        seq_parts = []
        for chrom, start, end, strand in parts:
            seq = genome[chrom][start-1:end].seq
            seq_parts.append(seq)
        full_seq = "".join(seq_parts)
        if feature_attrs[gene_id]['strand'] == '-':
            full_seq = reverse_complement(full_seq)
        gene_data[gene_id] = {
            'seq': full_seq,
            'length': len(full_seq),
            'tx_id': best_tx,
            'name': feature_attrs[gene_id].get('name', gene_id)
        }
        transcript_data[best_tx] = {
            'seq': full_seq,
            'length': len(full_seq),
            'gene_id': gene_id,
            'name': feature_attrs[best_tx].get('name', best_tx)
        }

    return gene_data, transcript_data


def read_counts_file(counts_file):
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
                count = float(parts[1])  # allow non-integer
            except ValueError:
                continue
            if count > 0:
                counts[gene_id] = count
    return counts


def strip_version(id_str):
    return id_str.split('.')[0] if '.' in id_str else id_str


def strip_prefix(id_str):
    for prefix in ['gene:', 'transcript:', 'locus:']:
        if id_str.startswith(prefix):
            return id_str[len(prefix):]
    return id_str


def match_ids(count_ids, gff_data, use_name=False):
    matched = {}
    unmatched = []
    strategies = {}
    id_to_key = {k: k for k in gff_data}
    name_to_key = {}
    if use_name:
        for key, info in gff_data.items():
            name = info.get('name')
            if name and name not in name_to_key:
                name_to_key[name] = key

    for cid in count_ids:
        found = False
        # exact
        if cid in id_to_key:
            matched[cid] = id_to_key[cid]
            strategies[cid] = 'exact'
            found = True
        # strip version
        elif not found:
            stripped = strip_version(cid)
            if stripped in id_to_key:
                matched[cid] = id_to_key[stripped]
                strategies[cid] = 'strip_version'
                found = True
        # strip prefix
        elif not found:
            no_prefix = strip_prefix(cid)
            if no_prefix in id_to_key:
                matched[cid] = id_to_key[no_prefix]
                strategies[cid] = 'strip_prefix'
                found = True
        # strip version + prefix
        elif not found:
            no_prefix = strip_prefix(cid)
            stripped = strip_version(no_prefix)
            if stripped in id_to_key:
                matched[cid] = id_to_key[stripped]
                strategies[cid] = 'strip_version+prefix'
                found = True
        # name match
        elif not found and use_name:
            if cid in name_to_key:
                matched[cid] = name_to_key[cid]
                strategies[cid] = 'name'
                found = True
            else:
                stripped = strip_version(cid)
                if stripped in name_to_key:
                    matched[cid] = name_to_key[stripped]
                    strategies[cid] = 'name_strip_version'
                    found = True
        if not found:
            unmatched.append(cid)

    strategy_counts = defaultdict(int)
    for v in strategies.values():
        strategy_counts[v] += 1
    if strategy_counts:
        print("   Matching strategy statistics:")
        for s, cnt in sorted(strategy_counts.items()):
            print(f"      {s}: {cnt}")
    if unmatched:
        print(f"   Unmatched count IDs (showing up to 10):")
        for uid in unmatched[:10]:
            print(f"      {uid}")
        if len(unmatched) > 10:
            print(f"      ... and {len(unmatched)-10} more")
    return matched


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
    parser.add_argument("--scale", type=float, default=None, help="Multiply all counts by this factor.")
    parser.add_argument("--target_reads", type=int, default=None, help="Target total read pairs (scales proportionally).")
    parser.add_argument("--gene_type", default="gene", help="Feature type for gene")
    parser.add_argument("--transcript_type", default="mRNA", help="Feature type for transcript")
    parser.add_argument("--feature_type", default="exon", help="Feature type for transcript parts")
    parser.add_argument("--id_type", choices=["gene", "transcript"], default="gene",
                        help="Whether counts are gene-level or transcript-level")
    parser.add_argument("--output_all_genes", action="store_true", default=False,
                        help="Also output a table with all GFF genes (including zero counts)")
    args = parser.parse_args()

    random.seed(args.seed)
    np.random.seed(args.seed)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    print("[1] Loading genome index...", flush=True)
    genome = Fasta(args.genome)

    print("[2] Reading counts file...", flush=True)
    counts_dict = read_counts_file(args.counts)
    print(f"   Total genes with positive counts: {len(counts_dict)}", flush=True)

    print("[3] Parsing GFF annotation file...", flush=True)
    gene_data, transcript_data = parse_gff(
        args.annotation, genome,
        gene_type=args.gene_type,
        transcript_type=args.transcript_type,
        feature_type=args.feature_type
    )
    print(f"   Extracted {len(gene_data)} genes and {len(transcript_data)} transcripts from GFF.", flush=True)

    if args.id_type == "gene":
        gff_data = gene_data
        use_name = True
    else:
        gff_data = transcript_data
        use_name = True

    print(f"[4] Matching count IDs to GFF {args.id_type} IDs...", flush=True)
    matched = match_ids(list(counts_dict.keys()), gff_data, use_name=use_name)
    print(f"   Successfully matched {len(matched)} out of {len(counts_dict)} IDs.", flush=True)

    if len(matched) == 0:
        print("Error: No IDs matched. Please check ID formats and consider using --id_type transcript.", flush=True)
        sys.exit(1)

    # Build valid genes with sequences; round counts to integers
    valid_genes = {}
    for count_id, gff_key in matched.items():
        info = gff_data[gff_key]
        cnt = int(round(counts_dict[count_id]))
        if cnt == 0:
            continue
        valid_genes[count_id] = {
            'seq': info['seq'],
            'length': info['length'],
            'count': cnt,
            'gff_id': gff_key
        }

    # Scaling
    if args.target_reads is not None:
        original_total = sum([v['count'] for v in valid_genes.values()])
        if original_total == 0:
            print("Error: Original total reads is zero.", flush=True)
            sys.exit(1)
        scale_factor = args.target_reads / original_total
        print(f"   Scaling all counts by factor {scale_factor:.6f} to reach target total reads {args.target_reads:,}", flush=True)
        for gid in valid_genes:
            valid_genes[gid]['count'] = int(round(valid_genes[gid]['count'] * scale_factor))
    elif args.scale is not None:
        print(f"   Scaling all counts by factor {args.scale}", flush=True)
        for gid in valid_genes:
            valid_genes[gid]['count'] = int(round(valid_genes[gid]['count'] * args.scale))

    valid_genes = {gid: info for gid, info in valid_genes.items() if info['count'] > 0}
    print(f"   After scaling, genes with positive counts: {len(valid_genes)}", flush=True)

    gene_ids = list(valid_genes.keys())
    counts_list = np.array([valid_genes[gid]['count'] for gid in gene_ids])
    seqs = [valid_genes[gid]['seq'] for gid in gene_ids]
    lengths = np.array([valid_genes[gid]['length'] for gid in gene_ids])

    total_reads = counts_list.sum()
    print(f"   Total read pairs to simulate: {total_reads:,}", flush=True)

    # Write gene info
    with open(outdir / "gene_info.tsv", "w") as f:
        f.write("count_id\tgff_id\tlength\tassigned_count\n")
        for gid in gene_ids:
            f.write(f"{gid}\t{valid_genes[gid]['gff_id']}\t{valid_genes[gid]['length']}\t{valid_genes[gid]['count']}\n")
    print(f"   Gene info written to {outdir / 'gene_info.tsv'}", flush=True)

    # ------------------------------------------------------------
    # Simulation
    # ------------------------------------------------------------
    print("[5] Simulating reads (this may take a while)...", flush=True)
    fq1 = gzip.open(outdir / "sample_01_1.fastq.gz", "wt")
    fq2 = gzip.open(outdir / "sample_01_2.fastq.gz", "wt")
    qual = "I" * args.readlen

    read_id = 0
    total_bp_written = 0
    short_total = 0

    for i, gid in enumerate(gene_ids):
        count = counts_list[i]
        if count <= 0:
            continue
        seq = seqs[i]
        seq_len = lengths[i]

        # For transcripts shorter than readlen: use full sequence + N padding
        if seq_len < args.readlen:
            short_total += count
            pad_len = args.readlen - seq_len
            seq_padded = seq + 'N' * pad_len
            seq_rc = reverse_complement(seq)
            seq2_padded = seq_rc + 'N' * pad_len
            for _ in range(count):
                write_fastq_pair(fq1, fq2, f"read_{read_id}", seq_padded, seq2_padded, qual)
                read_id += 1
                total_bp_written += 2 * args.readlen
                if read_id % 100000 == 0:
                    print(f"   Generated {read_id:,} read pairs", flush=True)
            continue

        # Normal case: seq_len >= readlen
        pos = np.random.randint(0, max(1, seq_len - args.frag_mean), size=count)
        frag_len = np.random.normal(args.frag_mean, args.frag_sd, size=count).astype(int)

        for j in range(count):
            max_len = min(seq_len - pos[j], 600)
            if frag_len[j] < args.readlen:
                frag_len[j] = args.readlen
            if frag_len[j] > max_len:
                frag_len[j] = max_len if max_len > args.readlen else args.readlen

            start1 = pos[j]
            end1 = start1 + args.readlen
            if end1 > seq_len:
                start1 = seq_len - args.readlen
                end1 = seq_len
            seq1 = seq[start1:end1]
            if len(seq1) != args.readlen:
                seq1 = seq1.ljust(args.readlen, 'N')[:args.readlen]

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
    print(f"   Total read pairs from short transcripts (padded with N): {short_total:,}", flush=True)
    print(f"   Total bases written: {total_bp_written:,} bp ({total_bp_written/1e9:.2f} Gb)", flush=True)

    # ------------------------------------------------------------
    # TPM calculation and output
    # ------------------------------------------------------------
    print("[6] Computing true TPM and outputting table...", flush=True)
    rpk = np.array([counts_list[i] / lengths[i] if lengths[i] > 0 else 0 for i in range(len(gene_ids))])
    tpm_real = rpk / np.sum(rpk) * 1e6 if np.sum(rpk) > 0 else np.zeros(len(gene_ids))
    tpm_dict = {gid: tpm_real[i] for i, gid in enumerate(gene_ids)}

    with open(outdir / "simulated_counts_tpm.tsv", "w") as f:
        f.write("gene_id\tlength\tcounts\tTPM\n")
        for gid in gene_ids:
            f.write(f"{gid}\t{valid_genes[gid]['length']}\t{valid_genes[gid]['count']}\t{tpm_dict[gid]:.2f}\n")

    if args.output_all_genes:
        print("   Generating full table with all GFF genes (including zero counts)...", flush=True)
        all_ids = list(gene_data.keys()) if args.id_type == "gene" else list(transcript_data.keys())
        full_counts = {gid: 0 for gid in all_ids}
        full_lengths = {gid: (gene_data[gid]['length'] if args.id_type == "gene" else transcript_data[gid]['length']) for gid in all_ids}
        for gid in gene_ids:
            gff_key = valid_genes[gid]['gff_id']
            if gff_key in full_counts:
                full_counts[gff_key] = valid_genes[gid]['count']
        rpk_all = np.array([full_counts[gid] / full_lengths[gid] if full_counts[gid] > 0 else 0 for gid in all_ids])
        sum_rpk = np.sum(rpk_all)
        tpm_all = rpk_all / sum_rpk * 1e6 if sum_rpk > 0 else np.zeros(len(all_ids))
        with open(outdir / "simulated_counts_tpm_full.tsv", "w") as f:
            f.write("gene_id\tlength\tcounts\tTPM\n")
            for idx, gid in enumerate(all_ids):
                f.write(f"{gid}\t{full_lengths[gid]}\t{full_counts[gid]}\t{tpm_all[idx]:.2f}\n")
        print(f"   Full table written to {outdir / 'simulated_counts_tpm_full.tsv'}", flush=True)

    print("All tasks completed. Output files are located in:", outdir, flush=True)


if __name__ == "__main__":
    main()

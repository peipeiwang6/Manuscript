#!/usr/bin/env python3
"""
07_count_overlap_genes.py
Count overlapping genes annotated in a GFF file.
Overlap is defined as two genes on the same chromosome having partially or completely overlapping intervals (strand ignored).
"""

import sys
import argparse
from collections import defaultdict

def parse_attributes(attr_str):
    """Parse GFF attribute column, return dict."""
    attrs = {}
    for item in attr_str.split(';'):
        if '=' in item:
            k, v = item.split('=', 1)
            attrs[k.strip()] = v.strip()
    return attrs

def main():
    parser = argparse.ArgumentParser(
        description='Count overlapping genes in GFF (based on position overlap)'
    )
    parser.add_argument('gff', help='Input GFF file path')
    parser.add_argument('--type', default='gene',
                        help='Feature type to count (default: gene)')
    parser.add_argument('--id-tag', default='ID',
                        help='Tag for gene identifier in attributes (default: ID)')
    parser.add_argument('--verbose', action='store_true',
                        help='Print overlapping gene groups to stdout (one group per line, tab-separated)')
    parser.add_argument('-o', '--output', metavar='FILE',
                        help='Write overlapping gene groups to FILE (one group per line, tab-separated)')
    args = parser.parse_args()

    # Store genes by chromosome: seqid -> list of (start, end, gene_id)
    genes_by_chr = defaultdict(list)

    with open(args.gff, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
            seqid, source, ftype, start, end, score, strand, phase, attrs = fields
            if ftype != args.type:
                continue
            start = int(start)
            end = int(end)
            attr_dict = parse_attributes(attrs)
            gene_id = attr_dict.get(args.id_tag)
            if gene_id is None:
                # If no ID, use "seqid:start-end" as unique identifier
                gene_id = f"{seqid}:{start}-{end}"
            genes_by_chr[seqid].append((start, end, gene_id))

    total_genes = sum(len(v) for v in genes_by_chr.values())
    overlap_ids = set()          # All overlapping gene IDs (for counting)
    overlap_groups = []          # List of groups, each group is a list of gene IDs

    for seqid, genes in genes_by_chr.items():
        # Sort by start position, then by end position
        genes.sort(key=lambda x: (x[0], x[1]))

        current_group = []      # Gene IDs in current connected component
        current_end = None      # Rightmost end of current connected component

        for start, end, gene_id in genes:
            if not current_group or start > current_end:
                # Finalize previous group if it had more than one gene
                if len(current_group) > 1:
                    overlap_ids.update(current_group)
                    overlap_groups.append(current_group)
                # Start new group
                current_group = [gene_id]
                current_end = end
            else:
                # Overlaps with current group, add to group and update right end
                current_group.append(gene_id)
                if end > current_end:
                    current_end = end

        # Handle last group
        if len(current_group) > 1:
            overlap_ids.update(current_group)
            overlap_groups.append(current_group)

    overlap_count = len(overlap_ids)

    # Print statistics to stdout
    print(f"Total genes: {total_genes}")
    print(f"Overlapping genes: {overlap_count}")
    if total_genes > 0:
        print(f"Overlap ratio: {overlap_count / total_genes * 100:.2f}%")

    # Output groups to file if specified
    if args.output:
        with open(args.output, 'w') as out_f:
            for group in overlap_groups:
                out_f.write('\t'.join(group) + '\n')

    # Output groups to stdout if verbose
    if args.verbose:
        for group in overlap_groups:
            print('\t'.join(group))

if __name__ == '__main__':
    main()
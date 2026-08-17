#!/usr/bin/env python3
"""
04_calculate_benchmark_accuracy.py - Calculate Pearson correlation coefficient (PCC), p-value, and MSE between two count files.

Usage:
    python 04_calculate_benchmark_accuracy.py --real real_counts.tsv --sim simulated_counts.tsv
    [--output plot.png] [--log] [--delimiter TAB]

Input files must be two-column TSV: gene_id (str), count (numeric).
The script computes:
    - PCC and p-value using scipy.stats.pearsonr
    - Mean Squared Error (MSE) = mean((real - sim)^2)
on the intersection of gene IDs.
Options:
    --log: Apply log2(count+1) transformation before correlation and MSE.
    --output: Save scatter plot to file (requires matplotlib).
"""

import argparse
import sys
import numpy as np

try:
    import scipy.stats as stats
    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False

try:
    import matplotlib.pyplot as plt
    HAS_MPL = True
except ImportError:
    HAS_MPL = False


def read_counts(filename, delimiter='\t'):
    """Read two-column count file, return dict {gene_id: count}."""
    counts = {}
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split(delimiter)
            if len(parts) < 2:
                continue
            gene = parts[0].strip()
            try:
                cnt = float(parts[1])
            except ValueError:
                continue
            if cnt > 0:
                counts[gene] = cnt
    return counts


def main():
    parser = argparse.ArgumentParser(description="Compute PCC, p-value, and MSE between two read count files.")
    parser.add_argument("--real", required=True, help="Real counts file (gene_id\\tcount)")
    parser.add_argument("--sim", required=True, help="Simulated counts file (gene_id\\tcount)")
    parser.add_argument("--delimiter", default="\t", help="Delimiter (default: tab)")
    parser.add_argument("--log", action="store_true", help="Apply log2(count+1) transformation before calculations")
    parser.add_argument("--output", help="Save scatter plot to this file (requires matplotlib)")
    args = parser.parse_args()

    if not HAS_SCIPY:
        print("Error: scipy is required for p-value calculation. Please install: pip install scipy", flush=True)
        sys.exit(1)

    print("[1] Reading count files...", flush=True)
    real = read_counts(args.real, args.delimiter)
    sim = read_counts(args.sim, args.delimiter)
    print(f"   Real genes: {len(real)}, Sim genes: {len(sim)}", flush=True)

    # Intersection
    common = set(real.keys()) & set(sim.keys())
    print(f"   Common genes: {len(common)}", flush=True)
    if len(common) == 0:
        print("Error: No common genes found. Please check gene IDs.", flush=True)
        sys.exit(1)

    # Extract values in consistent order
    genes = sorted(common)  # deterministic
    x = np.array([real[g] for g in genes], dtype=np.float64)
    y = np.array([sim[g] for g in genes], dtype=np.float64)

    # Optional log transformation
    transform_label = ""
    if args.log:
        x = np.log2(x + 1)
        y = np.log2(y + 1)
        transform_label = " (log2)"
        print("   Applied log2(count+1) transformation.", flush=True)

    # Compute PCC and p-value using scipy
    pcc, pval = stats.pearsonr(x, y)
    n = len(common)

    # Compute MSE
    mse = np.mean((x - y) ** 2)

    print(f"\n[2] Pearson correlation coefficient: {pcc:.6f}", flush=True)
    print(f"    p-value (two-sided): {pval:.6e}", flush=True)
    print(f"    Number of genes (n): {n}", flush=True)
    print(f"    Mean Squared Error (MSE){transform_label}: {mse:.6f}", flush=True)

    # Scatter plot if requested
    if args.output:
        if not HAS_MPL:
            print("Warning: matplotlib not installed. Cannot generate plot.", flush=True)
        else:
            print(f"\n[3] Generating scatter plot: {args.output}", flush=True)
            plt.figure(figsize=(6, 6))
            plt.scatter(x, y, s=5, alpha=0.5, edgecolors='none')
            plt.xlabel(f'Real counts{transform_label}')
            plt.ylabel(f'Simulated counts{transform_label}')
            plt.title(f'PCC = {pcc:.4f}, p = {pval:.2e}, MSE = {mse:.4f} (n={n})')
            plt.grid(True, linestyle='--', alpha=0.5)
            plt.tight_layout()
            plt.savefig(args.output, dpi=150)
            print(f"   Plot saved to {args.output}", flush=True)

    print("\nAll tasks completed.", flush=True)


if __name__ == "__main__":
    main()
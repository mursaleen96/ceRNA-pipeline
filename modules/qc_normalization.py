# modules/qc_normalization.py (simplified version)

import pandas as pd
import numpy as np
import yaml
import os
import sys

def filter_low_expression(counts_df, min_counts=5, min_samples_frac=0.8):
    """Filter out genes with fewer than min_counts in less than min_samples_frac fraction of samples"""
    n_samples = counts_df.shape[1]
    pass_filter = (counts_df >= min_counts).sum(axis=1) >= (n_samples * min_samples_frac)
    filtered_df = counts_df[pass_filter]
    return filtered_df

def normalize_counts_cpm(counts_df):
    """Simple CPM (Counts Per Million) normalization"""
    total_counts = counts_df.sum(axis=0)
    normalized_df = counts_df.div(total_counts, axis=1) * 1e6
    return normalized_df

def main():
    # Load config parameters
    cfg_path = "config/config.yaml"
    if os.path.exists(cfg_path):
        with open(cfg_path, "r") as f:
            cfg = yaml.safe_load(f)
    else:
        cfg = {}

    # Input counts
    counts_path = "data/input_counts.csv"
    counts_df = pd.read_csv(counts_path, index_col=0)

    print(f"Loaded {counts_df.shape[0]} genes across {counts_df.shape[1]} samples")

    # Basic sanity check
    if (counts_df < 0).any().any():
        raise ValueError("Input counts contain negative values. Check input file.")

    # Fill NA with zeros
    counts_df = counts_df.fillna(0)

    # Filter low-expression genes
    filtered_counts = filter_low_expression(
        counts_df,
        cfg.get("low_count_threshold", 5),
        cfg.get("sample_frac_threshold", 0.8)
    )
    print(f"Filtered out {counts_df.shape[0] - filtered_counts.shape[0]} genes due to low expression")

    # Simple CPM normalization (instead of edgeR TMM)
    print("Using CPM normalization (simplified method)")
    normalized_df = normalize_counts_cpm(filtered_counts)

    # Log2 transform with pseudocount
    normalized_df = np.log2(normalized_df + 1)
    
    # Save normalized counts and sample metadata
    os.makedirs("results", exist_ok=True)
    normalized_df.to_csv("results/norm_counts.csv")
    
    # Create simple sample metadata
    metadata = pd.DataFrame({
        'sample_id': counts_df.columns,
        'batch': ['batch1'] * len(counts_df.columns)
    })
    metadata.to_csv("results/sample_metadata.csv", index=False)

    print("QC and normalization completed successfully.")

if __name__ == "__main__":
    main()
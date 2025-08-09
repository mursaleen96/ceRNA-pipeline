import pandas as pd

def standardize_miRNA_counts(in_path, out_path):
    df = pd.read_csv(in_path, index_col=0)  # Assuming first column is gene/miRNA ID
    
    # Standardize miRNA rows (assuming they start with "MIR" or similar)
    def standardize(name):
        if name.upper().startswith("MIR"):
            # Add "hsa-miR-" prefix, lowercase, and assume no arm for simplicity
            return f"hsa-miR-{name[3:].lower()}"  # e.g., MIR3935 → hsa-miR-3935
        return name  # Leave non-miRNA rows unchanged
    
    df.index = [standardize(idx) for idx in df.index]
    
    # Drop duplicates if any after mapping
    df = df[~df.index.duplicated(keep='first')]
    
    df.to_csv(out_path)
    print(f"Standardized miRNA names and saved to {out_path}")

# Run with your file
standardize_miRNA_counts("GSE40419_RPKM.csv", "GSE40419_standardized.csv")
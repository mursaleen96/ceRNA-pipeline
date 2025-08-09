# modules/download_databases.py

import os
import pandas as pd
import shutil

def process_mirtarbase(src, dest):
    """Process miRTarBase txt file ('miRNA', 'Target'), rename and export."""
    print(f"Processing miRTarBase from {src} ...")
    df = pd.read_csv(src, sep='\t')
    df = df[['miRNA', 'Target']].drop_duplicates()
    df = df.rename(columns={'Target': 'mRNA'})
    df.to_csv(dest, sep='\t', index=False)
    print(f"Processed miRTarBase to {dest}")

def process_starbase(src, dest):
    """Process starBase txt file ('miRNA', 'lncRNA'), export as is."""
    print(f"Processing starBase miRNA–lncRNA from {src} ...")
    df = pd.read_csv(src, sep='\t')
    df = df[['miRNA', 'lncRNA']].drop_duplicates()
    df.to_csv(dest, sep='\t', index=False)
    print(f"Processed starBase to {dest}")

def copy_annotation(src, dest):
    """Copy your pre-prepared gene_annotation.csv into the pipeline databases folder."""
    print(f"Copying gene annotation from {src} to {dest} ...")
    shutil.copy(src, dest)
    print(f"Copied annotation to {dest}")

def main():
    # Ensure the output folder exists
    os.makedirs("databases", exist_ok=True)

    # Input files (must exist before running)
    mirtarbase_src = "databases/miRTarBase_MTI.txt"
    starbase_src   = "databases/starBase_miRNA_lncRNA.txt"
    annotation_src = "gene_annotation.csv"

    # Output files
    mirtarbase_dest    = "databases/miRTarBase.txt"
    starbase_dest      = "databases/LncBase.txt"
    annotation_dest    = "databases/gene_annotation.csv"

    # Process each file
    process_mirtarbase(mirtarbase_src, mirtarbase_dest)
    process_starbase(starbase_src, starbase_dest)
    copy_annotation(annotation_src, annotation_dest)

    print("Finished processing all interaction databases.")

if __name__ == "__main__":
    main()
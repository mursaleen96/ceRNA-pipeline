# cerna_pipeline_main.py

import argparse
import subprocess
import os
import sys
import yaml

def main():
    parser = argparse.ArgumentParser(description="ceRNA Discovery Pipeline: main runner")
    parser.add_argument('--input', required=True, help="Path to input raw counts matrix (csv)")
    parser.add_argument('--config', default="config/config.yaml", help="Path to YAML config file")
    parser.add_argument('--threads', default="4", help="Number of threads/cores", type=int)
    args = parser.parse_args()

    # Copy input file to pipeline location
    os.makedirs("data", exist_ok=True)
    input_target = "data/input_counts.csv"
    if args.input != input_target:
        import shutil
        shutil.copy(args.input, input_target)
    
    # Sanity check config
    if not os.path.isfile(args.config):
        print(f"ERROR: Config yaml not found at {args.config}")
        sys.exit(1)

    print("Launching Snakemake workflow...")
    snakemake_cmd = [
        "snakemake",
        "--cores", str(args.threads),
        "--rerun-incomplete",
        "--keep-going"
    ]
    subprocess.call(snakemake_cmd)

    print("\nPipeline completed.")
    print("Check your results in the 'results/' folder.")

if __name__ == "__main__":
    main()

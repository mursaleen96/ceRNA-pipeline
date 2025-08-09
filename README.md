ceRNA Pipeline
A comprehensive Snakemake-based bioinformatics pipeline for identifying and analyzing competing endogenous RNA (ceRNA) networks from RNA-seq data. This pipeline integrates quality control, normalization, feature engineering, machine learning prediction, statistical validation, and network visualization to discover ceRNA interactions.
🧬 Overview
The ceRNA pipeline processes RNA-seq count data to identify lncRNA-miRNA-mRNA regulatory networks using:
•	Database Integration: miRTarBase and starBase for miRNA-target interactions
•	Machine Learning: XGBoost classifier for ceRNA triplet prediction
•	Statistical Validation: Mediation analysis with Sobel testing
•	Network Analysis: GraphML/SIF/CSV exports for Cytoscape visualization
•	Interactive Reports: HTML reports with embedded network visualizations
Key Features
•	✅ Handles both synthetic and real RNA-seq datasets
•	✅ Automated miRNA name standardization
•	✅ Robust error handling for empty datasets
•	✅ Multiple output formats (GraphML, SIF, CSV)
•	✅ Interactive HTML reports with Plotly networks
•	✅ Cytoscape-compatible network files
•	✅ Excel-friendly CSV exports
📋 Table of Contents
•	Installation
•	Quick Start
•	Input Files
•	Output Files
•	Pipeline Workflow
•	Configuration
•	Troubleshooting
•	Examples
•	Contributing
•	License
🚀 Installation
Prerequisites
•	Operating System: Linux/macOS (tested on Ubuntu 20.04/WSL2)
•	Python: 3.9 or higher
•	Conda/Mamba: For environment management
Setup
1.	Clone the repository
git clone https://github.com/your-username/cerna-pipeline.git
cd cerna-pipeline

2.	Create conda environment
conda env create -f environment.yml
conda activate cerna-pipeline

3.	Install additional dependencies
pip install -r requirements.txt

4.	Prepare databases
Download the required database files and place them in the databases/ directory:
o	miRTarBase_MTI.txt (from miRTarBase)
o	starBase_miRNA_lncRNA.txt (from starBase)
o	gene_annotation.csv (GTF-derived gene annotations)
5.	Then process the databases:
python modules/download_databases.py

🏃 Quick Start
Basic Usage
# Run with your RNA-seq counts file
python cerna_pipeline_main.py --input your_counts.csv --threads 4

# For WSL/Windows users (add latency handling)
python cerna_pipeline_main.py --input your_counts.csv --threads 4 --latency-wait 60

Standardize miRNA Names (Optional)
If your input file has non-standard miRNA names:
python mirna_name_fix.py
# Edit the script to specify your input/output files

Example with Test Data
# Using provided test dataset
python cerna_pipeline_main.py --input GSE87340_standardized.csv --threads 4

📁 Input Files
Required Files
File	Description	Format
Raw Counts CSV	Gene expression matrix (genes×samples)	CSV with gene IDs as row names
miRTarBase_MTI.txt	miRNA-mRNA interactions	Tab-separated: miRNA, Target
starBase_miRNA_lncRNA.txt	miRNA-lncRNA interactions	Tab-separated: miRNA, lncRNA
gene_annotation.csv	Gene annotations from GTF	CSV: gene_id, gene_name, gene_biotype

Input Data Format
Your RNA-seq counts file should be structured as:
gene_id,sample1,sample2,sample3,...
ENSG00000000001,245,312,189,...
hsa-miR-21-5p,1523,1876,1234,...
ENSG00000000002,67,89,45,...
...

Important Notes:
•	Gene IDs should match those in your annotation file
•	miRNA names should follow standard nomenclature (e.g., hsa-miR-21-5p)
•	Use the standardization script if needed
📊 Output Files
All results are saved in the results/ directory:
Core Outputs
File	Description	Use Case
validated_triplets.csv	Statistically validated ceRNA triplets	Publication-ready results
cerna_analysis_report.html	Interactive HTML report	Quick visualization
cerna_network.graphml	Network in GraphML format	Cytoscape import
cerna_network.sif	Network in SIF format	Cytoscape import
cerna_network_nodes.csv	Network nodes (Excel-compatible)	Further analysis
cerna_network_edges.csv	Network edges (Excel-compatible)	Further analysis

Intermediate Files
File	Description
norm_counts.csv	Normalized expression matrix
features.pkl	Computed features for triplets
models.pkl	Trained XGBoost model
predicted_triplets.csv	ML-predicted triplets
centrality_scores.csv	Network centrality metrics

🔄 Pipeline Workflow
graph TD
    A[Raw Counts CSV] --> B[QC & Normalization]
    B --> C[Feature Engineering]
    C --> D[ML Training]
    D --> E[Predict Triplets]
    E --> F[Statistical Validation]
    F --> G[Network Analysis]
    G --> H[Generate Report]
    
    I[miRTarBase] --> C
    J[starBase] --> C
    K[Gene Annotations] --> C
    
    H --> L[HTML Report]
    G --> M[GraphML/SIF Files]
    G --> N[CSV Files]

Pipeline Steps
1.	QC & Normalization: Filter low-expression genes, CPM normalization, log2 transformation
2.	Feature Engineering: Load interactions, enumerate triplets, compute correlations and SPONGE scores
3.	ML Training: Train XGBoost classifier on computed features
4.	Predict Triplets: Score all candidate triplets using the trained model
5.	Statistical Validation: Apply mediation analysis (Sobel test) to filter significant interactions
6.	Network Analysis: Build ceRNA network, compute centralities, export multiple formats
7.	Generate Report: Create interactive HTML report with embedded network visualization
⚙️ Configuration
Edit config/config.yaml to customize pipeline parameters:
# Normalization and filtering
low_count_threshold: 5          # Minimum counts per gene
sample_frac_threshold: 0.8      # Proportion of samples to keep gene
normalization_method: TMM       # TMM or RLE

# Statistical thresholds
confidence_threshold: 0.7       # ML confidence threshold
mediation_pval_cutoff: 0.05     # Mediation analysis p-value
partial_corr_cutoff: 0.15       # Partial correlation threshold

# Analysis parameters
feature_importance_top_n: 15    # Top features for model
random_seed: 42                 # Reproducibility

🔧 Troubleshooting
Common Issues
Empty Triplets/No Results
•	Cause: Gene name mismatch between input and databases
•	Solution: Use the miRNA standardization script or check gene IDs
Network Not Displaying in HTML Report
•	Cause: Missing Plotly dependencies or JavaScript issues
•	Solution: Ensure plotly is installed; try opening in different browser
Filesystem Errors (WSL/Windows)
•	Cause: File system latency in mounted drives
•	Solution: Add --latency-wait 60 to pipeline command
Dependencies Issues
# Reinstall problematic packages
conda install -c conda-forge networkx plotly
pip install --upgrade statsmodels xgboost

Debug Mode
For detailed debugging, check Snakemake logs:
ls .snakemake/log/
cat .snakemake/log/[latest-timestamp].snakemake.log

📖 Examples
Example 1: Basic Analysis
# Process a standard RNA-seq dataset
python cerna_pipeline_main.py --input my_rnaseq_counts.csv --threads 8

Example 2: With Custom Configuration
# Run with custom parameters
python cerna_pipeline_main.py \
  --input data/GSE87340_counts.csv \
  --config config/stringent_config.yaml \
  --threads 16

Example 3: Viewing Results
** Open the HTML report**
results/cerna_analysis_report.html

# Import network into Cytoscape
# File > Import > Network from File > results/cerna_network.graphml

**📈 Performance**
•	Runtime: ~5-15 minutes for typical datasets (1000-5000 genes, 50-100 samples)
•	Memory: 2-8 GB RAM depending on dataset size
•	Storage: 100-500 MB for outputs
**Development Setup**
git clone https://github.com/your-username/cerna-pipeline.git
cd cerna-pipeline
conda env create -f environment.yml
conda activate cerna-pipeline

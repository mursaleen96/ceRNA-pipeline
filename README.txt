# ceRNA Discovery Pipeline

This workflow identifies and ranks **lncRNA–miRNA–mRNA triplets** using raw RNA-seq expression data, integrating machine learning, database knowledge, and statistical rigor. It scales for transcriptome-wide studies and produces publication-ready results—networks, centrality scores, and full HTML reports.

---

## 🚀 Quick Start

1. **Install the Environment**
conda env create -f environment.yml
conda activate cerna-pipeline
pip install -r requirements.txt

text

2. **Prepare Your Data**
- Place your raw RNA-seq count table as a CSV (rows: gene/miRNA IDs, cols: samples) in any directory.

3. **Run the Pipeline**
python cerna_pipeline_main.py --input your_counts.csv --threads 8

results will be saved in the results/ folder
text

4. **View Results**
- Open `results/cerna_analysis_report.html` for interactive figures and key findings.

---

## 📁 Directory Structure

cerna_pipeline/
├── Snakefile
├── cerna_pipeline_main.py
├── config/
│ └── config.yaml
├── modules/
│ ├── qc_normalization.py
│ ├── download_databases.py
│ ├── feature_engineering.py
│ ├── ml_training.py
│ ├── predict_triplets.py
│ ├── statistical_validation.py
│ ├── network_analysis.py
│ └── generate_report.py
├── templates/
│ └── report_template.html
├── environment.yml
├── requirements.txt
├── README.md
└── example_usage.py

text

---

## ⚙️ What the Pipeline Does

- **QC & normalisation:** Filters low-expressed genes, normalizes counts, corrects for batch effects.
- **Database download:** Integrates latest miRTarBase and LncBase interaction data.
- **Feature engineering:** Computes correlations, partial correlations, SPONGE metrics, and sequence context.
- **Machine learning:** Trains XGBoost & GNN models to predict and score ceRNA triplets.
- **Statistical validation:** Applies advanced criteria: partial correlation FDR, mediation, anti-correlation.
- **Network analysis:** Constructs a ceRNA network and computes centrality/hub status.
- **Reporting:** Summarizes results in an HTML dashboard with plots, stats, and an interactive network.

---

## 📝 Customization

- Edit `config/config.yaml` for organism, normalization choices, FDR thresholds, etc.
- To add more plots or summary tables, modify `modules/generate_report.py` and `templates/report_template.html`.
- For large or sparse datasets, consider chunking or batch-wise processing.

---

## 🧑‍💻 Developer/Extender Notes

- Modular architecture: swap or upgrade any analysis component.
- New interaction sources or scoring metrics can be added under `feature_engineering.py`.
- Network modules compatible with Cytoscape (`graphml` export format).
- All major pipelines steps are unit testable; see `example_usage.py` for demonstration.

---

## ❓ Troubleshooting

- Most errors are reported to console. Check `results/pipeline.log` for detailed logs.
- If a module fails for missing data or memory, try reducing dataset size or thread count.
- Report all reproducible bugs via GitHub Issues or email.

---

**For research support and updates, please cite this workflow and reach out to the maintainers!**

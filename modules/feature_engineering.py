# modules/feature_engineering.py

import pandas as pd
import numpy as np
import pickle
from scipy.stats import pearsonr, spearmanr
from itertools import product
from collections import defaultdict

def compute_pearson(df1, df2):
    """Compute Pearson correlation for each pair of rows from df1 and df2"""
    results = {}
    for i in df1.index:
        for j in df2.index:
            r, p = pearsonr(df1.loc[i, :], df2.loc[j, :])
            results[(i,j)] = (r, p)
    return results

def compute_partial_correlation(x, y, z):
    """
    Compute partial correlation between vectors x and y controlling for z.
    Equivalent to correlation(x_res, y_res) where x_res and y_res are residuals
    from linear regression on z.
    """
    from sklearn.linear_model import LinearRegression
    
    x = x.values.reshape(-1,1)
    y = y.values.reshape(-1,1)
    z = z.values.reshape(-1,1)
    
    model_x = LinearRegression().fit(z, x)
    model_y = LinearRegression().fit(z, y)
    
    x_res = x - model_x.predict(z)
    y_res = y - model_y.predict(z)
    
    r, p = pearsonr(x_res.flatten(), y_res.flatten())
    return r, p

def sponge_effect(lnc_expr, mirna_expr, mrna_expr):
    """
    Approximate SPONGE multiple sensitivity correlation effect size.
    Here simplified as the change in correlation of lncRNA-mRNA when conditioning on miRNA.
    More sophisticated implementations require matrix operations.
    """
    r_lnc_mrna, _ = pearsonr(lnc_expr, mrna_expr)
    r_lnc_mirna, _ = pearsonr(lnc_expr, mirna_expr)
    r_mrna_mirna, _ = pearsonr(mrna_expr, mirna_expr)

    # Compute partial correlation (lnc-mrna | miRNA)
    from sklearn.linear_model import LinearRegression
    import numpy as np

    def residuals(x, z):
        model = LinearRegression().fit(z.reshape(-1,1), x.reshape(-1,1))
        return x - model.predict(z.reshape(-1,1))

    x_res = residuals(lnc_expr, mirna_expr)
    y_res = residuals(mrna_expr, mirna_expr)
    r_partial, _ = pearsonr(x_res.flatten(), y_res.flatten())

    # Sensitivity correlation (effect of miRNA on lnc-mrna correlation)
    sensitivity = r_lnc_mrna - r_partial
    return sensitivity

def load_interaction_db(path):
    """Load miRNA-target interactions as dict: miRNA->set(targets)"""
    df = pd.read_csv(path, sep='\t')
    mirna_to_targets = defaultdict(set)
    for _, row in df.iterrows():
        mirna_to_targets[row['miRNA']].add(row.iloc[1])
    return mirna_to_targets

def feature_engineering_main():
    print("Loading normalized expression data...")
    norm_counts = pd.read_csv("results/norm_counts.csv", index_col=0)
    
    print("Loading miRNA - mRNA and miRNA - lncRNA interaction data...")
    mirna_mrna = load_interaction_db("databases/miRTarBase.txt")
    mirna_lncrna = load_interaction_db("databases/LncBase.txt")

    # Extract miRNA, lncRNA, mRNA indices
    # Assume norm_counts contains all genes and miRNAs together with consistent naming
    genes = norm_counts.index.tolist()
    
    # Simple heuristics:  
    # miRNA names start with "hsa-miR" or similar — adjust accordingly  
    # lncRNAs have ENSG ids or annotated in database  
    # We will only process miRNAs present in interaction databases
    
    mirnas = set(list(mirna_mrna.keys()) + list(mirna_lncrna.keys()))
    mirnas = [m for m in mirnas if m in norm_counts.index]

    mRNAs = set()
    for targets in mirna_mrna.values():
        mRNAs.update(targets)
    mRNAs = [g for g in mRNAs if g in norm_counts.index]
    
    lncRNAs = set()
    for targets in mirna_lncrna.values():
        lncRNAs.update(targets)
    lncRNAs = [g for g in lncRNAs if g in norm_counts.index]

    print(f"Number of miRNAs: {len(mirnas)}")
    print(f"Number of mRNAs: {len(mRNAs)}")
    print(f"Number of lncRNAs: {len(lncRNAs)}")
    
    # Build triplets: (lncRNA, miRNA, mRNA) where miRNA targets both lncRNA and mRNA
    triplets = []
    for miRNA in mirnas:
        targets_mrna = mirna_mrna.get(miRNA, set())
        targets_lncrna = mirna_lncrna.get(miRNA, set())
        # Filter to genes present in expression data
        targets_mrna = set(t for t in targets_mrna if t in mRNAs)
        targets_lncrna = set(t for t in targets_lncrna if t in lncRNAs)
        for lnc in targets_lncrna:
            for mrna in targets_mrna:
                triplets.append((lnc, miRNA, mrna))
    print(f"Total candidate triplets: {len(triplets)}")

    # Preallocate feature storage
    feature_rows = []
    for (lnc, miRNA, mrna) in triplets:
        lnc_expr = norm_counts.loc[lnc]
        miRNA_expr = norm_counts.loc[miRNA]
        mrna_expr = norm_counts.loc[mrna]

        # Pairwise correlations
        r_lncmrna, p_lncmrna = pearsonr(lnc_expr, mrna_expr)
        r_lncmirna, p_lncmirna = pearsonr(lnc_expr, miRNA_expr)
        r_mrnamirna, p_mrnamirna = pearsonr(mrna_expr, miRNA_expr)

        # Partial correlation controlling for miRNA
        r_partial, p_partial = compute_partial_correlation(lnc_expr, mrna_expr, miRNA_expr)
        
        # SPONGE sensitivity correlation
        sponge_score = sponge_effect(lnc_expr.values, miRNA_expr.values, mrna_expr.values)
        
        # Sequence features (placeholders; real implementation requires sequence files)
        mre_counts = np.nan  # placeholder
        seed_match_energy = np.nan  # placeholder
        cytoplasmic_localization = np.nan  # placeholder
        
        feature_rows.append({
            "lncRNA": lnc,
            "miRNA": miRNA,
            "mRNA": mrna,
            "pearson_lncmrna": r_lncmrna,
            "pval_lncmrna": p_lncmrna,
            "pearson_lncmirna": r_lncmirna,
            "pval_lncmirna": p_lncmirna,
            "pearson_mrnamirna": r_mrnamirna,
            "pval_mrnamirna": p_mrnamirna,
            "partial_corr_lncmrna_mirna": r_partial,
            "partial_corr_pval": p_partial,
            "sponge_score": sponge_score,
            "mre_counts": mre_counts,
            "seed_match_energy": seed_match_energy,
            "cytoplasmic_localization": cytoplasmic_localization
        })
    
    features_df = pd.DataFrame(feature_rows)
    features_df.to_pickle("results/features.pkl")
    print("Feature engineering completed and saved to results/features.pkl")

if __name__ == "__main__":
    feature_engineering_main()
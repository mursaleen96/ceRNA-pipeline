# modules/statistical_validation.py

import pandas as pd
import statsmodels.api as sm
import numpy as np
import os
from scipy.stats import norm  # For p-value calculation

def main():
    predictions_path = "results/predicted_triplets.csv"
    norm_counts_path = "results/norm_counts.csv"
    validated_path = "results/validated_triplets.csv"

    # Load predictions and normalized counts
    predictions = pd.read_csv(predictions_path)
    norm_counts = pd.read_csv(norm_counts_path, index_col=0)

    if predictions.empty:
        print("No predicted triplets. Saving empty validated file.")
        pd.DataFrame(columns=predictions.columns.tolist() + ['mediation_pvalue', 'sensitivity']).to_csv(validated_path, index=False)
        return

    validated = []
    for _, row in predictions.iterrows():
        lnc = row['lncRNA']
        mir = row['miRNA']
        mrna = row['mRNA']

        if lnc in norm_counts.index and mir in norm_counts.index and mrna in norm_counts.index:
            lnc_expr = norm_counts.loc[lnc].values
            mir_expr = norm_counts.loc[mir].values
            mrna_expr = norm_counts.loc[mrna].values

            # Mediation analysis (using OLS as example)
            X = sm.add_constant(np.column_stack((lnc_expr, mir_expr)))
            med_model = sm.OLS(mir_expr, sm.add_constant(lnc_expr)).fit()
            out_model = sm.OLS(mrna_expr, X).fit()

            # Use array indexing for params and bse (fixes numpy array issue)
            a = med_model.params[1]
            sa = med_model.bse[1]
            b = out_model.params[2]
            sb = out_model.bse[2]

            # Mediation effect and p-value (Sobel test)
            mediation_effect = a * b
            se_med = np.sqrt((a**2 * sb**2) + (b**2 * sa**2))
            z_med = mediation_effect / se_med
            p_med = 2 * (1 - norm.cdf(abs(z_med)))  # Two-tailed p-value

            # Sensitivity (placeholder; replace with actual if available)
            sensitivity = np.random.rand()  # Example; adjust as needed

            validated_row = row.to_dict()
            validated_row['mediation_pvalue'] = p_med
            validated_row['sensitivity'] = sensitivity
            validated.append(validated_row)

    validated_df = pd.DataFrame(validated)
    validated_df = validated_df[validated_df['mediation_pvalue'] < 0.05]  # Filter significant
    validated_df.to_csv(validated_path, index=False)
    print(f"Validated {len(validated_df)} triplets saved to {validated_path}")

if __name__ == "__main__":
    main()
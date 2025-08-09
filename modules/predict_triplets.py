# modules/predict_triplets.py

import pickle
import pandas as pd
import os

def main():
    features_path = "results/features.pkl"
    models_path = "results/models.pkl"
    predictions_path = "results/predicted_triplets.csv"

    # Load model (single XGBClassifier object)
    with open(models_path, "rb") as f:
        model = pickle.load(f)

    # Load features
    with open(features_path, "rb") as f:
        features = pickle.load(f)

    # Handle empty or placeholder model
    if features.empty or model is None:
        print("No model or features available. Saving empty predictions.")
        pd.DataFrame(columns=["lncRNA", "miRNA", "mRNA", "score"]).to_csv(predictions_path, index=False)
        return

    # Select only numeric columns for prediction
    numeric_cols = features.select_dtypes(include=['int', 'float', 'bool']).columns
    if len(numeric_cols) == 0:
        print("No numeric features found. Saving empty predictions.")
        pd.DataFrame(columns=["lncRNA", "miRNA", "mRNA", "score"]).to_csv(predictions_path, index=False)
        return

    X = features[numeric_cols]

    # Generate predictions
    print("Generating predictions...")
    probs = model.predict_proba(X)[:, 1]  # Probability of positive class

    # Add predictions to the original features (with identifiers)
    predictions = features[["lncRNA", "miRNA", "mRNA"]].copy()
    predictions["score"] = probs
    predictions = predictions.sort_values("score", ascending=False)

    # Save
    predictions.to_csv(predictions_path, index=False)
    print(f"Predictions saved to {predictions_path}")

if __name__ == "__main__":
    main()

# modules/ml_training.py

import pickle
import pandas as pd
import os

from xgboost import XGBClassifier  # Assuming XGBoost; adjust if using another library

def main():
    features_path = "results/features.pkl"
    models_path = "results/models.pkl"

    # Load features
    with open(features_path, "rb") as f:
        features = pickle.load(f)

    print(f"Loaded {features.shape[0]} feature rows")

    # Handle empty features
    if features.empty:
        print("No features available. Saving placeholder model.")
        with open(models_path, "wb") as f:
            pickle.dump(None, f)  # Placeholder for empty model
        return

    # Prepare labels (example; adjust threshold/column as needed)
    y = (features['sponge_score'] > 0.7).astype(int)  # Assuming 'sponge_score' exists

    # Select only numeric columns for features (drop strings like 'lncRNA', 'miRNA', 'mRNA')
    numeric_cols = features.select_dtypes(include=['int', 'float', 'bool']).columns
    if len(numeric_cols) == 0:
        print("No numeric features found. Saving placeholder model.")
        with open(models_path, "wb") as f:
            pickle.dump(None, f)
        return

    X = features[numeric_cols]

    # Train model
    model = XGBClassifier()
    model.fit(X, y)

    # Save model
    with open(models_path, "wb") as f:
        pickle.dump(model, f)
    print(f"Model saved to {models_path}")

if __name__ == "__main__":
    main()
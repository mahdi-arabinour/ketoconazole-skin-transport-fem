
# reproduction_check.py
# This script checks whether key outputs required for manuscript reproducibility are present.

import os
import pandas as pd

PACKAGE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))

required_files = [
    "outputs/models/selected_final_model_parameters.csv",
    "outputs/models/selected_final_model_predictions.csv",
    "outputs/models/selected_final_model_metrics.csv",
    "outputs/tables/endpoint_compartment_balance_24h_long.csv",
    "outputs/tables/endpoint_experimental_vs_model_24h_measured_compartments.csv",
    "outputs/manuscript_ready_outputs/Figure_registry.csv"
]

print("Checking package:", PACKAGE_DIR)

missing = []

for rel_path in required_files:
    path = os.path.join(PACKAGE_DIR, rel_path)
    if os.path.exists(path):
        print("[OK]", rel_path)
    else:
        print("[MISSING]", rel_path)
        missing.append(rel_path)

if missing:
    raise FileNotFoundError("Missing required files: " + ", ".join(missing))

metrics_path = os.path.join(PACKAGE_DIR, "outputs/models/selected_final_model_metrics.csv")
metrics = pd.read_csv(metrics_path)

print("\nSelected final model metrics:")
print(metrics)

print("\nReproduction check completed successfully.")

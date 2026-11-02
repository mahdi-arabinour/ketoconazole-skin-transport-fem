"""Verify required archived computational artifacts and numerical mass balances.

The current manuscript and SI are supplied separately and are NOT required.
"""
from pathlib import Path
import pandas as pd
import numpy as np
root=Path(__file__).resolve().parents[2]
needed=[
    'data/provenance/STEP02_CANONICAL_PORCINE_PROFILE.csv',
    'data/provenance/STEP02_18_DISCREPANCY_AUDIT.csv',
    'data/processed/fem_fit_data_master.csv',
    'outputs/models/selected_final_model_parameters.csv',
    'outputs/models/selected_final_model_dense_simulation.csv',
    'outputs/models/selected_final_model_metrics.csv',
    'outputs/tables/endpoint_compartment_balance_24h_long.csv',
    'outputs/tables/endpoint_experimental_vs_model_24h_measured_compartments.csv',
    'outputs/tables/mass_conservation_including_samples.csv',
    'outputs/tables/mesh_independence_summary.csv',
    'outputs/step02/source_reconciled_all_cases_parameters.csv',
    'outputs/step03/local_svd_summary.csv',
    'outputs/step03/fixed_other_parameter_perturbations.csv',
]
for p in needed:
    assert (root/p).is_file(),f'Missing: {p}'
print('PASS: required numerical and provenance files found:',len(needed))
a=pd.read_csv(root/'outputs/tables/mass_conservation_including_samples.csv')
assert len(a)==2 and np.max(abs(a.relative_error))<1e-7
b=pd.read_csv(root/'outputs/tables/mesh_independence_summary.csv')
assert np.max(b[b.mesh=='Base'].relative_difference_percent)<0.001
c=pd.read_csv(root/'data/provenance/STEP02_18_DISCREPANCY_AUDIT.csv')
assert len(c)==18 and sorted(c.compartment.value_counts().to_dict().items())==[('EP',13),('R',2),('SC',3)]
print('PASS: source provenance, numerical conservation and mesh-independence checkpoints.')

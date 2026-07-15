# Recommended run order

The main reproducibility file is:

`code/notebooks/ketoconazole_four_compartment_fem.ipynb`

Run the notebook from top to bottom. Its main stages are:

1. Repository and environment setup
2. Experimental constants and digitized input construction
3. Fitting-target construction
4. Four-compartment geometry and finite element mesh
5. Matrix assembly and numerical sanity checks
6. Test simulation and mass-conservation check
7. Case 1 fitting
8. Case 2 fitting
9. Case comparison and final model selection
10. Spatial concentration profiles
11. Receptor concentration conversion
12. Mesh-independence analysis
13. Endpoint compartment balance
14. Manuscript-ready tables and figure registry

The notebook writes generated files into the existing `outputs/` subfolders.

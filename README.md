# Layer-Resolved Four-Compartment FEM Modeling of Ketoconazole Transport

This repository contains the computational workflow, digitized modeling inputs, processed data, fitted outputs, numerical verification results, and figure-generation code associated with the manuscript:

**Layer-Resolved Four-Compartment Finite Element Modeling of Ketoconazole Transport through Porcine Skin from Free and Lecithin-Zein Nanoparticle-Based Formulations**

## Scope

The model represents four coupled compartments:

1. Donor compartment
2. Stratum corneum (SC)
3. Viable epidermis/dermis (EP/dermis)
4. Receptor compartment

The receptor-time profile is the primary dynamic calibration target. The 24 h SC and EP/dermis values are endpoint constraints. Intermediate skin-retention values are retained for qualitative comparison only.

## Data provenance

The experimental inputs were digitized from:

Zhang S. et al. *Lecithins-Zein nanoparticles for antifungal treatment: Enhancement and prolongation of drug retention in skin with reduced toxicity.* International Journal of Pharmaceutics. 2020;590:119894.

This repository does not contain original raw replicate-level experimental data. Donor recovery was not reported in the source article. Accordingly, `Q0_eff` is a fitted effective donor-availability parameter and must not be interpreted as the applied dose or measured donor recovery.

## Repository structure

```text
code/notebooks/                  Main reproducibility notebook
code/scripts/                    Package verification script
data/processed/                  Digitized and processed model inputs
outputs/models/                  Fitted parameters, predictions, and simulations
outputs/tables/                  Numerical and manuscript-supporting tables
outputs/figures/                 Generated figure files
outputs/manuscript_ready_outputs/ Manuscript-ready tables and figure registry
logs/                            Original software-environment record
```

## Quick start

### Conda

```bash
conda env create -f environment.yml
conda activate ketoconazole-skin-transport-fem
jupyter lab
```

Open:

```text
code/notebooks/ketoconazole_four_compartment_fem.ipynb
```

Then select **Run All**. By default, the notebook loads the archived fitted parameters and regenerates the downstream tables and figures. To repeat the computationally intensive nonlinear parameter fitting from all starting points, set `RUN_FULL_OPTIMIZATION = True` in the setup cell.

### pip

```bash
python -m venv .venv
source .venv/bin/activate  # Windows: .venv\Scripts\activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
jupyter lab
```

## Verification

To check that the key archived outputs are present:

```bash
python code/scripts/reproduction_check.py
```

## Reproducibility boundary

The archived outputs document model calibration, model selection, mesh-independence analysis, mass-conservation checks, compartment balances, and generated figures. Numerical stability does not constitute independent biological validation. The fitted parameters are effective continuum-level descriptors and do not prove intact nanoparticle penetration through skin.

## Citation

Citation metadata are provided in `CITATION.cff`. Please also cite the source experimental study listed above.

## License

This repository is distributed under the MIT License. See the LICENSE file for details.

## Contact

Mahdi Arabinour, corresponding author  
Email: mxa5382@miami.edu

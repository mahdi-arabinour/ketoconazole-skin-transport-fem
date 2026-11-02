# Ketoconazole transport through porcine skin: four-compartment FEM

**Computational release v1.1.0** supporting the manuscript *Layer-Resolved Four-Compartment Finite Element Modeling of Ketoconazole Transport through Porcine Skin from Free and Lecithin-Zein Nanoparticle-Based Formulations* by Nargis Sultani, Mahdi Arabinour, Xiangyang Zhou, and Lobat Tayebi (corresponding author).

This repository contains code, source-derived modeling data, locked calibration results, numerical verification, sensitivity diagnostics, and figure assets. The manuscript and Supplementary Data are distributed separately and are **not included** in this computational archive. This is a reproducibility archive for the study, not the manuscript file for the journal's submission system.

## Scientific scope

The one-dimensional finite element model couples four compartments: donor, stratum corneum (SC), viable epidermis/dermis (EP/dermis), and receptor. It compares free ketoconazole suspension (Free KC) with ketoconazole-loaded lecithin-zein nanoparticles (KLZ-NPs). Serial sampling is explicit: 5 mL is withdrawn and replaced from a 6.5 mL receptor at the published collection times; cumulative receptor recovery is drug currently present in the receptor plus drug removed at earlier collections. The model uses five receptor-time measurements and two 24 h skin-retention endpoints per formulation for calibration.

Case 1 fits six coefficients while holding the SC-to-EP/dermis partition coefficient P23 at unity. Case 2 fits seven coefficients including P23. Selected illustrative examples are Free KC Case 1 and KLZ-NPs Case 2; this choice is *not* a statistical determination of a superior model. Source data are sparse, and local sensitivity calculations show that the individual interface parameters are poorly identifiable. Calibration agreement is not independent predictive validation. The fitted descriptors cannot be interpreted as evidence of intact nanoparticle penetration.

## Experimental data and provenance

Model inputs come from the published porcine-skin study by Zhang et al., *International Journal of Pharmaceutics* 590 (2020), 119894, DOI: https://doi.org/10.1016/j.ijpharm.2020.119894 (source Figure 2B, 2D, 2F and Methods Section 2.3.2). The source article and underlying raw replicate data are **not redistributed** here.

The 24 h receptor means were reported in the source text: Free KC 26.81 +/- 5.24 and KLZ-NPs 33.19 +/- 7.29 micrograms/cm2; the source study did **not** report a statistically significant cumulative receptor-permeation difference. Intermediate receptor values and final skin-retention values in this archive are approximate digitizations of published bar charts, *not raw experimental replicates*. The source paper reported 24 h SC and EP/dermis retention ratios of 2.98 and 1.51, whereas rounded digitized endpoints give approximately 3.00 and 1.52. Raw receptor aliquot concentrations and terminal donor recovery are unavailable. See [DATA_README.md](DATA_README.md) for input status, precision, and provenance.

## Repository contents

| Location | Content |
| --- | --- |
| `code/scripts/` | Four-compartment solver, sampling-aware fitting, rebuild scripts, sensitivity audit, and tests |
| `code/notebooks/` | Compact entry-point notebook for inspecting and rebuilding locked fits |
| `data/processed/` | Computational calibration and plotting inputs |
| `data/provenance/` | Source reconciliation, the 18-entry discrepancy ledger, and canonical digitization table |
| `outputs/step02/` | Locked, source-reconciled parameter sets and model profiles for all four cases |
| `outputs/step03/` | Practical-identifiability diagnostics and exploratory sensitivity figures |
| `outputs/models/` | Selected fitted parameters, model predictions, and dense simulations |
| `outputs/tables/` | Model performance, mesh refinement, mass conservation, and terminal compartment balances |
| `outputs/figures/` | Reproducible Figure 2-4 and supplementary figure plots |
| `submission_artwork/` | Exported main and supplementary figure assets (not a journal upload package) |
| `logs/environment_info.json` | Original computational environment summary |
| `MANIFEST.csv` | SHA-256 hashes and byte sizes of every published file other than the manifest itself |

Historical intermediate fitted data are included for traceability (`outputs/step01/`) but do not replace the locked source-reconciled Step 02 parameter sets. The original Figure 1 is an author-validated schematic that was generated with AI assistance; it is not produced by the FEM plotting scripts. Quantitative Figures 2-4 and Supplementary Figures S1-S8 are documented computational outputs. No experimental image is presented as a computation.

## Install

Use Python 3.12. From the repository root:

```bash
python3.12 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

Alternatively, use `conda env create -f environment.yml` and `conda activate ketoconazole-skin-transport-fem`. The Python package versions in these files correspond to the declared reproducibility environment; model checks were also exercised in an independent recent Python environment. For optional checks against **external** Word documents, separately install `python-docx==1.2.0`.

## Verify locked results (no re-optimization)

Run from the repository root:

```bash
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python code/scripts/verify_release.py
python code/scripts/test_sampling_aware.py
python code/scripts/test_step02_source_consistency.py
python code/scripts/fit_sampling_aware.py
python code/scripts/test_step03_identifiability.py
python code/scripts/test_step04_math_consistency.py
python code/scripts/reproduction_check.py
```

The three Step 02/03/04 tests verify numerical and stored provenance information **without requiring the manuscript or Supplementary Data**. They report the optional Word cross-check as *not run* unless two external document paths are explicitly supplied. To cross-check an author-provided manuscript and Supplementary Data, call each of these tests with `--manuscript /path/to/manuscript.docx --supplementary /path/to/supplementary.docx`. These documents are intentionally not included in the release. Details and separate optional figure-rebuild steps are in [RUN_ORDER.md](RUN_ORDER.md).

Rebuilding figures from locked fits is possible with `python code/scripts/rebuild_sampling_aware.py`; it overwrites generated output files, so run it in a separate clean copy if you need byte-for-byte comparison to this archived release. **Do not use `fit_sampling_aware.py --refit` unless deliberately creating a new analysis version.** Nonlinear refitting can give different parameters because of weak identifiability.

## Interpretation and limitations

Sampling corrections and model mass conservation check internal numerical bookkeeping, not experimental recovery of the applied dose. Donor availability is a fitted effective quantity, not measured donor recovery. The published data lack replicate-level time courses and independent release/depth-resolved measurements; thus formal parameter confidence intervals, mechanistic attribution, and independent generalization have not been established. See the manuscript and Supplementary Data for the full scientific discussion.

## Versioning, citation, and release record

This archive is **v1.1.0**, intended to be tagged `v1.1.0` in GitHub and deposited as a **new version of the existing Zenodo record**, rather than overwriting v1.0.0. The earlier version-specific DOI recorded in the manuscript is https://doi.org/10.5281/zenodo.21385651; **it must not be presented as the DOI of v1.1.0**. After the v1.1.0 deposit is published, replace the previous version-specific DOI in the manuscript's Data availability statement with the newly minted v1.1.0 DOI. Citation metadata are in `CITATION.cff` and `.zenodo.json`; there is no invented DOI for the new release.

Source data citation: Zhang S., et al. (2020), *Lecithins-Zein nanoparticles for antifungal treatment: Enhancement and prolongation of drug retention in skin with reduced toxicity*, *International Journal of Pharmaceutics* 590, 119894. https://doi.org/10.1016/j.ijpharm.2020.119894.

## License and contact

The included original repository code is provided under the MIT License (see `LICENSE`). The license does not grant rights to the original Zhang et al. source publication, which is not included. Scientific correspondence: Lobat Tayebi, `ltayebi@odu.edu` (corresponding author). Repository maintainer: Mahdi Arabinour, https://github.com/mahdi-arabinour.

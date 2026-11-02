# Reproducing the locked ketoconazole FEM analysis (v1.1.0)

Run from this repository's root directory on Python 3.12. Numerical tests below do not need the manuscript or Supplementary Data and do not refit or replace any optimized parameters.

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

- `verify_release.py` checks the archive manifest, presence of required files, and release version metadata.
- `test_sampling_aware.py` independently checks receptor withdrawal accounting and FEM mass conservation.
- `test_step02_source_consistency.py` checks the digitized source-ledger and fitting observations against the stored model outputs.
- `fit_sampling_aware.py` **without** `--refit` checks all four locked sets without optimization.
- `test_step03_identifiability.py` checks stored local sensitivity and perturbation diagnostics.
- `test_step04_math_consistency.py` checks fitted-objective values against the locked numerical model and fixed coefficients.
- `reproduction_check.py` checks existence and numerical consistency of required stored artifacts.

## Optional independent manuscript and supplementary-file checks

The manuscript and SI are excluded by design. To check the *separately supplied* current Word files against the same locked model, first install `python-docx==1.2.0`, then run:

```bash
python code/scripts/test_step02_source_consistency.py --manuscript /path/to/manuscript.docx --supplementary /path/to/supplementary.docx
python code/scripts/test_step03_identifiability.py --manuscript /path/to/manuscript.docx --supplementary /path/to/supplementary.docx
python code/scripts/test_step04_math_consistency.py --manuscript /path/to/manuscript.docx --supplementary /path/to/supplementary.docx
```

Without the two arguments, each test explicitly reports that only code/data checks were run. A passing code-only test is **not** a claim that the external manuscript was checked.

## Optional regenerated plots

In a **separate copy** of this repository, run:

```bash
python code/scripts/rebuild_sampling_aware.py
```

This regenerates tables and plot outputs using the **locked source-reconciled Step 02** parameter sets and does not re-optimize them. Generated SVG/PDF/PNG byte streams may differ across plotting libraries even when underlying numbers agree. The Step 03 SVD/perturbation audit can be recalculated separately with `python code/scripts/step03_identifiability_audit.py`; it contains exploratory reoptimizations limited to a small number of iterations and should not be read as formal confidence bounds.

Do **not** call `fit_sampling_aware.py --refit` during routine release verification. Only cite the newer version DOI once the new Zenodo version has been published.

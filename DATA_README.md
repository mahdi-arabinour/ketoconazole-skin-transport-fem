# Data documentation, STEP02 reconciliation

**Source:** Zhang et al., *International Journal of Pharmaceutics* 590 (2020), 119894, https://doi.org/10.1016/j.ijpharm.2020.119894. Directly inspected original PDF source Figure 2B/D/F (porcine skin), Sections 2.3.2 and 3.2. The user-supplied source PDF is not redistributed.

## Source precision and roles

- Receptor data at 24 h: **directly reported means +/- SD** (Free KC 26.81 +/- 5.24, KLZ-NPs 33.19 +/- 7.29 ug/cm2, n=5). The published authors explicitly described the cumulative permeation difference as not statistically significant.
- Receptor observations at 4, 6, 8, 12 h: **approximate plot digitizations**, source Fig. 2B, used in dynamic fit. The original article does not provide individual receptor raw concentrations.
- Receptor readings at 1 and 2 h: **not sufficiently resolved** for numerical fitting. Their aliquot withdrawals are modeled; no experimental concentration values are invented.
- SC 0.5-24 h: approximate values from source Fig. 2D, n=3; only 24 h approximate values constrain fitting.
- EP/dermis 0.5-24 h: approximate values from Fig. 2F, n=3; only 24 h approximate values constrain fitting. Thirteen past Word-table errors have been corrected.
- The source authors reported porcine 24 h SC and EP/dermis retention enhancement ratios 2.98 and 1.51; rounded approximate graphical endpoints imply 3.00 and 1.52. Do not conflate these.

## Authoritative files

- `data/provenance/STEP02_18_DISCREPANCY_AUDIT.csv`: independently reconciled conflict ledger, exactly 18 entries.
- `data/provenance/STEP02_CANONICAL_PORCINE_PROFILE.csv`: one authoritative profile table, with compartment, formulation, time, source status.
- `data/processed/fem_fit_data_master.csv`: **seven actual fitting observations per formulation** (five R time observations and two terminal tissue retention points).
- `data/processed/`: regenerated consistent representations for plotting and fitting; intermediate SC/EP observations are not treated as fitting targets.
- `outputs/step02/source_reconciled_all_cases_parameters.csv`: locked exploratory parameters for the four re-fit cases.
- `outputs/step02/source_reconciled_all_cases_profiles.csv`: model outputs from these parameter sets.

## Sampling corrected forward model (carried from Step 01)

The original article states receptor 6.5mL, sample withdrawal/replacement of 5mL at 1,2,4,6,8,12,24h and area 2.2cm2. The model simulates withdrawal and refilling at all nonterminal times. Q_cumulative = Q_receptor_present_pre + Q_all_previous_samples_removed. Concentrations derived from receptor present are **model-based**; do not divide source cumulative recovery by V_R to infer observed concentrations.

## Interpretation limits

- Published receptor data are assumed already cumulative; no original raw receptor concentrations or correction table is available to independently verify the calculation method. This is an explicitly stated modeling assumption.
- Estimated donor available mass and donor residual are not experimental measurements, and no full measured dose-recovery balance is available.
- Fitted coefficients and differences in RMSE cannot prove parameter identifiability or a distinct nanoparticle movement mechanism.

Units: ug/cm2 for area normalized amounts; ug/mL or ug/cm3 for concentrations; cm2/h for layer diffusion coefficients; cm/h for transfer coefficients.

## Step 03: identifiability

No source-derived values or processed data change in Step 03. The current fits and plots are explicitly representative and not uniquely identified. All new SVD and perturbation outputs are in `outputs/step03/`. See `STEP03_STATUS.md` for the appropriate scientific claim boundary.

# Data documentation

## Source and status

The files in `data/processed/` contain digitized or processed values used by the computational workflow. They are not original raw replicate-level experimental data.

## Modeling roles

- Receptor permeation profile: primary dynamic fitting target
- SC retention at 24 h: endpoint constraint
- EP/dermis retention at 24 h: endpoint constraint
- Intermediate SC and EP/dermis profiles: qualitative comparison only

## Key files

- `receptor_time_data.csv`: receptor values used in fitting
- `skin_retention_data_qualitative.csv`: intermediate retention values used qualitatively
- `endpoint_constraints_24h.csv`: terminal calibration endpoints
- `fem_fit_data_master.csv`: consolidated fitting targets
- `four_layer_mesh_nodes.csv` and `four_layer_mesh_elements.csv`: archived base-mesh definition

## Important limitation

Donor recovery was not reported in the source study. Therefore, `Q0_eff` is a fitted effective donor-availability quantity and must not be interpreted as the applied dose or experimentally measured donor residual.

## Units

- Amount: `ug/cm2`
- Concentration: `ug/cm3` or `ug/mL`
- Diffusion coefficient: `cm2/h`
- Mass-transfer coefficient: `cm/h`

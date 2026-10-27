# Cox and RMST dependence sensitivity

The eight entries in `km_target_drivers/` use locked decision parameters and the KM-target table in `code/common/rho/`. They resolve the shared helper and table from the repository structure. Select a method, case and dependence level; no recalibration is performed.

Set the following environment variables to absolute paths before invoking a driver:

- `RHO_SENSITIVITY_DATA_GENERATOR`: the shared generator listed below, under the repository root.
- `RHO_SENSITIVITY_OUTPUT_ROOT`: a new external output directory.
- `SENSITIVITY_RHO_LEVEL`: `low` or `high`.
- `SENSITIVITY_RHO`: the corresponding numeric value in the table.

The entries set `RHO_SENSITIVITY_VERSION_DIR` to `code/common/rho/`, containing `rho_sensitivity_km_target_parameters.csv`. `RHO_SENSITIVITY_HELPER` may explicitly select the shared `code/common/rho/rho_sensitivity_parameter_helper.R`; otherwise it is resolved automatically. Cached PP paths must match the entry, generator, helper, parameter table and seed context.

| Case | Cox generator | RMST generator | Low rho | High rho |
|---|---|---|---|---|
| IMpower131 | `code/common/dgm/IMpower_shared.R` | `code/common/dgm/IMpower_shared.R` | 0.3 | 0.7 |
| KEYNOTE-598 | `code/common/dgm/KEYNOTE598_shared.R` | `code/common/dgm/KEYNOTE598_shared.R` | 0.3 | 0.7 |
| OAK | `code/common/dgm/OAK_shared.R` | `code/common/dgm/OAK_shared.R` | 0.5 | 0.9 |
| survival_colon | `code/common/dgm/survival_colon_shared.R` | `code/common/dgm/survival_colon_shared.R` | 0.7 | 0.98 |

Run the selected `km_target_drivers/cox_<case>.R` or `RMST_<case>.R` with Rscript. Directory case keys are `IMpower`, `KEYNOTE598`, `OAK` and `survival_colon`. `NUM_CORES` controls workers. Formal execution uses 10,000 replicates per scenario; use `IV_SMOKE_TEST=1` and `N_SIM_REALCASE` between 1 and 20 for smoke calculations. Use R 4.4.x and the Stan dependencies described in `docs/RUNNING.md`.

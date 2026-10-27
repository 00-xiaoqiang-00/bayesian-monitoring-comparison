# Entry-point navigation

The tables below identify the analysis entries and case arguments. Read [RUNNING.md](RUNNING.md) and the retained instructions for the selected analysis before choosing an output directory.

## Primary evaluation and calibration

| Method | Case | Independent-evaluation entry | Calibration-path entry |
|---|---|---|---|
| RMUT | IMpower, KEYNOTE598, OAK, survival_colon | `code/RMUT/validation/<case>/run_independent_validation_10000.R` | `code/RMUT/calibration/<case>/build_calibration_cache_10000.R` |
| Cox | IMpower, KEYNOTE598, OAK | `code/method_sources/cox_rmst/cox_<case>/independent_validation/run_independent_validation_10000.R` | `code/sensitivity_sources/calibration_cox_rmst/cox_<case>/calibration/build_calibration_cache_10000.R` |
| Cox | survival_colon | `code/method_sources/cox_rmst/cox_survival_colon/independent_validation/run_independent_validation_10000.R` | `code/sensitivity_entries/cox_survival_colon/realcase_calibrate_lambda_mcmc.R` |
| RMST | IMpower, KEYNOTE598, OAK, survival_colon | `code/method_sources/cox_rmst/RMST_<case>/independent_validation/run_independent_validation_10000.R` | `code/sensitivity_sources/calibration_cox_rmst/RMST_<case>/calibration/build_calibration_cache_10000.R` |
| Win Ratio | All four | `code/vendor/bayesian_win_ratio/validation/validate_locked_winratio.R --case=CASE_KEY` | `code/vendor/bayesian_win_ratio/calibration/calibrate_LG.R --case=CASE_KEY` |

Both Win Ratio drivers accept `CASE_KEY=IMpower131`, `KEYNOTE598`, `OAK` or `survival_colon`. Their case configurations are in `code/vendor/bayesian_win_ratio/cases/`. Other method directories use `IMpower` for the IMpower131 setting.

Calibration-path entries generate PP paths and candidate summaries. Formal RMUT selection uses `code/parameter_selection/rmut/select_rmut_parameters.R`; independent evaluation uses the supplied locked parameters.

## Sensitivity entries

| Analysis | Methods | Entry location |
|---|---|---|
| Interim 50/80 sensitivity | RMUT | `code/RMUT/S8/<case>/run_independent_validation_10000.R` |
| Interim 50/80 sensitivity | RMST | `code/interim_sensitivity/RMST_<case>_sensitivity.R` |
| Interim 50/80 sensitivity | Cox | `code/sensitivity_entries/cox_interim/cox_<case>/interim_50_80/run_independent_validation_10000.R` |
| Interim-time sensitivity | Win Ratio | `code/sensitivity_entries/win_ratio_interim/run_interim_time_sensitivity.R` |
| Dependence sensitivity | RMUT | `code/RMUT/rho/<case>/run_independent_validation_10000.R` |
| Dependence sensitivity | Cox/RMST | `code/sensitivity_sources/dependence_cox_rmst/` |
| Dependence sensitivity | Win Ratio | `code/sensitivity_sources/dependence_win_ratio/validate_winratio_rho.R` in a prepared external runtime |

## Shared sources

- `code/common/dgm/` contains the case data generators.
- `code/common/rmst/` contains RMST fitting and validation helpers.
- `code/RMUT/common/` contains RMUT fitting, pairing and runtime helpers.
- `code/vendor/bayesian_win_ratio/` contains one R/C++ engine, four case configurations, source licenses and shared calibration and independent-evaluation entries.
- Win Ratio rho runtime preparation combines the registered C++ source in the formal package, its `shared_model/` adapter and case configurations. See the local runtime README.

## Population effects

| Analysis | Entry | Fixed inputs |
|---|---|---|
| H1-standard final effects | `code/estimands/compute_dgm_estimands.py --recompute` | `data/configuration/` |
| H1-observed final effects | `code/estimands/compute_h1_observed_estimands.py` | `data/configuration/` |
| Analysis-time RMST/RMUT effects | `code/estimands/compute_interim_estimands.py` | `data/configuration/` |

These entries require a new output directory outside the repository. Their local READMEs give complete commands and numerical definitions.

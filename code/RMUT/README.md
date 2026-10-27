# RMUT source package

The package covers IMpower, KEYNOTE598, OAK, and survival_colon. Case entry points use the fitting, posterior-pairing and runtime helpers in `common/`, together with `code/common/analysis_helpers.R` and the generators in `code/common/dgm/`. Dependence parameters are shared in `code/common/rho/`.

| Directory | Entry point | Analysis |
|---|---|---|
| `calibration/<case>` | `build_calibration_cache_10000.R` | Calibration |
| `validation/<case>` | `run_independent_validation_10000.R` | Independent validation |
| `S8/<case>` | `run_independent_validation_10000.R` | Interim analyses at 50% and 80% |
| `rho/<case>` | `run_independent_validation_10000.R` | Dependence sensitivity |

Use R 4.4.x with `brms`, `rstan`, `posterior`, `survival`, and `MASS`; Stan entries require a C++ toolchain.

`CURRENT_PARAMETER_DEFAULTS.tsv` gives the default RMUT thresholds. A frozen four-row TSV can be supplied through `RMUT_LOCKED_PARAMETER_MANIFEST` with columns `case`, `candidate_id`, `c`, and `gamma`. Thresholds classify posterior-predictive paths without changing data generation or model fitting.

Death and composite-event posterior rows are paired before integration. Equal-length samples retain their order; unequal samples are subsampled without replacement using an explicit seed. The helper restores the RNG kind and state. Failed fits remain in the output; unconditional summaries count them as no-go and use `Nmax` for ESS.

Use `code/parameter_selection/rmut/select_rmut_parameters.R` for formal parameter selection from the calibration PP paths. Calibration builders produce PP paths and candidate summaries.

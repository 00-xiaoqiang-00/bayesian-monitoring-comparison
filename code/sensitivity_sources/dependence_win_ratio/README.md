# Win Ratio dependence sensitivity

The entries evaluate low and high dependence for the four cases using the supplied KM-target parameters. Decision thresholds, prior settings, trial schedules and seed rules use the formal case configuration.

Use R 4.4.x with `Rcpp`, `RcppArmadillo`, `RcppDist`, `MASS`, `survival`, `digest`, and `jsonlite`, plus a working C++ toolchain. On Windows, set `RTOOLS44_HOME` if the compiler and make tools are not on PATH.

Keep the repository structure intact. From its root, prepare a new external runtime:

```text
Rscript code/sensitivity_sources/dependence_win_ratio/prepare_runtime.R --output=/absolute/path/to/new_winratio_runtime
```

The command assembles each case's runtime from the registered C++ sources, `RcppExports.R`, author metadata and license notices in `code/vendor/bayesian_win_ratio/engine/`; the adapter and locked parameters in `shared_model/`; the selected configuration in `code/vendor/bayesian_win_ratio/cases/`; the generators in `code/common/dgm/`; and the dependence helper and parameter table in `code/common/rho/`. The assembled runtime uses `vendor_snapshot/<case>/` for case sources and local paths for the shared helper and parameter table. Source identities are recorded in `SOURCE_PROVENANCE.json`; the build checks routine names and arities and records the built runtime in `RUNTIME_SHA256SUMS.tsv` and `RUNTIME_BUILD_RECEIPT.json`.

In the prepared runtime, run:

```text
bash run_winratio_rho.sh smoke
bash run_winratio_rho.sh formal
```

The runner accepts `smoke`, `pilot` or `formal`. Formal mode generates 10,000 replicates per scenario and requires completed smoke calculations. `RSCRIPT` selects Rscript, `NUM_CORES` selects workers, and `PILOT_N` controls pilot replicates. The per-case entry is `validate_winratio_rho.R`, with `--case`, `--level`, `--mode` and `--cores`. Supported cases are `IMpower131`, `KEYNOTE598`, `OAK` and `survival_colon`; levels are `low` and `high`.

Generated results and process logs use the runtime's `outputs/` and `logs/` directories. The loader uses R's platform-specific `.so` or `.dll` suffix. On Windows, invoke the Rscript entry directly with `--cores=1` because the analysis uses fork parallelism. The KM-target acceptance input is `rho_sensitivity_km_target_acceptance.csv`.

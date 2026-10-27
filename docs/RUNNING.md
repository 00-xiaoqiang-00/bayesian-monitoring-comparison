# Running the source package

Use R 4.4.x and install the packages required by the selected entry point, including `brms`, `rstan`, `posterior`, `survival`, and `MASS` where applicable. Stan-based entries require a working C++ toolchain. Win Ratio entries require `pkgload` to load the supplied engine and compile its C++ interface in the target environment.

Shell runners under `code/` use POSIX utilities. On Windows, use Git Bash or WSL for those wrappers and configure R/Rtools for the selected R and C++ entry.

Run commands from the repository root or pass explicit paths through the documented arguments. Use a new output directory outside the repository where the entry supports an output argument. Other entries use their documented output directories.

## Trial simulation

The method-by-case entries are listed in `docs/ENTRYPOINTS.md`. Formal independent evaluation and sensitivity entries generate 10,000 trial replicates per scenario unless the entry specifies a different count. Monte Carlo workers use the entry's `NUM_CORES` setting; posterior fits use the sampler settings in the source.

Data generators are under `code/common/dgm/`; shared RMST helpers are under `code/common/rmst/`. RMUT entries use `code/RMUT/common/` for fitting, posterior pairing and runtime helpers. A four-row TSV supplied through `RMUT_LOCKED_PARAMETER_MANIFEST` must contain `case`, `candidate_id`, `c`, and `gamma`.

Calibration and parameter selection use H0 and H1-standard PP paths and candidate summaries. The Cox and RMST IMpower131 builders additionally report H1-observed descriptively; those outputs do not enter parameter selection. Use `code/parameter_selection/rmut/select_rmut_parameters.R` for formal RMUT parameter selection; the fixed inputs and saved selections are documented in its README. Formal independent evaluation, including H1-observed, uses locked parameters through the separate validation entries listed in `docs/ENTRYPOINTS.md`.

## Win Ratio entries

From the repository root, the shared independent-evaluation entry selects a case explicitly:

```text
Rscript code/vendor/bayesian_win_ratio/validation/validate_locked_winratio.R --case=IMpower131
```

The supported case keys are `IMpower131`, `KEYNOTE598`, `OAK` and `survival_colon`. The common R/C++ package is in `code/vendor/bayesian_win_ratio/engine/`; case configurations are in `cases/` alongside it. Parameters and seeds come from `validation/`. Results and PP caches use `output/validation/<case>/<mode>/` under the Win Ratio directory.

Calibration uses `code/vendor/bayesian_win_ratio/calibration/calibrate_LG.R` with the same `--case` keys and writes to `output/cases/<configuration>/` under the Win Ratio directory.

Dependence sensitivity uses a separate external runtime. Follow `code/sensitivity_sources/dependence_win_ratio/README.md` to prepare it from the repository's registered C++ source, shared model adapter and case configurations before running the rho entries.

## Fixed-input calculations

Eight fixed calibration PP inputs support RMUT selection. Their paths and identities are listed in `data/parameter_selection/rmut/fixed_input_manifest.csv`.

Population-effect Python entries read the DGM and use NumPy/SciPy as listed in `code/estimands/requirements.txt`. H1-standard, H1-observed and analysis-time calculations are documented in `code/estimands/README.md`.

Method-specific input and runtime checks are part of the corresponding analysis entries.

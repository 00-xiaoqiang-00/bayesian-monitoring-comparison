# Win Ratio independent evaluation

The shared entry uses the locked parameter table and seed registry in this directory. The formal design evaluates 40%, 60% and 80% interim looks plus the final analysis, with 10,000 replicates per scenario.

From the repository root:

```text
Rscript code/vendor/bayesian_win_ratio/validation/validate_locked_winratio.R --case=IMpower131
```

Supported case keys are `IMpower131`, `KEYNOTE598`, `OAK` and `survival_colon`. Use R 4.4.x and a working C++ toolchain. The driver loads the shared `engine/` package and the selected configuration from `cases/`. It resolves the repository from its script location; `PROJECT_ROOT` can explicitly identify the repository root. The shell runner uses its own directory by default; `CONTROL_ROOT` can override that location and `RUN_R` selects Rscript.

Each case writes results and PP caches under `code/vendor/bayesian_win_ratio/output/validation/<case_key>/<mode>/`. Calibration output uses `output/cases/<config_id>/` and is separate from these independent-validation caches.

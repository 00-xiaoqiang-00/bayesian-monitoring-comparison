# Population effects

These calculations evaluate population effects under the locked DGM. WR is the death-prioritized IPCW pairwise target under the final follow-up process; Cox HR is the equal-allocation partial-likelihood projection for min(R,D) under that process. Delta RMST and delta RMUT use their fixed method-specific horizons and utility weight.

## Saved results

- `data/estimands/H1_standard_true_estimands.csv`: H1-standard final effects, with `computation_config.json`.
- `data/estimands/H1_observed_true_estimands.csv`: H1-observed final effects, with `H1_observed_computation_config.json`.
- `data/estimands/interim_true_estimands.csv`: H1-standard and H1-observed RMST/RMUT differences at the reported analysis times.

The final-effect files contain population-effect approximations and their numerical settings. MCSE and numerical intervals describe numerical Monte Carlo uncertainty, not uncertainty in the fitted DGM or original trial population. Cox MCSE does not measure finite-sample approximation bias. Quadrature error estimates and analytic normal-tail bounds are distinct quantities.

## Running

From the repository root, install the dependencies listed in this directory's `requirements.txt`. Choose a new output directory outside the repository for each calculation.

```sh
python code/estimands/compute_dgm_estimands.py --recompute --output-dir /path/to/new_standard_results
python code/estimands/compute_h1_observed_estimands.py --output-dir /path/to/new_observed_results
python code/estimands/compute_interim_estimands.py --output-dir /path/to/new_interim_results
```

The H1-standard entry, when called without `--recompute`, reads and checks its saved archive and locked configuration without numerical integration. H1-observed integrates independent censoring analytically conditional on latent events for WR, reducing numerical variance while preserving the target. Its Cox calculation uses the locked treatment distribution for H1-observed. Deterministic conditional restricted means are used for its RMST/RMUT calculation and for the analysis-time effects.

The analysis-time schedule supplies separate `RMST_tau_months` and `RMUT_tau_months`. In particular, the final KEYNOTE-598 horizons are 12 and 18 months respectively. The interim calculator's deterministic integration may differ in unreported digits from the saved grid approximation; both evaluate the same DGM quantities.

# RMUT interim-time sensitivity analysis

The four case entry points evaluate interim looks at 50% and 80%. They load the generators in `code/common/dgm/` and the fitting, posterior-pairing and runtime helpers in `code/RMUT/common/`.

Run a preflight without generating data or fitting models:

```text
Rscript CASE/run_independent_validation_10000.R --preflight-only
```

Formal mode generates 10,000 datasets per scenario with the included generator and explicit posterior-pairing seeds. Set `RMUT_LOCKED_PARAMETER_MANIFEST` to a four-row TSV with columns `case`, `candidate_id`, `c`, and `gamma` when a frozen parameter table is required.

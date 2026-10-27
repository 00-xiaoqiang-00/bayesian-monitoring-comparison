# RMUT dependence sensitivity analysis

The four case entry points load the generators in `code/common/dgm/` and the fitting, posterior-pairing and runtime helpers in `code/RMUT/common/`. The death marginal, mixture weights, and distribution families are held fixed; recurrence parameters are refitted at each rho level against the same PFS or recurrence cumulative-incidence targets.

Set `SENSITIVITY_RHO_LEVEL` to `low` or `high` and set `SENSITIVITY_RHO` to the matching value. Formal mode generates 10,000 datasets per scenario with explicit seeds for data generation, model fitting, and posterior pairing.

The frozen KM-target parameter table and helper are shared in `code/common/rho/`. Entry points resolve that directory directly and include its files in preflight and cache checks. Set `RMUT_LOCKED_PARAMETER_MANIFEST` when a frozen threshold manifest is required.

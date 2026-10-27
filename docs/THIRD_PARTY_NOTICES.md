# Source notices

The Win Ratio engine in `code/vendor/bayesian_win_ratio/engine/` retains the upstream `BayesianWinRatio` package metadata and author attribution. It includes the upstream `LICENSE` copyright notice and `LICENSE.md` MIT permission text, retained byte for byte from upstream commit `4c9e5264ba3799c1e18fefe1f7ffffcd263012b0`. Runtime preparation copies this metadata and both license notices into each case directory of the external dependence-analysis runtime. The upstream files are [LICENSE](https://github.com/xinran-h/BayesianWinRatio/blob/4c9e5264ba3799c1e18fefe1f7ffffcd263012b0/LICENSE) and [LICENSE.md](https://github.com/xinran-h/BayesianWinRatio/blob/4c9e5264ba3799c1e18fefe1f7ffffcd263012b0/LICENSE.md). These notices apply to the included third-party components; project analysis scripts and configurations are separate work.

Rcpp, RcppArmadillo, RcppDist, brms, rstan, survival, and other packages used by the scripts are external dependencies and are not vendored here.

The reconstructed survival inputs are derived from publicly available sources. The survival_colon data are taken from the public `survival::colon` dataset; no restricted patient-level clinical-trial data are distributed.

No project-wide license is asserted for the complete source distribution.

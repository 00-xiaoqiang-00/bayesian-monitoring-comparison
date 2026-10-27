# Four Bayesian monitoring strategies

Source code and fixed inputs for comparing Bayesian Win Ratio, Cox, restricted mean survival time (RMST), and restricted mean utility time (RMUT) monitoring in the IMpower131, KEYNOTE-598, OAK, and survival_colon settings.

The repository contains the four method implementations; calibration, parameter selection and independent simulation evaluation; the required data-generating mechanisms, parameters, seeds and fixed inputs; population-effect calculations; and running instructions and source notices. Reported interim-time and dependence analyses are included.

## Start here

1. Select a method, case and analysis in [docs/ENTRYPOINTS.md](docs/ENTRYPOINTS.md).
2. Read [docs/RUNNING.md](docs/RUNNING.md) and the retained instructions for the selected analysis.
3. Use the locked inputs and parameters in `data/configuration/`; direct generated files to a new output directory outside the repository where the entry supports it.

Use R 4.4.x. Stan and Win Ratio C++ entries need the applicable compiler toolchain. Python population-effect calculations use the dependencies in `code/estimands/requirements.txt`.

## Contents

| Location | Contents |
|---|---|
| `code/common/` | Shared case generators, monitoring helpers, RMST helpers and dependence parameters |
| `code/RMUT/` | RMUT calibration, independent simulation evaluation and sensitivity analyses |
| `code/method_sources/cox_rmst/` | Cox and RMST analysis implementations |
| `code/vendor/bayesian_win_ratio/` | One Win Ratio R/C++ engine, four case configurations and shared analysis entries |
| `code/parameter_selection/` | Parameter-selection sources, including the formal RMUT fixed-input rule |
| `code/interim_sensitivity/`, `code/sensitivity_entries/`, `code/sensitivity_sources/` | Reported interim-time and dependence analyses |
| `code/estimands/` | H1-standard, H1-observed and analysis-time population-effect calculations |
| `data/configuration/` | DGM, trial process, analysis schedule, locked parameters and seeds |
| `data/reconstruction/` | Public marginal inputs and source/reconstruction provenance |
| `data/parameter_selection/` | Required fixed RMUT calibration PP inputs and saved selections |
| `data/estimands/` | Saved population effects and numerical settings |
| `docs/` | Running instructions, source layout, analysis navigation and third-party notices |

## Analysis definitions

H0 and H1-standard are used in calibration under the documented error, power and efficiency requirements. H1-observed is evaluated after parameters are locked and is not used for parameter selection. Decision parameters and seeds are recorded in `data/configuration/`.

Independent simulation evaluation generates new trial replicates using the locked parameters. Each entry checks its required inputs and writes the requested analysis outputs.

The analyses start from the supplied marginal PFS and OS inputs. These records do not identify matched patient pairs; joint dependence is specified by the DGM. The survival_colon input is derived from the public `survival` package.

## Population effects

See [code/estimands/README.md](code/estimands/README.md) for H1-standard/H1-observed WR, Cox projection HR, delta RMST and delta RMUT, including numerical uncertainty. Time-specific calculations preserve the separate RMST and RMUT horizons.


## Attribution

Third-party attribution and component license notices are retained in [docs/THIRD_PARTY_NOTICES.md](docs/THIRD_PARTY_NOTICES.md) and the included source packages. No repository-wide license is asserted beyond those notices.

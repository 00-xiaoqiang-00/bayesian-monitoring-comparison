# RMUT parameter selection

`select_rmut_parameters.R` is the formal RMUT selection entry. It reads the eight H0/H1-standard calibration PP inputs in `data/parameter_selection/rmut/fixed_pp_inputs/`. `fixed_input_manifest.csv` records their paths and identities; `formal_selection/` contains the four eligible sets and selected parameter rows. H1-observed does not enter selection.

From the repository root:

```text
Rscript code/parameter_selection/rmut/select_rmut_parameters.R
Rscript code/parameter_selection/rmut/select_rmut_parameters.R verify /path/to/new_selection_summary
Rscript code/parameter_selection/rmut/select_rmut_parameters.R rescreen /path/to/new_selection_results
```

The default checks the saved candidates against the fixed PP paths and applies the paired-MCSE selection rule. The optional `verify` output directory must be new and outside the repository.

`rescreen NEW_OUTPUT_DIRECTORY [CASE ...]` evaluates the complete parameter grid using the same saved PP inputs. Both modes classify existing PP paths without generating patients or fitting models. `RMUT_FIXED_INPUT_MANIFEST` can identify another complete eight-row input manifest; relative input paths resolve from that manifest.

Calibration-path builders in `code/RMUT/` generate inputs and candidate summaries; formal parameter selection uses this entry.

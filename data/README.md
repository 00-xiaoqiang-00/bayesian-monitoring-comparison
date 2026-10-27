# Fixed inputs and reported numeric results

| Directory | Role |
|---|---|
| `configuration/` | Locked data-generating parameters, trial process, analysis schedule, decision parameters and seeds |
| `reconstruction/` | Supplied public marginal inputs and source provenance |
| `parameter_selection/` | Eight fixed RMUT calibration PP inputs, saved eligible sets and selections |
| `estimands/` | Final and analysis-time population effects and numerical settings |

The fixed PP files contain numeric decision probabilities, not fitted model objects or posterior parameter draws. Their manifests record source identities and file hashes. Retain these files when reproducing parameter selection.

The reconstructed PFS and OS inputs do not establish matched patient pairs. Joint dependence is specified by the DGM.


The seed registry's analysis_entry field links each seed specification to its analysis source. The scenario_label field is the exact text used in deterministic seed derivation. Its H1_sens2 label identifies the locked H1-standard distribution for IMpower131 RMST; it does not define an additional evaluation scenario. Keep this label unchanged when reproducing the seed sequence.

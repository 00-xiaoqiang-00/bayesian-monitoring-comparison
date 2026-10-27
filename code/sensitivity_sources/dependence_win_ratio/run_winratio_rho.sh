#!/usr/bin/env bash
set -euo pipefail
VERSION_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$VERSION_DIR"
mode="${1:?Pass smoke, pilot or formal}"
case "$mode" in smoke|pilot|formal) ;; *) exit 2;; esac
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
RS="${RSCRIPT:-Rscript}"
command -v "$RS" >/dev/null || { echo "Rscript not found: $RS" >&2; exit 2; }
exec 9>"$VERSION_DIR/.runner.lock"
flock -n 9 || { echo 'Another rho runner is active' >&2; exit 2; }
"$RS" --vanilla verify_runtime_sources.R "$VERSION_DIR"
mkdir -p "logs/$mode"
"$RS" --vanilla preflight_winratio_rho.R >"logs/$mode/preflight.log" 2>&1
if [[ "$mode" == formal ]]; then
  # Verify smoke outputs and their source signatures before a formal run.
  "$RS" --vanilla audit_winratio_rho.R smoke >"logs/$mode/smoke_gate.log" 2>&1
fi
cores="${NUM_CORES:-1}"
if [[ "$mode" == smoke ]]; then cores=1; fi
for case_key in IMpower131 KEYNOTE598 OAK survival_colon; do
  for level in low high; do
    echo "[START] $mode $case_key $level cores=$cores"
    extra=()
    if [[ "$mode" == pilot ]]; then extra+=("--n=${PILOT_N:-20}"); fi
    "$RS" --vanilla validate_winratio_rho.R "--mode=$mode" "--case=$case_key" "--level=$level" "--cores=$cores" "${extra[@]}" >"logs/$mode/${case_key}_${level}.log" 2>&1
    echo "[COMPLETE] $mode $case_key $level"
  done
done
"$RS" --vanilla audit_winratio_rho.R "$mode" >"logs/$mode/audit.log" 2>&1
echo "RUN_ALL_COMPLETE $mode"

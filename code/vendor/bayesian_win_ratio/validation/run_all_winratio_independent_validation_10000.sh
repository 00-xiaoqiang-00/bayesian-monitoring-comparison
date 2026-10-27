#!/usr/bin/env bash
set -euo pipefail

CONTROL_ROOT="${CONTROL_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
SCRIPT="$CONTROL_ROOT/validate_locked_winratio.R"
RUN_R="${RUN_R:-Rscript}"
RUN_ID="$(date +%Y%m%d_%H%M%S)"
LOG_DIR="$CONTROL_ROOT/logs/run_$RUN_ID"
STATUS_FILE="$LOG_DIR/STATUS.txt"
mkdir -p "$LOG_DIR"

printf 'status=RUNNING\nrun_id=%s\nstarted_at=%s\n' "$RUN_ID" "$(date -Is)" > "$STATUS_FILE"

for CASE_KEY in IMpower131 KEYNOTE598 OAK survival_colon; do
  LOG_FILE="$LOG_DIR/${CASE_KEY}.log"
  printf '[%s] START %s\n' "$(date -Is)" "$CASE_KEY" | tee -a "$LOG_FILE"
  "$RUN_R" "$SCRIPT" "--case=$CASE_KEY" 2>&1 | tee -a "$LOG_FILE"
  printf '[%s] COMPLETE %s\n' "$(date -Is)" "$CASE_KEY" | tee -a "$LOG_FILE"
done

printf 'status=COMPLETE\nrun_id=%s\nstarted_at=%s\ncompleted_at=%s\n' \
  "$RUN_ID" "$(sed -n 's/^started_at=//p' "$STATUS_FILE")" "$(date -Is)" > "$STATUS_FILE"

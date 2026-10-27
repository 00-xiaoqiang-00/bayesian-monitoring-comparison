#!/usr/bin/env bash
set -euo pipefail

CONTROL_ROOT="${CONTROL_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
SCRIPT="$CONTROL_ROOT/run_interim_time_sensitivity.R"
RUN_R="${RUN_R:-Rscript}"
FORMAL_PID_FILE="${FORMAL_PID_FILE:-$CONTROL_ROOT/formal_run.pid}"

MODE="formal"
EXTRA_ARGS=()
for arg in "$@"; do
  case "$arg" in
    --dry-run|--smoke|--replay-only|--rebuild-look50|--force)
      if [[ "$arg" == "--force" ]]; then
        FORCE=1
      else
        EXTRA_ARGS+=("$arg")
        if [[ "$arg" == "--smoke" ]]; then MODE="smoke"; fi
        if [[ "$arg" == "--dry-run" ]]; then MODE="dry-run"; fi
        if [[ "$arg" == "--replay-only" ]]; then MODE="replay-only"; fi
      fi
      ;;
    *)
      echo "Unknown argument: $arg" >&2
      exit 1
      ;;
  esac
done
FORCE="${FORCE:-0}"

if [[ "$MODE" == "formal" && "$FORCE" != "1" && -f "$FORMAL_PID_FILE" ]]; then
  FORMAL_PID="$(tr -d '[:space:]' < "$FORMAL_PID_FILE")"
  if [[ -n "$FORMAL_PID" ]] && kill -0 "$FORMAL_PID" 2>/dev/null; then
    echo "Formal independent validation still running (pid $FORMAL_PID)." >&2
    echo "Refusing to start sensitivity while independent validation is running. Wait, or pass --force / --replay-only / --smoke / --dry-run." >&2
    exit 1
  fi
fi

RUN_ID="$(date +%Y%m%d_%H%M%S)"
LOG_DIR="$CONTROL_ROOT/logs/run_$RUN_ID"
STATUS_FILE="$LOG_DIR/STATUS.txt"
mkdir -p "$LOG_DIR"
echo $$ > "$CONTROL_ROOT/ACTIVE_RUN.pid"
printf 'status=RUNNING\nrun_id=%s\nstarted_at=%s\nmode=%s\n' \
  "$RUN_ID" "$(date -Is)" "$MODE" > "$STATUS_FILE"

for CASE_KEY in IMpower131 KEYNOTE598 OAK survival_colon; do
  LOG_FILE="$LOG_DIR/${CASE_KEY}.log"
  printf '[%s] START %s\n' "$(date -Is)" "$CASE_KEY" | tee -a "$LOG_FILE"
  "$RUN_R" "$SCRIPT" "--case=$CASE_KEY" "${EXTRA_ARGS[@]+"${EXTRA_ARGS[@]}"}" 2>&1 | tee -a "$LOG_FILE"
  printf '[%s] COMPLETE %s\n' "$(date -Is)" "$CASE_KEY" | tee -a "$LOG_FILE"
done

"$RUN_R" "$CONTROL_ROOT/summarize_interim_time_sensitivity.R" "--mode=$( [[ "$MODE" == "smoke" ]] && echo smoke || echo formal_10000 )" 2>&1 | tee -a "$LOG_DIR/summarize.log" || true

printf 'status=COMPLETE\nrun_id=%s\nstarted_at=%s\ncompleted_at=%s\nmode=%s\n' \
  "$RUN_ID" "$(sed -n 's/^started_at=//p' "$STATUS_FILE")" "$(date -Is)" "$MODE" > "$STATUS_FILE"
rm -f "$CONTROL_ROOT/ACTIVE_RUN.pid"

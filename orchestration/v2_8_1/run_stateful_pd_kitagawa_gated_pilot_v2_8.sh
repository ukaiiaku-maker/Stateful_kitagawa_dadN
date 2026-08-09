#!/usr/bin/env bash
set -euo pipefail
HERE=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
STANDALONE_ROOT=${STANDALONE_ROOT:?set STANDALONE_ROOT}
OUTROOT=${OUTROOT:-$HERE/runs/stateful_pd_kitagawa_v2_8_gated_pilot}
GATE=${GATE:-failure}
MAX_BLOCKS=${MAX_BLOCKS:-20000}
CHECKPOINT_EVERY_BLOCKS=${CHECKPOINT_EVERY_BLOCKS:-10}
PRINT_EVERY=${PRINT_EVERY:-100}

exec python -B "$HERE/run_stateful_pd_kitagawa_gated_pilot_v2_8.py" \
  --standalone-root "$STANDALONE_ROOT" \
  --outroot "$OUTROOT" \
  --gate "$GATE" \
  --max-blocks "$MAX_BLOCKS" \
  --checkpoint-every-blocks "$CHECKPOINT_EVERY_BLOCKS" \
  --print-every "$PRINT_EVERY" \
  "$@"

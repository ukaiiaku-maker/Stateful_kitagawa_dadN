#!/usr/bin/env bash
set -euo pipefail
HERE=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
python -B "$HERE/verify_stateful_pd_kitagawa_gated_pilot_v2_8_1.py"
echo "STATEFUL_PD_KITAGAWA_GATED_PILOT_V2_8_1 verification OK"

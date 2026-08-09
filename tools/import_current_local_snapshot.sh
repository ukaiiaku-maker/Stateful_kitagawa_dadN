#!/usr/bin/env bash
set -euo pipefail

# Import the authoritative local Stateful-PD source tree into this Git repository.
# This intentionally copies source/tests only; runs and large generated outputs are excluded.

SRC="${STATEFUL_SOURCE_ROOT:-/Volumes/Data/Data/Nanopillar_calculation/Fatigue-PF/dist/stateful_pd_kitagawa_production_v2_standalone}"
REGISTRY="${STATEFUL_REGISTRY_ROOT:-/Volumes/Data/Data/Nanopillar_calculation/stateful_pd_kitagawa_production_v2_persistent_registry}"
REPO="$(git rev-parse --show-toplevel)"

cd "$REPO"

case "$(git branch --show-current)" in
  codex/sn-endurance-limit|codex/sn-endurance-limit-*) ;;
  *)
    echo "ERROR: run this only on codex/sn-endurance-limit or a child branch." >&2
    echo "Current branch: $(git branch --show-current)" >&2
    exit 2
    ;;
esac

if [[ ! -d "$SRC/arrhenius_fracture" ]]; then
  echo "ERROR: source tree not found: $SRC/arrhenius_fracture" >&2
  exit 2
fi

mkdir -p arrhenius_fracture tests campaigns/local_registry reference

# Core package and tests from the actual local standalone installation.
rsync -a --exclude '__pycache__/' --exclude '*.pyc' \
  "$SRC/arrhenius_fracture/" arrhenius_fracture/

if [[ -d "$SRC/tests" ]]; then
  rsync -a --exclude '__pycache__/' --exclude '*.pyc' \
    "$SRC/tests/" tests/
fi

# Top-level source helpers/runners from the standalone, but no run directories/data.
shopt -s nullglob
for f in \
  "$SRC"/run_*.py "$SRC"/run_*.sh \
  "$SRC"/summarize_*.py "$SRC"/analyze_*.py \
  "$SRC"/verify_*.sh "$SRC"/environment*.yml "$SRC"/requirements*.txt; do
  cp -p "$f" ./
done

# Preserve the latest registry-side campaign/orchestration code separately from core source.
if [[ -d "$REGISTRY" ]]; then
  for f in \
    "$REGISTRY"/run_stateful_pd_kitagawa*.py \
    "$REGISTRY"/run_stateful_pd_kitagawa*.sh \
    "$REGISTRY"/analyze_stateful_pd_kitagawa*.py \
    "$REGISTRY"/summarize_stateful_pd_kitagawa*.py; do
    [[ -f "$f" ]] && cp -p "$f" campaigns/local_registry/
  done
fi

# Record hashes without embedding absolute local filesystem paths in the public repo.
python - <<'PY'
from pathlib import Path
from datetime import datetime, timezone
import hashlib, json, subprocess

root = Path.cwd()
interesting = [
    "arrhenius_fracture/sn_pd2d_stateful_v8_7_generalized_features.py",
    "arrhenius_fracture/stateful_peridynamics_v8_7_local_front_spacing.py",
    "arrhenius_fracture/sn_feature_geometry_v8_7.py",
    "arrhenius_fracture/sn_pd2d_stateful_v8_3.py",
    "arrhenius_fracture/stateful_peridynamics_v8_3.py",
]
hashes = {}
for rel in interesting:
    p = root / rel
    if p.is_file():
        hashes[rel] = hashlib.sha256(p.read_bytes()).hexdigest()
manifest = {
    "manifest_version": "LOCAL_STATEFUL_SOURCE_IMPORT_V1",
    "created_utc": datetime.now(timezone.utc).isoformat(),
    "git_branch": subprocess.check_output(["git", "branch", "--show-current"], text=True).strip(),
    "source_hashes": hashes,
}
(root / "reference" / "LOCAL_IMPORT_MANIFEST.json").write_text(
    json.dumps(manifest, indent=2, sort_keys=True) + "\n"
)
print(json.dumps(manifest, indent=2, sort_keys=True))
PY

python -m compileall -q arrhenius_fracture tests campaigns

# Run tests that are already in the imported snapshot. Do not launch physical campaigns.
if compgen -G 'tests/test_stateful_pd_v8_3_*.py' >/dev/null; then
  PYTHONPATH="$REPO" python -m unittest discover -s tests -p 'test_stateful_pd_v8_3_*.py' -v
fi

echo
echo "Import complete. No physical simulation was launched."
echo "Review before committing:"
git status --short

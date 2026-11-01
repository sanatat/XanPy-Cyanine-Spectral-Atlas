#!/usr/bin/env bash
set -euo pipefail
PYTHON_BIN="${XANPY_PYTHON:-$HOME/Downloads/.xanpy-venv/bin/python}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
"$PYTHON_BIN" "$HERE/stage2.py"
"$PYTHON_BIN" "$HERE/stage3.py" emission_nm
"$PYTHON_BIN" "$HERE/stage3.py" absorption_nm
"$PYTHON_BIN" "$HERE/stage4.py"
"$PYTHON_BIN" "$HERE/stage5.py" emission_nm
for stage in 6 7 8 9 10 11 12 13; do
    "$PYTHON_BIN" "$HERE/stage${stage}.py"
done

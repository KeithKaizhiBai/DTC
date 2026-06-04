#!/bin/bash
set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PYTHON_BIN="${PYTHON_BIN:-python}"

"${PYTHON_BIN}" "${PROJECT_DIR}/reproduce_interacting.py" \
  --params "${PROJECT_DIR}/params/interacting_params_local.json" \
  --output-root "${PROJECT_DIR}" \
  --device cuda

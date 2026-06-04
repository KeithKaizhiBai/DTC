#!/bin/bash
#SBATCH --job-name=dtc_interacting
#SBATCH --partition=NV_H100
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --exclude=gpuh01

set -euo pipefail

PROJECT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
PYTHON_BIN="${PYTHON_BIN:-python}"
PARAMS_FILE="${PARAMS_FILE:-${PROJECT_DIR}/params/interacting_params_hpc.json}"
LOG_DIR="${PROJECT_DIR}/logs"
RESULT_DIR="${PROJECT_DIR}/results"
mkdir -p "${LOG_DIR}" "${RESULT_DIR}"

echo "timestamp=$(date -Iseconds)" > "${LOG_DIR}/env_snapshot.txt"
echo "hostname=$(hostname)" >> "${LOG_DIR}/env_snapshot.txt"
echo "project_dir=${PROJECT_DIR}" >> "${LOG_DIR}/env_snapshot.txt"
echo "params_file=${PARAMS_FILE}" >> "${LOG_DIR}/env_snapshot.txt"
echo "exclude_nodes=gpuh01" >> "${LOG_DIR}/env_snapshot.txt"
echo "cuda_visible_devices=${CUDA_VISIBLE_DEVICES:-}" >> "${LOG_DIR}/env_snapshot.txt"
command -v nvidia-smi >/dev/null 2>&1 && nvidia-smi >> "${LOG_DIR}/env_snapshot.txt" || true

"${PYTHON_BIN}" "${PROJECT_DIR}/reproduce_interacting.py" \
  --params "${PARAMS_FILE}" \
  --output-root "${RESULT_DIR}" \
  --device cuda \
  > "${LOG_DIR}/interacting_run.log" 2>&1

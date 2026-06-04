# DTC Reproduction Project

This repository contains reading notes and numerical reproductions for
"Solvable model for discrete time crystal enforced by nonsymmorphic dynamical
symmetry", Physical Review Research 5, L032024 (2023).

## Contents

- `reading_report/`: English reading report in Markdown, LaTeX, and PDF.
- `reproduction_noninteracting/`: noninteracting two-level reproduction, data,
  figures, and LaTeX/PDF reports.
- `reproduction_interacting/`: interacting spin-chain reproduction, local/HPC
  parameters, Slurm launchers, data, figures, and LaTeX/PDF report.
- `docs/PROJECT_SUMMARY.md`: project memory and current status.

The project uses the dimensionless convention `hbar = 1`, `omega = 1`,
`T = 2*pi`, and `Omega = alpha/8`.

## Main Results

- The noninteracting reproduction compares fourth-order Magnus evolution with
  RK4 and uses both refined `alpha_13` and detuned `alpha = 81.60`.
- The interacting reproduction uses `J = 0.2`, open boundaries, and an initial
  all-plus `x`-polarized product state.
- The latest H100 run focuses on `L = 8`; refined `alpha_13` keeps the
  stroboscopic order parameter above `Z(n)=0.8` until `n=82839`, while the
  detuned `alpha=81.60` run crosses that threshold at `n=12`.

## Running

Local short interacting check:

```bash
python reproduction_interacting/reproduce_interacting.py \
  --params reproduction_interacting/params/interacting_params_local.json \
  --output-root reproduction_interacting \
  --device cuda
```

HPC run:

```bash
cd reproduction_interacting
PYTHON_BIN=$(which python) sbatch run_interacting_hpc.sh
```

The HPC launcher excludes `gpuh01`.

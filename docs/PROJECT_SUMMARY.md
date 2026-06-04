# Project Summary

## 1. Project Goal

Reproduce the paper "Solvable model for discrete time crystal enforced by nonsymmorphic dynamical symmetry" (Physical Review Research 5, L032024, 2023). The immediate goals are an English reading report, a Python/CUDA reproduction of the noninteracting two-level dynamics, and a first interacting spin-chain reproduction with an HPC-ready package for longer runs.

## 2. Current Architecture

The project is organized around the source PDFs and two output folders:

- `reading_report/`: English literature notes, reading-stage plan, LaTeX reading report, and compiled reading-report PDF.
- `reproduction_noninteracting/`: Python/CUDA scripts, generated figures, data, reproduction report, and numerical-methods supplement for the noninteracting model.
- `reproduction_interacting/`: exact state-vector CUDA/PyTorch script, local and HPC parameter files, generated small-chain data/figures, LaTeX interacting report, and the compressed HPC upload package.
- `output/pdf/`: Copies of final compiled report PDFs.
- `tmp/pdfs/`: temporary extracted PDF text and rendered page images used for local verification.

The numerical workflows use explicit functions in single scripts. They set `hbar = 1` and `omega = 1`, so the drive period is `T = 2*pi` and `alpha = 8*Omega/omega`.

## 3. Important Files

- `Repeat.md`: User-provided task definition.
- `readme.md`: GitHub-facing project overview and run instructions.
- `.gitignore`: Excludes LaTeX intermediates, temporary files, logs, and bulky copied output folders.
- `.gitattributes`: Forces LF line endings for source, scripts, data, and LaTeX files.
- `requirements.txt`: Minimal Python package list for local reproduction scripts.
- `prr_2023_Solvable model for discrete time crystal enforced by nonsymmorphic dynamical symmetry.pdf`: Main paper.
- `prr_2023_SI_Solvable model for discrete time crystal enforced by nonsymmorphic dynamical symmetry.pdf`: Supplemental material.
- `reading_report/reading_report.md`: Markdown source notes for the English reading report.
- `reading_report/reading_report.tex`: LaTeX source for the English reading report.
- `reading_report/reading_report.pdf`: Compiled reading-report PDF.
- `reproduction_noninteracting/reproduce_noninteracting.py`: Main reproduction script for the noninteracting dynamics.
- `reproduction_noninteracting/reproduction_report.md`: Markdown source notes for the noninteracting reproduction report.
- `reproduction_noninteracting/reproduction_report.tex`: LaTeX source for the noninteracting reproduction report.
- `reproduction_noninteracting/reproduction_report.pdf`: Compiled noninteracting reproduction-report PDF.
- `reproduction_noninteracting/numerical_methods_supplement.tex`: LaTeX supplement explaining fourth-order Magnus and RK4 theory.
- `reproduction_noninteracting/numerical_methods_supplement.pdf`: Compiled numerical-methods supplement.
- `reproduction_interacting/reproduce_interacting.py`: Exact state-vector CUDA/PyTorch implementation of the interacting chain.
- `reproduction_interacting/params/interacting_params_local.json`: Short local validation parameters.
- `reproduction_interacting/params/interacting_params_hpc.json`: Current cluster-run parameters focused on `L = 8`, with Fig. 4-scale long stroboscopic sampling.
- `reproduction_interacting/run_interacting_hpc.sh`: Slurm launcher with `#SBATCH --exclude=gpuh01`.
- `reproduction_interacting/reproduction_interacting_report.tex`: LaTeX source for the interacting reproduction report.
- `reproduction_interacting/reproduction_interacting_report.pdf`: Compiled interacting reproduction report.
- `output/pdf/reading_report.pdf`: Final PDF copy for the reading report.
- `output/pdf/reproduction_report.pdf`: Final PDF copy for the noninteracting reproduction report.
- `output/pdf/numerical_methods_supplement.pdf`: Final PDF copy for the numerical-methods supplement.
- `output/pdf/interacting_reproduction_report.pdf`: Final PDF copy for the interacting reproduction report.
- `output/interacting_reproduction_hpc_package.zip`: Compressed HPC upload package.

## 4. Major Decisions

- Date: 2026-06-04
- Decision: Use a single function-based Python script with explicit constants near the top.
- Reason: The project is a research reproduction task where transparent equations, parameters, and saved outputs are more valuable than a general framework.
- Consequences: Future changes should extend the existing functions instead of adding classes or configuration managers unless stateful many-body simulations require them.

- Date: 2026-06-04
- Decision: Use a fourth-order Gauss-Legendre Magnus step for the unitary solver and classical RK4 as the comparison solver.
- Reason: The user requested a fourth-order Magnus method that preserves unitarity and RK4 as an accuracy comparison.
- Consequences: Magnus evolution should show negligible norm drift, while RK4 can be monitored through explicit norm diagnostics.

- Date: 2026-06-04
- Decision: Use the dimensionless convention `omega = 1`, `hbar = 1`, and `Omega = alpha/8`.
- Reason: The paper defines `alpha = 8*Omega/omega`, and this convention keeps all numerical times in units of the drive angular frequency.
- Consequences: Plots label time as `t/T`, and spectra label angular frequency as `omega_tilde/omega`.

- Date: 2026-06-04
- Decision: Implement the interacting chain as an exact state-vector simulation using bit-flip tensor indexing instead of dense many-body matrices.
- Reason: The local and cluster targets need transparent Hamiltonian action while avoiding dense `2^L x 2^L` memory growth.
- Consequences: Local validation is limited to small chains, while the same script can push longer stroboscopic runs on the GPU cluster.

- Date: 2026-06-04
- Decision: Use RK4 with per-step normalization for the interacting state evolution.
- Reason: It gives a simple, inspectable time-domain implementation for this first interacting-stage check, and normalization keeps local finite-step norm drift controlled.
- Consequences: Diagnostics record maximum norm error, but high-precision lifetime scaling should still be rerun with tighter step counts on the cluster if needed.

- Date: 2026-06-04
- Decision: Package cluster scripts with `gpuh01` excluded.
- Reason: The available cluster guidance indicates that `gpuh01` has an older CUDA environment.
- Consequences: `run_interacting_hpc.sh` writes an environment snapshot and uses `#SBATCH --exclude=gpuh01`.

- Date: 2026-06-04
- Decision: Change the interacting HPC stroboscopic job to `L = 8` and use one-period Floquet-operator logarithmic sampling.
- Reason: Direct step-by-step RK4 to the Fig. 4 lifetime scales would require too many periods and would write impractically large CSV files.
- Consequences: The package now samples `alpha=81.60` to `n = 100000` and refined `alpha_13` to `n = 100000000`, while keeping the time trace at `L = 8` for the first 40 periods.

- Date: 2026-06-04
- Decision: Treat `Z(n)<0.8` as the first operational lifetime threshold for the current single-size HPC check.
- Reason: The paper describes the stroboscopic signal dropping from the high plateau before later collapse; a fixed threshold gives a concrete comparison between alpha values before full finite-size fitting.
- Consequences: For the `L=8` H100 run, the detuned case crosses below `0.8` at `n=12`, while refined `alpha_13` first crosses below `0.8` at `n=82839`.

## 5. Completed Milestones

- Date: 2026-06-04
- Milestone: Parsed the user task and extracted readable text from the main paper and supplemental material.
- Evidence: `tmp/pdfs/paper.txt`, `tmp/pdfs/si.txt`, and rendered page checks in `tmp/pdfs/`.

- Date: 2026-06-04
- Milestone: Verified the key noninteracting Hamiltonian from rendered PDF pages.
- Evidence: The reproduction code uses `H0(t) = 0.5*Omega*sin(t)*sigma_x + Omega*sin(t/2)^2*sigma_y`.

- Date: 2026-06-04
- Milestone: Completed the first noninteracting reproduction run on CUDA.
- Evidence: `reproduction_noninteracting/data/diagnostics.json` reports `device: cuda`, RTX 4090, `alpha_13_refined = 80.04624211689759`, and generated PNG/CSV outputs.

- Date: 2026-06-04
- Milestone: Rebuilt both reports as LaTeX sources and compiled PDFs.
- Evidence: `reading_report/reading_report.pdf`, `reproduction_noninteracting/reproduction_report.pdf`, and copied final PDFs under `output/pdf/`.

- Date: 2026-06-04
- Milestone: Implemented `Repeat.md` section `I. 修正指示---无相互作用`.
- Evidence: `reproduction_noninteracting/reproduce_noninteracting.py` now computes solver traces for both refined `alpha_13` and detuned `alpha=81.60`, generates fractional-period projection data at `0,T/4,T/2,3T/4`, and the corrected LaTeX/PDF report plus numerical-methods supplement were compiled.

- Date: 2026-06-04
- Milestone: Completed the local first pass of `Repeat.md` stage 3, interacting reproduction.
- Evidence: `reproduction_interacting/data/diagnostics.json` reports CUDA execution, `J = 0.2`, open boundaries, maximum normalized-step norm error near `1.2e-7`, and dominant Fourier peaks near normalized frequency `0.5`.

- Date: 2026-06-04
- Milestone: Prepared the interacting-stage report and HPC upload package.
- Evidence: `reproduction_interacting/reproduction_interacting_report.pdf`, `output/pdf/interacting_reproduction_report.pdf`, and `output/interacting_reproduction_hpc_package.zip`.

- Date: 2026-06-04
- Milestone: Uploaded the screened project files to GitHub.
- Evidence: Repository `https://github.com/KeithKaizhiBai/DTC`, branch `main`, commit `392796e`.

## 6. Recent Major Changes

- Date: 2026-06-04
- Files changed: `docs/PROJECT_SUMMARY.md`
- Summary: Created the project memory file required by the project rules.
- Why it matters: Future sessions can recover the reproduction scope, conventions, and current state without rereading all context.

- Date: 2026-06-04
- Files changed: `reading_report/stage_plan.md`, `reading_report/reading_report.md`
- Summary: Added the reading-stage plan and English report covering the DTC background, nonsymmorphic dynamical symmetry mechanism, model, parameters, and observables.
- Why it matters: This completes the first literature-reading deliverable requested in `Repeat.md`.

- Date: 2026-06-04
- Files changed: `reproduction_noninteracting/stage_plan.md`, `reproduction_noninteracting/reproduce_noninteracting.py`, `reproduction_noninteracting/reproduction_report.md`
- Summary: Added a function-based CUDA/Python reproduction script with fourth-order Magnus evolution, RK4 comparison, `alpha_13` numerical refinement, CSV data, diagnostics, and PNG figures.
- Why it matters: This completes a verified first pass at the noninteracting dynamics reproduction.

- Date: 2026-06-04
- Files changed: `reading_report/reading_report.tex`, `reading_report/reading_report.pdf`, `reproduction_noninteracting/reproduction_report.tex`, `reproduction_noninteracting/reproduction_report.pdf`, `output/pdf/reading_report.pdf`, `output/pdf/reproduction_report.pdf`
- Summary: Rebuilt the reading and reproduction reports in LaTeX and compiled them to PDF with `pdflatex`.
- Why it matters: The report deliverables now match the requested LaTeX plus compiled PDF format.

- Date: 2026-06-04
- Files changed: `reproduction_noninteracting/reproduce_noninteracting.py`, `reproduction_noninteracting/reproduction_report.tex`, `reproduction_noninteracting/reproduction_report.pdf`, `reproduction_noninteracting/numerical_methods_supplement.tex`, `reproduction_noninteracting/numerical_methods_supplement.pdf`, `reproduction_noninteracting/data/*alpha_13*`, `reproduction_noninteracting/data/fractional_projections_*`, `reproduction_noninteracting/figures/*alpha_13*`, `reproduction_noninteracting/figures/fractional_projections_*`, `output/pdf/reproduction_report.pdf`, `output/pdf/numerical_methods_supplement.pdf`
- Summary: Executed the no-interaction correction instructions: both solver and projection sections now compare the accurate and detuned parameters; projection sampling includes `0,T/4,T/2,3T/4`; numerical-method theory is moved into a standalone LaTeX/PDF supplement.
- Why it matters: The no-interaction reproduction now addresses the identified gaps in parameter comparison, projection sampling, and numerical-method exposition.

- Date: 2026-06-04
- Files changed: `reproduction_interacting/stage_plan.md`, `reproduction_interacting/reproduce_interacting.py`, `reproduction_interacting/params/interacting_params_local.json`, `reproduction_interacting/params/interacting_params_hpc.json`, `reproduction_interacting/run_interacting_hpc.sh`, `reproduction_interacting/run_interacting_local.sh`, `reproduction_interacting/data/*`, `reproduction_interacting/figures/*`
- Summary: Added and ran a first interacting exact-state-vector CUDA reproduction for `J = 0.2`, open boundaries, `alpha=81.60`, and refined `alpha_13`.
- Why it matters: This verifies that the interacting implementation produces the expected period-doubled response before spending cluster time on larger and longer simulations.

- Date: 2026-06-04
- Files changed: `reproduction_interacting/reproduction_interacting_report.tex`, `reproduction_interacting/reproduction_interacting_report.pdf`, `output/pdf/interacting_reproduction_report.pdf`, `reproduction_interacting/hpc_package/interacting_reproduction_package/*`, `output/interacting_reproduction_hpc_package.zip`
- Summary: Created the LaTeX/PDF interacting report, visually checked the rendered PDF, and built the self-contained HPC upload package.
- Why it matters: The interacting-stage deliverables now match the requested report format and are ready for cluster upload.

- Date: 2026-06-04
- Files changed: `reproduction_interacting/reproduce_interacting.py`, `reproduction_interacting/params/interacting_params_hpc.json`, `reproduction_interacting/hpc_package/interacting_reproduction_package/*`, `output/interacting_reproduction_hpc_package.zip`
- Summary: Reworked the HPC package to focus on `L = 8`; the stroboscopic run now uses a dense one-period Floquet operator and logarithmic sampling up to `n = 100000` for `alpha=81.60` and `n = 100000000` for refined `alpha_13`.
- Why it matters: This matches the period range needed for a Fig. 4-style `L=8` check without attempting infeasible step-by-step evolution through every period.

- Date: 2026-06-04
- Files changed: `output/interacting_reproduction_package/results/*`, `reproduction_interacting/figures/hpc_z_l8_longtime.png`, `reproduction_interacting/data/hpc_diagnostics_l8.json`, `reproduction_interacting/reproduction_interacting_report.tex`, `reproduction_interacting/reproduction_interacting_report.pdf`, `output/pdf/interacting_reproduction_report.pdf`
- Summary: Ingested the completed H100 `L=8` cluster run, generated a full-range long-time `Z(n)` figure, and updated the interacting LaTeX/PDF report with HPC diagnostics and larger-system parameter guidance.
- Why it matters: The interacting report now contains actual cluster data rather than only local validation and package preparation.

- Date: 2026-06-04
- Files changed: `readme.md`, `.gitignore`, `.gitattributes`, `requirements.txt`, selected report/source/data/figure files in `reading_report/`, `reproduction_noninteracting/`, `reproduction_interacting/`, and `docs/PROJECT_SUMMARY.md`
- Summary: Screened upload-worthy files and pushed them to `https://github.com/KeithKaizhiBai/DTC` on branch `main`.
- Why it matters: The GitHub repository now contains the reproducible code, LaTeX/PDF reports, parameters, core data, and figures while excluding original paper PDFs, cluster manuals, temporary files, logs, and bulky copied output folders.

## 7. Validation and Tests

Planned validation commands:

```powershell
python reproduction_noninteracting\reproduce_noninteracting.py
```

Latest validation checks:

- CUDA device was used: NVIDIA GeForce RTX 4090.
- Refined `alpha_13`: Magnus4 norm drift `4.02e-05`, RK4 norm drift `1.59e-03`, maximum `sigma_x(t)` solver difference `8.35e-03`, dominant normalized spectrum frequency `0.4998`.
- Detuned `alpha=81.60`: Magnus4 norm drift `3.92e-05`, RK4 norm drift `1.79e-03`, maximum `sigma_x(t)` solver difference `6.50e-03`, dominant normalized spectrum frequency `0.5498`.
- Refined `alpha_13`: `80.04624211689759`; one-period return probability: `3.57e-18`.
- Fractional-period projections are sampled at `0,T/4,T/2,3T/4` for both refined `alpha_13` and detuned `alpha=81.60`.

Report build checks:

```powershell
pdflatex -interaction=nonstopmode -halt-on-error reading_report.tex
pdflatex -interaction=nonstopmode -halt-on-error reproduction_report.tex
pdflatex -interaction=nonstopmode -halt-on-error numerical_methods_supplement.tex
pdftoppm -png -r 120 reading_report\reading_report.pdf tmp\pdfs\report_checks\reading_report
pdftoppm -png -r 120 reproduction_noninteracting\reproduction_report.pdf tmp\pdfs\report_checks\reproduction_report
pdftoppm -png -r 120 reproduction_noninteracting\numerical_methods_supplement.pdf tmp\pdfs\report_checks\numerical_methods_supplement
pdftoppm -png -r 120 reproduction_interacting\reproduction_interacting_report.pdf tmp\pdfs\report_checks\reproduction_interacting_report
```

Rendered PDF pages were visually checked for clipping, blank pages, missing figures, and overlapping text.

Interacting-stage validation command:

```powershell
python reproduction_interacting\reproduce_interacting.py --params reproduction_interacting\params\interacting_params_local.json --output-root reproduction_interacting --device cuda
```

Latest interacting validation checks:

- CUDA device was used: NVIDIA GeForce RTX 4090.
- Model: open-boundary spin chain with `J = 0.2` and `sigma_x sigma_x` interaction.
- Time-trace run: `L = 10`, 40 periods, `alpha=81.60` and refined `alpha_13 = 80.04624211689759`.
- Dominant normalized spectrum frequency: `0.4998` for both tested alpha values.
- Maximum normalized-step norm error: `1.19e-7`.
- Stroboscopic check over 80 periods: larger local chains `L = 6, 8, 10` keep positive `Z(n)` around `0.95` to `0.98`, while `L = 4` shows strong finite-size oscillations.
- `output/interacting_reproduction_hpc_package.zip` contains only the source script, parameter files, launchers, readme, and empty `logs/` and `results/` directories.
- Latest package smoke test: a reduced `log_floquet` run with `L=4`, tiny period counts, and `steps_per_period=4` completed successfully in `tmp/floquet_smoke`.

Latest H100 `L=8` cluster-result checks:

- Run location copied locally under `output/interacting_reproduction_package/`.
- GPU/node: NVIDIA H100 80GB HBM3 on `gpuh14`; `gpuh01` was excluded.
- Time trace: `L=8`, 40 periods, 120 steps per period; dominant normalized frequency `0.4998958550` for both alpha values.
- Stroboscopic `alpha=81.60`: sampled to `n=100000`, 1793 samples, max norm error `2.39e-5`, first `Z(n)<0.8` at `n=12`, first `Z(n)<0.5` at `n=2018`, first `Z(n)<0` at `n=3642`.
- Stroboscopic refined `alpha_13`: sampled to `n=100000000`, 1801 samples, max norm error `1.64e-5`, first `Z(n)<0.8` at `n=82839`, first `Z(n)<0.5` at `n=138896`, first `Z(n)<0` at `n=212790`.
- Updated PDF report was compiled with `pdflatex` and rendered via `pdftoppm` as `tmp/pdfs/report_checks/reproduction_interacting_report_updated-*.png`.

## 8. Open Problems

- The local interacting run is a first small-chain validation; full long-time finite-size lifetime scaling for Fig. 4(c,d) still needs the prepared HPC run.
- The current completed HPC run targets only `L = 8`; the next larger-system step should be `L = 10` before attempting dense-Floquet `L = 12`.
- Exact Heun-function evaluation is not implemented; special `alpha_n` values are refined numerically through high-accuracy period evolution.
- The current reading report is concise; broader many-body localization and experimental DTC background can be expanded later if needed.
- CUDA is used, but the two-level problem is too small to benefit from GPU acceleration. Increasing `STEPS_PER_PERIOD` gives higher resolution but longer runtime.

## 9. Next Steps

- If higher publication-quality accuracy is needed, rerun the noninteracting script with larger `STEPS_PER_PERIOD` and compare diagnostics.
- For the next cluster run, edit `params/interacting_params_hpc.json` to set `stroboscopic.lengths` to `[10]` and set `periods_by_alpha.alpha_13` to about `10000000000` if reproducing the Fig. 4 scale for `L=10`.
- After the `L=10` cluster run finishes, compare threshold periods such as first `Z(n)<0.8` against the `L=8` result before attempting dense-Floquet `L=12`.
- Expand the reading report's broader literature background if the final deliverable requires a more formal review section.

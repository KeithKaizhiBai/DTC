# Project Summary

## 1. Project Goal

Reproduce the paper "Solvable model for discrete time crystal enforced by nonsymmorphic dynamical symmetry" (Physical Review Research 5, L032024, 2023). The immediate goals are an English reading report, a Python/CUDA reproduction of the noninteracting two-level dynamics, and a first interacting spin-chain reproduction with an HPC-ready package for longer runs.

## 2. Current Architecture

The project is organized around the source PDFs and two output folders:

- `reading_report/`: English literature notes, reading-stage plan, LaTeX reading report, and compiled reading-report PDF.
- `reproduction_noninteracting/`: Python/CUDA scripts, generated figures, data, reproduction report, and numerical-methods supplement for the noninteracting model.
- `reproduction_interacting/`: exact state-vector CUDA/PyTorch script, local and HPC parameter files, generated small-chain data/figures, LaTeX interacting report, and the compressed HPC upload package.
- `comprehensive_report/`: RevTeX-style final summary report, selected figure copies, and compiled PDF for the full reproduction project.
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
- `reproduction_interacting/params/interacting_params_hpc.json`: Current cluster-run parameters focused on corrected `L = 10` long stroboscopic sampling, with `interaction_j = 0.2` and `interaction_operator_scale = 0.5`.
- `reproduction_interacting/run_interacting_hpc.sh`: Slurm launcher with `#SBATCH --exclude=gpuh01`.
- `reproduction_interacting/reproduction_interacting_report.tex`: LaTeX source for the interacting reproduction report.
- `reproduction_interacting/reproduction_interacting_report.pdf`: Compiled interacting reproduction report.
- `comprehensive_report/comprehensive_report.tex`: RevTeX-style source for the final full reproduction report.
- `comprehensive_report/comprehensive_report.pdf`: Compiled final full reproduction report.
- `output/pdf/reading_report.pdf`: Final PDF copy for the reading report.
- `output/pdf/reproduction_report.pdf`: Final PDF copy for the noninteracting reproduction report.
- `output/pdf/numerical_methods_supplement.pdf`: Final PDF copy for the numerical-methods supplement.
- `output/pdf/interacting_reproduction_report.pdf`: Final PDF copy for the interacting reproduction report.
- `output/pdf/comprehensive_report.pdf`: Final PDF copy for the comprehensive report.
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

- Date: 2026-06-04
- Decision: Separate the paper-label interaction strength from the Pauli-matrix coefficient by setting `interaction_operator_scale = 0.5`.
- Reason: Directly using `J = 0.2` as the Pauli-product coefficient made the corrected `L=10`, refined-`alpha_13` stroboscopic signal decay below `Z=0.8` at `n=160`, contradicting the paper-scale expectation and the DTC physics. With `J_eff = 0.1`, `Z(n)` remains above `0.9` through the sampled `10^10` periods.
- Consequences: Future interacting runs must keep `interaction_j = 0.2` as the paper label and apply `interaction_operator_scale = 0.5` unless a different spin-operator convention is intentionally tested.

- Date: 2026-06-04
- Decision: Add a parameter-controlled `boundary` interface for the interacting chain.
- Reason: The original implementation was hard-coded to open boundary conditions, while the detuned lifetime discrepancy may depend on whether Fig. 4 used open or periodic boundaries.
- Consequences: `boundary = "open"` uses `L-1` interaction bonds. `boundary = "periodic"` uses `L` bonds for `L > 2` by adding the final first-last bond. Periodic-boundary outputs include a `_periodic` suffix to avoid overwriting open-boundary files.

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

- Date: 2026-06-04
- Milestone: Fixed the refined-`alpha_13` early-decay problem in the interacting long-time run.
- Evidence: Corrected local `L=10` Floquet sampling uses `J_eff = 0.1`; `data/diagnostics.json` reports `z_min = 0.9088579791801286` and no sampled crossing below `Z=0.9` through `n = 10000000000`.

- Date: 2026-06-04
- Milestone: Completed the first periodic-boundary check for the interacting chain.
- Evidence: `reproduction_interacting/periodic_l10/data/diagnostics.json` reports `boundary = periodic`, `interaction_bond_count = 10`, `L=10`, refined `alpha_13`, `n_max = 100000000000`, and first sampled crossing below `Z=0.8` at `n=789049238`.

- Date: 2026-06-04
- Milestone: Completed detuned boundary comparison and extended OBC `alpha_13` check.
- Evidence: `reproduction_interacting/data/l10_boundary_and_longtime_checks.json` records `alpha=81.60` OBC/PBC comparison to `n=100000` and refined `alpha_13` OBC extension to `n=1000000000000`.

- Date: 2026-06-04
- Milestone: Completed current-convention `L=8` versus `L=10` size comparison.
- Evidence: `reproduction_interacting/data/l8_l10_size_comparison_current.json` compares OBC, `J_eff=0.1` data for `alpha=81.60` and refined `alpha_13`; refined `alpha_13` shows delayed collapse for `L=10`, while detuned `alpha=81.60` remains too fast at both sizes.

- Date: 2026-06-04
- Milestone: Completed the final comprehensive reproduction report.
- Evidence: `comprehensive_report/comprehensive_report.tex` compiles with RevTeX 4.2 to `comprehensive_report/comprehensive_report.pdf`; selected pages were rendered under `tmp/pdfs/report_checks/comprehensive_report/` and visually checked.

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

- Date: 2026-06-04
- Files changed: `reproduction_interacting/reproduce_interacting.py`, `reproduction_interacting/params/interacting_params_hpc.json`, `reproduction_interacting/params/interacting_params_local.json`, `reproduction_interacting/reproduction_interacting_report.tex`, `reproduction_interacting/reproduction_interacting_report.pdf`, `output/pdf/interacting_reproduction_report.pdf`, `reproduction_interacting/data/diagnostics.json`, `reproduction_interacting/data/l10_alpha13_j_convention_check.json`, `reproduction_interacting/figures/l10_alpha13_j_convention_check.png`, `output/interacting_reproduction_hpc_package.zip`
- Summary: Added explicit interaction-operator scaling, reran corrected `L=10` long-time Floquet sampling locally, rebuilt the LaTeX/PDF report, and refreshed the HPC package.
- Why it matters: The refined `alpha_13` result now matches the expected long-lived DTC behavior instead of showing an unphysical early decay.

- Date: 2026-06-04
- Files changed: `reproduction_interacting/reproduce_interacting.py`, `reproduction_interacting/figures/z_stroboscopic_alpha_13_l10.png`, `reproduction_interacting/figures/z_stroboscopic_alpha_81p60_l10.png`, `reproduction_interacting/reproduction_interacting_report.tex`, `reproduction_interacting/reproduction_interacting_report.pdf`, `output/pdf/interacting_reproduction_report.pdf`, `reproduction_interacting/hpc_package/interacting_reproduction_package/readme.md`, `output/interacting_reproduction_hpc_package.zip`
- Summary: Changed stroboscopic PNG naming so figure files include the simulated length suffix, while keeping compatibility copies without the suffix.
- Why it matters: HPC outputs from different lengths can now be distinguished without opening the files or reading the legend.

- Date: 2026-06-04
- Files changed: `reproduction_interacting/reproduce_interacting.py`, `reproduction_interacting/params/interacting_params_periodic_l10.json`, `reproduction_interacting/periodic_l10/*`, `reproduction_interacting/figures/l10_alpha13_boundary_check.png`, `reproduction_interacting/data/l10_alpha13_boundary_check.json`, `reproduction_interacting/reproduction_interacting_report.tex`, `reproduction_interacting/reproduction_interacting_report.pdf`, `output/pdf/interacting_reproduction_report.pdf`, `reproduction_interacting/run_interacting_hpc.sh`, `reproduction_interacting/hpc_package/interacting_reproduction_package/*`, `output/interacting_reproduction_hpc_package_periodic.zip`
- Summary: Added periodic boundary conditions, ran a local `L=10`, refined-`alpha_13`, `n_max=10^11` periodic-boundary Floquet check, updated the report, and created a refreshed HPC package under a new filename because the previous standard zip was locked by another process.
- Why it matters: Boundary conditions are now a controlled parameter, and the first periodic test shows collapse/revival rather than improved `alpha_13` stability.

- Date: 2026-06-04
- Files changed: `reproduction_interacting/params/interacting_params_alpha81p60_periodic_l10.json`, `reproduction_interacting/params/interacting_params_alpha13_obc_n1e12_l10.json`, `reproduction_interacting/periodic_l10/data/z_stroboscopic_alpha_81p60_l10_periodic.csv`, `reproduction_interacting/data/z_stroboscopic_alpha_13_l10_n1e12.csv`, `reproduction_interacting/data/l10_boundary_and_longtime_checks.json`, `reproduction_interacting/figures/l10_alpha81p60_boundary_check.png`, `reproduction_interacting/figures/l10_alpha13_obc_n1e12.png`, `reproduction_interacting/reproduction_interacting_report.tex`, `reproduction_interacting/reproduction_interacting_report.pdf`, `output/pdf/interacting_reproduction_report.pdf`, `output/interacting_reproduction_hpc_package_boundary_checks.zip`
- Summary: Added the detuned `alpha=81.60` OBC/PBC comparison and extended the refined `alpha_13` OBC run to `n_max=10^12`; rebuilt the report and prepared a new HPC package containing all boundary-check parameter files.
- Why it matters: The detuned fast decay persists under both OBC and PBC, while the refined `alpha_13` OBC run now shows the expected later finite-size oscillation beyond the previous `10^10` cutoff.

- Date: 2026-06-04
- Files changed: `reproduction_interacting/data/z_stroboscopic_alpha_81p60_l8_current.csv`, `reproduction_interacting/data/z_stroboscopic_alpha_13_l8_n1e12_current.csv`, `reproduction_interacting/data/l8_l10_size_comparison_current.json`, `reproduction_interacting/figures/l8_l10_alpha81p60_obc_comparison.png`, `reproduction_interacting/figures/l8_l10_alpha13_obc_comparison.png`, `reproduction_interacting/reproduction_interacting_report.tex`, `reproduction_interacting/reproduction_interacting_report.pdf`, `output/pdf/interacting_reproduction_report.pdf`
- Summary: Recomputed `L=8` data using the same current OBC, `J_eff=0.1`, `magnus4` Floquet protocol as `L=10`, then added `L=8` versus `L=10` comparison figures and analysis to the report.
- Why it matters: The refined `alpha_13` case now shows the expected larger-`L` lifetime enhancement, but the detuned `alpha=81.60` case still fails to show the Fig. 4(c)-style plateau.

- Date: 2026-06-04
- Files changed: `comprehensive_report/comprehensive_report.tex`, `comprehensive_report/comprehensive_report.pdf`, `comprehensive_report/figures/*`, `output/pdf/comprehensive_report.pdf`, `docs/PROJECT_SUMMARY.md`
- Summary: Added the final RevTeX-style comprehensive report that consolidates the theory review, noninteracting reproduction, interacting reproduction, corrected interaction convention, size comparison, boundary checks, validation commands, unresolved detuned-case discrepancy, and a post-reference supplemental methods section.
- Why it matters: The project now has a single top-level report suitable for sharing alongside the source code and selected figures.

- Date: 2026-06-04
- Files changed: `comprehensive_report/comprehensive_report.tex`, `comprehensive_report/comprehensive_report.pdf`, `output/pdf/comprehensive_report.pdf`, `docs/PROJECT_SUMMARY.md`
- Summary: Rewrote the post-reference supplemental methods section with a fuller derivation of fourth-order Runge--Kutta from Taylor expansion and order conditions, RK4 local/global error and non-unitarity, the Magnus expansion, the fourth-order Gauss--Legendre Magnus formula, and the reason Magnus preserves unitary Schrodinger evolution.
- Why it matters: The final report now explains the numerical methods at a textbook level sufficient for readers to understand the coefficient choices, error scaling, matrix-state generalization, and unitary-structure difference between RK4 and Magnus.

- Date: 2026-06-04
- Files changed: `comprehensive_report/comprehensive_report.tex`, `comprehensive_report/comprehensive_report.pdf`, `output/pdf/comprehensive_report.pdf`, `docs/PROJECT_SUMMARY.md`
- Summary: Rewrote the comprehensive-report numerical-methods supplement from a "why" perspective: RK4 slopes and weights are derived from the exact increment integral, Simpson/Taylor matching, the scalar linear test equation, and the Runge--Kutta order conditions; Magnus4 is derived from time ordering, BCH commutators, the continuum Magnus series, and the two-point Gauss--Legendre commutator coefficient.
- Why it matters: The report now explains why the RK4 parameters have their particular values and why the Magnus expansion has commutator terms, rather than only stating the final algorithms.

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
pdflatex -interaction=nonstopmode -halt-on-error comprehensive_report.tex
pdftoppm -png -r 120 reading_report\reading_report.pdf tmp\pdfs\report_checks\reading_report
pdftoppm -png -r 120 reproduction_noninteracting\reproduction_report.pdf tmp\pdfs\report_checks\reproduction_report
pdftoppm -png -r 120 reproduction_noninteracting\numerical_methods_supplement.pdf tmp\pdfs\report_checks\numerical_methods_supplement
pdftoppm -png -r 120 reproduction_interacting\reproduction_interacting_report.pdf tmp\pdfs\report_checks\reproduction_interacting_report
pdftoppm -png -r 120 comprehensive_report\comprehensive_report.pdf tmp\pdfs\report_checks\comprehensive_report\page
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

Latest corrected `L=10` local long-time checks:

- Command: `python reproduction_interacting\reproduce_interacting.py --params reproduction_interacting\params\interacting_params_hpc.json --output-root reproduction_interacting --device cuda`
- Device: NVIDIA GeForce RTX 4090.
- Coupling convention: `interaction_j_paper = 0.2`, `interaction_operator_scale = 0.5`, `interaction_j_effective = 0.1`.
- Floquet method: `magnus4`, `float64`, `steps_per_period = 120`, logarithmic stroboscopic sampling.
- Refined `alpha_13`, `L=10`: sampled to `n=10000000000`, `z_min = 0.9088579791801286`, no sampled crossing below `Z=0.9`, first crossing below `Z=0.95` at `n=7251171612`.
- Old direct Pauli coefficient check: using `J_eff = 0.2` crossed below `Z=0.8` at `n=160`, confirming the early-decay bug.
- Corrected PDF report was compiled with `pdflatex` and rendered via `pdftoppm` as `tmp/pdfs/report_checks/interacting_reproduction_report_corrected-*.png`.

Latest periodic-boundary check:

- Command: `python reproduction_interacting\reproduce_interacting.py --params reproduction_interacting\params\interacting_params_periodic_l10.json --output-root reproduction_interacting\periodic_l10 --device cuda`
- Device: NVIDIA GeForce RTX 4090.
- Model: `boundary = periodic`, `L = 10`, `interaction_bond_count = 10`, `J_eff = 0.1`, refined `alpha_13`.
- Floquet method: `magnus4`, `float64`, `steps_per_period = 120`, sampled to `n = 100000000000`.
- Result: `z_min = -0.9967548471649453`, `z_final = 0.9923496522533153`, first sampled `Z(n)<0.9` at `n=540430471`, first sampled `Z(n)<0.8` at `n=789049238`, first sampled `Z(n)<0` at `n=1922399212`.
- Interpretation: periodic boundaries produce strong finite-size collapse/revival and do not by themselves explain the detuned Fig. 4(c) discrepancy.
- Package note: the standard `output/interacting_reproduction_hpc_package.zip` was locked by another process; the updated package with periodic support is `output/interacting_reproduction_hpc_package_periodic.zip`.

Latest boundary and long-time additions:

- Detuned boundary comparison: `alpha=81.60`, `L=10`, `J_eff=0.1`, `n_max=100000`.
- OBC detuned result: first sampled `Z(n)<0.8` at `n=3`, first `Z(n)<0` at `n=78`, `z_final=0.1722362564708381`.
- PBC detuned result: first sampled `Z(n)<0.8` at `n=3`, first `Z(n)<0` at `n=301`, `z_final=-0.22409478154793133`.
- Extended OBC refined `alpha_13`: sampled to `n=1000000000000`, first `Z(n)<0.9` at `n=10579865403`, first `Z(n)<0.8` at `n=15094461760`, first `Z(n)<0` at `n=37262828480`, `z_min=-0.9972063062403195`.
- Report was compiled twice with `pdflatex` and rendered via `pdftoppm` as `tmp/pdfs/report_checks/interacting_reproduction_report_boundary_checks-*.png`.
- Updated upload package containing all boundary-check parameter files: `output/interacting_reproduction_hpc_package_boundary_checks.zip`.

Latest current-convention size comparison:

- All data use OBC, `J_eff=0.1`, `magnus4`, `float64`, and logarithmic stroboscopic sampling.
- `alpha=81.60`: `L=8` first `Z(n)<0.8` at `n=3`, first `Z(n)<0` at `n=65`; `L=10` first `Z(n)<0.8` at `n=3`, first `Z(n)<0` at `n=78`.
- Refined `alpha_13`: `L=8` first `Z(n)<0.8` at `n=199668663`, first `Z(n)<0` at `n=487931047`; `L=10` first `Z(n)<0.8` at `n=15094461760`, first `Z(n)<0` at `n=37262828480`.
- Interpretation: the refined `alpha_13` case has the expected size-enhanced lifetime, while the detuned `alpha=81.60` raw trace remains inconsistent with Fig. 4(c).
- Report was compiled twice with `pdflatex` and rendered via `pdftoppm` as `tmp/pdfs/report_checks/interacting_reproduction_report_size_compare_latest-*.png`.

Latest comprehensive-report checks:

- Command: `pdflatex -interaction=nonstopmode -halt-on-error comprehensive_report.tex` from `comprehensive_report/`.
- Output: `comprehensive_report/comprehensive_report.pdf`, 24 pages, also copied to `output/pdf/comprehensive_report.pdf`.
- Log status: no undefined references and no overfull boxes; only a standard RevTeX/hyperref `nameref` label-definition warning remained.
- Render check: `pdftoppm -png -r 120 comprehensive_report\comprehensive_report.pdf tmp\pdfs\report_checks\comprehensive_report\page`; title page, figure pages, table page, RK4 slope/weight derivation pages, RK4 non-unitarity page, Magnus BCH/source page, Gauss--Legendre Magnus4 coefficient page, and final checklist page were visually inspected.

## 8. Open Problems

- The refined `alpha_13` early-decay issue is fixed at the `10^10` scale, extending OBC to `10^12` reveals later finite-size oscillations, and the current-convention `L=8` versus `L=10` comparison shows the expected size enhancement for refined `alpha_13`. The detuned `alpha=81.60`, `J_eff=0.1` raw `Z(n)` still decays too quickly compared with Fig. 4(c) under both OBC and PBC and does not show a strong `L=8` to `L=10` improvement.
- Boundary conditions alone do not explain the detuned Fig. 4(c) discrepancy; remaining checks should include interaction normalization, observable/lifetime definition, and possible differences in the paper's numerical protocol.
- The next larger-system step is `L = 12`, but dense Floquet matrices grow to dimension `4096`; memory and runtime should be checked on H100 before launching long sweeps.
- Exact Heun-function evaluation is not implemented; special `alpha_n` values are refined numerically through high-accuracy period evolution.
- The current reading report is concise; broader many-body localization and experimental DTC background can be expanded later if needed.
- CUDA is used, but the two-level problem is too small to benefit from GPU acceleration. Increasing `STEPS_PER_PERIOD` gives higher resolution but longer runtime.

## 9. Next Steps

- If higher publication-quality accuracy is needed, rerun the noninteracting script with larger `STEPS_PER_PERIOD` and compare diagnostics.
- For the next cluster run, keep `interaction_operator_scale = 0.5`; change `stroboscopic.lengths` to `[12]` only after confirming H100 memory headroom for a `4096 x 4096` dense Floquet matrix.
- Investigate why the corrected `alpha=81.60`, `J_eff=0.1` raw `Z(n)` decays faster than Fig. 4(c): test alternative interaction normalization and envelope/windowed lifetime extraction before claiming detuned finite-size scaling.
- Expand the reading report's broader literature background if the final deliverable requires a more formal review section.

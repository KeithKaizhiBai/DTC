# Noninteracting Reproduction Stage Plan

## Objective

Reproduce the noninteracting two-level dynamics of the paper with Python/CUDA. The main target is the `sigma_x(t)` trace and spectrum analogous to Fig. 4(a), plus solver validation between a fourth-order unitary Magnus method and RK4.

## Numerical Steps

1. Use the dimensionless convention `hbar = 1`, `omega = 1`, `T = 2*pi`, and `Omega = alpha/8`.
2. Implement the two-level Hamiltonian

   ```text
   H0(t) = 0.5*Omega*sin(t)*sigma_x + Omega*sin(t/2)^2*sigma_y.
   ```

3. Use the initial state `|+x> = (1, 1)/sqrt(2)`.
4. Implement a fourth-order Magnus step with two Gauss-Legendre nodes and matrix exponential propagation. This keeps the propagator unitary.
5. Implement a fourth-order Runge-Kutta step for comparison.
6. Refine the special value near `alpha_13` by minimizing the one-period return probability using high-accuracy SciPy integration.
7. Run 40 periods for `alpha = 81.60`, the detuned value used in the paper's noninteracting Fig. 4(a).
8. Save time traces, Fourier spectra, stroboscopic projections, diagnostic JSON, and publication-style PNG figures.
9. Compare Magnus and RK4 traces through maximum norm drift and maximum observable difference.

## Deliverables

- `reproduction_noninteracting/reproduce_noninteracting.py`
- `reproduction_noninteracting/data/*.csv`
- `reproduction_noninteracting/data/diagnostics.json`
- `reproduction_noninteracting/figures/*.png`
- `reproduction_noninteracting/reproduction_report.md`

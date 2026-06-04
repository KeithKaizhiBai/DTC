# Interacting Reproduction Stage Plan

## Objective

Reproduce the interacting spin-chain dynamics in the nonsymmorphic DTC paper for small system sizes on the local RTX 4090 GPU, then prepare a self-contained HPC calculation package for larger simulations.

## Model

Use the paper's interacting Hamiltonian

```text
H(t) = sum_i H0_i(t) + J sum_i sigma_x^i sigma_x^{i+1}
```

with

```text
H0_i(t) = 0.5*Omega*sin(t)*sigma_x^i + Omega*sin(t/2)^2*sigma_y^i
alpha = 8*Omega
hbar = omega = 1
J = 0.2
```

The main paper target is `mu=x`, initial state polarized in the `+x` direction, `alpha=81.60`, and first 40 periods. The refined special value `alpha_13=80.04624211689759` is also included for comparison.

## Steps

1. Implement exact state-vector dynamics for small chains using CUDA tensors.
2. Use explicit Hamiltonian-action functions instead of building dense matrices.
3. Compute average magnetization `m_x(t) = (1/L) sum_i <sigma_x^i(t)>`.
4. Compute Fourier spectra of `m_x(t)`.
5. Compute stroboscopic `Z(n)=(-1)^n m_x(nT)` for several small chain lengths.
6. Generate figures for the first 40 periods and stroboscopic finite-size comparison.
7. Write a LaTeX/PDF reproduction report.
8. Prepare a compressed HPC package containing scripts, parameters, bash launcher, logs and results folders. The launcher must exclude `gpuh01`.

## Local Validation Targets

- Confirm CUDA is used on the RTX 4090.
- Confirm `alpha=81.60`, `L=10`, `J=0.2` produces a smoother period-doubled `m_x(t)` response than the noninteracting detuned result.
- Confirm the spectrum has a strong feature near half the drive frequency.
- Confirm small-length stroboscopic `Z(n)` remains positive for an initial time window.

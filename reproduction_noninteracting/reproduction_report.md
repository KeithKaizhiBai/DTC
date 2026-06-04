# Noninteracting Reproduction Report

## Scope

This report documents the first reproduction stage for the noninteracting two-level model in Hu, Fu, Li, and Shen, Physical Review Research 5, L032024 (2023). The target is the single-spin dynamics behind the paper's Fig. 4(a), not the interacting spin-chain results.

## Model and Units

The reproduced Hamiltonian is

```text
H0(t) = (1/2) Omega sin(t) sigma_x + Omega sin^2(t/2) sigma_y,
```

with `hbar = 1` and drive angular frequency `omega = 1`. Therefore,

```text
T = 2*pi
alpha = 8*Omega
Omega = alpha/8.
```

The initial state is

```text
|+x> = (1, 1) / sqrt(2).
```

The main observable is

```text
<sigma_x(t)> = <psi(t)|sigma_x|psi(t)> / <psi(t)|psi(t)>.
```

For the chosen initial state `|+x>`, the `x`-spin autocorrelation

```text
C_x(t) = <+x|U(t)^dagger sigma_x U(t) sigma_x|+x>
```

equals `<sigma_x(t)>`, because `sigma_x|+x> = |+x>`. Thus the plotted trace also serves as the requested single-spin `x` autocorrelation for this initial condition.

## Numerical Methods

Two fourth-order propagators are implemented in `reproduce_noninteracting.py`.

The primary solver is a fourth-order Magnus method using two Gauss-Legendre nodes in each time step:

```text
A(t) = -i H(t)
Omega_M = dt/2*(A1 + A2) - sqrt(3)*dt^2/12*[A1, A2]
psi(t + dt) = exp(Omega_M) psi(t).
```

Because `Omega_M` is anti-Hermitian for a Hermitian Hamiltonian, this step is unitary up to floating-point error.

The comparison solver is standard RK4 applied directly to the Schrodinger equation. RK4 is fourth order but not exactly unitary, so the script reports its norm drift.

## Parameters

The paper lists `alpha_13 = 80.07` in Table I. A high-accuracy SciPy period-evolution minimization in the local convention refines the nearby zero of the one-period return probability to approximately

```text
alpha_13 = 80.0462421169.
```

The detuned noninteracting Fig. 4(a) target uses

```text
alpha = 81.60.
```

## Generated Outputs

Running

```powershell
python reproduction_noninteracting\reproduce_noninteracting.py
```

generates:

- `data/sigma_x_alpha_81p60.csv`
- `data/spectrum_alpha_81p60.csv`
- `data/stroboscopic_alpha_13.csv`
- `data/diagnostics.json`
- `figures/sigma_x_alpha_81p60.png`
- `figures/magnus_rk4_comparison_alpha_81p60.png`
- `figures/stroboscopic_alpha_13.png`

## Validation Run

The latest run used the local RTX 4090 CUDA device:

```text
device: cuda
cuda_device_name: NVIDIA GeForce RTX 4090
alpha_13_refined: 80.04624211689759
alpha_13_one_period_return_probability: 3.57e-18
alpha_detuned: 81.60
periods_to_plot: 40
steps_per_period: 80
rk4_effective_steps_per_period: 320
magnus4_max_norm_drift: 3.92e-05
rk4_max_norm_drift: 1.79e-03
max_abs_sigma_x_difference_magnus4_vs_rk4: 6.50e-03
dominant_spectrum_omega_over_drive: 0.5498
```

The RK4 solver uses four substeps per plotted Magnus step. This keeps the plotted sampling grid fixed while making RK4 a meaningful comparison for this rapidly oscillating detuned case.

## Expected Physical Behavior

For `alpha = alpha_13`, stroboscopic projections should collapse to two alternating points, expressing the `2T` period of the state. For `alpha = 81.60`, which is detuned from `alpha_13`, the noninteracting `sigma_x(t)` trace should show rapid irregular oscillations and multiple spectral peaks. This is the behavior emphasized in the paper before interactions stabilize the smoother subharmonic response.

The generated figures match these expectations: the detuned `sigma_x(t)` curve has rapid non-sinusoidal oscillations with several Fourier features, while the refined `alpha_13` stroboscopic projection collapses to the two expected points.

## Current Limitations

The script does not evaluate confluent Heun functions directly. Instead, it refines `alpha_13` through high-accuracy numerical period evolution. This is sufficient for the present reproduction stage because the final dynamics are generated from the Hamiltonian itself.

The interacting many-body Hamiltonian and finite-size lifetime analysis are not implemented yet.

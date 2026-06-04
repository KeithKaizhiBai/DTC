# Reading Report

## Paper

Zi-Ang Hu, Bo Fu, Xiao Li, and Shun-Qing Shen, "Solvable model for discrete time crystal enforced by nonsymmorphic dynamical symmetry," Physical Review Research 5, L032024 (2023).

## Background

A time crystal is a phase or dynamical regime where time-translation symmetry is broken in the observed motion. For continuously time-independent equilibrium systems this is forbidden, but periodically driven Floquet systems can show discrete time-translation symmetry breaking. In a discrete time crystal (DTC), the Hamiltonian has period `T`, while robust observables respond with a larger period, commonly `2T`.

Earlier DTC work often relied on disorder, many-body localization, or prethermal regimes to prevent generic driven systems from heating to featureless states. The present paper has a different starting point. It asks whether a nonsymmorphic dynamical symmetry can enforce period extension already at the single-particle or single-spin level. Interactions then play a stabilizing role by suppressing fast off-diagonal oscillations and producing a cleaner prethermal DTC response.

Many-body localization (MBL) is important in the DTC literature because it can prevent Floquet heating in disordered interacting systems. In an MBL Floquet DTC, local memory survives for long times, and a subharmonic response can remain rigid against small drive imperfections. However, strict MBL is sensitive to dimensionality, disorder assumptions, and finite-size limitations. This motivates prethermal DTCs, where a high-frequency or otherwise constrained drive produces an exponentially long-lived regime before eventual heating. The Hu-Fu-Li-Shen model is not primarily an MBL construction. Its distinctive point is that a nonsymmorphic dynamical symmetry creates a single-spin period extension mechanism first, and interactions are then used to stabilize observable subharmonic motion rather than to create the basic Mobius twist.

## Central Model

The main two-level Hamiltonian is

```text
H0(t) = (1/2) hbar Omega sin(omega t) sigma_x
      + hbar Omega sin^2(omega t / 2) sigma_y.
```

The drive period is

```text
T = 2*pi/omega.
```

The paper defines the dimensionless bandwidth-to-frequency ratio

```text
alpha = 8*Omega/omega.
```

The reproduction uses `hbar = 1` and `omega = 1`, hence `Omega = alpha/8` and `T = 2*pi`.

The instantaneous eigenstates are

```text
phi_+(t) = ( exp(-i omega t/2), 1 ) / sqrt(2)
phi_-(t) = ( -exp(-i omega t/2), 1 ) / sqrt(2)
```

up to the sign convention used in the paper. At `t = 0`, `phi_+(0)` is the `+x` polarized state.

## Nonsymmorphic Dynamical Symmetry

The Hamiltonian has a dynamical glide symmetry. The symmetry operator has the same period as the Hamiltonian, but its eigenvalues and eigenstates interchange after one drive period and return only after two periods. This is the Mobius twist of the instantaneous eigenstate bundle.

The key physical point is that the Hamiltonian itself remains `T` periodic, while the instantaneous eigenstates carry a representation with period `2T`. This mismatch can enforce a subharmonic response without requiring an ordinary avoided-crossing adiabatic picture.

## Exact Solution and Special Parameters

In the instantaneous eigenbasis, the coefficients obey a second-order differential equation whose solution is written using confluent Heun functions. The probability for an initial `phi_+(0)` state to remain in `phi_+(0)` after one period is controlled by a Heun function value at `x = pi/2`.

There is a discrete sequence `alpha_n` where this one-period return probability vanishes. At those points, the state after one period is mapped to the partner state up to a phase, and after two periods it returns to the initial state with a Berry phase of `pi`.

The paper lists approximate values:

```text
n:       1     2      3      4      5      6      7
alpha:  4.21  10.73  17.11  23.44  29.75  36.07  42.36

n:       8      9      10     11     12     13     14
alpha:  48.66  54.95  61.24  67.46  73.76  80.07  86.34
```

For large `n`, the asymptotic estimate is approximately

```text
alpha_n ~= 2*pi*n - pi/2.
```

For the reproduction target near Fig. 4(a), the important values are `alpha_13` near `80.07` and the detuned value `alpha = 81.60`.

## Observables

The noninteracting reproduction focuses on

```text
<sigma_x(t)> = <psi(t)|sigma_x|psi(t)>
```

for the initial state polarized in the `+x` direction. In the noninteracting model, when `alpha` is detuned from an `alpha_n`, the off-diagonal term in the instantaneous basis produces rapid oscillations through the dynamical phase factor

```text
exp[-i alpha sin^2(omega t/4)].
```

This causes a complicated `sigma_x(t)` trace and multiple Fourier peaks, as shown in Fig. 4(a). In the interacting case, the paper argues that interaction-induced decoherence suppresses these rapid oscillations and stabilizes a smoother period-doubled average magnetization.

## Reproduction Targets

The first numerical stage should reproduce the noninteracting dynamics:

1. Implement the two-level Hamiltonian directly.
2. Use the `+x` polarized initial state.
3. Evolve over the first 40 drive periods.
4. Compute `sigma_x(t)` and its Fourier spectrum for `alpha = 81.60`.
5. Compare a unitary fourth-order Magnus solver with a fourth-order Runge-Kutta solver.
6. Verify that the calibrated `alpha_13` produces near-perfect two-period stroboscopic behavior.

The interacting many-body figures are not part of this first reproduction stage.

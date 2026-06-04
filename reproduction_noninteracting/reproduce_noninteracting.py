import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import torch
from scipy.integrate import solve_ivp
from scipy.optimize import minimize_scalar


HBAR = 1.0
DRIVE_OMEGA = 1.0
PERIOD = 2.0 * np.pi / DRIVE_OMEGA
ALPHA_TABLE_13 = 80.07
ALPHA_DETUNED = 81.60
PERIODS_TO_PLOT = 40
STEPS_PER_PERIOD = 80
RK4_SUBSTEPS = 4
STROBOSCOPIC_PERIODS = 60
FRACTIONAL_SAMPLES_PER_PERIOD = 4

ROOT_SEARCH_BOUNDS = (79.7, 80.4)

SCRIPT_DIR = Path(__file__).resolve().parent
DATA_DIR = SCRIPT_DIR / "data"
FIGURE_DIR = SCRIPT_DIR / "figures"


def torch_device():
    return torch.device("cuda" if torch.cuda.is_available() else "cpu")


def torch_dtypes(device):
    if device.type == "cuda":
        return torch.float32, torch.complex64
    return torch.float64, torch.complex128


def pauli_matrices(device):
    _, complex_dtype = torch_dtypes(device)
    sx = torch.tensor([[0.0, 1.0], [1.0, 0.0]], dtype=complex_dtype, device=device)
    sy = torch.tensor([[0.0, -1.0j], [1.0j, 0.0]], dtype=complex_dtype, device=device)
    sz = torch.tensor([[1.0, 0.0], [0.0, -1.0]], dtype=complex_dtype, device=device)
    return sx, sy, sz


def hamiltonian_torch(t, alpha, sx, sy):
    real_dtype, _ = torch_dtypes(sx.device)
    omega_amp = alpha / 8.0
    t_tensor = torch.as_tensor(t, dtype=real_dtype, device=sx.device)
    return (
        0.5 * omega_amp * torch.sin(t_tensor) * sx
        + omega_amp * torch.sin(t_tensor / 2.0) ** 2 * sy
    )


def generator_torch(t, alpha, sx, sy):
    return -1.0j * hamiltonian_torch(t, alpha, sx, sy) / HBAR


def magnus4_step(state, t, dt, alpha, sx, sy):
    c1 = 0.5 - np.sqrt(3.0) / 6.0
    c2 = 0.5 + np.sqrt(3.0) / 6.0
    a1 = generator_torch(t + c1 * dt, alpha, sx, sy)
    a2 = generator_torch(t + c2 * dt, alpha, sx, sy)
    commutator = a1 @ a2 - a2 @ a1
    omega_m = 0.5 * dt * (a1 + a2) - np.sqrt(3.0) * dt * dt * commutator / 12.0
    return apply_traceless_antihermitian_exponential(omega_m, state)


def apply_traceless_antihermitian_exponential(omega_m, state):
    theta_squared = torch.real(-0.5 * torch.trace(omega_m @ omega_m))
    theta = torch.sqrt(torch.clamp(theta_squared, min=0.0))
    identity = torch.eye(2, dtype=state.dtype, device=state.device)
    small = theta < 1e-12
    scale = torch.where(small, 1.0 - theta_squared / 6.0, torch.sin(theta) / theta)
    unitary = torch.cos(theta) * identity + scale * omega_m
    return unitary @ state


def rk4_step(state, t, dt, alpha, sx, sy):
    k1 = generator_torch(t, alpha, sx, sy) @ state
    k2 = generator_torch(t + 0.5 * dt, alpha, sx, sy) @ (state + 0.5 * dt * k1)
    k3 = generator_torch(t + 0.5 * dt, alpha, sx, sy) @ (state + 0.5 * dt * k2)
    k4 = generator_torch(t + dt, alpha, sx, sy) @ (state + dt * k3)
    return state + dt * (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0


def normalized_expectation(state, operator):
    numerator = torch.vdot(state, operator @ state)
    denominator = torch.vdot(state, state)
    return torch.real(numerator / denominator)


def evolve_time_trace(alpha, periods, steps_per_period, method, device):
    sx, sy, _ = pauli_matrices(device)
    dt = PERIOD / steps_per_period
    steps = int(periods * steps_per_period)
    _, complex_dtype = torch_dtypes(device)
    state = torch.tensor([1.0, 1.0], dtype=complex_dtype, device=device) / np.sqrt(2.0)
    times = np.empty(steps + 1)
    sigma_x = np.empty(steps + 1)
    norms = np.empty(steps + 1)

    times[0] = 0.0
    sigma_x[0] = normalized_expectation(state, sx).detach().cpu().item()
    norms[0] = torch.linalg.norm(state).detach().cpu().item()

    stepper = magnus4_step if method == "magnus4" else rk4_step
    for step in range(1, steps + 1):
        t0 = (step - 1) * dt
        state = stepper(state, t0, dt, alpha, sx, sy)
        times[step] = step * dt
        sigma_x[step] = normalized_expectation(state, sx).detach().cpu().item()
        norms[step] = torch.linalg.norm(state).detach().cpu().item()

    return times, sigma_x, norms


def evolve_rk4_time_trace_with_substeps(alpha, periods, output_steps_per_period, substeps, device):
    sx, sy, _ = pauli_matrices(device)
    output_dt = PERIOD / output_steps_per_period
    dt = output_dt / substeps
    output_steps = int(periods * output_steps_per_period)
    state = torch.tensor([1.0, 1.0], dtype=torch_dtypes(device)[1], device=device) / np.sqrt(2.0)
    times = np.empty(output_steps + 1)
    sigma_x = np.empty(output_steps + 1)
    norms = np.empty(output_steps + 1)

    times[0] = 0.0
    sigma_x[0] = normalized_expectation(state, sx).detach().cpu().item()
    norms[0] = torch.linalg.norm(state).detach().cpu().item()

    for output_step in range(1, output_steps + 1):
        for substep in range(substeps):
            t0 = ((output_step - 1) * substeps + substep) * dt
            state = rk4_step(state, t0, dt, alpha, sx, sy)
        times[output_step] = output_step * output_dt
        sigma_x[output_step] = normalized_expectation(state, sx).detach().cpu().item()
        norms[output_step] = torch.linalg.norm(state).detach().cpu().item()

    return times, sigma_x, norms


def evolve_stroboscopic_projections(alpha, periods, steps_per_period, device):
    sx, sy, _ = pauli_matrices(device)
    dt = PERIOD / steps_per_period
    _, complex_dtype = torch_dtypes(device)
    state = torch.tensor([1.0, 1.0], dtype=complex_dtype, device=device) / np.sqrt(2.0)
    plus_x = torch.tensor([1.0, 1.0], dtype=complex_dtype, device=device) / np.sqrt(2.0)
    minus_x = torch.tensor([1.0, -1.0], dtype=complex_dtype, device=device) / np.sqrt(2.0)

    a_plus = []
    a_minus = []
    for n in range(periods + 1):
        a_plus.append(torch.abs(torch.vdot(plus_x, state)).detach().cpu().item())
        a_minus.append(torch.abs(torch.vdot(minus_x, state)).detach().cpu().item())
        if n == periods:
            break
        for local_step in range(steps_per_period):
            t0 = (n * steps_per_period + local_step) * dt
            state = magnus4_step(state, t0, dt, alpha, sx, sy)

    return np.arange(periods + 1), np.array(a_plus), np.array(a_minus)


def evolve_fractional_period_projections(alpha, periods, steps_per_period, samples_per_period, device):
    if steps_per_period % samples_per_period != 0:
        raise ValueError("steps_per_period must be divisible by samples_per_period")

    sx, sy, _ = pauli_matrices(device)
    sample_steps = steps_per_period // samples_per_period
    dt = PERIOD / steps_per_period
    _, complex_dtype = torch_dtypes(device)
    state = torch.tensor([1.0, 1.0], dtype=complex_dtype, device=device) / np.sqrt(2.0)
    plus_x = torch.tensor([1.0, 1.0], dtype=complex_dtype, device=device) / np.sqrt(2.0)
    minus_x = torch.tensor([1.0, -1.0], dtype=complex_dtype, device=device) / np.sqrt(2.0)

    total_samples = periods * samples_per_period
    sample_index = []
    t_over_t = []
    phase_fraction = []
    a_plus = []
    a_minus = []

    for sample in range(total_samples + 1):
        sample_index.append(sample)
        t_over_t.append(sample / samples_per_period)
        phase_fraction.append((sample % samples_per_period) / samples_per_period)
        a_plus.append(torch.abs(torch.vdot(plus_x, state)).detach().cpu().item())
        a_minus.append(torch.abs(torch.vdot(minus_x, state)).detach().cpu().item())

        if sample == total_samples:
            break
        for local_step in range(sample_steps):
            global_step = sample * sample_steps + local_step
            t0 = global_step * dt
            state = magnus4_step(state, t0, dt, alpha, sx, sy)

    return (
        np.array(sample_index),
        np.array(t_over_t),
        np.array(phase_fraction),
        np.array(a_plus),
        np.array(a_minus),
    )


def spectrum(times, values):
    dt = times[1] - times[0]
    centered = values - np.mean(values)
    window = np.hanning(len(centered))
    amplitudes = np.abs(np.fft.rfft(centered * window))
    amplitudes /= np.max(amplitudes) if np.max(amplitudes) > 0 else 1.0
    angular_frequency = 2.0 * np.pi * np.fft.rfftfreq(len(centered), d=dt)
    return angular_frequency / DRIVE_OMEGA, amplitudes


def hamiltonian_numpy(t, alpha):
    sx = np.array([[0.0, 1.0], [1.0, 0.0]], dtype=complex)
    sy = np.array([[0.0, -1.0j], [1.0j, 0.0]], dtype=complex)
    omega_amp = alpha / 8.0
    return 0.5 * omega_amp * np.sin(t) * sx + omega_amp * np.sin(t / 2.0) ** 2 * sy


def one_period_return_probability(alpha):
    plus_x = np.array([1.0, 1.0], dtype=complex) / np.sqrt(2.0)

    def rhs(t, y):
        state = y[:2] + 1.0j * y[2:]
        derivative = -1.0j * hamiltonian_numpy(t, alpha).dot(state)
        return np.r_[derivative.real, derivative.imag]

    y0 = np.r_[plus_x.real, plus_x.imag]
    solution = solve_ivp(rhs, (0.0, PERIOD), y0, method="DOP853", rtol=1e-10, atol=1e-12)
    final_state = solution.y[:2, -1] + 1.0j * solution.y[2:, -1]
    return float(abs(np.vdot(plus_x, final_state)) ** 2)


def refine_alpha_13():
    result = minimize_scalar(
        one_period_return_probability,
        bounds=ROOT_SEARCH_BOUNDS,
        method="bounded",
        options={"xatol": 1e-8},
    )
    return float(result.x), float(result.fun)


def save_trace_csv(path, times, magnus_values, rk4_values, magnus_norms, rk4_norms):
    data = np.column_stack([times / PERIOD, times, magnus_values, rk4_values, magnus_norms, rk4_norms])
    header = "t_over_T,t,sigma_x_magnus4,sigma_x_rk4,norm_magnus4,norm_rk4"
    np.savetxt(path, data, delimiter=",", header=header, comments="")


def save_spectrum_csv(path, frequency, amplitudes):
    data = np.column_stack([frequency, amplitudes])
    np.savetxt(path, data, delimiter=",", header="omega_tilde_over_omega,normalized_amplitude", comments="")


def save_fractional_projection_csv(path, sample_index, t_over_t, phase_fraction, a_plus, a_minus):
    data = np.column_stack([sample_index, t_over_t, phase_fraction, a_plus, a_minus])
    header = "sample_index,t_over_T,phase_fraction,a_plus_abs,a_minus_abs"
    np.savetxt(path, data, delimiter=",", header=header, comments="")


def plot_noninteracting_trace(times, sigma_x, frequency, amplitudes, alpha, path):
    fig, axes = plt.subplots(2, 1, figsize=(8.0, 5.8), constrained_layout=True)
    axes[0].plot(times / PERIOD, sigma_x, lw=0.8)
    axes[0].set_xlim(0, PERIODS_TO_PLOT)
    axes[0].set_ylim(-1.1, 1.1)
    axes[0].set_xlabel("t / T")
    axes[0].set_ylabel("<sigma_x>")
    axes[0].set_title(f"Noninteracting dynamics, alpha = {alpha:.2f}")

    axes[1].plot(frequency, amplitudes, lw=0.8)
    axes[1].set_xlim(0, 3)
    axes[1].set_xlabel("omega_tilde / omega")
    axes[1].set_ylabel("normalized amplitude")
    axes[1].set_title("Fourier spectrum")
    fig.savefig(path, dpi=220)
    plt.close(fig)


def plot_trace_comparison(case_results, path):
    fig, axes = plt.subplots(2, 1, figsize=(8.2, 5.8), constrained_layout=True)
    for result in case_results:
        axes[0].plot(
            result["times"] / PERIOD,
            result["sigma_x_magnus"],
            lw=0.8,
            label=result["display_name"],
        )
        axes[1].plot(
            result["frequency"],
            result["amplitudes"],
            lw=0.8,
            label=result["display_name"],
        )
    axes[0].set_xlim(0, PERIODS_TO_PLOT)
    axes[0].set_ylim(-1.1, 1.1)
    axes[0].set_xlabel("t / T")
    axes[0].set_ylabel("<sigma_x>")
    axes[0].set_title("Noninteracting dynamics for exact and detuned parameters")
    axes[0].legend(frameon=False)

    axes[1].set_xlim(0, 3)
    axes[1].set_xlabel("omega_tilde / omega")
    axes[1].set_ylabel("normalized amplitude")
    axes[1].set_title("Fourier spectra")
    axes[1].legend(frameon=False)
    fig.savefig(path, dpi=220)
    plt.close(fig)


def plot_solver_comparison(times, magnus_values, rk4_values, magnus_norms, rk4_norms, alpha, path):
    fig, axes = plt.subplots(2, 1, figsize=(8.0, 5.8), constrained_layout=True)
    axes[0].plot(times / PERIOD, magnus_values, label="Magnus4", lw=0.9)
    axes[0].plot(times / PERIOD, rk4_values, "--", label="RK4", lw=0.9)
    axes[0].set_xlim(0, PERIODS_TO_PLOT)
    axes[0].set_ylim(-1.1, 1.1)
    axes[0].set_xlabel("t / T")
    axes[0].set_ylabel("<sigma_x>")
    axes[0].set_title(f"Solver comparison, alpha = {alpha:.2f}")
    axes[0].legend(frameon=False)

    axes[1].plot(times / PERIOD, np.abs(magnus_norms - 1.0), label="Magnus4")
    axes[1].plot(times / PERIOD, np.abs(rk4_norms - 1.0), label="RK4")
    axes[1].set_yscale("log")
    axes[1].set_xlabel("t / T")
    axes[1].set_ylabel("|norm - 1|")
    axes[1].set_title("Norm drift")
    axes[1].legend(frameon=False)
    fig.savefig(path, dpi=220)
    plt.close(fig)


def plot_integer_stroboscopic(a_plus, a_minus, alpha, path):
    fig, ax = plt.subplots(figsize=(4.8, 4.8), constrained_layout=True)
    ax.scatter(a_minus, a_plus, s=14, alpha=0.75)
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-0.02, 1.02)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("|a_-|")
    ax.set_ylabel("|a_+|")
    ax.set_title(f"Stroboscopic projections, alpha = {alpha:.8f}")
    fig.savefig(path, dpi=220)
    plt.close(fig)


def plot_fractional_stroboscopic(t_over_t, phase_fraction, a_plus, a_minus, alpha, path):
    fig, ax = plt.subplots(figsize=(5.3, 4.9), constrained_layout=True)
    scatter = ax.scatter(
        a_minus,
        a_plus,
        c=phase_fraction,
        cmap="viridis",
        s=18,
        alpha=0.85,
        edgecolors="none",
    )
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-0.02, 1.02)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("|a_-|")
    ax.set_ylabel("|a_+|")
    ax.set_title(f"Fractional-period projections, alpha = {alpha:.8f}")
    cbar = fig.colorbar(scatter, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("sample phase within T")
    cbar.set_ticks([0.0, 0.25, 0.5, 0.75])
    cbar.set_ticklabels(["0", "T/4", "T/2", "3T/4"])
    fig.savefig(path, dpi=220)
    plt.close(fig)


def plot_fractional_stroboscopic_comparison(fractional_results, path):
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.2), constrained_layout=True)
    for ax, result in zip(axes, fractional_results):
        scatter = ax.scatter(
            result["a_minus"],
            result["a_plus"],
            c=result["phase_fraction"],
            cmap="viridis",
            s=16,
            alpha=0.85,
            edgecolors="none",
        )
        ax.set_xlim(-0.02, 1.02)
        ax.set_ylim(-0.02, 1.02)
        ax.set_aspect("equal", adjustable="box")
        ax.set_xlabel("|a_-|")
        ax.set_ylabel("|a_+|")
        ax.set_title(result["display_name"])
    cbar = fig.colorbar(scatter, ax=axes.ravel().tolist(), fraction=0.046, pad=0.04)
    cbar.set_label("sample phase within T")
    cbar.set_ticks([0.0, 0.25, 0.5, 0.75])
    cbar.set_ticklabels(["0", "T/4", "T/2", "3T/4"])
    fig.savefig(path, dpi=220)
    plt.close(fig)


def run_parameter_case(label, display_name, alpha, device):
    times, sigma_x_magnus, norms_magnus = evolve_time_trace(
        alpha, PERIODS_TO_PLOT, STEPS_PER_PERIOD, "magnus4", device
    )
    _, sigma_x_rk4, norms_rk4 = evolve_rk4_time_trace_with_substeps(
        alpha, PERIODS_TO_PLOT, STEPS_PER_PERIOD, RK4_SUBSTEPS, device
    )
    frequency, amplitudes = spectrum(times, sigma_x_magnus)
    save_trace_csv(
        DATA_DIR / f"sigma_x_{label}.csv",
        times,
        sigma_x_magnus,
        sigma_x_rk4,
        norms_magnus,
        norms_rk4,
    )
    save_spectrum_csv(DATA_DIR / f"spectrum_{label}.csv", frequency, amplitudes)
    plot_noninteracting_trace(
        times,
        sigma_x_magnus,
        frequency,
        amplitudes,
        alpha,
        FIGURE_DIR / f"sigma_x_{label}.png",
    )
    plot_solver_comparison(
        times,
        sigma_x_magnus,
        sigma_x_rk4,
        norms_magnus,
        norms_rk4,
        alpha,
        FIGURE_DIR / f"magnus_rk4_comparison_{label}.png",
    )
    return {
        "label": label,
        "display_name": display_name,
        "alpha": alpha,
        "times": times,
        "sigma_x_magnus": sigma_x_magnus,
        "sigma_x_rk4": sigma_x_rk4,
        "norms_magnus": norms_magnus,
        "norms_rk4": norms_rk4,
        "frequency": frequency,
        "amplitudes": amplitudes,
        "magnus4_max_norm_drift": float(np.max(np.abs(norms_magnus - 1.0))),
        "rk4_max_norm_drift": float(np.max(np.abs(norms_rk4 - 1.0))),
        "max_abs_sigma_x_difference_magnus4_vs_rk4": float(
            np.max(np.abs(sigma_x_magnus - sigma_x_rk4))
        ),
        "dominant_spectrum_omega_over_drive": float(frequency[np.argmax(amplitudes[1:]) + 1]),
    }


def run_fractional_projection_case(label, display_name, alpha, device):
    sample_index, t_over_t, phase_fraction, a_plus, a_minus = evolve_fractional_period_projections(
        alpha,
        STROBOSCOPIC_PERIODS,
        STEPS_PER_PERIOD,
        FRACTIONAL_SAMPLES_PER_PERIOD,
        device,
    )
    save_fractional_projection_csv(
        DATA_DIR / f"fractional_projections_{label}.csv",
        sample_index,
        t_over_t,
        phase_fraction,
        a_plus,
        a_minus,
    )
    plot_fractional_stroboscopic(
        t_over_t,
        phase_fraction,
        a_plus,
        a_minus,
        alpha,
        FIGURE_DIR / f"fractional_projections_{label}.png",
    )
    return {
        "label": label,
        "display_name": display_name,
        "alpha": alpha,
        "sample_index": sample_index,
        "t_over_t": t_over_t,
        "phase_fraction": phase_fraction,
        "a_plus": a_plus,
        "a_minus": a_minus,
        "unique_phase_fractions": sorted(set(np.round(phase_fraction, 8).tolist())),
    }


def main():
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    device = torch_device()

    alpha_13_refined, alpha_13_probability = refine_alpha_13()

    parameter_cases = [
        ("alpha_13", "refined alpha_13", alpha_13_refined),
        ("alpha_81p60", "detuned alpha = 81.60", ALPHA_DETUNED),
    ]
    case_results = [run_parameter_case(label, display_name, alpha, device) for label, display_name, alpha in parameter_cases]
    fractional_results = [
        run_fractional_projection_case(label, display_name, alpha, device)
        for label, display_name, alpha in parameter_cases
    ]

    _, a_plus_integer, a_minus_integer = evolve_stroboscopic_projections(
        alpha_13_refined, STROBOSCOPIC_PERIODS, STEPS_PER_PERIOD, device
    )
    np.savetxt(
        DATA_DIR / "stroboscopic_alpha_13.csv",
        np.column_stack([np.arange(len(a_plus_integer)), a_plus_integer, a_minus_integer]),
        delimiter=",",
        header="n,a_plus_abs,a_minus_abs",
        comments="",
    )
    plot_integer_stroboscopic(
        a_plus_integer,
        a_minus_integer,
        alpha_13_refined,
        FIGURE_DIR / "stroboscopic_alpha_13.png",
    )
    plot_trace_comparison(case_results, FIGURE_DIR / "sigma_x_parameter_comparison.png")
    plot_fractional_stroboscopic_comparison(
        fractional_results, FIGURE_DIR / "fractional_projections_parameter_comparison.png"
    )

    diagnostics = {
        "device": str(device),
        "cuda_device_name": torch.cuda.get_device_name(0) if torch.cuda.is_available() else None,
        "alpha_13_table": ALPHA_TABLE_13,
        "alpha_13_refined": alpha_13_refined,
        "alpha_13_one_period_return_probability": alpha_13_probability,
        "alpha_detuned": ALPHA_DETUNED,
        "periods_to_plot": PERIODS_TO_PLOT,
        "steps_per_period": STEPS_PER_PERIOD,
        "rk4_substeps_per_output_step": RK4_SUBSTEPS,
        "rk4_effective_steps_per_period": STEPS_PER_PERIOD * RK4_SUBSTEPS,
        "fractional_samples_per_period": FRACTIONAL_SAMPLES_PER_PERIOD,
        "fractional_sample_phases": ["0", "T/4", "T/2", "3T/4"],
        "cases": {
            result["label"]: {
                "alpha": result["alpha"],
                "magnus4_max_norm_drift": result["magnus4_max_norm_drift"],
                "rk4_max_norm_drift": result["rk4_max_norm_drift"],
                "max_abs_sigma_x_difference_magnus4_vs_rk4": result[
                    "max_abs_sigma_x_difference_magnus4_vs_rk4"
                ],
                "dominant_spectrum_omega_over_drive": result["dominant_spectrum_omega_over_drive"],
            }
            for result in case_results
        },
    }
    with open(DATA_DIR / "diagnostics.json", "w", encoding="utf-8") as handle:
        json.dump(diagnostics, handle, indent=2)

    print(json.dumps(diagnostics, indent=2))


if __name__ == "__main__":
    main()

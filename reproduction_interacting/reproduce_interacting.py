import argparse
import json
import time
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import torch


HBAR = 1.0
DRIVE_OMEGA = 1.0
PERIOD = 2.0 * np.pi / DRIVE_OMEGA

SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_PARAMS = SCRIPT_DIR / "params" / "interacting_params_local.json"


def parse_args():
    parser = argparse.ArgumentParser(description="Interacting spin-chain reproduction for the nonsymmorphic DTC model.")
    parser.add_argument("--params", default=str(DEFAULT_PARAMS), help="Path to a JSON parameter file.")
    parser.add_argument("--output-root", default=str(SCRIPT_DIR), help="Directory where data, figures, and diagnostics are saved.")
    parser.add_argument("--device", default="cuda", choices=["cuda", "cpu"], help="Torch device preference.")
    return parser.parse_args()


def load_params(path):
    with open(path, "r", encoding="utf-8") as handle:
        return json.load(handle)


def select_device(preference):
    if preference == "cuda" and torch.cuda.is_available():
        return torch.device("cuda")
    return torch.device("cpu")


def torch_dtypes(device):
    if device.type == "cuda":
        return torch.float32, torch.complex64
    return torch.float64, torch.complex128


def make_output_dirs(output_root):
    output_root = Path(output_root)
    data_dir = output_root / "data"
    figure_dir = output_root / "figures"
    log_dir = output_root / "logs"
    data_dir.mkdir(parents=True, exist_ok=True)
    figure_dir.mkdir(parents=True, exist_ok=True)
    log_dir.mkdir(parents=True, exist_ok=True)
    return output_root, data_dir, figure_dir, log_dir


def make_spin_ops(length, device):
    dim = 1 << length
    basis = torch.arange(dim, dtype=torch.long, device=device)
    site_masks = torch.tensor([1 << site for site in range(length)], dtype=torch.long, device=device)
    x_flip = torch.bitwise_xor(basis.unsqueeze(0), site_masks.unsqueeze(1))

    if length > 1:
        xx_masks = torch.tensor(
            [(1 << site) ^ (1 << (site + 1)) for site in range(length - 1)],
            dtype=torch.long,
            device=device,
        )
        xx_flip = torch.bitwise_xor(basis.unsqueeze(0), xx_masks.unsqueeze(1))
    else:
        xx_flip = torch.empty((0, dim), dtype=torch.long, device=device)

    bits = torch.bitwise_and(basis.unsqueeze(0), site_masks.unsqueeze(1)) != 0
    _, complex_dtype = torch_dtypes(device)
    y_phase = torch.where(
        bits,
        torch.tensor(1.0j, dtype=complex_dtype, device=device),
        torch.tensor(-1.0j, dtype=complex_dtype, device=device),
    )

    return {
        "length": length,
        "dim": dim,
        "basis": basis,
        "x_flip": x_flip,
        "xx_flip": xx_flip,
        "y_phase": y_phase,
    }


def initial_x_product_state(length, signs, device):
    dim = 1 << length
    _, complex_dtype = torch_dtypes(device)
    if signs == "all_plus":
        sign_values = [1] * length
    else:
        sign_values = [int(value) for value in signs]
        if len(sign_values) != length:
            raise ValueError("The initial-state sign list must match the chain length.")

    amplitudes = torch.ones(dim, dtype=complex_dtype, device=device)
    basis = torch.arange(dim, dtype=torch.long, device=device)
    for site, sign in enumerate(sign_values):
        if sign not in (-1, 1):
            raise ValueError("Initial x-product signs must be +1 or -1.")
        if sign == -1:
            bit_is_one = (torch.bitwise_and(basis, 1 << site) != 0)
            amplitudes = torch.where(bit_is_one, -amplitudes, amplitudes)
    return amplitudes / np.sqrt(dim)


def hamiltonian_action(state, t, alpha, interaction_j, ops):
    omega_amp = alpha / 8.0
    sx_coeff = 0.5 * omega_amp * np.sin(t)
    sy_coeff = omega_amp * np.sin(t / 2.0) ** 2

    x_terms = state[ops["x_flip"]]
    sx_sum = x_terms.sum(dim=0)
    y_phase = ops["y_phase"] if state.ndim == 1 else ops["y_phase"].unsqueeze(-1)
    sy_sum = (y_phase * x_terms).sum(dim=0)
    if ops["xx_flip"].shape[0] > 0:
        xx_sum = state[ops["xx_flip"]].sum(dim=0)
    else:
        xx_sum = torch.zeros_like(state)
    return sx_coeff * sx_sum + sy_coeff * sy_sum + interaction_j * xx_sum


def rhs(state, t, alpha, interaction_j, ops):
    return -1.0j * hamiltonian_action(state, t, alpha, interaction_j, ops) / HBAR


def normalize_state(state):
    return state / torch.linalg.norm(state)


def rk4_step(state, t, dt, alpha, interaction_j, ops, normalize_each_step):
    k1 = rhs(state, t, alpha, interaction_j, ops)
    k2 = rhs(state + 0.5 * dt * k1, t + 0.5 * dt, alpha, interaction_j, ops)
    k3 = rhs(state + 0.5 * dt * k2, t + 0.5 * dt, alpha, interaction_j, ops)
    k4 = rhs(state + dt * k3, t + dt, alpha, interaction_j, ops)
    next_state = state + dt * (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0
    if normalize_each_step:
        next_state = normalize_state(next_state)
    return next_state


def sigma_x_sites(state, ops):
    norm = torch.real(torch.vdot(state, state))
    values = (torch.conj(state).unsqueeze(0) * state[ops["x_flip"]]).sum(dim=1)
    return torch.real(values / norm)


def average_mx(state, ops):
    return torch.mean(sigma_x_sites(state, ops)).detach().cpu().item()


def dense_hamiltonian(t, alpha, interaction_j, ops):
    dim = int(ops["dim"])
    _, complex_dtype = torch_dtypes(ops["basis"].device)
    basis_vectors = torch.eye(dim, dtype=complex_dtype, device=ops["basis"].device)
    return hamiltonian_action(basis_vectors, t, alpha, interaction_j, ops)


def one_period_floquet_operator(length, alpha, interaction_j, params, device):
    ops = make_spin_ops(length, device)
    strobe_params = params["stroboscopic"]
    steps_per_period = int(strobe_params["steps_per_period"])
    dt = PERIOD / steps_per_period
    dim = int(ops["dim"])
    _, complex_dtype = torch_dtypes(device)
    floquet = torch.eye(dim, dtype=complex_dtype, device=device)

    for step in range(steps_per_period):
        midpoint_t = (step + 0.5) * dt
        h_mid = dense_hamiltonian(midpoint_t, alpha, interaction_j, ops)
        step_unitary = torch.linalg.matrix_exp((-1.0j * dt / HBAR) * h_mid)
        floquet = step_unitary @ floquet

    return floquet, ops


def stroboscopic_sample_points(max_period, strobe_params):
    linear_prefix = int(strobe_params.get("linear_prefix_periods", 200))
    log_count = int(strobe_params.get("log_sample_count", 1200))
    linear_stop = min(max_period, linear_prefix)
    points = set(range(linear_stop + 1))
    if max_period > linear_stop:
        log_start = max(1, linear_stop + 1)
        log_points = np.unique(np.geomspace(log_start, max_period, log_count).astype(np.int64))
        points.update(int(value) for value in log_points)
        points.add(int(max_period))
    return np.array(sorted(points), dtype=np.int64)


def average_mx_numpy(state, x_flip):
    norm = np.vdot(state, state).real
    site_values = np.sum(np.conj(state)[None, :] * state[x_flip], axis=1) / norm
    return float(np.real(np.mean(site_values)))


def evolve_stroboscopic_floquet(length, alpha, interaction_j, params, device, max_period):
    floquet, ops = one_period_floquet_operator(length, alpha, interaction_j, params, device)
    state0 = initial_x_product_state(length, params["initial_state"]["signs"], device)
    sample_points = stroboscopic_sample_points(max_period, params["stroboscopic"])

    floquet_np = floquet.detach().cpu().numpy().astype(np.complex128)
    state0_np = state0.detach().cpu().numpy().astype(np.complex128)
    x_flip_np = ops["x_flip"].detach().cpu().numpy()

    eigenvalues, eigenvectors = np.linalg.eig(floquet_np)
    eigenvalues = eigenvalues / np.maximum(np.abs(eigenvalues), 1.0e-300)
    coefficients = np.linalg.solve(eigenvectors, state0_np)

    rows = []
    max_norm_error = 0.0
    for period_index in sample_points:
        phase = np.exp(period_index * np.log(eigenvalues))
        state = eigenvectors @ (phase * coefficients)
        norm = float(np.linalg.norm(state))
        max_norm_error = max(max_norm_error, abs(norm - 1.0))
        state = state / norm
        mx = average_mx_numpy(state, x_flip_np)
        z_value = (1.0 if period_index % 2 == 0 else -1.0) * mx
        rows.append((period_index, period_index * PERIOD, mx, z_value, norm))

    return np.array(rows, dtype=np.float64), max_norm_error


def evolve_trace(length, alpha, interaction_j, params, device):
    ops = make_spin_ops(length, device)
    state = initial_x_product_state(length, params["initial_state"]["signs"], device)
    trace_params = params["time_trace"]
    steps_per_period = int(trace_params["steps_per_period"])
    periods = int(trace_params["periods"])
    record_stride = int(trace_params.get("record_stride", 1))
    normalize_each_step = bool(trace_params.get("normalize_each_step", True))
    dt = PERIOD / steps_per_period
    total_steps = periods * steps_per_period

    rows = []
    max_norm_error = 0.0
    for step in range(total_steps + 1):
        if step % record_stride == 0:
            norm = torch.linalg.norm(state).detach().cpu().item()
            max_norm_error = max(max_norm_error, abs(norm - 1.0))
            rows.append((step * dt / PERIOD, step * dt, average_mx(state, ops), norm))
        if step == total_steps:
            break
        state = rk4_step(state, step * dt, dt, alpha, interaction_j, ops, normalize_each_step)

    return np.array(rows), max_norm_error


def evolve_stroboscopic(length, alpha, interaction_j, params, device, max_period=None):
    ops = make_spin_ops(length, device)
    state = initial_x_product_state(length, params["initial_state"]["signs"], device)
    strobe_params = params["stroboscopic"]
    steps_per_period = int(strobe_params["steps_per_period"])
    periods = int(max_period if max_period is not None else strobe_params["periods"])
    normalize_each_step = bool(strobe_params.get("normalize_each_step", True))
    dt = PERIOD / steps_per_period

    rows = []
    max_norm_error = 0.0
    for period_index in range(periods + 1):
        mx = average_mx(state, ops)
        z_value = ((-1) ** period_index) * mx
        norm = torch.linalg.norm(state).detach().cpu().item()
        max_norm_error = max(max_norm_error, abs(norm - 1.0))
        rows.append((period_index, period_index * PERIOD, mx, z_value, norm))
        if period_index == periods:
            break
        for local_step in range(steps_per_period):
            t = (period_index * steps_per_period + local_step) * dt
            state = rk4_step(state, t, dt, alpha, interaction_j, ops, normalize_each_step)

    return np.array(rows), max_norm_error


def spectrum(times, values):
    dt = times[1] - times[0]
    centered = values - np.mean(values)
    window = np.hanning(len(centered))
    amplitudes = np.abs(np.fft.rfft(centered * window))
    amplitudes /= np.max(amplitudes) if np.max(amplitudes) > 0 else 1.0
    angular_frequency = 2.0 * np.pi * np.fft.rfftfreq(len(centered), d=dt)
    return angular_frequency / DRIVE_OMEGA, amplitudes


def save_csv(path, data, header):
    np.savetxt(path, data, delimiter=",", header=header, comments="")


def plot_trace_and_spectrum(trace, frequency, amplitudes, title, path):
    fig, axes = plt.subplots(2, 1, figsize=(8.0, 5.8), constrained_layout=True)
    axes[0].plot(trace[:, 0], trace[:, 2], lw=0.9)
    axes[0].set_xlim(0, trace[-1, 0])
    axes[0].set_ylim(-1.1, 1.1)
    axes[0].set_xlabel("t / T")
    axes[0].set_ylabel("<m_x>")
    axes[0].set_title(title)

    axes[1].plot(frequency, amplitudes, lw=0.9)
    axes[1].set_xlim(0, 3)
    axes[1].set_xlabel("omega_tilde / omega")
    axes[1].set_ylabel("normalized amplitude")
    axes[1].set_title("Fourier spectrum")
    fig.savefig(path, dpi=220)
    plt.close(fig)


def plot_trace_comparison(trace_results, path):
    fig, axes = plt.subplots(2, 1, figsize=(8.2, 5.8), constrained_layout=True)
    for result in trace_results:
        axes[0].plot(result["trace"][:, 0], result["trace"][:, 2], lw=0.9, label=result["display_name"])
        axes[1].plot(result["frequency"], result["amplitudes"], lw=0.9, label=result["display_name"])
    axes[0].set_xlim(0, max(item["trace"][-1, 0] for item in trace_results))
    axes[0].set_ylim(-1.1, 1.1)
    axes[0].set_xlabel("t / T")
    axes[0].set_ylabel("<m_x>")
    axes[0].set_title("Interacting dynamics")
    axes[0].legend(frameon=False)
    axes[1].set_xlim(0, 3)
    axes[1].set_xlabel("omega_tilde / omega")
    axes[1].set_ylabel("normalized amplitude")
    axes[1].set_title("Fourier spectra")
    axes[1].legend(frameon=False)
    fig.savefig(path, dpi=220)
    plt.close(fig)


def plot_z_by_length(strobe_results, title, path, use_log_x=False):
    fig, ax = plt.subplots(figsize=(7.2, 4.8), constrained_layout=True)
    for result in strobe_results:
        rows = result["rows"]
        if use_log_x:
            positive = rows[:, 0] > 0
            ax.plot(rows[positive, 0], rows[positive, 3], lw=1.0, label=f"L={result['length']}")
        else:
            ax.plot(rows[:, 0], rows[:, 3], lw=1.0, label=f"L={result['length']}")
    if use_log_x:
        ax.set_xscale("log")
    ax.set_xlabel("n")
    ax.set_ylabel("Z(n) = (-1)^n m_x(nT)")
    ax.set_ylim(-0.1, 1.05)
    ax.set_title(title)
    ax.legend(frameon=False, ncol=2)
    fig.savefig(path, dpi=220)
    plt.close(fig)


def safe_label(label):
    return label.replace(".", "p").replace("=", "").replace(" ", "_")


def run(params, output_root, device):
    output_root, data_dir, figure_dir, log_dir = make_output_dirs(output_root)
    interaction_j = float(params["interaction_j"])
    if params.get("interaction_axis", "x") != "x":
        raise ValueError("This reproduction script currently implements the paper's mu=x interaction.")
    if params.get("boundary", "open") != "open":
        raise ValueError("This reproduction script currently uses open boundary conditions.")

    diagnostics = {
        "timestamp_unix": time.time(),
        "device": str(device),
        "cuda_device_name": torch.cuda.get_device_name(0) if device.type == "cuda" else None,
        "interaction_j": interaction_j,
        "boundary": params.get("boundary", "open"),
        "trace": {},
        "stroboscopic": {},
    }

    trace_results = []
    trace_length = int(params["time_trace"]["length"])
    for alpha_label in params["time_trace"]["alpha_labels"]:
        alpha = float(params["alpha_values"][alpha_label])
        label = safe_label(alpha_label)
        trace, max_norm_error = evolve_trace(trace_length, alpha, interaction_j, params, device)
        frequency, amplitudes = spectrum(trace[:, 1], trace[:, 2])
        dominant_frequency = float(frequency[np.argmax(amplitudes[1:]) + 1])

        save_csv(data_dir / f"mx_trace_{label}_l{trace_length}.csv", trace, "t_over_T,t,mx,norm")
        save_csv(
            data_dir / f"mx_spectrum_{label}_l{trace_length}.csv",
            np.column_stack([frequency, amplitudes]),
            "omega_tilde_over_omega,normalized_amplitude",
        )
        plot_trace_and_spectrum(
            trace,
            frequency,
            amplitudes,
            f"Interacting mx(t), {alpha_label}, L={trace_length}, J={interaction_j}",
            figure_dir / f"mx_trace_spectrum_{label}_l{trace_length}.png",
        )

        trace_results.append(
            {
                "alpha_label": alpha_label,
                "display_name": f"{alpha_label}, L={trace_length}",
                "alpha": alpha,
                "trace": trace,
                "frequency": frequency,
                "amplitudes": amplitudes,
            }
        )
        diagnostics["trace"][alpha_label] = {
            "alpha": alpha,
            "length": trace_length,
            "max_norm_error": max_norm_error,
            "dominant_spectrum_omega_over_drive": dominant_frequency,
            "mx_min": float(np.min(trace[:, 2])),
            "mx_max": float(np.max(trace[:, 2])),
        }

    plot_trace_comparison(trace_results, figure_dir / f"mx_trace_comparison_l{trace_length}.png")

    for alpha_label in params["stroboscopic"]["alpha_labels"]:
        alpha = float(params["alpha_values"][alpha_label])
        label = safe_label(alpha_label)
        strobe_results = []
        strobe_params = params["stroboscopic"]
        periods_by_alpha = strobe_params.get("periods_by_alpha", {})
        max_period = int(periods_by_alpha.get(alpha_label, strobe_params["periods"]))
        sampling_mode = strobe_params.get("sampling_mode", "all_periods")
        diagnostics["stroboscopic"][alpha_label] = {
            "alpha": alpha,
            "max_period": max_period,
            "sampling_mode": sampling_mode,
            "lengths": {},
        }
        for length in params["stroboscopic"]["lengths"]:
            length = int(length)
            if sampling_mode == "log_floquet":
                rows, max_norm_error = evolve_stroboscopic_floquet(length, alpha, interaction_j, params, device, max_period)
            else:
                rows, max_norm_error = evolve_stroboscopic(length, alpha, interaction_j, params, device, max_period)
            save_csv(data_dir / f"z_stroboscopic_{label}_l{length}.csv", rows, "n,t,mx,z,norm")
            strobe_results.append({"length": length, "rows": rows})
            diagnostics["stroboscopic"][alpha_label]["lengths"][str(length)] = {
                "max_norm_error": max_norm_error,
                "sample_count": int(len(rows)),
                "z_min": float(np.min(rows[:, 3])),
                "z_mean_first_20": float(np.mean(rows[: min(21, len(rows)), 3])),
                "z_final": float(rows[-1, 3]),
            }
        plot_z_by_length(
            strobe_results,
            f"Stroboscopic Z(n), {alpha_label}, J={interaction_j}",
            figure_dir / f"z_stroboscopic_{label}.png",
            use_log_x=(sampling_mode == "log_floquet"),
        )

    with open(data_dir / "diagnostics.json", "w", encoding="utf-8") as handle:
        json.dump(diagnostics, handle, indent=2)

    with open(log_dir / "run_summary.json", "w", encoding="utf-8") as handle:
        json.dump(diagnostics, handle, indent=2)
    return diagnostics


def main():
    args = parse_args()
    params = load_params(args.params)
    device = select_device(args.device)
    diagnostics = run(params, args.output_root, device)
    print(json.dumps(diagnostics, indent=2))


if __name__ == "__main__":
    main()

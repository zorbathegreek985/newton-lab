"""Phase 48: pinned path relaxation and continuum heat-mode comparison.

The agent states are dimensionless estimates, not temperatures.  A uniform
path with two fixed zero-valued leaders obeys the same second-difference
operator as a finite-difference heat equation on a unit interval.
"""

from __future__ import annotations

import csv
import json
import math
import platform
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scipy  # type: ignore[import-untyped]
from scipy.integrate import solve_ivp  # type: ignore[import-untyped]

from newton_lab.exceptions import ScientificValidationError
from newton_lab.modal_research import heat_mode_decay_rate

GRID_SIZES = (7, 15, 31, 63, 127)
MODE_NUMBERS = (1, 2, 3)
OUTPUT_POINTS = (501, 2001)
SOLVER_SETTINGS = ((1e-9, 1e-11), (1e-12, 1e-14))
SCALED_TIME_END = 5.0


@dataclass(frozen=True)
class ModalRun:
    """One deterministic numerical integration and its exact references."""

    agent_count: int
    mode_number: int
    h: float
    output_points: int
    relative_tolerance: float
    absolute_tolerance: float
    discrete_rate_per_tau: float
    continuum_rate_per_tau: float
    numerical_rate_per_tau: float
    numerical_rate_relative_error: float
    discrete_rate_relative_error: float
    maximum_absolute_state_error: float
    maximum_state_error_relative_to_initial_peak: float
    solver_success: bool
    solver_message: str


def _positive_integer(value: int, name: str) -> None:
    if isinstance(value, bool) or not isinstance(value, int) or value < 1:
        raise ScientificValidationError(f"{name} must be a positive integer")


def validate_mode(agent_count: int, mode_number: int) -> None:
    """Validate an interior-agent count and a resolvable sine mode."""
    _positive_integer(agent_count, "agent_count")
    _positive_integer(mode_number, "mode_number")
    if mode_number > agent_count:
        raise ScientificValidationError("mode_number must be between 1 and agent_count")


def path_laplacian_operator(agent_count: int) -> np.ndarray:
    """Return the positive Dirichlet path Laplacian divided by ``h**2``.

    Boundary leaders are fixed at zero and are not part of this matrix.
    Consequently ``dz/dtau = -A @ z``.
    """
    _positive_integer(agent_count, "agent_count")
    h = 1.0 / (agent_count + 1)
    operator = np.diag(np.full(agent_count, 2.0 / h**2, dtype=np.float64))
    if agent_count > 1:
        off_diagonal = np.full(agent_count - 1, -1.0 / h**2, dtype=np.float64)
        operator += np.diag(off_diagonal, 1) + np.diag(off_diagonal, -1)
    return operator


def initial_sine_mode(agent_count: int, mode_number: int) -> np.ndarray:
    """Return samples ``sin(n*pi*i*h)`` at the interior agent locations."""
    validate_mode(agent_count, mode_number)
    h = 1.0 / (agent_count + 1)
    indices = np.arange(1, agent_count + 1, dtype=np.float64)
    return np.sin(mode_number * math.pi * indices * h)


def discrete_mode_decay_rate(agent_count: int, mode_number: int) -> float:
    """Return the exact eigenvalue ``4 h^-2 sin²(n*pi*h/2)``."""
    validate_mode(agent_count, mode_number)
    h = 1.0 / (agent_count + 1)
    return 4.0 / h**2 * math.sin(mode_number * math.pi * h / 2.0) ** 2


def continuum_mode_decay_rate(mode_number: int) -> float:
    """Return the unit-interval, unit-diffusivity heat rate ``(n*pi)**2``."""
    _positive_integer(mode_number, "mode_number")
    return heat_mode_decay_rate(1.0, mode_number, 1.0)


def modal_projection(states: np.ndarray, eigenvector: np.ndarray) -> np.ndarray:
    """Project one state vector or a row-wise state matrix onto a mode."""
    values = np.asarray(states, dtype=np.float64)
    vector = np.asarray(eigenvector, dtype=np.float64)
    if vector.ndim != 1 or len(vector) == 0 or not np.isfinite(vector).all():
        raise ScientificValidationError("eigenvector must be a nonempty finite vector")
    if values.ndim not in (1, 2) or values.shape[-1] != len(vector):
        raise ScientificValidationError("states must end in the eigenvector dimension")
    if not np.isfinite(values).all():
        raise ScientificValidationError("states must contain only finite values")
    denominator = float(vector @ vector)
    if denominator == 0.0:
        raise ScientificValidationError("eigenvector must have nonzero norm")
    return values @ vector / denominator


def estimate_modal_decay_rate(times: np.ndarray, amplitudes: np.ndarray) -> float:
    """Estimate positive decay rate by least-squares fitting log amplitude."""
    t = np.asarray(times, dtype=np.float64)
    a = np.asarray(amplitudes, dtype=np.float64)
    if (
        t.ndim != 1
        or a.ndim != 1
        or len(t) != len(a)
        or len(t) < 3
        or not np.isfinite(t).all()
        or not np.isfinite(a).all()
        or np.any(np.diff(t) <= 0.0)
        or np.any(a == 0.0)
    ):
        raise ScientificValidationError("rate fit requires finite, ordered samples")
    slope = float(np.polyfit(t, np.log(np.abs(a / a[0])), 1)[0])
    if not math.isfinite(slope) or slope >= 0.0:
        raise ScientificValidationError("projected mode did not show positive decay")
    return -slope


def simulate_pinned_path_mode(
    agent_count: int,
    mode_number: int,
    *,
    output_points: int = 501,
    relative_tolerance: float = 1e-9,
    absolute_tolerance: float = 1e-11,
) -> tuple[ModalRun, np.ndarray, np.ndarray, np.ndarray]:
    """Integrate one pinned-path eigenmode with DOP853.

    Samples span dimensionless continuum time ``s=lambda_cont*tau`` from 0 to
    5. The returned arrays are times, numerical states, and exact discrete
    states; the result reports independent numerical and spatial errors.
    """
    validate_mode(agent_count, mode_number)
    _positive_integer(output_points, "output_points")
    if output_points < 2:
        raise ScientificValidationError("output_points must be at least two")
    tolerances = (relative_tolerance, absolute_tolerance)
    if not np.isfinite(tolerances).all() or min(tolerances) <= 0.0:
        raise ScientificValidationError("solver tolerances must be finite and positive")
    continuum_rate = continuum_mode_decay_rate(mode_number)
    discrete_rate = discrete_mode_decay_rate(agent_count, mode_number)
    initial = initial_sine_mode(agent_count, mode_number)
    times = np.linspace(0.0, SCALED_TIME_END / continuum_rate, output_points)
    operator = path_laplacian_operator(agent_count)
    solution = solve_ivp(
        lambda _time, state: -(operator @ state),
        (float(times[0]), float(times[-1])),
        initial,
        method="DOP853",
        t_eval=times,
        rtol=relative_tolerance,
        atol=absolute_tolerance,
    )
    if (
        not solution.success
        or not isinstance(solution.y, np.ndarray)
        or solution.y.shape != (agent_count, output_points)
    ):
        raise RuntimeError(f"pinned-path integration failed: {solution.message}")
    numerical_states = solution.y.T.copy()
    exact_states = np.exp(-discrete_rate * times[:, None]) * initial[None, :]
    eigenvector = initial_sine_mode(agent_count, mode_number)
    projected = modal_projection(numerical_states, eigenvector)
    numerical_rate = estimate_modal_decay_rate(times, projected)
    state_error = np.abs(numerical_states - exact_states)
    run = ModalRun(
        agent_count=agent_count,
        mode_number=mode_number,
        h=1.0 / (agent_count + 1),
        output_points=output_points,
        relative_tolerance=relative_tolerance,
        absolute_tolerance=absolute_tolerance,
        discrete_rate_per_tau=discrete_rate,
        continuum_rate_per_tau=continuum_rate,
        numerical_rate_per_tau=numerical_rate,
        numerical_rate_relative_error=abs(numerical_rate - discrete_rate)
        / discrete_rate,
        discrete_rate_relative_error=abs(discrete_rate - continuum_rate)
        / continuum_rate,
        maximum_absolute_state_error=float(np.max(state_error)),
        maximum_state_error_relative_to_initial_peak=float(
            np.max(state_error) / np.max(np.abs(initial))
        ),
        solver_success=bool(solution.success),
        solver_message=str(solution.message),
    )
    return run, times, numerical_states, exact_states


def _fitted_order(h: np.ndarray, error: np.ndarray) -> float:
    slope = np.polyfit(np.log(h), np.log(error), 1)[0]
    return float(slope)


def _write_csv(path: Path, records: list[dict[str, Any]]) -> None:
    if not records:
        raise ValueError("cannot write an empty study table")
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)


def run_pinned_path_study(output_directory: Path) -> dict[str, Any]:
    """Run the frozen Phase 48 grid and write tables, metadata, and figures."""
    output_directory = Path(output_directory)
    output_directory.mkdir(parents=True, exist_ok=True)
    per_condition: list[dict[str, Any]] = []
    runs: dict[
        tuple[int, int, int, float], tuple[ModalRun, np.ndarray, np.ndarray, np.ndarray]
    ] = {}
    for count in GRID_SIZES:
        for mode in MODE_NUMBERS:
            for points in OUTPUT_POINTS:
                for relative_tolerance, absolute_tolerance in SOLVER_SETTINGS:
                    run_data = simulate_pinned_path_mode(
                        count,
                        mode,
                        output_points=points,
                        relative_tolerance=relative_tolerance,
                        absolute_tolerance=absolute_tolerance,
                    )
                    run, times, numerical, exact = run_data
                    runs[(count, mode, points, relative_tolerance)] = run_data
                    row = asdict(run)
                    row["scaled_time_start"] = 0.0
                    row["scaled_time_end"] = SCALED_TIME_END
                    row["scaled_time_points"] = points
                    row["maximum_absolute_state_error"] = (
                        run.maximum_absolute_state_error
                    )
                    per_condition.append(row)

    summaries: list[dict[str, Any]] = []
    for mode in MODE_NUMBERS:
        selected = [
            min(
                (
                    runs[(count, mode, points, rtol)][0]
                    for points in OUTPUT_POINTS
                    for rtol, _atol in SOLVER_SETTINGS
                ),
                key=lambda result: result.numerical_rate_relative_error,
            )
            for count in GRID_SIZES
        ]
        h_values = np.array([result.h for result in selected[-4:]])
        errors = np.array(
            [result.discrete_rate_relative_error for result in selected[-4:]]
        )
        ratios = np.array(
            [
                result.discrete_rate_per_tau / result.continuum_rate_per_tau
                for result in selected
            ]
        )
        refinement_errors = np.array(
            [result.discrete_rate_relative_error for result in selected]
        )
        order = _fitted_order(h_values, errors)
        decreasing = bool(np.all(np.diff(refinement_errors) < 0.0))
        below = bool(np.all(ratios < 1.0))
        order_valid = bool(1.7 <= order <= 2.3)
        summaries.append(
            {
                "mode_number": mode,
                "continuum_rate_per_tau": continuum_mode_decay_rate(mode),
                "grid_sizes": ";".join(str(value) for value in GRID_SIZES),
                "discrete_rate_ratios": ";".join(f"{value:.12g}" for value in ratios),
                "relative_spatial_errors": ";".join(
                    f"{value:.12g}" for value in refinement_errors
                ),
                "fitted_order_finest_four": order,
                "all_ratios_below_one": below,
                "error_decreases_each_refinement": decreasing,
                "order_within_1_7_to_2_3": order_valid,
                "prediction_supported_for_mode": below and decreasing and order_valid,
            }
        )

    figures_directory = output_directory
    _plot_rate_convergence(figures_directory / "modal_rate_convergence.png")
    _plot_trajectory_comparison(
        runs, figures_directory / "modal_trajectory_comparison.png"
    )
    _write_csv(output_directory / "per_condition_results.csv", per_condition)
    _write_csv(output_directory / "convergence_summary.csv", summaries)
    metadata = {
        "study": "Phase 48 pinned-path modal decay and convergence validation",
        "model": "dz_i/dtau = (z_(i-1)-2*z_i+z_(i+1))/h^2, z_0=z_(N+1)=0",
        "grid_sizes": list(GRID_SIZES),
        "mode_numbers": list(MODE_NUMBERS),
        "output_points": list(OUTPUT_POINTS),
        "solver": "scipy.integrate.solve_ivp(method='DOP853')",
        "solver_settings": [{"rtol": r, "atol": a} for r, a in SOLVER_SETTINGS],
        "scaled_time_interval": [0.0, SCALED_TIME_END],
        "spatial_rate_formula": "4*h^-2*sin(n*pi*h/2)^2",
        "continuum_rate_formula": "(n*pi)^2; alpha=1, L=1",
        "continuum_reference_source": (
            "newton_lab.modal_research.heat_mode_decay_rate(1, n, 1)"
        ),
        "fitted_order_rule": (
            "OLS slope log(relative spatial error) vs log(h), finest four N"
        ),
        "acceptance_rule": (
            "all ratios < 1, strictly decreasing errors, and fitted order "
            "in [1.7, 2.3], for every mode"
        ),
        "environment": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "matplotlib": matplotlib.__version__,
        },
        "row_count_per_condition_results": len(per_condition),
        "row_count_convergence_summary": len(summaries),
        "prediction_supported_all_modes": all(
            row["prediction_supported_for_mode"] for row in summaries
        ),
    }
    (output_directory / "metadata.json").write_text(
        json.dumps(metadata, indent=2) + "\n", encoding="utf-8"
    )
    return metadata


def _plot_rate_convergence(path: Path) -> None:
    fig, ax = plt.subplots(figsize=(8, 5))
    h_values = np.array([1.0 / (count + 1) for count in GRID_SIZES])
    for mode in MODE_NUMBERS:
        error = np.array(
            [
                1.0
                - discrete_mode_decay_rate(count, mode)
                / continuum_mode_decay_rate(mode)
                for count in GRID_SIZES
            ]
        )
        ax.loglog(h_values, error, marker="o", label=f"mode n={mode}")
    reference = (math.pi * h_values) ** 2 / 12.0
    ax.loglog(h_values, reference, "k--", label="O(h²) reference, n=1")
    ax.set_xlabel("Grid spacing h")
    ax.set_ylabel("Relative rate error 1 - λ(N,n)/λ(n)")
    ax.set_title("Pinned-path modal rate spatial convergence")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def _plot_trajectory_comparison(
    runs: dict[
        tuple[int, int, int, float], tuple[ModalRun, np.ndarray, np.ndarray, np.ndarray]
    ],
    path: Path,
) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), sharey=True)
    for ax, count in zip(axes, (7, 31), strict=True):
        mode = 1
        continuum = continuum_mode_decay_rate(mode)
        for points, rtol, label in (
            (501, 1e-9, "numerical, 501 points"),
            (2001, 1e-12, "numerical, 2001 points"),
        ):
            _run, times, numerical, _exact = runs[(count, mode, points, rtol)]
            ax.plot(continuum * times, numerical[:, count // 2], label=label)
        _run, times, _numerical, exact = runs[(count, mode, 2001, 1e-12)]
        ax.plot(
            continuum * times,
            exact[:, count // 2],
            "k--",
            label="exact discrete trajectory",
        )
        ax.set_title(f"N={count}, mode n=1")
        ax.set_xlabel("Scaled time s = λcontinuum τ")
        ax.grid(True, alpha=0.25)
        ax.legend(fontsize=8)
    axes[0].set_ylabel("Midpoint agent state")
    fig.suptitle("Numerical pinned-path trajectory vs exact discrete mode")
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)

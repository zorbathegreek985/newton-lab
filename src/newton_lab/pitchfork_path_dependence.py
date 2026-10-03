"""Bounded Phase 62 parameter-path study for the supercritical pitchfork."""

from __future__ import annotations

import csv
import hashlib
import json
import platform
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scipy  # type: ignore[import-untyped]
from numpy.typing import NDArray
from scipy.integrate import solve_ivp  # type: ignore[import-untyped]

from newton_lab.exceptions import IntegrationError, ScientificValidationError
from newton_lab.simulation import ODESolverConfiguration

PROTOCOL_FILENAME = "protocol.json"
FROZEN_PROTOCOL_SHA256 = (
    "A2D35E4AD3B756AA6071AC3C7B4F0141ACEC05BC9DDE9BFE09981997A66351E6"
)
MODEL_ID = "supercritical_pitchfork"


@dataclass(frozen=True, slots=True)
class RampLeg:
    """One continuous linear control ramp and its scalar state trajectory."""

    direction: str
    mu: NDArray[np.float64]
    time: NDArray[np.float64]
    state: NDArray[np.float64]
    function_evaluations: int


@dataclass(frozen=True, slots=True)
class PitchforkCycle:
    """Upward and downward ramps with state carried through the turning point."""

    leg_duration: float
    initial_state: float
    upward: RampLeg
    downward: RampLeg


def pitchfork_equilibrium_magnitude(
    mu: NDArray[np.float64] | float,
) -> NDArray[np.float64] | float:
    """Return the magnitude of the stable equilibrium branch of dx=mu*x-x^3."""
    values = np.asarray(mu, dtype=np.float64)
    if not np.isfinite(values).all():
        raise ScientificValidationError("mu must contain only finite values")
    result = np.sqrt(np.maximum(values, 0.0))
    if np.ndim(mu) == 0:
        return float(result)
    return np.asarray(result, dtype=np.float64)


def _simulate_leg(
    *,
    mu_start: float,
    mu_end: float,
    initial_state: float,
    duration: float,
    output_interval: float,
    solver: ODESolverConfiguration,
    direction: str,
) -> RampLeg:
    if (
        not np.isfinite(
            (mu_start, mu_end, initial_state, duration, output_interval)
        ).all()
        or duration <= 0.0
        or output_interval <= 0.0
        or direction not in {"up", "down"}
    ):
        raise ScientificValidationError("invalid pitchfork ramp inputs")
    point_count = int(round(duration / output_interval)) + 1
    if point_count < 2:
        raise ScientificValidationError("a ramp requires at least two output points")
    time = np.linspace(0.0, duration, point_count, dtype=np.float64)
    mu = np.linspace(mu_start, mu_end, point_count, dtype=np.float64)
    ramp_rate = (mu_end - mu_start) / duration

    def derivative(
        time_value: float, state: NDArray[np.float64]
    ) -> NDArray[np.float64]:
        mu_value = mu_start + ramp_rate * time_value
        value = float(state[0])
        return np.asarray([mu_value * value - value**3], dtype=np.float64)

    solution = solve_ivp(
        derivative,
        (0.0, duration),
        np.asarray([initial_state], dtype=np.float64),
        method=solver.method,
        t_eval=time,
        rtol=solver.relative_tolerance,
        atol=solver.absolute_tolerance,
        max_step=np.inf if solver.maximum_step is None else solver.maximum_step,
    )
    if not solution.success:
        raise IntegrationError(f"pitchfork {direction} ramp failed: {solution.message}")
    if solution.y.shape != (1, point_count) or not np.isfinite(solution.y).all():
        raise IntegrationError(f"pitchfork {direction} ramp returned invalid output")
    return RampLeg(
        direction=direction,
        mu=mu,
        time=np.asarray(solution.t, dtype=np.float64).copy(),
        state=np.asarray(solution.y[0], dtype=np.float64).copy(),
        function_evaluations=int(solution.nfev),
    )


def simulate_pitchfork_cycle(
    *,
    mu_min: float,
    mu_max: float,
    initial_state: float,
    leg_duration: float,
    output_interval: float,
    solver: ODESolverConfiguration,
) -> PitchforkCycle:
    """Integrate one up/down pitchfork sweep, carrying state between legs.

    This model-specific continuation helper uses the project's explicit SciPy
    ODE settings. It does not introduce a general solver or infer equilibrium
    from solver success.
    """
    if not mu_min < mu_max:
        raise ScientificValidationError("mu_min must be less than mu_max")
    upward = _simulate_leg(
        mu_start=mu_min,
        mu_end=mu_max,
        initial_state=initial_state,
        duration=leg_duration,
        output_interval=output_interval,
        solver=solver,
        direction="up",
    )
    downward = _simulate_leg(
        mu_start=mu_max,
        mu_end=mu_min,
        initial_state=float(upward.state[-1]),
        duration=leg_duration,
        output_interval=output_interval,
        solver=solver,
        direction="down",
    )
    return PitchforkCycle(
        leg_duration=leg_duration,
        initial_state=initial_state,
        upward=upward,
        downward=downward,
    )


def _loop_area(cycle: PitchforkCycle) -> float:
    down_magnitude_increasing_mu = np.abs(cycle.downward.state[::-1])
    difference = np.abs(np.abs(cycle.upward.state) - down_magnitude_increasing_mu)
    return float(np.trapezoid(difference, cycle.upward.mu))


def _equilibrium_rmse(leg: RampLeg) -> float:
    equilibrium = np.asarray(pitchfork_equilibrium_magnitude(leg.mu))
    return float(np.sqrt(np.mean((np.abs(leg.state) - equilibrium) ** 2)))


def _cycle_metrics(cycle: PitchforkCycle, reference_area: float) -> dict[str, Any]:
    return {
        "leg_duration": cycle.leg_duration,
        "initial_state": cycle.initial_state,
        "initial_sign": int(np.sign(cycle.initial_state)),
        "loop_area": _loop_area(cycle),
        "normalized_loop_area": _loop_area(cycle) / reference_area,
        "up_equilibrium_rmse": _equilibrium_rmse(cycle.upward),
        "down_equilibrium_rmse": _equilibrium_rmse(cycle.downward),
        "upper_turning_point_state": float(cycle.upward.state[-1]),
        "lower_turning_point_state": float(cycle.downward.state[-1]),
        "function_evaluations": (
            cycle.upward.function_evaluations + cycle.downward.function_evaluations
        ),
        "maximum_absolute_state": float(
            max(
                np.max(np.abs(cycle.upward.state)),
                np.max(np.abs(cycle.downward.state)),
            )
        ),
    }


def _cycle_rows(cycle: PitchforkCycle) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for leg in (cycle.upward, cycle.downward):
        equilibrium = np.asarray(pitchfork_equilibrium_magnitude(leg.mu))
        for index, (time_value, mu_value, state_value, eq_value) in enumerate(
            zip(leg.time, leg.mu, leg.state, equilibrium, strict=True)
        ):
            rows.append(
                {
                    "leg_duration": cycle.leg_duration,
                    "initial_state": cycle.initial_state,
                    "direction": leg.direction,
                    "sample_index": index,
                    "time_in_leg": float(time_value),
                    "mu": float(mu_value),
                    "state": float(state_value),
                    "absolute_state": abs(float(state_value)),
                    "stable_equilibrium_magnitude": float(eq_value),
                }
            )
    return rows


def _interpolated_state_difference(
    base: PitchforkCycle, refined: PitchforkCycle
) -> float:
    differences: list[float] = []
    for base_leg, refined_leg in (
        (base.upward, refined.upward),
        (base.downward, refined.downward),
    ):
        if base_leg.direction == "down":
            comparison_mu = base_leg.mu[::-1]
            comparison_state = base_leg.state[::-1]
            target_mu = refined_leg.mu[::-1]
            target_state = refined_leg.state[::-1]
        else:
            comparison_mu = base_leg.mu
            comparison_state = base_leg.state
            target_mu = refined_leg.mu
            target_state = refined_leg.state
        interpolated = np.interp(target_mu, comparison_mu, comparison_state)
        differences.append(float(np.max(np.abs(interpolated - target_state))))
    return max(differences)


def _write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise ScientificValidationError(f"refusing to write empty result file: {path}")
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def _save_figures(
    output: Path,
    cycles: list[PitchforkCycle],
    metric_rows: list[dict[str, Any]],
    numerical_rows: list[dict[str, Any]],
    reference_area: float,
) -> None:
    figure, axes = plt.subplots(1, 2, figsize=(12, 4.8), sharey=True)
    for axis, duration in zip(axes, (20.0, 80.0), strict=True):
        selected_cycles = [cycle for cycle in cycles if cycle.leg_duration == duration]
        for cycle in selected_cycles:
            if cycle.initial_state == 0.0:
                continue
            label = f"x0={cycle.initial_state:+g}"
            axis.plot(cycle.upward.mu, np.abs(cycle.upward.state), label=f"up {label}")
            axis.plot(
                cycle.downward.mu,
                np.abs(cycle.downward.state),
                linestyle="--",
                label=f"down {label}",
            )
        mu_reference = np.linspace(-0.16, 0.16, 300)
        axis.plot(
            mu_reference,
            pitchfork_equilibrium_magnitude(mu_reference),
            color="black",
            linewidth=1.5,
            label="stable equilibrium |x*|",
        )
        axis.axvline(0.0, color="gray", linestyle=":")
        axis.set_title(f"Ramp leg duration {duration:g}")
        axis.set_xlabel("Control parameter mu")
        axis.grid(alpha=0.25)
        axis.legend(fontsize=7)
    axes[0].set_ylabel("Absolute state |x|")
    figure.suptitle("Finite-rate pitchfork path response (dimensionless)")
    figure.tight_layout()
    figure.savefig(
        output / "pitchfork_path_branches.png",
        dpi=150,
        metadata={"Software": "Newton Lab Phase 62"},
    )
    plt.close(figure)

    duration_rows = [
        row for row in metric_rows if row["initial_state"] in (0.01, -0.01)
    ]
    figure, axis = plt.subplots(figsize=(7.5, 4.5))
    for sign, marker, label in ((1, "o", "+epsilon"), (-1, "s", "-epsilon")):
        selected = sorted(
            (row for row in duration_rows if row["initial_sign"] == sign),
            key=lambda row: row["leg_duration"],
        )
        axis.plot(
            [row["leg_duration"] for row in selected],
            [row["normalized_loop_area"] for row in selected],
            marker=marker,
            label=label,
        )
    for row in numerical_rows:
        axis.errorbar(
            row["leg_duration"],
            row["primary_normalized_loop_area"],
            yerr=row["loop_area_absolute_difference"] / reference_area,
            fmt="none",
            color="gray",
            capsize=3,
        )
    axis.set_xlabel("Ramp leg duration (dimensionless time)")
    axis.set_ylabel("Loop area / equilibrium response-area reference")
    axis.set_title("Path-loop metric and numerical sensitivity")
    axis.grid(alpha=0.25)
    axis.legend()
    figure.tight_layout()
    figure.savefig(
        output / "loop_area_by_ramp_duration.png",
        dpi=150,
        metadata={"Software": "Newton Lab Phase 62"},
    )
    plt.close(figure)


def _classify_hypotheses(
    metric_rows: list[dict[str, Any]],
    numerical_rows: list[dict[str, Any]],
    cycles: list[PitchforkCycle],
    reference_area: float,
    protocol: dict[str, Any],
) -> dict[str, Any]:
    relative_acceptance = protocol["metrics"][
        "solver_acceptance_relative_to_reference_area"
    ]
    resolution_by_duration: dict[float, list[dict[str, Any]]] = {}
    for row in numerical_rows:
        resolution_by_duration.setdefault(row["leg_duration"], []).append(row)
    resolution_errors = {
        duration: max(row["loop_area_absolute_difference"] for row in rows)
        for duration, rows in resolution_by_duration.items()
    }
    resolution_pass = all(
        row["loop_area_absolute_difference"] <= relative_acceptance * reference_area
        for row in numerical_rows
    )
    numeric_floor = 1e-12 * reference_area
    threshold_by_duration = {
        duration: 10.0 * max(error, numeric_floor)
        for duration, error in resolution_errors.items()
    }

    nonzero_rows = [row for row in metric_rows if row["initial_state"] != 0.0]
    detected = [
        row["loop_area"] > threshold_by_duration[row["leg_duration"]]
        for row in nonzero_rows
    ]
    per_duration_signs = {
        duration: all(
            row["loop_area"] > threshold_by_duration[duration]
            for row in nonzero_rows
            if row["leg_duration"] == duration
        )
        for duration in threshold_by_duration
    }
    if not resolution_pass:
        h1 = "Inconclusive"
    elif any(per_duration_signs.values()):
        h1 = "Supported within the tested conditions"
    elif not any(detected):
        h1 = "Not supported within the tested conditions"
    else:
        h1 = "Inconclusive"

    fast_duration, slow_duration = sorted(threshold_by_duration)
    fast_area = float(
        np.mean(
            [
                row["loop_area"]
                for row in nonzero_rows
                if row["leg_duration"] == fast_duration
            ]
        )
    )
    slow_area = float(
        np.mean(
            [
                row["loop_area"]
                for row in nonzero_rows
                if row["leg_duration"] == slow_duration
            ]
        )
    )
    rate_difference = fast_area - slow_area
    rate_threshold = 10.0 * max(
        resolution_errors[fast_duration],
        resolution_errors[slow_duration],
        numeric_floor,
    )
    if not resolution_pass:
        h2 = "Inconclusive"
    elif rate_difference > rate_threshold:
        h2 = "Supported within the tested conditions"
    elif abs(rate_difference) <= rate_threshold:
        h2 = "Not supported within the tested conditions"
    else:
        h2 = "Inconclusive"

    signed_pairs = {
        duration: (
            next(
                row
                for row in nonzero_rows
                if row["leg_duration"] == duration and row["initial_state"] > 0.0
            ),
            next(
                row
                for row in nonzero_rows
                if row["leg_duration"] == duration and row["initial_state"] < 0.0
            ),
        )
        for duration in threshold_by_duration
    }
    branch_signs_match = all(
        np.sign(positive["upper_turning_point_state"]) > 0.0
        and np.sign(negative["upper_turning_point_state"]) < 0.0
        and np.sign(positive["lower_turning_point_state"])
        == np.sign(positive["upper_turning_point_state"])
        and np.sign(negative["lower_turning_point_state"])
        == np.sign(negative["upper_turning_point_state"])
        for positive, negative in signed_pairs.values()
    )
    state_error_by_duration = {
        duration: max(
            row["maximum_state_difference"]
            for row in numerical_rows
            if row["leg_duration"] == duration
            and row["comparison"] == "tighter_solver_settings"
        )
        for duration in threshold_by_duration
    }
    symmetry_difference_by_duration: dict[float, float] = {}
    for duration in threshold_by_duration:
        positive_cycle = next(
            cycle
            for cycle in cycles
            if cycle.leg_duration == duration and cycle.initial_state > 0.0
        )
        negative_cycle = next(
            cycle
            for cycle in cycles
            if cycle.leg_duration == duration and cycle.initial_state < 0.0
        )
        symmetry_difference_by_duration[duration] = max(
            float(
                np.max(np.abs(np.abs(positive_leg.state) - np.abs(negative_leg.state)))
            )
            for positive_leg, negative_leg in (
                (positive_cycle.upward, negative_cycle.upward),
                (positive_cycle.downward, negative_cycle.downward),
            )
        )
    symmetry_matches = all(
        symmetry_difference_by_duration[duration]
        <= 10.0 * max(state_error_by_duration[duration], 1e-12)
        for duration in threshold_by_duration
    )
    if resolution_pass and branch_signs_match and symmetry_matches:
        h3 = "Supported within the tested conditions"
    elif resolution_pass and branch_signs_match is False:
        h3 = "Not supported within the tested conditions"
    else:
        h3 = "Inconclusive"

    zero_rows = [row for row in metric_rows if row["initial_state"] == 0.0]
    zero_invariant = all(row["maximum_absolute_state"] <= 1e-14 for row in zero_rows)
    return {
        "numerical_resolution_accepted": resolution_pass,
        "maximum_loop_area_change_by_leg_duration": {
            str(duration): error for duration, error in resolution_errors.items()
        },
        "H1_finite_rate_path_difference": {
            "classification": h1,
            "evidence": {
                "loop_area_exceeds_10x_numerical_change_by_duration_and_sign": {
                    str(duration): passed
                    for duration, passed in per_duration_signs.items()
                },
                "null_equilibrium_loop_area": 0.0,
            },
        },
        "H2_slower_ramp_reduces_loop": {
            "classification": h2,
            "evidence": {
                "fast_leg_duration": fast_duration,
                "slow_leg_duration": slow_duration,
                "fast_mean_loop_area": fast_area,
                "slow_mean_loop_area": slow_area,
                "fast_minus_slow_area": rate_difference,
                "frozen_resolution_threshold": rate_threshold,
            },
        },
        "H3_initial_sign_selects_pitchfork_branch": {
            "classification": h3,
            "evidence": {
                "branch_signs_match_initial_signs": branch_signs_match,
                "absolute_response_symmetry_matches_numerical_tolerance": (
                    symmetry_matches
                ),
                "maximum_absolute_response_symmetry_difference": {
                    str(duration): value
                    for duration, value in symmetry_difference_by_duration.items()
                },
                "upper_turning_point_states": {
                    str(duration): {
                        "positive_initial": positive["upper_turning_point_state"],
                        "negative_initial": negative["upper_turning_point_state"],
                    }
                    for duration, (positive, negative) in signed_pairs.items()
                },
            },
        },
        "zero_initial_state_invariant_control": {
            "observed_invariant": zero_invariant,
            "maximum_absolute_state": max(
                row["maximum_absolute_state"] for row in zero_rows
            ),
            "interpretation": (
                "The exactly zero deterministic state remains invariant even after "
                "it becomes unstable; it is not a noisy or robust physical "
                "initialization."
            ),
        },
        "quasi_static_equilibrium_hysteresis": {
            "classification": (
                "Not supported by the declared pitchfork equilibrium structure"
            ),
            "reason": (
                "The stable equilibrium magnitude is single-valued across mu and "
                "stable branches do not coexist across a control-parameter interval."
            ),
            "not_a_claim": (
                "This does not rule out rate-dependent loops or hysteresis in "
                "other models."
            ),
        },
    }


def _read_protocol(protocol_path: Path) -> tuple[dict[str, Any], str]:
    raw = protocol_path.read_bytes()
    digest = hashlib.sha256(raw).hexdigest().upper()
    if digest != FROZEN_PROTOCOL_SHA256:
        raise ScientificValidationError("Phase 62 frozen protocol hash mismatch")
    protocol = json.loads(raw)
    if (
        protocol.get("phase") != 62
        or protocol.get("protocol_status") != "frozen_before_primary_evaluation"
        or protocol.get("phase32_status") != "BLOCKED_AUTHORIZATION"
        or protocol.get("phase52_status") != "NO-GO"
        or protocol.get("kuramoto_thread_status")
        != "CLOSED_INCONCLUSIVE_UNDER_DECLARED_PROTOCOL"
    ):
        raise ScientificValidationError("Phase 62 protocol status is invalid")
    return protocol, digest


def run_pitchfork_path_study(output_directory: Path) -> dict[str, Any]:
    """Run the frozen path-dependence benchmark and write manifested outputs."""
    output = Path(output_directory)
    protocol_path = output / PROTOCOL_FILENAME
    if not protocol_path.is_file():
        raise FileNotFoundError(f"frozen Phase 62 protocol not found: {protocol_path}")
    unexpected = [
        path.name for path in output.iterdir() if path.name != PROTOCOL_FILENAME
    ]
    if unexpected:
        raise FileExistsError(
            "refusing to overwrite Phase 62 outputs: " + ", ".join(sorted(unexpected))
        )
    protocol, protocol_hash = _read_protocol(protocol_path)
    sweep = protocol["sweep"]
    solver_data = protocol["solver"]["primary"]
    solver = ODESolverConfiguration(
        method=solver_data["method"],
        relative_tolerance=solver_data["relative_tolerance"],
        absolute_tolerance=solver_data["absolute_tolerance"],
        maximum_step=solver_data["maximum_step"],
    )
    mu_min = float(sweep["mu_min"])
    mu_max = float(sweep["mu_max"])
    reference_area = (2.0 / 3.0) * mu_max**1.5
    durations = tuple(float(value) for value in sweep["ramp_leg_durations"])
    initial_states = tuple(float(value) for value in sweep["initial_states"])
    base_interval = float(sweep["base_output_interval"])

    cycles: list[PitchforkCycle] = []
    trajectory_rows: list[dict[str, Any]] = []
    metric_rows: list[dict[str, Any]] = []
    for duration in durations:
        for initial_state in initial_states:
            cycle = simulate_pitchfork_cycle(
                mu_min=mu_min,
                mu_max=mu_max,
                initial_state=initial_state,
                leg_duration=duration,
                output_interval=base_interval,
                solver=solver,
            )
            cycles.append(cycle)
            metric_rows.append(_cycle_metrics(cycle, reference_area))
            trajectory_rows.extend(_cycle_rows(cycle))

    numerical_rows: list[dict[str, Any]] = []
    numerical_output_interval = float(
        protocol["solver"]["output_resolution_check"]["output_interval"]
    )
    tight_data = protocol["solver"]["tolerance_check"]
    tight_solver = ODESolverConfiguration(
        method=tight_data["method"],
        relative_tolerance=tight_data["relative_tolerance"],
        absolute_tolerance=tight_data["absolute_tolerance"],
        maximum_step=tight_data["maximum_step"],
    )
    resolution_count = 0
    for duration in durations:
        primary = next(
            cycle
            for cycle in cycles
            if cycle.leg_duration == duration
            and cycle.initial_state == initial_states[0]
        )
        output_refined = simulate_pitchfork_cycle(
            mu_min=mu_min,
            mu_max=mu_max,
            initial_state=initial_states[0],
            leg_duration=duration,
            output_interval=numerical_output_interval,
            solver=solver,
        )
        tight = simulate_pitchfork_cycle(
            mu_min=mu_min,
            mu_max=mu_max,
            initial_state=initial_states[0],
            leg_duration=duration,
            output_interval=base_interval,
            solver=tight_solver,
        )
        for comparison, candidate in (
            ("output_grid_refinement", output_refined),
            ("tighter_solver_settings", tight),
        ):
            numerical_rows.append(
                {
                    "leg_duration": duration,
                    "comparison": comparison,
                    "primary_loop_area": _loop_area(primary),
                    "refined_loop_area": _loop_area(candidate),
                    "loop_area_absolute_difference": abs(
                        _loop_area(primary) - _loop_area(candidate)
                    ),
                    "maximum_state_difference": _interpolated_state_difference(
                        primary, candidate
                    ),
                    "primary_normalized_loop_area": _loop_area(primary)
                    / reference_area,
                    "refined_normalized_loop_area": _loop_area(candidate)
                    / reference_area,
                    "analytical_reference_area": reference_area,
                }
            )
            resolution_count += 2

    hypotheses = _classify_hypotheses(
        metric_rows, numerical_rows, cycles, reference_area, protocol
    )
    _write_csv(output / "ramp_trajectories.csv", trajectory_rows)
    _write_csv(output / "path_metrics.csv", metric_rows)
    _write_csv(output / "numerical_checks.csv", numerical_rows)
    (output / "hypothesis_results.json").write_text(
        json.dumps(hypotheses, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    metadata = {
        "phase": 62,
        "study": "supercritical_pitchfork_finite_rate_path_dependence",
        "protocol_sha256": protocol_hash,
        "primary_integrations": len(durations) * len(initial_states) * 2,
        "numerical_resolution_integrations": resolution_count,
        "total_successful_integrations": len(durations) * len(initial_states) * 2
        + resolution_count,
        "trajectory_rows": len(trajectory_rows),
        "path_metric_rows": len(metric_rows),
        "numerical_check_rows": len(numerical_rows),
        "analytical_reference_area": reference_area,
        "hypothesis_classifications": {
            key: value["classification"]
            for key, value in hypotheses.items()
            if isinstance(value, dict) and "classification" in value
        },
        "numerical_resolution_accepted": hypotheses["numerical_resolution_accepted"],
        "python_version": platform.python_version(),
        "platform": platform.platform(),
        "numpy_version": np.__version__,
        "scipy_version": scipy.__version__,
        "matplotlib_version": matplotlib.__version__,
        "phase32_status": protocol["phase32_status"],
        "phase52_status": protocol["phase52_status"],
        "kuramoto_thread_status": protocol["kuramoto_thread_status"],
        "external_data_or_services": False,
    }
    (output / "metadata.json").write_text(
        json.dumps(metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    _save_figures(output, cycles, metric_rows, numerical_rows, reference_area)
    manifest = {
        path.name: hashlib.sha256(path.read_bytes()).hexdigest().upper()
        for path in sorted(output.iterdir())
        if path.is_file() and path.name != "sha256_manifest.json"
    }
    (output / "sha256_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return metadata

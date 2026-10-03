"""Bounded matched-spectrum transient-growth study for a 2x2 matrix family."""

from __future__ import annotations

import csv
import hashlib
import json
import platform
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import scipy  # type: ignore[import-untyped]
from numpy.typing import NDArray
from scipy.integrate import solve_ivp  # type: ignore[import-untyped]
from scipy.linalg import expm  # type: ignore[import-untyped]

from newton_lab.exceptions import ScientificValidationError
from newton_lab.simulation import ODESolverConfiguration

PROTOCOL_FILENAME = "protocol.json"
EXPLORATION_FILENAME = "design_exploration.json"
FROZEN_PROTOCOL_SHA256 = (
    "6AE93A4020FE12E19D1B26F5EF36DF7EE07F547ACBAB5C590A74B73B78B1D523"
)

FloatMatrix = NDArray[np.float64]


@dataclass(frozen=True, slots=True)
class MatrixPair:
    """Normal and non-normal matrices with matched stable eigenvalues."""

    normal: FloatMatrix
    non_normal: FloatMatrix
    eigenvalues: tuple[float, float]
    a: float
    b: float
    kappa: float


@dataclass(frozen=True, slots=True)
class GainCurve:
    """Induced 2-norm gain samples over a finite time grid."""

    times: NDArray[np.float64]
    gain: NDArray[np.float64]


def build_matrix_pair(a: float, b: float, kappa: float) -> MatrixPair:
    """Build the declared two-state family after validating its parameters."""
    values = np.asarray((a, b, kappa), dtype=np.float64)
    if not np.isfinite(values).all():
        raise ScientificValidationError("a, b, and kappa must be finite")
    if a <= 0.0 or b <= 0.0:
        raise ScientificValidationError("a and b must be positive for stability")
    if a == b:
        raise ScientificValidationError("a and b must differ in this protocol")
    normal = np.asarray([[-a, 0.0], [0.0, -b]], dtype=np.float64)
    non_normal = np.asarray([[-a, kappa], [0.0, -b]], dtype=np.float64)
    return MatrixPair(
        normal=normal,
        non_normal=non_normal,
        eigenvalues=(-a, -b),
        a=float(a),
        b=float(b),
        kappa=float(kappa),
    )


def analytic_propagator(a: float, b: float, kappa: float, time: float) -> FloatMatrix:
    """Return exp(A_kappa*time) for a != b and nonnegative finite time."""
    values = np.asarray((a, b, kappa, time), dtype=np.float64)
    if not np.isfinite(values).all():
        raise ScientificValidationError("propagator inputs must be finite")
    if a <= 0.0 or b <= 0.0 or a == b:
        raise ScientificValidationError("require positive, distinct a and b")
    if time < 0.0:
        raise ScientificValidationError("time must be non-negative")
    decay_a = float(np.exp(-a * time))
    decay_b = float(np.exp(-b * time))
    coupling = kappa * (decay_a - decay_b) / (b - a)
    return np.asarray([[decay_a, coupling], [0.0, decay_b]], dtype=np.float64)


def calculate_gain_curve(
    a: float, b: float, kappa: float, times: NDArray[np.float64]
) -> GainCurve:
    """Calculate sampled induced 2-norm gains for a finite increasing grid."""
    samples = np.asarray(times, dtype=np.float64)
    if (
        samples.ndim != 1
        or len(samples) < 2
        or not np.isfinite(samples).all()
        or samples[0] < 0.0
        or not np.all(np.diff(samples) > 0.0)
    ):
        raise ScientificValidationError(
            "times must be a finite, increasing grid of nonnegative values"
        )
    gains = np.empty(len(samples), dtype=np.float64)
    for index, time in enumerate(samples):
        propagator = analytic_propagator(a, b, kappa, float(time))
        gains[index] = np.linalg.svd(propagator, compute_uv=False)[0]
    return GainCurve(times=samples.copy(), gain=gains)


def _leading_right_singular_vector(matrix: FloatMatrix) -> NDArray[np.float64]:
    _, _, right_vectors = np.linalg.svd(matrix)
    vector = np.asarray(right_vectors[0], dtype=np.float64).copy()
    pivot = int(np.argmax(np.abs(vector)))
    if vector[pivot] < 0.0:
        vector *= -1.0
    return vector


def _write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise ScientificValidationError(f"refusing to write empty result table: {path}")
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def _read_frozen_inputs(output: Path) -> tuple[dict[str, Any], dict[str, Any]]:
    protocol_path = output / PROTOCOL_FILENAME
    exploration_path = output / EXPLORATION_FILENAME
    raw = protocol_path.read_bytes()
    digest = hashlib.sha256(raw).hexdigest().upper()
    if digest != FROZEN_PROTOCOL_SHA256:
        raise ScientificValidationError("Phase 64 frozen protocol hash mismatch")
    protocol = json.loads(raw)
    exploration = json.loads(exploration_path.read_text(encoding="utf-8"))
    if (
        protocol.get("phase") != 64
        or protocol.get("protocol_status") != "frozen_before_primary_evaluation"
        or protocol.get("design_gate", {}).get("status")
        != "PASS_SCOPED_ABSTRACT_STATE_SPACE"
        or protocol.get("preserved_status", {}).get("phase32")
        != "BLOCKED_AUTHORIZATION"
        or protocol.get("preserved_status", {}).get("phase52") != "NO-GO"
        or exploration.get("stage") != "exploratory_design_only_not_primary_evaluation"
    ):
        raise ScientificValidationError("Phase 64 protocol or design gate is invalid")
    return protocol, exploration


def _integrate_condition(
    *,
    condition_id: str,
    kappa: float,
    initial_state: NDArray[np.float64],
    times: NDArray[np.float64],
    solver_data: dict[str, Any],
    trajectory_acceptance_tolerance: float,
    tight: bool,
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    pair = build_matrix_pair(1.0, 2.0, kappa)
    matrix = pair.non_normal
    solver = ODESolverConfiguration(
        method=solver_data["method"],
        relative_tolerance=solver_data["relative_tolerance"],
        absolute_tolerance=solver_data["absolute_tolerance"],
        maximum_step=solver_data["maximum_step"],
    )
    try:
        solution = solve_ivp(
            lambda _time, state: matrix @ state,
            (float(times[0]), float(times[-1])),
            initial_state.copy(),
            method=solver.method,
            t_eval=times,
            rtol=solver.relative_tolerance,
            atol=solver.absolute_tolerance,
            max_step=(np.inf if solver.maximum_step is None else solver.maximum_step),
        )
        solver_success = bool(solution.success)
        message = str(solution.message)
        returned_times = np.asarray(solution.t, dtype=np.float64)
        returned_states = np.asarray(solution.y, dtype=np.float64)
        function_evaluations = int(solution.nfev)
        if (
            not solver_success
            or returned_states.shape != (2, len(returned_times))
            or len(returned_times) != len(times)
            or not np.isfinite(returned_times).all()
            or not np.isfinite(returned_states).all()
        ):
            row = {
                "condition_id": condition_id,
                "kappa": kappa,
                "tight_solver": tight,
                "solver_status": "failed",
                "solver_success": solver_success,
                "message": message or "solver returned incomplete or nonfinite data",
                "function_evaluations": function_evaluations,
                "max_absolute_state_error": None,
                "numerical_accepted": False,
                "initial_x1": float(initial_state[0]),
                "initial_x2": float(initial_state[1]),
            }
            return row, []
        exact_states = np.column_stack(
            [
                analytic_propagator(1.0, 2.0, kappa, float(time)) @ initial_state
                for time in times
            ]
        )
        state_errors = np.linalg.norm(returned_states - exact_states, axis=0)
        maximum_error = float(np.max(state_errors))
        accepted = maximum_error <= trajectory_acceptance_tolerance
        row = {
            "condition_id": condition_id,
            "kappa": kappa,
            "tight_solver": tight,
            "solver_status": "success" if accepted else "error_above_tolerance",
            "solver_success": solver_success,
            "message": message,
            "function_evaluations": function_evaluations,
            "max_absolute_state_error": maximum_error,
            "numerical_accepted": accepted,
            "initial_x1": float(initial_state[0]),
            "initial_x2": float(initial_state[1]),
        }
        samples = [
            {
                "condition_id": condition_id,
                "kappa": kappa,
                "tight_solver": tight,
                "time": float(time),
                "x1": float(returned_states[0, index]),
                "x2": float(returned_states[1, index]),
                "state_norm": float(np.linalg.norm(returned_states[:, index])),
                "analytic_state_error": float(state_errors[index]),
            }
            for index, time in enumerate(returned_times)
        ]
        return row, samples
    except Exception as error:  # Preserve solver failures as explicit evidence.
        return (
            {
                "condition_id": condition_id,
                "kappa": kappa,
                "tight_solver": tight,
                "solver_status": "exception",
                "solver_success": False,
                "message": f"{type(error).__name__}: {error}",
                "function_evaluations": None,
                "max_absolute_state_error": None,
                "numerical_accepted": False,
                "initial_x1": float(initial_state[0]),
                "initial_x2": float(initial_state[1]),
            },
            [],
        )


def run_non_normal_transient_growth_study(output_directory: Path) -> dict[str, Any]:
    """Run the frozen, model-specific analytic/numerical Phase 64 study."""
    output = Path(output_directory)
    expected_inputs = {PROTOCOL_FILENAME, EXPLORATION_FILENAME}
    unexpected = [
        path.name for path in output.iterdir() if path.name not in expected_inputs
    ]
    if unexpected:
        raise FileExistsError(
            "refusing to overwrite Phase 64 outputs: " + ", ".join(sorted(unexpected))
        )
    protocol, exploration = _read_frozen_inputs(output)
    model = protocol["model"]
    params = model["parameters"]
    a = float(params["a"])
    b = float(params["b"])
    kappas = tuple(float(value) for value in params["kappa_values"])
    protocol_section = protocol["protocol"]
    horizon_start = float(protocol_section["horizon"][0])
    horizon_end = float(protocol_section["horizon"][1])
    primary_times = np.linspace(
        horizon_start,
        horizon_end,
        int(protocol_section["primary_time_grid"]["points"]),
        dtype=np.float64,
    )
    refined_times = np.linspace(
        horizon_start,
        horizon_end,
        int(protocol_section["resolution_grid"]["points"]),
        dtype=np.float64,
    )

    gain_rows: list[dict[str, Any]] = []
    summary_rows: list[dict[str, Any]] = []
    all_curves: dict[float, tuple[GainCurve, GainCurve]] = {}
    numerical_resolution_pass = True
    for kappa in kappas:
        primary_curve = calculate_gain_curve(a, b, kappa, primary_times)
        refined_curve = calculate_gain_curve(a, b, kappa, refined_times)
        all_curves[kappa] = (primary_curve, refined_curve)
        primary_peak_index = int(np.argmax(primary_curve.gain))
        refined_peak_index = int(np.argmax(refined_curve.gain))
        primary_gain = float(primary_curve.gain[primary_peak_index])
        refined_gain = float(refined_curve.gain[refined_peak_index])
        primary_peak_time = float(primary_curve.times[primary_peak_index])
        refined_peak_time = float(refined_curve.times[refined_peak_index])
        positive_peak = float(np.max(primary_curve.gain[1:]))
        positive_peak_time = float(
            primary_curve.times[1 + np.argmax(primary_curve.gain[1:])]
        )
        resolution_gain_difference = abs(primary_gain - refined_gain)
        resolution_time_difference = abs(primary_peak_time - refined_peak_time)
        row_resolution_pass = resolution_gain_difference <= float(
            protocol_section["resolution_grid"]["acceptance_peak_gain_change_at_most"]
        ) and resolution_time_difference <= float(
            protocol_section["resolution_grid"]["acceptance_peak_time_change_at_most"]
        )
        numerical_resolution_pass &= row_resolution_pass
        normal_baseline = float(np.max(all_curves[0.0][0].gain))
        summary_rows.append(
            {
                "kappa": kappa,
                "normality": "normal" if kappa == 0.0 else "non_normal",
                "eigenvalue_1": -a,
                "eigenvalue_2": -b,
                "asymptotically_stable": bool(a > 0.0 and b > 0.0),
                "sampled_max_gain_primary": primary_gain,
                "sampled_peak_time_primary": primary_peak_time,
                "sampled_max_positive_time_gain_primary": positive_peak,
                "sampled_positive_peak_time_primary": positive_peak_time,
                "sampled_max_gain_refined": refined_gain,
                "sampled_peak_time_refined": refined_peak_time,
                "matched_normal_primary_max_gain": normal_baseline,
                "gain_difference_from_normal": primary_gain - normal_baseline,
                "grid_peak_gain_absolute_change": resolution_gain_difference,
                "grid_peak_time_absolute_change": resolution_time_difference,
                "resolution_accepted": row_resolution_pass,
            }
        )
        for index, time in enumerate(primary_curve.times):
            gain_rows.append(
                {
                    "kappa": kappa,
                    "time": float(time),
                    "gain_2_norm": float(primary_curve.gain[index]),
                    "matched_normal_gain": float(all_curves[0.0][0].gain[index])
                    if 0.0 in all_curves
                    else float(np.exp(-a * time)),
                }
            )

    normal_maximum = float(np.max(all_curves[0.0][0].gain))
    positive_threshold = 1.0 + float(protocol_section["gain_effect_threshold"])
    h1_conditions = [
        kappa
        for kappa in kappas
        if kappa != 0.0
        and float(np.max(all_curves[kappa][0].gain[1:])) > positive_threshold
        and float(np.max(all_curves[kappa][1].gain[1:])) > positive_threshold
    ]
    h2_conditions = [
        kappa
        for kappa in kappas
        if kappa != 0.0
        and float(np.max(all_curves[kappa][0].gain))
        >= normal_maximum + float(protocol_section["gain_effect_threshold"])
        and float(np.max(all_curves[kappa][1].gain))
        >= float(np.max(all_curves[0.0][1].gain))
        + float(protocol_section["gain_effect_threshold"])
    ]

    analytic_check_rows: list[dict[str, Any]] = []
    fixed_times = (0.0, 0.37, 1.2, 4.0)
    for kappa in kappas:
        pair = build_matrix_pair(a, b, kappa)
        for time in fixed_times:
            closed_form = analytic_propagator(a, b, kappa, time)
            library_reference = expm(pair.non_normal * time)
            difference = float(np.max(np.abs(closed_form - library_reference)))
            analytic_check_rows.append(
                {
                    "check": "closed_form_vs_scipy_expm",
                    "kappa": kappa,
                    "time": time,
                    "maximum_absolute_matrix_difference": difference,
                    "accepted": difference <= 1e-12,
                }
            )
    semigroup_pairs = ((0.2, 0.7), (0.8, 1.1))
    for kappa in kappas:
        for first, second in semigroup_pairs:
            left = analytic_propagator(a, b, kappa, first + second)
            right = analytic_propagator(a, b, kappa, first) @ analytic_propagator(
                a, b, kappa, second
            )
            difference = float(np.max(np.abs(left - right)))
            analytic_check_rows.append(
                {
                    "check": "semigroup_property",
                    "kappa": kappa,
                    "time": first + second,
                    "maximum_absolute_matrix_difference": difference,
                    "accepted": difference <= 1e-12,
                }
            )

    structure_rows: list[dict[str, Any]] = []
    for kappa in kappas:
        pair = build_matrix_pair(a, b, kappa)
        normal_eigenvalues = np.sort_complex(np.linalg.eigvals(pair.normal))
        non_normal_eigenvalues = np.sort_complex(np.linalg.eigvals(pair.non_normal))
        eigenvalue_difference = float(
            np.max(np.abs(normal_eigenvalues - non_normal_eigenvalues))
        )
        normality_residual = float(
            np.linalg.norm(
                pair.non_normal.T @ pair.non_normal
                - pair.non_normal @ pair.non_normal.T,
                ord=2,
            )
        )
        normal_matrix_residual = float(
            np.linalg.norm(
                pair.normal.T @ pair.normal - pair.normal @ pair.normal.T,
                ord=2,
            )
        )
        structure_rows.append(
            {
                "kappa": kappa,
                "matched_eigenvalues": eigenvalue_difference <= 1e-14,
                "eigenvalue_max_difference": eigenvalue_difference,
                "strictly_stable": bool(np.max(np.real(non_normal_eigenvalues)) < 0.0),
                "normal_matrix_normality_residual": normal_matrix_residual,
                "normal_matrix_is_normal": normal_matrix_residual <= 1e-14,
                "normality_residual": normality_residual,
                "normality_matches_declared_coupling": bool(
                    (kappa == 0.0 and normality_residual == 0.0)
                    or (kappa != 0.0 and normality_residual > 0.0)
                ),
            }
        )

    trajectory_rows: list[dict[str, Any]] = []
    trajectory_samples: list[dict[str, Any]] = []
    conditions: list[tuple[str, float, NDArray[np.float64], bool]] = []
    for item in protocol_section["numerical_comparisons"]:
        kappa = float(item["kappa"])
        curve = all_curves[kappa][0]
        peak_index = int(np.argmax(curve.gain))
        if "leading right singular vector" in item["initial_condition"]:
            p = analytic_propagator(a, b, kappa, float(curve.times[peak_index]))
            initial_state = _leading_right_singular_vector(p)
            condition_id = f"primary_kappa_{kappa:g}_operator_gain_direction"
        else:
            initial_state = np.asarray([0.0, 1.0], dtype=np.float64)
            condition_id = f"primary_kappa_{kappa:g}_coordinate_e2"
        conditions.append((condition_id, kappa, initial_state, False))
    for kappa in (4.0, 8.0):
        curve = all_curves[kappa][0]
        peak_index = int(np.argmax(curve.gain))
        p = analytic_propagator(a, b, kappa, float(curve.times[peak_index]))
        conditions.append(
            (
                f"tight_kappa_{kappa:g}_operator_gain_direction",
                kappa,
                _leading_right_singular_vector(p),
                True,
            )
        )
    for condition_id, kappa, initial_state, tight in conditions:
        solver_data = protocol_section["tight_solver" if tight else "solver"]
        output_times = np.linspace(
            horizon_start,
            horizon_end,
            int(solver_data["output_points"]),
            dtype=np.float64,
        )
        row, samples = _integrate_condition(
            condition_id=condition_id,
            kappa=kappa,
            initial_state=initial_state,
            times=output_times,
            solver_data=solver_data,
            trajectory_acceptance_tolerance=float(
                protocol_section["trajectory_acceptance_tolerance"]
            ),
            tight=tight,
        )
        trajectory_rows.append(row)
        trajectory_samples.extend(samples)

    all_structure_checks = all(
        row["matched_eigenvalues"]
        and row["strictly_stable"]
        and row["normal_matrix_is_normal"]
        and row["normality_matches_declared_coupling"]
        for row in structure_rows
    )
    normal_reference_curve = np.maximum(
        np.exp(-a * primary_times), np.exp(-b * primary_times)
    )
    zero_coupling_difference = float(
        np.max(np.abs(all_curves[0.0][0].gain - normal_reference_curve))
    )
    analytic_check_rows.append(
        {
            "check": "zero_coupling_matches_normal_gain",
            "kappa": 0.0,
            "time": None,
            "maximum_absolute_matrix_difference": zero_coupling_difference,
            "accepted": zero_coupling_difference <= 1e-14,
        }
    )
    all_analytic_checks = all(row["accepted"] for row in analytic_check_rows)
    all_trajectories_accepted = all(
        row["numerical_accepted"] for row in trajectory_rows
    )
    h3_supported = (
        all_analytic_checks
        and all_structure_checks
        and numerical_resolution_pass
        and all_trajectories_accepted
        and normal_maximum <= 1.0 + 1e-12
    )
    unresolved = not (
        all_analytic_checks and all_structure_checks and numerical_resolution_pass
    )
    hypotheses = {
        "H1_finite_time_amplification": {
            "classification": (
                "INCONCLUSIVE"
                if unresolved
                else "SUPPORTED_WITHIN_FROZEN_FAMILY"
                if h1_conditions
                else "NOT_SUPPORTED_ON_FROZEN_GRID"
            ),
            "effect_threshold": positive_threshold,
            "kappa_values_exceeding_threshold_on_both_grids": h1_conditions,
        },
        "H2_exceeds_matched_normal": {
            "classification": (
                "INCONCLUSIVE"
                if unresolved
                else "SUPPORTED_WITHIN_FROZEN_FAMILY"
                if h2_conditions
                else "NOT_SUPPORTED_ON_FROZEN_GRID"
            ),
            "absolute_gain_difference_threshold": float(
                protocol_section["gain_effect_threshold"]
            ),
            "normal_sampled_max_gain": normal_maximum,
            "kappa_values_exceeding_threshold_on_both_grids": h2_conditions,
        },
        "H3_analytic_numerical_agreement": {
            "classification": (
                "SUPPORTED_WITHIN_FROZEN_TOLERANCES"
                if h3_supported
                else "NOT_SUPPORTED_OR_UNRESOLVED"
            ),
            "all_analytic_propagator_checks_accepted": all_analytic_checks,
            "all_structure_checks_accepted": all_structure_checks,
            "time_grid_resolution_accepted": numerical_resolution_pass,
            "all_trajectory_checks_accepted": all_trajectories_accepted,
        },
        "normal_baseline_non_amplifying": {
            "classification": (
                "SUPPORTED" if normal_maximum <= 1.0 + 1e-12 else "FAILED_CHECK"
            ),
            "sampled_max_gain": normal_maximum,
        },
        "zero_coupling_matches_normal": {
            "classification": (
                "SUPPORTED" if zero_coupling_difference <= 1e-14 else "FAILED_CHECK"
            ),
            "maximum_gain_curve_difference": zero_coupling_difference,
        },
    }

    _write_csv(output / "gain_curves.csv", gain_rows)
    _write_csv(output / "gain_summary.csv", summary_rows)
    _write_csv(output / "analytic_reference_checks.csv", analytic_check_rows)
    _write_csv(output / "matrix_structure_checks.csv", structure_rows)
    _write_csv(output / "trajectory_checks.csv", trajectory_rows)
    if trajectory_samples:
        _write_csv(output / "trajectory_samples.csv", trajectory_samples)
    (output / "hypothesis_results.json").write_text(
        json.dumps(hypotheses, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    metadata = {
        "phase": 64,
        "study": "matched_spectrum_non_normal_finite_time_gain",
        "protocol_sha256": FROZEN_PROTOCOL_SHA256,
        "design_exploration_stage": exploration["stage"],
        "state_space": "abstract_dimensionless_R2",
        "norm": "Euclidean_induced_matrix_2_norm",
        "kappa_values": list(kappas),
        "primary_grid_points_per_kappa": len(primary_times),
        "resolution_grid_points_per_kappa": len(refined_times),
        "gain_curve_rows": len(gain_rows),
        "analytic_reference_check_rows": len(analytic_check_rows),
        "matrix_structure_check_rows": len(structure_rows),
        "trajectory_integrations_attempted": len(trajectory_rows),
        "trajectory_integrations_successful": sum(
            row["solver_success"] for row in trajectory_rows
        ),
        "trajectory_sample_rows": len(trajectory_samples),
        "python_version": platform.python_version(),
        "platform": platform.platform(),
        "numpy_version": np.__version__,
        "scipy_version": scipy.__version__,
        "phase32_status": protocol["preserved_status"]["phase32"],
        "phase52_status": protocol["preserved_status"]["phase52"],
        "kuramoto_thread_status": protocol["preserved_status"]["kuramoto"],
        "phase62_status": protocol["preserved_status"]["phase62"],
        "external_data_or_services": False,
    }
    (output / "metadata.json").write_text(
        json.dumps(metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    manifest = {
        path.name: hashlib.sha256(path.read_bytes()).hexdigest().upper()
        for path in sorted(output.iterdir())
        if path.is_file() and path.name != "sha256_manifest.json"
    }
    (output / "sha256_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return metadata

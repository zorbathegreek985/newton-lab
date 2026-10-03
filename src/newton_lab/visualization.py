"""Headless scientific plots and Markdown reporting for Newton Lab results.

Matplotlib is optional. Install ``newton-lab[visualization]`` to import and use
this module. Functions consume existing results or analysis projections and
never rerun a simulation.
"""

from __future__ import annotations

import os
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, TypeAlias
from urllib.parse import quote

import numpy as np

if TYPE_CHECKING:
    from matplotlib.axes import Axes
    from matplotlib.figure import Figure

from newton_lab.boundary_value import BoundaryValueResult
from newton_lab.dynamics import DampedOscillatorResult
from newton_lab.exceptions import NewtonLabError, ScientificValidationError
from newton_lab.experiment_analysis import (
    ExperimentCollectionSummary,
    MetricComparison,
    MetricObservation,
    MetricStatistics,
    ParameterSweepAnalysis,
    RunAnalysisProjection,
    analyze_parameter_sweep,
    extract_metric_observations,
    project_experiment_runs,
    summarize_experiments,
    summarize_metric,
)
from newton_lab.experiments import (
    BVPReferenceComparison,
    ExperimentRecord,
    ExperimentRunStatus,
)
from newton_lab.pendulum import DampedPendulumResult

ODETrajectory: TypeAlias = DampedOscillatorResult | DampedPendulumResult
BVPProfileSource: TypeAlias = BoundaryValueResult | RunAnalysisProjection


class VisualizationError(NewtonLabError):
    """Raised when supplied numerical data cannot be plotted faithfully."""


@dataclass(slots=True)
class VisualizationArtifact:
    """Headless Matplotlib figure and provenance for the plotted source data.

    The caller owns the figure lifecycle. Call ``close()`` after display/export
    or when the figure is no longer needed.
    """

    figure: Figure
    saved_path: Path | None
    source_ids: tuple[str, ...]
    notes: tuple[str, ...] = ()

    def close(self) -> None:
        """Clear artists and release the figure's in-memory plot data."""
        self.figure.clear()


@dataclass(frozen=True, slots=True)
class ReportPlotReference:
    """A saved figure included in a Markdown experiment report."""

    path: Path
    label: str


@dataclass(frozen=True, slots=True)
class ExperimentReport:
    """Deterministic report content and validated references to saved plots."""

    body: str
    source_experiment_ids: tuple[str, ...]
    summary: ExperimentCollectionSummary
    plots: tuple[ReportPlotReference, ...]

    def render_markdown(self, *, relative_to: Path | None = None) -> str:
        """Render report text, using plot links relative to a destination folder."""
        if not self.plots:
            return self.body
        lines = [self.body.rstrip(), "", "## Generated plots", ""]
        for plot in self.plots:
            if not plot.path.is_file():
                raise FileNotFoundError(f"report plot does not exist: {plot.path}")
            shown_path = (
                Path(os.path.relpath(plot.path, relative_to))
                if relative_to is not None
                else plot.path
            )
            link_path = quote(shown_path.as_posix(), safe="/._-")
            lines.append(f"- [{plot.label}]({link_path})")
        return "\n".join(lines) + "\n"


def _new_figure(
    rows: int = 1, *, height: float = 4.0
) -> tuple[Figure, tuple[Axes, ...]]:
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.figure import Figure

    figure = Figure(figsize=(8.0, height * rows), layout="constrained")
    FigureCanvasAgg(figure)
    axes = tuple(figure.add_subplot(rows, 1, index + 1) for index in range(rows))
    return figure, axes


def _format_label(name: str, unit: str | None) -> str:
    return f"{name} ({unit})" if unit else f"{name} (unit unknown)"


def _with_context(kind: str, source: str, context_title: str | None) -> str:
    if context_title:
        return f"{context_title}: {kind}: {source}"
    return f"{kind}: {source}"


def _finish(
    figure: Figure,
    *,
    source_ids: tuple[str, ...],
    notes: tuple[str, ...] = (),
    save_path: str | Path | None = None,
    overwrite: bool = False,
    create_parent_dirs: bool = False,
) -> VisualizationArtifact:
    target: Path | None = None
    if save_path is not None:
        target = Path(save_path)
        suffix = target.suffix.lower()
        if suffix not in {".png", ".svg"}:
            raise ScientificValidationError("plot output path must end in .png or .svg")
        if target.exists() and not overwrite:
            raise FileExistsError(f"plot output already exists: {target}")
        if target.exists() and target.is_dir():
            raise IsADirectoryError(f"plot output path is a directory: {target}")
        if not target.parent.exists():
            if create_parent_dirs:
                target.parent.mkdir(parents=True, exist_ok=True)
            else:
                raise FileNotFoundError(
                    f"plot parent directory does not exist: {target.parent}"
                )
        metadata = (
            {"Software": "Newton Lab"}
            if suffix == ".png"
            else {"Creator": "Newton Lab", "Date": None}
        )
        try:
            import matplotlib as mpl

            with mpl.rc_context({"svg.hashsalt": "newton-lab"}):
                figure.savefig(
                    target,
                    format=suffix.removeprefix("."),
                    dpi=120,
                    bbox_inches="tight",
                    metadata=metadata,
                )
        except Exception as exc:
            raise VisualizationError(
                f"could not save plot to {target}: {type(exc).__name__}: {exc}"
            ) from exc
    return VisualizationArtifact(
        figure=figure,
        saved_path=target,
        source_ids=source_ids,
        notes=notes,
    )


def _vector(values: object, *, name: str, minimum_size: int = 2) -> np.ndarray:
    try:
        vector = np.asarray(values, dtype=np.float64)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ScientificValidationError(f"{name} must be a numeric vector") from exc
    if vector.ndim != 1 or len(vector) < minimum_size:
        raise ScientificValidationError(
            f"{name} must be a one-dimensional vector with at least "
            f"{minimum_size} values"
        )
    if not np.isfinite(vector).all():
        raise ScientificValidationError(f"{name} must contain only finite values")
    return vector


def plot_ode_trajectory(
    result: ODETrajectory,
    *,
    components: tuple[str, ...],
    source_id: str | None = None,
    title: str | None = None,
    save_path: str | Path | None = None,
    overwrite: bool = False,
    create_parent_dirs: bool = False,
) -> VisualizationArtifact:
    """Plot selected ODE state components at their actual returned time samples."""
    if isinstance(result, DampedOscillatorResult):
        source = source_id or "damped_harmonic_oscillator"
        time_values = result.time_s
        available = {
            "displacement_m": (result.displacement_m, "m"),
            "velocity_m_per_s": (result.velocity_m_per_s, "m/s"),
        }
    else:
        source = source_id or "damped_nonlinear_pendulum"
        time_values = result.time_s
        available = {
            "angle_rad": (result.angle_rad, "rad"),
            "angular_velocity_rad_per_s": (
                result.angular_velocity_rad_per_s,
                "rad/s",
            ),
        }
    if not components or len(components) != len(set(components)):
        raise ScientificValidationError(
            "select one or more unique ODE state components"
        )
    unknown = set(components) - set(available)
    if unknown:
        raise ScientificValidationError(f"unknown ODE components: {sorted(unknown)}")
    time = _vector(time_values, name="time samples")
    if not np.all(np.diff(time) > 0.0):
        raise ScientificValidationError("time samples must be strictly increasing")
    selected: list[tuple[str, np.ndarray, str]] = []
    for component_name in components:
        raw, unit = available[component_name]
        values = _vector(raw, name=component_name)
        if len(values) != len(time):
            raise ScientificValidationError(
                f"{component_name} length must match the returned time samples"
            )
        selected.append((component_name, values, unit))
    figure, axes = _new_figure(len(selected))
    figure.suptitle(_with_context("ODE trajectory", source, title))
    for axis, (component_name, values, unit) in zip(axes, selected, strict=True):
        axis.plot(time, values, label=_format_label(component_name, unit))
        axis.set_ylabel(_format_label(component_name, unit))
        axis.legend(loc="best")
        axis.grid(True, alpha=0.25)
    axes[-1].set_xlabel("Time (s)")
    return _finish(
        figure,
        source_ids=(source,),
        notes=(
            "Uses returned samples only; sampled extrema may miss "
            "between-sample behavior.",
        ),
        save_path=save_path,
        overwrite=overwrite,
        create_parent_dirs=create_parent_dirs,
    )


def _bvp_profile_data(
    source: BVPProfileSource,
) -> tuple[
    str,
    str,
    tuple[str, str],
    tuple[str | None, str | None],
    np.ndarray,
    np.ndarray,
    tuple[float, float],
    tuple[float, ...],
    tuple[str, ...],
    BVPReferenceComparison | None,
]:
    if isinstance(source, BoundaryValueResult):
        problem = source.problem
        source_id = problem.problem_id
        names: tuple[str, ...] = tuple(state.name for state in problem.states)
        units: tuple[str | None, ...] = tuple(state.unit for state in problem.states)
        domain = problem.domain
        domain_unit = problem.domain_unit
        mesh = np.asarray(source.mesh, dtype=np.float64)
        solution = np.asarray(source.solution, dtype=np.float64)
        parameters = tuple(item.value for item in problem.parameters)
        assumptions = problem.assumptions
        comparison = None
    else:
        if source.problem_family != "boundary_value" or source.bvp_diagnostics is None:
            raise ScientificValidationError(
                "run projection must contain BVP diagnostics"
            )
        diagnostics = source.bvp_diagnostics
        source_id = f"{source.experiment_id}#{source.run_index}"
        names = diagnostics.state_names
        units = diagnostics.state_units
        domain = diagnostics.domain
        domain_unit = diagnostics.domain_unit
        mesh = np.asarray(diagnostics.mesh, dtype=np.float64)
        solution = np.asarray(diagnostics.solution, dtype=np.float64)
        repro = source.reproducibility
        parameters = (
            () if repro is None else tuple(item.value for item in repro.parameters)
        )
        assumptions = () if repro is None else repro.assumptions
        comparison = source.bvp_reference_comparison
    if len(names) != 2 or len(units) != 2:
        raise ScientificValidationError(
            "BVP profile must describe exactly two state components"
        )
    if (
        mesh.ndim != 1
        or len(mesh) < 2
        or solution.shape != (2, len(mesh))
        or not np.isfinite(mesh).all()
        or not np.isfinite(solution).all()
        or not np.all(np.diff(mesh) > 0.0)
        or mesh[0] != domain[0]
        or mesh[-1] != domain[1]
    ):
        raise ScientificValidationError(
            "BVP mesh and profile must be finite, aligned, and domain preserving"
        )
    return (
        source_id,
        domain_unit,
        names,
        units,
        mesh,
        solution,
        domain,
        parameters,
        assumptions,
        comparison,
    )


def plot_bvp_profile(
    source: BVPProfileSource,
    *,
    components: tuple[str, ...],
    include_reference: bool = False,
    reference_comparison: BVPReferenceComparison | None = None,
    title: str | None = None,
    save_path: str | Path | None = None,
    overwrite: bool = False,
    create_parent_dirs: bool = False,
) -> VisualizationArtifact:
    """Plot returned BVP state profiles without resampling or smoothing."""
    (
        source_id,
        domain_unit,
        state_names,
        state_units,
        mesh,
        solution,
        domain,
        parameter_values,
        assumptions,
        embedded_comparison,
    ) = _bvp_profile_data(source)
    if not components or len(components) != len(set(components)):
        raise ScientificValidationError("select one or more unique BVP states")
    unknown = set(components) - set(state_names)
    if unknown:
        raise ScientificValidationError(f"unknown BVP states: {sorted(unknown)}")
    comparison = reference_comparison or embedded_comparison
    reference: np.ndarray | None = None
    if include_reference:
        if comparison is None or comparison.status != "passed":
            raise ScientificValidationError(
                "a passed, compatible BVP reference comparison is required"
            )
        if comparison.state_name != "temperature" or comparison.unit != "K":
            raise ScientificValidationError(
                "reference comparison is incompatible with temperature"
            )
        required_assumptions = {
            "Steady state.",
            "One-dimensional homogeneous material with constant thermal conductivity.",
            "No internal heat generation.",
        }
        param_map = dict(
            zip(_bvp_parameter_names(source), parameter_values, strict=True)
        )
        if (
            _bvp_problem_id(source) != "steady_linear_heat_conduction"
            or not required_assumptions.issubset(assumptions)
            or "temperature" not in state_names
            or state_units[state_names.index("temperature")] != "K"
            or not {"left_temperature_k", "right_temperature_k"}.issubset(param_map)
        ):
            raise ScientificValidationError(
                "heat reference assumptions do not match this BVP"
            )
        left = param_map["left_temperature_k"]
        right = param_map["right_temperature_k"]
        expected = left + (right - left) * (mesh - domain[0]) / (domain[1] - domain[0])
        actual = solution[state_names.index("temperature")]
        maximum_error = float(np.max(np.abs(actual - expected)))
        rms_error = float(np.sqrt(np.mean((actual - expected) ** 2)))
        if (
            comparison.maximum_absolute_error is None
            or comparison.rms_error is None
            or not np.isclose(
                maximum_error, comparison.maximum_absolute_error, rtol=1e-10, atol=1e-14
            )
            or not np.isclose(rms_error, comparison.rms_error, rtol=1e-10, atol=1e-14)
        ):
            raise ScientificValidationError(
                "reference comparison metrics do not match this returned BVP profile"
            )
        reference = expected
    figure, axes = _new_figure(len(components))
    figure.suptitle(_with_context("BVP spatial profile", source_id, title))
    for axis, state_name in zip(axes, components, strict=True):
        index = state_names.index(state_name)
        unit = state_units[index] or "unit unknown"
        axis.plot(mesh, solution[index], marker="o", markersize=3, label=state_name)
        if reference is not None and state_name == "temperature":
            axis.plot(mesh, reference, linestyle="--", label="Analytical reference")
        axis.set_ylabel(_format_label(state_name, unit))
        axis.set_xlabel(_format_label("Spatial coordinate", domain_unit))
        axis.legend(loc="best")
        axis.grid(True, alpha=0.25)
    notes = [
        "Uses the returned spatial mesh and nodal values without smoothing "
        "or resampling.",
        "Sampled profile appearance does not establish interior accuracy "
        "or physical validity.",
    ]
    if include_reference:
        notes.append(
            "Reference overlay uses the passed linear heat-profile comparison."
        )
    return _finish(
        figure,
        source_ids=(source_id,),
        notes=tuple(notes),
        save_path=save_path,
        overwrite=overwrite,
        create_parent_dirs=create_parent_dirs,
    )


def _bvp_problem_id(source: BVPProfileSource) -> str:
    if isinstance(source, BoundaryValueResult):
        return source.problem.problem_id
    if source.bvp_diagnostics is None:
        raise ScientificValidationError("run projection must contain BVP diagnostics")
    return source.bvp_diagnostics.problem_id


def _bvp_parameter_names(source: BVPProfileSource) -> tuple[str, ...]:
    if isinstance(source, BoundaryValueResult):
        return tuple(item.name for item in source.problem.parameters)
    if source.reproducibility is None:
        return ()
    return tuple(item.name for item in source.reproducibility.parameters)


def _status_for_observation(observation: MetricObservation) -> str:
    if observation.execution_status == ExperimentRunStatus.FAILED:
        return "execution failed"
    if observation.metric_status == "valid":
        return "valid"
    return {
        "missing": "metric missing",
        "non_finite": "non-finite metric",
        "invalid": "invalid metric",
        "valid": "valid",
    }[observation.metric_status]


def _plot_sweep_points(
    analysis: ParameterSweepAnalysis,
    *,
    kind: str,
    title: str | None,
    allow_unreliable_metadata: bool,
    extra_notes: tuple[str, ...],
    save_path: str | Path | None,
    overwrite: bool,
    create_parent_dirs: bool,
) -> VisualizationArtifact:
    if not analysis.ordering_reliable and not allow_unreliable_metadata:
        raise ScientificValidationError(
            f"sweep ordering is unreliable: {analysis.message}"
        )
    observations = tuple(point.observation for point in analysis.points)
    available = tuple(
        item
        for item in observations
        if item.metric_status == "valid" and item.value is not None
    )
    if not available:
        metric_unit = next((item.unit for item in observations if item.unit), None)
        definition = next(
            (item.definition for item in observations if item.definition), None
        )
    else:
        origins = {
            (item.model_id, item.problem_id, item.metric_id) for item in available
        }
        metric_units = {item.unit for item in available}
        metric_definitions = {item.definition for item in available}
        if len(origins) != 1 or len(metric_units) != 1 or len(metric_definitions) != 1:
            raise ScientificValidationError(
                "sweep points contain incompatible metric definitions or units"
            )
        metric_unit = available[0].unit
        definition = available[0].definition
    parameter_units = {
        next(
            (
                item.unit
                for item in point.observation.parameters
                if item.name == analysis.parameter_name
            ),
            None,
        )
        for point in analysis.points
    }
    if len(parameter_units) > 1:
        raise ScientificValidationError("sweep parameter units differ across runs")
    parameter_unit = next(iter(parameter_units), None)
    source_ids = tuple(
        f"{point.observation.experiment_id}#{point.run_index}"
        for point in analysis.points
    )
    figure, axes = _new_figure(rows=2, height=3.0)
    axes[0].set_title(_with_context(kind, analysis.experiment_id, title))
    for point in analysis.points:
        x = point.parameter_value
        if x is None:
            continue
        status = _status_for_observation(point.observation)
        if status == "valid" and point.observation.value is not None:
            axes[0].scatter([x], [point.observation.value], color="#1f77b4", zorder=3)
            axes[0].annotate(
                f"r{point.run_index}",
                (x, point.observation.value),
                xytext=(3, 4),
                textcoords="offset points",
                fontsize=7,
            )
        axes[1].scatter(
            [x],
            [status],
            marker="o" if status == "valid" else "x",
            color="#1f77b4" if status == "valid" else "#d62728",
        )
        axes[1].annotate(
            f"r{point.run_index}",
            (x, status),
            xytext=(3, 2),
            textcoords="offset points",
            fontsize=7,
        )
    axes[0].set_ylabel(
        f"{analysis.metric_id} ({metric_unit or 'unit unknown'})"
        + (f"\n{definition}" if definition else "")
    )
    axes[0].grid(True, alpha=0.25)
    axes[0].text(
        0.01,
        0.99,
        "Points only; no interpolation or connecting line.",
        transform=axes[0].transAxes,
        va="top",
        fontsize=8,
    )
    axes[1].set_yticks(
        [
            "valid",
            "execution failed",
            "metric missing",
            "non-finite metric",
            "invalid metric",
        ]
    )
    axes[1].set_ylabel("Run / metric status")
    axes[1].set_xlabel(_format_label(analysis.parameter_name, parameter_unit))
    axes[1].grid(True, axis="x", alpha=0.25)
    notes = list(extra_notes)
    if not analysis.ordering_reliable:
        notes.append(f"Sweep metadata unreliable: {analysis.message}")
    if any(item.unit in {None, "unknown"} for item in observations):
        notes.append("Some metric units are unavailable or unknown.")
    if any(item.metric_status != "valid" for item in observations):
        notes.append(
            "Failed or unavailable metric points are marked in the status panel."
        )
    return _finish(
        figure,
        source_ids=source_ids,
        notes=tuple(notes),
        save_path=save_path,
        overwrite=overwrite,
        create_parent_dirs=create_parent_dirs,
    )


def plot_metric_sweep(
    analysis: ParameterSweepAnalysis,
    *,
    title: str | None = None,
    allow_unreliable_metadata: bool = False,
    save_path: str | Path | None = None,
    overwrite: bool = False,
    create_parent_dirs: bool = False,
) -> VisualizationArtifact:
    """Plot a selected analysis metric as unsorted points plus a status panel."""
    return _plot_sweep_points(
        analysis,
        kind="Parameter sweep",
        title=title,
        allow_unreliable_metadata=allow_unreliable_metadata,
        extra_notes=(),
        save_path=save_path,
        overwrite=overwrite,
        create_parent_dirs=create_parent_dirs,
    )


def plot_algebraic_sweep(
    record: ExperimentRecord,
    *,
    parameter_name: str,
    metric_id: str,
    title: str | None = None,
    allow_unreliable_metadata: bool = False,
    save_path: str | Path | None = None,
    overwrite: bool = False,
    create_parent_dirs: bool = False,
) -> VisualizationArtifact:
    """Plot algebraic solution or residual metrics as disconnected candidates."""
    if record.specification.model_id != "nonlinear_algebraic_system":
        raise ScientificValidationError(
            "record is not a nonlinear algebraic experiment"
        )
    if not (
        metric_id == "residual_norm"
        or metric_id.startswith("solution:")
        or metric_id.startswith("residual:")
    ):
        raise ScientificValidationError(
            "algebraic plots require a solution, residual, or residual_norm metric"
        )
    analysis = analyze_parameter_sweep(
        record, parameter_name=parameter_name, metric_id=metric_id
    )
    return _plot_sweep_points(
        analysis,
        kind="Algebraic candidate sweep",
        title=title,
        allow_unreliable_metadata=allow_unreliable_metadata,
        extra_notes=(
            "Candidate points are not connected: root-branch identity is not "
            "established.",
            "A local solve and plotted candidate do not establish uniqueness "
            "or root completeness.",
        ),
        save_path=save_path,
        overwrite=overwrite,
        create_parent_dirs=create_parent_dirs,
    )


def plot_reference_error_sweep(
    analysis: ParameterSweepAnalysis,
    *,
    title: str | None = None,
    save_path: str | Path | None = None,
    overwrite: bool = False,
    create_parent_dirs: bool = False,
) -> VisualizationArtifact:
    """Plot recorded scalar or profile reference errors over a reliable sweep."""
    if not analysis.ordering_reliable:
        raise ScientificValidationError(
            f"sweep ordering is unreliable: {analysis.message}"
        )
    points = []
    for point in analysis.points:
        observation = point.observation
        reference = observation.reference_comparison
        bvp_reference = observation.bvp_reference_comparison
        if reference is not None:
            value = reference.absolute_error
            unit = reference.unit
            description = f"Absolute error for {reference.metric_id}"
            status = reference.status
        elif bvp_reference is not None:
            value = bvp_reference.maximum_absolute_error
            unit = bvp_reference.unit
            description = f"Maximum profile error for {bvp_reference.state_name}"
            status = bvp_reference.status
        else:
            value = None
            unit = None
            description = "Reference error"
            status = "unavailable"
        replacement = observation.model_copy(
            update={
                "metric_id": "reference_absolute_error",
                "value": value,
                "unit": unit,
                "definition": description,
                "metric_status": "valid" if value is not None else "missing",
                "observation_status": "valid" if value is not None else "missing",
            }
        )
        if status == "unavailable":
            replacement = replacement.model_copy(update={"metric_status": "missing"})
        points.append(point.model_copy(update={"observation": replacement}))
    if not any(point.observation.value is not None for point in points):
        raise ScientificValidationError("sweep contains no available reference errors")
    error_analysis = analysis.model_copy(
        update={"metric_id": "reference_absolute_error", "points": tuple(points)}
    )
    return _plot_sweep_points(
        error_analysis,
        kind="Reference comparison sweep",
        title=title,
        allow_unreliable_metadata=False,
        extra_notes=(
            "Reference error is distinct from solver convergence and "
            "numerical acceptance.",
        ),
        save_path=save_path,
        overwrite=overwrite,
        create_parent_dirs=create_parent_dirs,
    )


def plot_collection_summary(
    summary: ExperimentCollectionSummary,
    *,
    title: str | None = None,
    save_path: str | Path | None = None,
    overwrite: bool = False,
    create_parent_dirs: bool = False,
) -> VisualizationArtifact:
    """Plot family and execution-status run counts from Phase 14 summaries."""
    figure, axes = _new_figure(rows=2, height=3.2)
    figure.suptitle(title or "Experiment run summary")
    family_names = [item.category for item in summary.counts_by_family]
    family_counts = [item.count for item in summary.counts_by_family]
    status_names = [item.category for item in summary.counts_by_execution_status]
    status_counts = [item.count for item in summary.counts_by_execution_status]
    axes[0].bar(family_names or ["no runs"], family_counts or [0], color="#4c78a8")
    axes[0].set_ylabel(f"Runs (n={summary.run_count})")
    axes[0].set_title("Runs by problem family")
    axes[1].bar(status_names or ["no runs"], status_counts or [0], color="#72b7b2")
    axes[1].set_ylabel(f"Runs (n={summary.run_count})")
    axes[1].set_title("Runs by execution status")
    for axis in axes:
        axis.grid(True, axis="y", alpha=0.25)
        axis.tick_params(axis="x", rotation=15)
    return _finish(
        figure,
        source_ids=(f"{summary.record_count} experiment records",),
        notes=(
            "Execution status counts come from stored run status, not metric presence.",
        ),
        save_path=save_path,
        overwrite=overwrite,
        create_parent_dirs=create_parent_dirs,
    )


def plot_metric_distribution(
    comparison: MetricComparison,
    *,
    statistics: MetricStatistics | None = None,
    title: str | None = None,
    save_path: str | Path | None = None,
    overwrite: bool = False,
    create_parent_dirs: bool = False,
) -> VisualizationArtifact:
    """Plot observed values of a compatible metric, with exclusions called out."""
    if comparison.status == "not_comparable":
        raise ScientificValidationError(
            f"metric observations are not comparable: {comparison.explanation}"
        )
    stats = statistics or summarize_metric(comparison)
    valid = [
        item.value
        for item in comparison.observations
        if item.execution_status == ExperimentRunStatus.SUCCEEDED
        and item.metric_status == "valid"
        and item.value is not None
    ]
    figure, axes = _new_figure(rows=2, height=3.0)
    source = f"{len(comparison.observations)} metric observations"
    figure.suptitle(_with_context("Observed metric values", source, title))
    if valid:
        axes[0].hist(valid, bins="auto", color="#4c78a8", edgecolor="white")
        axes[0].set_ylabel("Observation count")
    else:
        axes[0].text(0.5, 0.5, "No valid observations", ha="center", va="center")
        axes[0].set_yticks([])
    axes[0].set_xlabel(
        f"{comparison.shared_definition or 'Selected compatible metric'} "
        f"({comparison.unit or 'unit unknown'})"
    )
    axes[0].text(
        0.01,
        0.99,
        "Observed values; not a fitted probability distribution.",
        transform=axes[0].transAxes,
        va="top",
        fontsize=8,
    )
    exclusions = (
        ("valid", stats.valid_count),
        ("failed execution", stats.execution_failed_count),
        ("missing metric", stats.missing_metric_count),
        ("non-finite metric", stats.non_finite_metric_count),
        ("invalid metric", stats.invalid_metric_count),
    )
    axes[1].bar(
        [label for label, _count in exclusions],
        [count for _label, count in exclusions],
        color=["#59a14f", "#e15759", "#f28e2b", "#b07aa1", "#9c755f"],
    )
    axes[1].set_ylabel(f"Observation count (n={stats.observation_count})")
    axes[1].set_title(
        f"Included: {stats.valid_count}; excluded: {stats.excluded_count}"
    )
    axes[1].tick_params(axis="x", rotation=15)
    axes[1].grid(True, axis="y", alpha=0.25)
    source_ids = tuple(
        dict.fromkeys(
            f"{item.experiment_id}#{item.record_index}:{item.run_index}"
            for item in comparison.observations
        )
    )
    return _finish(
        figure,
        source_ids=source_ids,
        notes=(
            "Only finite metrics from successful runs enter the observed-value "
            "histogram.",
            "Exclusion reason counts can overlap; excluded_count counts each "
            "observation once.",
            "The plot does not assert a probability distribution or predictive "
            "validity.",
        ),
        save_path=save_path,
        overwrite=overwrite,
        create_parent_dirs=create_parent_dirs,
    )


def _markdown_escape(value: object) -> str:
    return str(value).replace("|", "\\|").replace("\n", " ")


def _format_metric(value: float | None) -> str:
    return "unavailable" if value is None else f"{value:.8g}"


def build_experiment_report(
    records: Sequence[ExperimentRecord],
    *,
    selected_metrics: tuple[str, ...] = (),
    plots: tuple[VisualizationArtifact, ...] = (),
    title: str = "Newton Lab experiment report",
) -> ExperimentReport:
    """Build deterministic Markdown from records and Phase 14 analysis outputs."""
    records = tuple(records)
    if len(selected_metrics) != len(set(selected_metrics)):
        raise ScientificValidationError("selected metric IDs must be unique")
    summary = summarize_experiments(records)
    projections = project_experiment_runs(records)
    families_by_record = {
        projection.record_index: projection.problem_family for projection in projections
    }
    lines = [
        f"# {_markdown_escape(title)}",
        "",
        "Generated from existing experiment records; no simulations were run.",
        "",
        "## Collection summary",
        "",
        f"- Records: {summary.record_count}",
        f"- Runs: {summary.run_count}",
        f"- Succeeded: {summary.successful_run_count}",
        f"- Failed: {summary.failed_run_count}",
        f"- Solver success / failure / unavailable: {summary.solver_success_count} / "
        f"{summary.solver_failure_count} / {summary.solver_success_unavailable_count}",
        "- Numerical acceptance / rejection / unavailable: "
        f"{summary.numerical_acceptance_count} / {summary.numerical_rejection_count} / "
        f"{summary.numerical_acceptance_unavailable_count}",
        "",
        "## Source records",
        "",
    ]
    for record_index, record in enumerate(records):
        spec = record.specification
        family = families_by_record.get(record_index, "unknown")
        lines.extend(
            [
                f"### Record {record_index}: `{_markdown_escape(spec.experiment_id)}`",
                "",
                f"- Model: `{_markdown_escape(spec.model_id)}`",
                f"- Record status: `{record.status.value}`",
                f"- Problem family: `{_markdown_escape(family)}`",
                "- Parameters and units are recorded for each run below.",
            ]
        )
        if spec.equation_record_ids:
            lines.append(
                "- Equation records: "
                + ", ".join(
                    f"`{_markdown_escape(value)}`" for value in spec.equation_record_ids
                )
            )
        if spec.source_reference_ids:
            lines.append(
                "- Source references: "
                + ", ".join(
                    f"`{_markdown_escape(value)}`"
                    for value in spec.source_reference_ids
                )
            )
        lines.append("")
        if spec.assumptions:
            lines.append("Assumptions:")
            lines.extend(f"- {_markdown_escape(item)}" for item in spec.assumptions)
            lines.append("")
        lines.extend(
            [
                "| Run | Outcome | Solver | Solver success | Numerical acceptance | "
                "Parameters | Failure |",
                "| ---: | --- | --- | --- | --- | --- | --- |",
            ]
        )
        for run in record.runs:
            solver_success = (
                run.solver_success if run.solver_success is not None else "unavailable"
            )
            numerical_acceptance = (
                run.numerical_acceptance
                if run.numerical_acceptance is not None
                else "unavailable"
            )
            params = "; ".join(
                f"{_markdown_escape(item.name)}={item.value:.8g} "
                f"{_markdown_escape(item.unit)}"
                for item in run.parameters
            )
            failure = (
                ""
                if run.status == ExperimentRunStatus.SUCCEEDED
                else f"{run.failure_kind}: {run.error_type}: {run.message}"
            )
            lines.append(
                f"| {run.run_index} | {run.status.value} | "
                f"{_markdown_escape(run.solver_id)} / "
                f"{_markdown_escape(run.solver_method)} | "
                f"{solver_success} | {numerical_acceptance} | "
                f"{_markdown_escape(params)} | {_markdown_escape(failure)} |"
            )
            provenance = run.reproducibility
            provenance_summary = (
                f"configuration={provenance.configuration_version}; "
                f"Python={provenance.python_version}; "
                f"NumPy={provenance.numpy_version}; "
                f"SciPy={provenance.scipy_version}; seed={provenance.random_seed}; "
                f"rtol={provenance.relative_tolerance}; "
                f"atol={provenance.absolute_tolerance}; "
                f"residual tolerance={provenance.residual_tolerance}; "
                f"BVP tolerance={provenance.bvp_solver_tolerance}"
            )
            lines.append(
                f"\nRun source `"
                f"{_markdown_escape(spec.experiment_id)}#{run.run_index}`; "
                f"solver metadata: `{_markdown_escape(provenance_summary)}`."
            )
            if provenance.equation_record_ids:
                lines.append(
                    " Equation records: "
                    + ", ".join(
                        f"`{_markdown_escape(value)}`"
                        for value in provenance.equation_record_ids
                    )
                    + "."
                )
            if provenance.source_reference_ids:
                lines.append(
                    " Source references: "
                    + ", ".join(
                        f"`{_markdown_escape(value)}`"
                        for value in provenance.source_reference_ids
                    )
                    + "."
                )
            if provenance.assumptions:
                lines.extend(
                    ["", "Run assumptions:"]
                    + [f"- {_markdown_escape(item)}" for item in provenance.assumptions]
                )
            if run.reference_comparison is not None:
                ref = run.reference_comparison
                lines.append(
                    f"\nReference `{_markdown_escape(ref.metric_id)}`: "
                    f"{ref.status}; absolute error "
                    f"{_format_metric(ref.absolute_error)} "
                    f"{_markdown_escape(ref.unit)}. {_markdown_escape(ref.message)}"
                )
            if run.bvp_reference_comparison is not None:
                ref_bvp = run.bvp_reference_comparison
                lines.append(
                    f"\nBVP reference `{_markdown_escape(ref_bvp.state_name)}`: "
                    f"{ref_bvp.status}; maximum profile error "
                    f"{_format_metric(ref_bvp.maximum_absolute_error)} "
                    f"{_markdown_escape(ref_bvp.unit)}. "
                    f"{_markdown_escape(ref_bvp.message)}"
                )
        lines.append("")
    if selected_metrics:
        lines.extend(["## Selected metrics", ""])
        for metric_id in selected_metrics:
            observations = extract_metric_observations(records, metric_id)
            lines.extend(
                [
                    f"### `{_markdown_escape(metric_id)}`",
                    "",
                    "| Record | Run | Outcome | Metric status | Value | Unit | "
                    "Definition |",
                    "| ---: | ---: | --- | --- | ---: | --- | --- |",
                ]
            )
            for item in observations:
                lines.append(
                    f"| {item.record_index} "
                    f"(`{_markdown_escape(item.experiment_id)}`) | "
                    f"{item.run_index} | {item.execution_status.value} | "
                    f"{item.metric_status} | {_format_metric(item.value)} | "
                    f"{_markdown_escape(item.unit or 'unit unavailable')} | "
                    f"{_markdown_escape(item.definition or 'definition unavailable')} |"
                )
            lines.append("")
    missing_provenance = sum(
        projection.reproducibility is None for projection in projections
    )
    lines.extend(
        [
            "## Provenance and limitations",
            "",
            f"- Runs without reproducibility metadata: {missing_provenance} "
            f"of {len(projections)}.",
            "- Run parameters, solver status, and family-specific diagnostics "
            "remain in the source records.",
            "- Numerical acceptance does not establish physical validity, "
            "uniqueness, or model accuracy.",
            "- Descriptive summaries and visual appearance do not establish "
            "causality or predictive power.",
            "- A smaller metric value is not inherently better; interpret "
            "each quantity using its definition and assumptions.",
            "",
        ]
    )
    plot_references: list[ReportPlotReference] = []
    for index, plot in enumerate(plots):
        if plot.saved_path is None:
            raise ScientificValidationError(
                "reports can include only saved plot artifacts"
            )
        resolved = plot.saved_path.resolve(strict=True)
        if not resolved.is_file():
            raise FileNotFoundError(f"report plot is not a file: {resolved}")
        label = (
            f"Plot {index + 1}: {', '.join(plot.source_ids) or 'unspecified source'}"
        )
        plot_references.append(ReportPlotReference(path=resolved, label=label))
    return ExperimentReport(
        body="\n".join(lines),
        source_experiment_ids=tuple(
            record.specification.experiment_id for record in records
        ),
        summary=summary,
        plots=tuple(plot_references),
    )


def write_experiment_report(
    report: ExperimentReport,
    path: str | Path,
    *,
    overwrite: bool = False,
    create_parent_dirs: bool = False,
) -> Path:
    """Write Markdown without overwriting or creating folders unless requested."""
    target = Path(path)
    if target.suffix.lower() != ".md":
        raise ScientificValidationError("experiment report path must end in .md")
    if target.exists() and not overwrite:
        raise FileExistsError(f"report output already exists: {target}")
    if target.exists() and target.is_dir():
        raise IsADirectoryError(f"report output path is a directory: {target}")
    if not target.parent.exists():
        if create_parent_dirs:
            target.parent.mkdir(parents=True, exist_ok=True)
        else:
            raise FileNotFoundError(
                f"report parent directory does not exist: {target.parent}"
            )
    text = report.render_markdown(relative_to=target.parent.resolve())
    try:
        target.write_text(text, encoding="utf-8", newline="\n")
    except OSError as exc:
        raise VisualizationError(
            f"could not write experiment report to {target}: {exc}"
        ) from exc
    return target

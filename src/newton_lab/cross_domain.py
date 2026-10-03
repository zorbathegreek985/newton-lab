"""A bounded oscillator-to-adjustment mathematical case study.

The target model is an abstract setpoint adjustment equation. It is not a model
of prices, returns, portfolios, or a complete financial control process.
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import Literal

import numpy as np
from numpy.typing import NDArray
from pydantic import Field, field_validator, model_validator

from newton_lab.discovery import (
    ApplicationHypothesis,
    DiscoveryEvidenceStatus,
    DiscoveryResult,
    ObservationAvailability,
    run_discovery,
)
from newton_lab.dynamics import simulate_damped_oscillator
from newton_lab.exceptions import IntegrationError, ScientificValidationError
from newton_lab.experiments import (
    ExperimentRecord,
    make_oscillator_experiment,
    run_experiment,
)
from newton_lab.simulation import ODESolverConfiguration, SimulationModel


class AdjustmentMetrics(SimulationModel):
    """Defined finite-horizon metrics for the abstract adjustment response."""

    peak_absolute_rate_per_s: float
    maximum_overshoot: float
    sampled_settling_time_s: float | None
    target_crossing_count: int = Field(ge=0)
    maximum_analytical_absolute_error: float

    @field_validator(
        "peak_absolute_rate_per_s",
        "maximum_overshoot",
        "sampled_settling_time_s",
        "maximum_analytical_absolute_error",
    )
    @classmethod
    def finite_metrics(cls, value: float | None) -> float | None:
        if value is not None and (not np.isfinite(value) or value < 0.0):
            raise ValueError("adjustment metrics must be finite and non-negative")
        return value


class AdjustmentSweepPoint(SimulationModel):
    """One ordered damping-ratio setting and its sampled target response."""

    damping_ratio: float = Field(ge=0.0)
    status: Literal["succeeded", "failed"]
    time_s: tuple[float, ...] = ()
    adjustment_state: tuple[float, ...] = ()
    adjustment_rate_per_s: tuple[float, ...] = ()
    metrics: AdjustmentMetrics | None = None
    solver_method: str
    failure_type: str | None = None
    failure_message: str | None = None

    @model_validator(mode="after")
    def validate_point_outcome(self) -> AdjustmentSweepPoint:
        if self.status == "succeeded":
            if self.metrics is None or self.failure_type is not None:
                raise ValueError(
                    "successful sweep points require metrics and no failure"
                )
            if not (
                len(self.time_s)
                == len(self.adjustment_state)
                == len(self.adjustment_rate_per_s)
                >= 2
            ):
                raise ValueError("successful target trajectories must be aligned")
            if (
                not np.isfinite(self.time_s).all()
                or not np.isfinite(self.adjustment_state).all()
                or not np.isfinite(self.adjustment_rate_per_s).all()
                or not np.all(np.diff(self.time_s) > 0.0)
            ):
                raise ValueError("successful target samples must be finite and ordered")
        elif (
            self.metrics is not None
            or self.failure_type is None
            or self.failure_message is None
        ):
            raise ValueError(
                "failed sweep points require failure details and no metrics"
            )
        return self


class DampingTransferCaseStudy(SimulationModel):
    """Reproducible comparison linked to the unmodified Phase 16 hypothesis."""

    case_study_id: str
    source_experiment_id: str
    source_experiment: ExperimentRecord
    source_hypothesis: ApplicationHypothesis
    source_peak_rate_endpoint_change: float
    source_metric_id: str
    source_metric_unit: str
    source_metric_parameter_values: tuple[float, ...] = Field(min_length=1)
    source_metric_values: tuple[float, ...] = Field(min_length=1)
    source_damping_ratio_values: tuple[float, ...] = Field(min_length=1)
    source_normalized_peak_rate_values_per_s: tuple[float, ...] = Field(min_length=1)
    source_target_peak_rate_absolute_differences_per_s: tuple[float | None, ...] = (
        Field(min_length=1)
    )
    amplitude_scale_m: float
    source_observation_ids: tuple[str, ...] = Field(min_length=1)
    target_equation: str
    target_damping_ratios: tuple[float, ...] = Field(min_length=3)
    baseline_damping_ratio: float
    time_constant_s: float
    initial_error: float
    initial_rate_per_s: float
    duration_s: float
    output_point_count: int = Field(ge=2)
    solver_configuration: ODESolverConfiguration
    sweep_points: tuple[AdjustmentSweepPoint, ...] = Field(min_length=1)
    peak_rate_trend_on_tested_grid: str
    overshoot_trend_on_tested_grid: str
    settling_time_trend_on_tested_grid: str
    mathematical_correspondence: Literal[
        "equivalent_under_declared_constant_target_transformation"
    ]
    target_domain_interpretation: Literal["structural_analogy_only"]
    real_world_application_status: Literal["unvalidated"]
    assumptions: tuple[str, ...] = Field(min_length=1)
    comparison_results: tuple[str, ...] = Field(min_length=1)
    limitations: tuple[str, ...] = Field(min_length=1)
    report_markdown: str

    @model_validator(mode="after")
    def validate_linked_evidence(self) -> DampingTransferCaseStudy:
        if self.source_hypothesis.validated_application:
            raise ValueError(
                "the source application hypothesis must remain unvalidated"
            )
        if self.source_hypothesis.evidence_status != (
            DiscoveryEvidenceStatus.PROPOSED_APPLICATION_HYPOTHESIS
        ):
            raise ValueError("the source must be the existing proposed hypothesis")
        if tuple(point.damping_ratio for point in self.sweep_points) != (
            self.target_damping_ratios
        ):
            raise ValueError("sweep points must preserve the declared damping order")
        if (
            self.source_experiment.specification.experiment_id
            != self.source_experiment_id
        ):
            raise ValueError(
                "embedded source record must match its source experiment ID"
            )
        if len(self.source_metric_values) != len(self.source_metric_parameter_values):
            raise ValueError("source metric values and parameter values must align")
        if not (
            len(self.source_metric_values)
            == len(self.source_damping_ratio_values)
            == len(self.source_normalized_peak_rate_values_per_s)
            == len(self.source_target_peak_rate_absolute_differences_per_s)
        ):
            raise ValueError("source and mapped target metric provenance must align")
        if (
            not np.isfinite(self.source_metric_values).all()
            or not np.isfinite(self.source_damping_ratio_values).all()
        ):
            raise ValueError("source sweep values and metrics must be finite")
        if self.amplitude_scale_m <= 0.0 or not np.isfinite(self.amplitude_scale_m):
            raise ValueError("amplitude_scale_m must be finite and positive")
        if self.baseline_damping_ratio not in self.target_damping_ratios:
            raise ValueError("baseline damping ratio must occur in the sweep")
        if self.real_world_application_status != "unvalidated":
            raise ValueError(
                "toy-model results cannot validate a real-world application"
            )
        return self


DEFAULT_DAMPING_RATIOS = (0.0, 0.025, 0.05, 0.7, 1.0, 2.0)
DEFAULT_BASELINE_DAMPING_RATIO = 0.025
DEFAULT_TIME_CONSTANT_S = 0.5
DEFAULT_DURATION_S = 60.0
DEFAULT_OUTPUT_POINT_COUNT = 6001
_DEFAULT_SOLVER_CONFIGURATION = ODESolverConfiguration(
    method="DOP853", relative_tolerance=1e-10, absolute_tolerance=1e-12
)


def analytical_step_response(
    time_s: Sequence[float] | NDArray[np.float64],
    *,
    damping_ratio: float,
    time_constant_s: float,
    initial_error: float = 1.0,
) -> tuple[float, ...]:
    """Return the exact constant-setpoint response for zero initial rate.

    The returned values are dimensionless adjustment state samples. This
    closed form is for the residual error after a unit setpoint change, with
    constant parameters and no later disturbance or saturation.
    """
    times = np.asarray(time_s, dtype=np.float64)
    damping = _finite_nonnegative(damping_ratio, "damping_ratio")
    tau = _finite_positive(time_constant_s, "time_constant_s")
    initial = _finite(initial_error, "initial_error")
    if times.ndim != 1 or len(times) < 1 or not np.isfinite(times).all():
        raise ScientificValidationError(
            "time_s must be a finite one-dimensional sequence"
        )
    if np.any(times < 0.0) or np.any(np.diff(times) < 0.0):
        raise ScientificValidationError("time_s must be non-negative and ordered")

    omega = 1.0 / tau
    error_0 = initial
    if damping < 1.0:
        damped_omega = omega * np.sqrt(1.0 - damping**2)
        error = (
            error_0
            * np.exp(-damping * omega * times)
            * (
                np.cos(damped_omega * times)
                + damping / np.sqrt(1.0 - damping**2) * np.sin(damped_omega * times)
            )
        )
    elif damping == 1.0:
        error = error_0 * (1.0 + omega * times) * np.exp(-omega * times)
    else:
        root = np.sqrt(damping**2 - 1.0)
        fast_root = -omega * (damping + root)
        slow_root = -omega * (damping - root)
        fast_coefficient = -slow_root / (fast_root - slow_root)
        slow_coefficient = fast_root / (fast_root - slow_root)
        error = error_0 * (
            fast_coefficient * np.exp(fast_root * times)
            + slow_coefficient * np.exp(slow_root * times)
        )
    return tuple(float(value) for value in error)


def _finite(value: float, name: str) -> float:
    if isinstance(value, bool) or not np.isfinite(value):
        raise ScientificValidationError(f"{name} must be finite")
    return float(value)


def _finite_positive(value: float, name: str) -> float:
    result = _finite(value, name)
    if result <= 0.0:
        raise ScientificValidationError(f"{name} must be positive")
    return result


def _finite_nonnegative(value: float, name: str) -> float:
    result = _finite(value, name)
    if result < 0.0:
        raise ScientificValidationError(f"{name} must be non-negative")
    return result


def _observed_metrics(
    *,
    time_s: np.ndarray,
    state: np.ndarray,
    rate: np.ndarray,
    initial_error: float,
    analytic: np.ndarray,
) -> AdjustmentMetrics:
    """Compute explicitly sampled metrics; no continuous extrema are claimed."""
    error = state
    crossing_tolerance = max(1e-10, 1e-8 * abs(initial_error))
    signs = np.where(
        error > crossing_tolerance,
        1,
        np.where(error < -crossing_tolerance, -1, 0),
    )
    nonzero_signs = signs[signs != 0]
    crosses = int(np.count_nonzero(nonzero_signs[:-1] != nonzero_signs[1:]))
    band = 0.02 * abs(initial_error)
    if band == 0.0:
        settling_time = None
    else:
        inside = np.abs(error) <= band
        remains_inside = np.logical_and.accumulate(inside[::-1])[::-1]
        in_band_indices = np.flatnonzero(remains_inside)
        settling_time = (
            None
            if len(in_band_indices) == 0
            else float(time_s[int(in_band_indices[0])])
        )
    return AdjustmentMetrics(
        peak_absolute_rate_per_s=float(np.max(np.abs(rate))),
        maximum_overshoot=float(max(0.0, float(-np.min(error)))),
        sampled_settling_time_s=settling_time,
        target_crossing_count=crosses,
        maximum_analytical_absolute_error=float(np.max(np.abs(state - analytic))),
    )


def _grid_trend(values: Sequence[float | None]) -> str:
    """Describe adjacent finite differences on exactly the supplied sweep grid."""
    if any(value is None for value in values):
        return "unavailable_incomplete_sweep"
    numeric = np.asarray(values, dtype=np.float64)
    if len(numeric) < 2:
        return "unavailable_insufficient_values"
    changes = np.diff(numeric)
    tolerance = 1e-9 * max(1.0, float(np.max(np.abs(numeric))))
    if np.all(changes <= tolerance):
        return "nonincreasing_on_tested_grid"
    if np.all(changes >= -tolerance):
        return "nondecreasing_on_tested_grid"
    return "nonmonotonic_on_tested_grid"


def _render_report(case: DampingTransferCaseStudy) -> str:
    lines = [
        "# Damping transfer case study",
        "",
        f"Source experiment: `{case.source_experiment_id}`",
        f"Source hypothesis: `{case.source_hypothesis.hypothesis_id}` "
        "(proposed; unvalidated)",
        "",
        "## Target model",
        "",
        f"`{case.target_equation}`",
        "",
        "The state q(t) is dimensionless residual error from a setpoint, and "
        "q'(t) has units s^-1 by convention. A unit target step establishes "
        "q(0)=1 and q'(0)=0; the equilibrium is q=0, q'=0. The damping ratio "
        "zeta is dimensionless and tau is a characteristic time in seconds. "
        "There are no later inputs, noise, constraints, or saturation.",
        "",
        "## Explicit correspondence",
        "",
        "After the initial step, q satisfies "
        "q'' + (2 zeta/tau)q' + q/tau^2 = 0. The physical equation is "
        "m x'' + c x' + k x = 0. Set x=Aq, "
        "then q'' + (c/m)q' + (k/m)q = 0. Choosing "
        "tau = sqrt(m/k), and zeta = c/(2 sqrt(m k)) makes both coefficients "
        "agree. For the source example A=0.1 m per state unit, "
        "m=1 kg, k=4 N/m, tau=0.5 s, and the physical damping sweep "
        "c=(0, 0.1, 0.2) kg/s maps to zeta=(0, 0.025, 0.05). Thus the "
        "initial conditions also agree (x(0)=0.1 m maps from q(0)=1). "
        "Conversely c=2 zeta m/tau and k=m/tau^2. This exact equation "
        "transformation does not identify a financial quantity with physical "
        "displacement or establish a market model.",
        "",
        "## Controlled sweep",
        "",
        "Metrics use returned samples: peak |q'| (s^-1); overshoot beyond "
        "equilibrium max(0, -min(q)) (dimensionless); zero crossings (count); "
        "sampled 2% settling time (s), defined by the first sample after which "
        "all remaining samples satisfy |q| <= 0.02; and maximum error against "
        "the closed-form response (dimensionless).",
        "",
    ]
    lines.extend(
        [
            "",
            "Source physical sweep (the endpoint change used by the original "
            "hypothesis is last minus first):",
            "",
            "| Physical damping c (kg/s) | Mapped zeta | Peak |x'| (m/s) | "
            "Normalized source peak |q'| (s^-1) | Target-source absolute "
            "difference (s^-1) |",
            "| ---: | ---: | ---: | ---: | ---: |",
        ]
    )
    lines.extend(
        f"| {damping:g} | {ratio:g} | {rate:.7g} | {normalized:.7g} | "
        f"{('unavailable' if difference is None else f'{difference:.3g}')} |"
        for damping, ratio, rate, normalized, difference in zip(
            case.source_metric_parameter_values,
            case.source_damping_ratio_values,
            case.source_metric_values,
            case.source_normalized_peak_rate_values_per_s,
            case.source_target_peak_rate_absolute_differences_per_s,
            strict=True,
        )
    )
    lines.extend(
        [
            "",
            f"Original Phase 16 hypothesis: {case.source_hypothesis.statement}",
            "Its source measurements are sampled physical velocity in m/s; "
            "normalization divides by A=0.1 m to express rate in q-units/s. "
            "Target-source differences are finite-grid discrepancies, not "
            "domain performance differences.",
            "",
            "Target sweep:",
            "",
            "| Damping ratio | Peak |q'| (s^-1) | Overshoot | Crossings | "
            "Settling time (s) | Max analytical error | Status |",
            "| ---: | ---: | ---: | ---: | ---: | ---: | --- |",
        ]
    )
    for point in case.sweep_points:
        if point.metrics is None:
            lines.append(
                f"| {point.damping_ratio:g} | unavailable | unavailable | "
                f"unavailable | unavailable | unavailable | failed: "
                f"{point.failure_type} |"
            )
            continue
        metrics = point.metrics
        settling = (
            "not reached"
            if metrics.sampled_settling_time_s is None
            else f"{metrics.sampled_settling_time_s:.6g}"
        )
        lines.append(
            f"| {point.damping_ratio:g} | {metrics.peak_absolute_rate_per_s:.7g} | "
            f"{metrics.maximum_overshoot:.7g} | {metrics.target_crossing_count} | "
            f"{settling} | {metrics.maximum_analytical_absolute_error:.3g} | "
            f"{point.status} |"
        )
    lines.extend(
        [
            "",
            "On this tested grid, peak rate was "
            f"**{case.peak_rate_trend_on_tested_grid}** "
            f"and overshoot was **{case.overshoot_trend_on_tested_grid}**; sampled "
            f"settling time was **{case.settling_time_trend_on_tested_grid}**. "
            "These are finite-grid descriptions, not general monotonicity claims.",
            "",
            "## Evidence classification",
            "",
            "- Mathematical result: correspondence established for the stated "
            "constant-target toy equation and parameter transformation.",
            "- Domain interpretation: structural analogy only; this model has no "
            "empirical financial semantics.",
            "- Real-world application: **unvalidated**. No price, return, portfolio, "
            "or transaction-cost data were used.",
            "",
            "## What challenges or limits the hypothesis",
            "",
            "Increasing damping across the tested grid reduces overshoot and "
            "peak rate in this step-response experiment, but the overdamped "
            "setting can take longer to settle than the critical setting. Thus "
            "lower oscillation is not equivalent to faster adjustment or better "
            "overall performance. The zero-damping point fails to settle in the "
            "observed interval, and damping ratios at or above one show no "
            "repeated crossings in this finite run.",
            "",
            "## Limitations and next evidence",
            "",
        ]
    )
    lines.extend(f"- {item}" for item in case.limitations)
    lines.extend(
        [
            "- A real-world study would need suitable empirical data, "
            "prespecified defensible baselines, transaction costs and market "
            "impact where relevant, risk measures, chronological out-of-sample "
            "evaluation, and explicit overfitting controls.",
            "- Solver success and agreement with this closed form validate only "
            "the numerical implementation of the toy equation, not its domain "
            "interpretation or usefulness.",
            "",
        ]
    )
    return "\n".join(lines)


def run_damping_transfer_case_study(
    *,
    source_record: ExperimentRecord | None = None,
    damping_ratios: tuple[float, ...] = DEFAULT_DAMPING_RATIOS,
    baseline_damping_ratio: float = DEFAULT_BASELINE_DAMPING_RATIO,
    time_constant_s: float = DEFAULT_TIME_CONSTANT_S,
    duration_s: float = DEFAULT_DURATION_S,
    output_point_count: int = DEFAULT_OUTPUT_POINT_COUNT,
    solver_configuration: ODESolverConfiguration = _DEFAULT_SOLVER_CONFIGURATION,
) -> DampingTransferCaseStudy:
    """Reproduce the oscillator hypothesis and test its bounded toy-model form.

    If no physical record is supplied, a deterministic Phase 16 oscillator
    sweep is generated. A supplied record is read and passed to discovery but
    never modified. Target-domain sweep order is retained exactly as given.
    """
    requested_tau = _finite_positive(time_constant_s, "time_constant_s")
    duration = _finite_positive(duration_s, "duration_s")
    baseline = _finite_nonnegative(baseline_damping_ratio, "baseline_damping_ratio")
    if (
        isinstance(output_point_count, bool)
        or not isinstance(output_point_count, int)
        or output_point_count < 2
    ):
        raise ScientificValidationError(
            "output_point_count must be an integer of at least 2"
        )
    ratios = tuple(
        _finite_nonnegative(value, "damping_ratio") for value in damping_ratios
    )
    if (
        len(ratios) < 3
        or baseline not in ratios
        or not any(value < baseline for value in ratios)
        or not any(value > baseline for value in ratios)
    ):
        raise ScientificValidationError(
            "damping sweep must include at least three settings, the baseline, "
            "and settings both below and above the baseline"
        )
    record = source_record or run_experiment(
        make_oscillator_experiment(
            experiment_id="phase18_source_oscillator_damping_sweep",
            mass_kg=1.0,
            damping_coefficient_kg_per_s=0.2,
            stiffness_n_per_m=4.0,
            initial_displacement_m=0.1,
            initial_velocity_m_per_s=0.0,
            duration_s=10.0,
            sweep_values=(0.0, 0.1, 0.2),
            metrics=("maximum_absolute_velocity_m_per_s",),
        )
    )
    source_before = record.model_dump(mode="python")
    spec = record.specification
    if spec.model_id != "damped_harmonic_oscillator":
        raise ScientificValidationError(
            "source record must use the damped oscillator model"
        )
    if (
        spec.sweep is None
        or spec.sweep.parameter_name != "damping_coefficient_kg_per_s"
    ):
        raise ScientificValidationError(
            "source record must sweep damping_coefficient_kg_per_s"
        )
    if "maximum_absolute_velocity_m_per_s" not in spec.metrics:
        raise ScientificValidationError(
            "source record must request maximum_absolute_velocity_m_per_s"
        )

    base_parameters = {item.name: item for item in spec.parameters}
    expected_parameters = {
        "mass_kg": (1.0, "kg"),
        "stiffness_n_per_m": (4.0, "N/m"),
        "initial_displacement_m": (0.1, "m"),
        "initial_velocity_m_per_s": (0.0, "m/s"),
    }
    for name, (expected_value, expected_unit) in expected_parameters.items():
        actual = base_parameters.get(name)
        if (
            actual is None
            or actual.value != expected_value
            or actual.unit != expected_unit
        ):
            raise ScientificValidationError(
                f"source parameter {name!r} must be {expected_value:g} {expected_unit}"
            )
    damping_parameter = base_parameters.get("damping_coefficient_kg_per_s")
    if damping_parameter is None or damping_parameter.unit != "kg/s":
        raise ScientificValidationError("source damping parameter must use kg/s")
    if "No external driving force." not in spec.assumptions:
        raise ScientificValidationError(
            "source must declare the unforced oscillator assumption"
        )
    if not {
        "Linear spring with constant stiffness.",
        "Linear viscous damping with constant coefficient.",
    }.issubset(spec.assumptions):
        raise ScientificValidationError(
            "source does not declare the linear constant-parameter assumptions"
        )

    mass = base_parameters["mass_kg"].value
    stiffness = base_parameters["stiffness_n_per_m"].value
    amplitude_scale_m = base_parameters["initial_displacement_m"].value
    tau = float(np.sqrt(mass / stiffness))
    if not np.isclose(requested_tau, tau, rtol=1e-12, atol=0.0):
        raise ScientificValidationError(
            f"time_constant_s must equal sqrt(m/k)={tau:g} s for this source"
        )
    if duration < 10.0 * tau:
        raise ScientificValidationError(
            "duration_s must cover at least ten time constants for this comparison"
        )
    source_damping_values = tuple(spec.sweep.values)
    if len(record.runs) != len(source_damping_values):
        raise ScientificValidationError(
            "source run count does not match its declared sweep"
        )
    if any(run.status.value != "succeeded" for run in record.runs):
        raise ScientificValidationError(
            "source failures cannot produce the existing complete Phase 16 hypothesis"
        )
    source_zeta_values = tuple(
        value / (2.0 * np.sqrt(mass * stiffness)) for value in source_damping_values
    )
    if not set(source_zeta_values).issubset(ratios):
        raise ScientificValidationError(
            "target damping ratios must include every mapped source sweep value"
        )
    source_metric_values: list[float] = []
    for index, (run, damping_value) in enumerate(
        zip(record.runs, source_damping_values, strict=True)
    ):
        if (
            run.run_index != index
            or run.swept_parameter_name != "damping_coefficient_kg_per_s"
        ):
            raise ScientificValidationError(
                "source run order/provenance does not match damping sweep"
            )
        if run.reproducibility.model_id != spec.model_id:
            raise ScientificValidationError(
                f"source run {index} has incompatible model provenance"
            )
        actual_damping = next(
            (
                item
                for item in run.parameters
                if item.name == "damping_coefficient_kg_per_s"
            ),
            None,
        )
        if (
            actual_damping is None
            or actual_damping.value != damping_value
            or actual_damping.unit != "kg/s"
        ):
            raise ScientificValidationError(
                f"source run {index} damping parameter provenance is inconsistent"
            )
        if run.swept_parameter_value != damping_value:
            raise ScientificValidationError(
                "source run damping value does not match sweep order"
            )
        run_parameters = {item.name: item for item in run.parameters}
        for name, (expected_value, expected_unit) in expected_parameters.items():
            item = run_parameters.get(name)
            if (
                item is None
                or item.value != expected_value
                or item.unit != expected_unit
            ):
                raise ScientificValidationError(
                    f"source run {index} has incompatible {name}"
                )
        metric = next(
            (
                item
                for item in run.metrics
                if item.metric_id == "maximum_absolute_velocity_m_per_s"
            ),
            None,
        )
        if metric is None or metric.unit != "m/s":
            raise ScientificValidationError(
                f"source run {index} is missing its finite maximum velocity "
                "metric in m/s"
            )
        source_metric_values.append(metric.value)

    discovery: DiscoveryResult = run_discovery((record,))
    matching_hypotheses = tuple(
        hypothesis
        for hypothesis in discovery.hypotheses
        if hypothesis.source_experiment_id == record.specification.experiment_id
        and hypothesis.target_domain == "quantitative_finance_risk_control_research"
        and hypothesis.evidence_status
        == DiscoveryEvidenceStatus.PROPOSED_APPLICATION_HYPOTHESIS
    )
    if len(matching_hypotheses) != 1:
        raise ScientificValidationError(
            "source record must produce exactly one existing oscillator "
            "risk-control hypothesis"
        )
    source_hypothesis = matching_hypotheses[0]
    if source_hypothesis.source_equation_id != "damped_harmonic_oscillator":
        raise ScientificValidationError(
            "source hypothesis is not linked to the oscillator equation"
        )
    source_change = next(
        (
            item
            for item in discovery.observations
            if item.observation_id in source_hypothesis.source_behavior_observation_ids
            and item.descriptor == "sweep_endpoint_change"
            and item.availability == ObservationAvailability.AVAILABLE
        ),
        None,
    )
    if source_change is None or source_change.value is None:
        raise ScientificValidationError(
            "source hypothesis has no available sweep evidence"
        )
    if not np.isclose(
        source_change.value,
        source_metric_values[-1] - source_metric_values[0],
        rtol=1e-12,
        atol=1e-15,
    ):
        raise ScientificValidationError(
            "source hypothesis endpoint change does not match its recorded run metrics"
        )
    if not np.isclose(
        source_change.value,
        source_metric_values[-1] - source_metric_values[0],
        rtol=1e-12,
        atol=1e-15,
    ):
        raise ScientificValidationError(
            "source hypothesis endpoint change does not match its recorded run metrics"
        )

    if record.model_dump(mode="python") != source_before:
        raise RuntimeError(
            "discovery unexpectedly changed the source experiment record"
        )

    target_equation = "q'' + (2*zeta/tau)q' + q/tau^2 = 0, q(0)=1, q'(0)=0"
    target_initial_error = 1.0
    target_initial_rate = 0.0
    output_times = np.linspace(0.0, duration, output_point_count, dtype=np.float64)
    points: list[AdjustmentSweepPoint] = []
    for damping_ratio in ratios:
        try:
            oscillator = simulate_damped_oscillator(
                mass_kg=mass,
                damping_coefficient_kg_per_s=2.0
                * damping_ratio
                * np.sqrt(mass * stiffness),
                stiffness_n_per_m=stiffness,
                initial_displacement_m=amplitude_scale_m * target_initial_error,
                initial_velocity_m_per_s=amplitude_scale_m * target_initial_rate,
                duration_s=duration,
                num_points=output_point_count,
                method=solver_configuration.method,
                relative_tolerance=solver_configuration.relative_tolerance,
                absolute_tolerance=solver_configuration.absolute_tolerance,
                maximum_step_s=solver_configuration.maximum_step,
            )
        except IntegrationError as exc:
            points.append(
                AdjustmentSweepPoint(
                    damping_ratio=damping_ratio,
                    status="failed",
                    solver_method=solver_configuration.method,
                    failure_type=type(exc).__name__,
                    failure_message=str(exc),
                )
            )
            continue
        analytic_state = np.asarray(
            analytical_step_response(
                output_times,
                damping_ratio=damping_ratio,
                time_constant_s=tau,
                initial_error=target_initial_error,
            ),
            dtype=np.float64,
        )
        adjustment_state = oscillator.displacement_m / amplitude_scale_m
        adjustment_rate = oscillator.velocity_m_per_s / amplitude_scale_m
        metrics = _observed_metrics(
            time_s=oscillator.time_s,
            state=adjustment_state,
            rate=adjustment_rate,
            initial_error=target_initial_error,
            analytic=analytic_state,
        )
        points.append(
            AdjustmentSweepPoint(
                damping_ratio=damping_ratio,
                status="succeeded",
                time_s=tuple(float(value) for value in oscillator.time_s),
                adjustment_state=tuple(float(value) for value in adjustment_state),
                adjustment_rate_per_s=tuple(float(value) for value in adjustment_rate),
                metrics=metrics,
                solver_method=solver_configuration.method,
            )
        )

    points_tuple = tuple(points)
    successful = tuple(point for point in points_tuple if point.metrics is not None)
    if not successful:
        raise IntegrationError("all target adjustment sweep points failed")
    source_normalized_rates = tuple(
        value / amplitude_scale_m for value in source_metric_values
    )
    source_target_differences: list[float | None] = []
    for zeta, normalized_rate in zip(
        source_zeta_values, source_normalized_rates, strict=True
    ):
        target_point = next(
            (point for point in points_tuple if point.damping_ratio == zeta), None
        )
        source_target_differences.append(
            None
            if target_point is None or target_point.metrics is None
            else abs(target_point.metrics.peak_absolute_rate_per_s - normalized_rate)
        )
    peak_trend = _grid_trend(
        tuple(
            None if point.metrics is None else point.metrics.peak_absolute_rate_per_s
            for point in points_tuple
        )
    )
    overshoot_trend = _grid_trend(
        tuple(
            None if point.metrics is None else point.metrics.maximum_overshoot
            for point in points_tuple
        )
    )
    settling_values = tuple(
        None if point.metrics is None else point.metrics.sampled_settling_time_s
        for point in points_tuple
    )
    settling_trend = (
        "unavailable_censored_no_settling_within_horizon"
        if any(
            point.metrics is not None and point.metrics.sampled_settling_time_s is None
            for point in points_tuple
        )
        else _grid_trend(settling_values)
    )
    if record.model_dump(mode="python") != source_before:
        raise RuntimeError(
            "case study unexpectedly changed the source experiment record"
        )

    succeeded_text = (
        "All target-model sweep points integrated successfully."
        if len(successful) == len(points_tuple)
        else f"{len(points_tuple) - len(successful)} target-model sweep points failed."
    )
    limitations = (
        "The target state is abstract and dimensionless; no financial quantity is "
        "identified with physical displacement.",
        "The setpoint is stepped initially and then held fixed. The model has no "
        "later commands, noise, delays, constraints, or saturation.",
        "The mathematical correspondence follows from a chosen equation "
        "transformation; it is not evidence that the equation describes finance.",
        "Metrics use a finite uniform grid and do not prove continuous extrema, "
        "asymptotic stability, or trends outside tested damping ratios.",
        "The source endpoint change is a single last-minus-first comparison; "
        "the target sweep reports adjacent changes on its explicit grid.",
        "No empirical data, trading rule, profitability test, or real-world "
        "financial application is included.",
    )
    comparison_results = (
        "At each tested damping ratio, target samples are compared with the "
        "closed-form response and oscillator integrator output.",
        "For damping ratio at or above one, the exact response is "
        "non-oscillatory; solver success alone is not used to infer stability.",
        "Peak rate, overshoot, and sampled settling time are separate metrics "
        "and can have different trends.",
        succeeded_text,
    )
    shell = DampingTransferCaseStudy(
        case_study_id="oscillator_to_second_order_adjustment",
        source_experiment_id=record.specification.experiment_id,
        source_experiment=record,
        source_hypothesis=source_hypothesis,
        source_peak_rate_endpoint_change=source_change.value,
        source_metric_id="maximum_absolute_velocity_m_per_s",
        source_metric_unit="m/s",
        source_metric_parameter_values=source_damping_values,
        source_metric_values=tuple(source_metric_values),
        source_damping_ratio_values=source_zeta_values,
        source_normalized_peak_rate_values_per_s=source_normalized_rates,
        source_target_peak_rate_absolute_differences_per_s=tuple(
            source_target_differences
        ),
        amplitude_scale_m=amplitude_scale_m,
        source_observation_ids=source_hypothesis.source_behavior_observation_ids,
        target_equation=target_equation,
        target_damping_ratios=ratios,
        baseline_damping_ratio=baseline,
        time_constant_s=tau,
        initial_error=target_initial_error,
        initial_rate_per_s=target_initial_rate,
        duration_s=duration,
        output_point_count=output_point_count,
        solver_configuration=solver_configuration,
        sweep_points=points_tuple,
        peak_rate_trend_on_tested_grid=peak_trend,
        overshoot_trend_on_tested_grid=overshoot_trend,
        settling_time_trend_on_tested_grid=settling_trend,
        mathematical_correspondence="equivalent_under_declared_constant_target_transformation",
        target_domain_interpretation="structural_analogy_only",
        real_world_application_status="unvalidated",
        assumptions=(
            "The target reference is constant after its initial step.",
            "The response is linear, time-invariant, deterministic, and unforced "
            "after that step.",
            "Adjustment state is dimensionless; adjustment rate uses s^-1 by "
            "convention.",
            "The initial target error is +1 with zero initial error rate.",
            "The state scale is x=A*q with A=0.1 m per dimensionless state unit.",
            "The mapping uses tau=sqrt(m/k) and zeta=c/(2*sqrt(m*k)); c is "
            "kg/s while zeta is dimensionless.",
        ),
        comparison_results=comparison_results,
        limitations=limitations,
        report_markdown="",
    )
    return shell.model_copy(update={"report_markdown": _render_report(shell)})

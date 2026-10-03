"""A bounded analytical heat-diffusion and synthetic-signal case study.

The existing heat experiment remains a steady spatial BVP. This module adds a
separate closed-form transient model for homogeneous Dirichlet perturbations;
it does not turn the BVP solver into a time-dependent solver.
"""

from __future__ import annotations

from collections.abc import Sequence
from math import exp, isfinite, pi, sqrt
from typing import Literal

import numpy as np
from pydantic import Field, field_validator, model_validator

from newton_lab.discovery import DiscoveryResult, run_discovery
from newton_lab.exceptions import ScientificValidationError
from newton_lab.experiments import (
    ExperimentRecord,
    ExperimentRunStatus,
    make_heat_conduction_experiment,
    run_experiment,
)
from newton_lab.simulation import SimulationModel
from newton_lab.visualization import build_experiment_report


class SpatialMode(SimulationModel):
    """One sine mode; amplitude has the field's units."""

    mode_number: int = Field(ge=1, strict=True)
    amplitude: float

    @field_validator("amplitude")
    @classmethod
    def finite_amplitude(cls, value: float) -> float:
        if isinstance(value, bool) or not isfinite(value):
            raise ValueError("mode amplitude must be finite")
        return value


class HeatDiffusionModel(SimulationModel):
    """Finite sine-series solution about a fixed-end steady temperature.

    ``T(x,t) = T_base(x) + sum(a_n sin(n*pi*x/L) exp(-alpha*(n*pi/L)^2*t))``.
    The endpoints remain at their fixed temperatures because the perturbation
    is zero there. The model assumes a homogeneous 1D medium, constant positive
    diffusivity, no internal source, and fixed endpoint temperatures.
    """

    length_m: float = Field(gt=0)
    left_temperature_k: float
    right_temperature_k: float
    diffusivity_m2_per_s: float = Field(gt=0)
    modes: tuple[SpatialMode, ...] = Field(min_length=1)

    @field_validator(
        "length_m", "left_temperature_k", "right_temperature_k", "diffusivity_m2_per_s"
    )
    @classmethod
    def finite_parameters(cls, value: float) -> float:
        if isinstance(value, bool) or not isfinite(value):
            raise ValueError("model parameters must be finite")
        return value

    @model_validator(mode="after")
    def unique_modes(self) -> HeatDiffusionModel:
        numbers = tuple(mode.mode_number for mode in self.modes)
        if len(numbers) != len(set(numbers)):
            raise ValueError("mode numbers must be unique")
        return self


class HeatDiffusionEvaluation(SimulationModel):
    """Analytical temperatures and modal evidence on requested finite grids."""

    positions_m: tuple[float, ...] = Field(min_length=2)
    times_s: tuple[float, ...] = Field(min_length=1)
    temperature_k: tuple[tuple[float, ...], ...]
    modal_amplitudes: tuple[tuple[float, ...], ...]
    attenuation_factors: tuple[tuple[float, ...], ...]
    perturbation_rms_gradient_k_per_m: tuple[float, ...]

    @model_validator(mode="after")
    def aligned_results(self) -> HeatDiffusionEvaluation:
        if len(self.temperature_k) != len(self.times_s):
            raise ValueError("temperature rows must match requested times")
        if any(len(row) != len(self.positions_m) for row in self.temperature_k):
            raise ValueError("temperature rows must match requested positions")
        if not (
            len(self.modal_amplitudes)
            == len(self.attenuation_factors)
            == len(self.perturbation_rms_gradient_k_per_m)
            == len(self.times_s)
        ):
            raise ValueError("modal results must match requested times")
        if any(
            len(amplitudes) != len(factors)
            for amplitudes, factors in zip(
                self.modal_amplitudes, self.attenuation_factors, strict=True
            )
        ):
            raise ValueError("amplitude and attenuation rows must align")
        return self


def evaluate_heat_diffusion(
    model: HeatDiffusionModel,
    positions_m: Sequence[float],
    times_s: Sequence[float],
) -> HeatDiffusionEvaluation:
    """Evaluate the exact finite modal solution without mutating inputs.

    The RMS gradient is the exact spatial RMS of the transient perturbation
    derivative over ``[0,L]`` (not a finite-difference estimate). Orthogonality
    of the cosine derivatives gives ``sqrt(sum((a_n*n*pi/L)^2*factor_n^2)/2)``.
    """
    positions = _finite_sequence(positions_m, "positions_m", minimum_length=2)
    times = _finite_sequence(times_s, "times_s", minimum_length=1)
    if any(x < 0.0 or x > model.length_m for x in positions):
        raise ScientificValidationError("positions_m must lie within [0, length_m]")
    if any(t < 0.0 for t in times):
        raise ScientificValidationError("times_s must be non-negative")

    x = np.asarray(positions, dtype=np.float64)
    base = np.asarray(steady_temperature_profile(model, positions), dtype=np.float64)
    temperatures: list[tuple[float, ...]] = []
    all_amplitudes: list[tuple[float, ...]] = []
    all_factors: list[tuple[float, ...]] = []
    variations: list[float] = []
    for time_s in times:
        factors = tuple(
            exp(
                -model.diffusivity_m2_per_s
                * (mode.mode_number * pi / model.length_m) ** 2
                * time_s
            )
            for mode in model.modes
        )
        amplitudes = tuple(
            mode.amplitude * factor
            for mode, factor in zip(model.modes, factors, strict=True)
        )
        row = base.copy()
        for mode, amplitude in zip(model.modes, amplitudes, strict=True):
            row += amplitude * np.sin(mode.mode_number * pi * x / model.length_m)
        variation = sqrt(
            sum(
                (mode.amplitude * mode.mode_number * pi / model.length_m * factor) ** 2
                for mode, factor in zip(model.modes, factors, strict=True)
            )
            / 2.0
        )
        if not np.isfinite(row).all() or not isfinite(variation):
            raise ScientificValidationError(
                "requested model values overflowed finite floating-point range"
            )
        temperatures.append(tuple(float(value) for value in row))
        all_amplitudes.append(amplitudes)
        all_factors.append(factors)
        variations.append(variation)
    return HeatDiffusionEvaluation(
        positions_m=positions,
        times_s=times,
        temperature_k=tuple(temperatures),
        modal_amplitudes=tuple(all_amplitudes),
        attenuation_factors=tuple(all_factors),
        perturbation_rms_gradient_k_per_m=tuple(variations),
    )


def steady_temperature_profile(
    model: HeatDiffusionModel, positions_m: Sequence[float]
) -> tuple[float, ...]:
    """Return the exact linear steady profile for fixed endpoint temperatures."""
    positions = _finite_sequence(positions_m, "positions_m", minimum_length=1)
    if any(x < 0.0 or x > model.length_m for x in positions):
        raise ScientificValidationError("positions_m must lie within [0, length_m]")
    result = tuple(
        model.left_temperature_k
        + (model.right_temperature_k - model.left_temperature_k) * x / model.length_m
        for x in positions
    )
    if not np.isfinite(result).all():
        raise ScientificValidationError("steady profile exceeds finite numeric range")
    return result


class SyntheticSignalComparison(SimulationModel):
    """Synthetic abstract signal filtered by the same modal factors."""

    unit: str = "arbitrary signal units"
    initial_mode_amplitudes: tuple[float, ...] = Field(min_length=2)
    mode_numbers: tuple[int, ...] = Field(min_length=2)
    times_s: tuple[float, ...] = Field(min_length=1)
    modal_amplitudes: tuple[tuple[float, ...], ...]
    sampled_values: tuple[tuple[float, ...], ...]

    @model_validator(mode="after")
    def validate_signal(self) -> SyntheticSignalComparison:
        if len(self.initial_mode_amplitudes) != len(self.mode_numbers):
            raise ValueError("signal amplitudes and mode numbers must align")
        if len(self.modal_amplitudes) != len(self.times_s) or len(
            self.sampled_values
        ) != len(self.times_s):
            raise ValueError("signal samples must match requested times")
        if any(len(row) != len(self.mode_numbers) for row in self.modal_amplitudes):
            raise ValueError("signal modal amplitudes must align with modes")
        if any(not np.isfinite(values).all() for values in self.sampled_values):
            raise ValueError("signal samples must be finite")
        return self


class HeatDiffusionSmoothingCaseStudy(SimulationModel):
    """Reproducible steady-source, transient-mode, and synthetic-signal report."""

    case_study_id: Literal["phase19_heat_diffusion_smoothing"] = (
        "phase19_heat_diffusion_smoothing"
    )
    source_experiment: ExperimentRecord
    source_bvp_run_index: int = Field(ge=0)
    source_bvp_numerically_accepted: bool
    source_bvp_max_profile_error_k: float
    source_equation_record_ids: tuple[str, ...]
    source_experiment_report: str
    discovery: DiscoveryResult
    diffusion_model: HeatDiffusionModel
    evaluation: HeatDiffusionEvaluation
    synthetic_signal: SyntheticSignalComparison
    mathematical_relationship: str
    evidence_classification: Literal["analytical_model_result"] = (
        "analytical_model_result"
    )
    cross_domain_status: Literal["synthetic_structural_analogy_only"] = (
        "synthetic_structural_analogy_only"
    )
    real_world_or_financial_validation: Literal[False] = False
    limitations: tuple[str, ...] = Field(min_length=1)
    falsification_tests: tuple[str, ...] = Field(min_length=1)
    report_markdown: str


def run_heat_diffusion_smoothing_case_study(
    *,
    source_record: ExperimentRecord | None = None,
    length_m: float = 1.0,
    left_temperature_k: float = 300.0,
    right_temperature_k: float = 400.0,
    diffusivity_m2_per_s: float = 1.0e-4,
    modes: tuple[SpatialMode, ...] = (
        SpatialMode(mode_number=1, amplitude=1.0),
        SpatialMode(mode_number=4, amplitude=0.5),
    ),
    times_s: tuple[float, ...] = (0.0, 60.0, 300.0, 1200.0),
    sample_count: int = 101,
) -> HeatDiffusionSmoothingCaseStudy:
    """Run the existing steady BVP and an independent exact transient example."""
    if (
        isinstance(sample_count, bool)
        or not isinstance(sample_count, int)
        or sample_count < 2
    ):
        raise ScientificValidationError("sample_count must be an integer >= 2")
    model = HeatDiffusionModel(
        length_m=length_m,
        left_temperature_k=left_temperature_k,
        right_temperature_k=right_temperature_k,
        diffusivity_m2_per_s=diffusivity_m2_per_s,
        modes=modes,
    )
    times = _finite_sequence(times_s, "times_s", minimum_length=1)
    if any(time < 0.0 for time in times):
        raise ScientificValidationError("times_s must be non-negative")

    record = source_record
    if record is None:
        specification = make_heat_conduction_experiment(
            experiment_id="phase19_steady_heat_source",
            length_m=length_m,
            left_temperature_k=left_temperature_k,
            right_temperature_k=right_temperature_k,
            mesh_points=21,
            sweep_values=None,
            metrics=("mesh_point_count", "boundary_residual_max_abs"),
        )
        record = run_experiment(specification)
    if record.specification.bvp_problem is None:
        raise ScientificValidationError("source record must contain a BVP problem")
    if record.specification.model_id != "boundary_value_problem":
        raise ScientificValidationError("source record must be a BVP experiment")
    if len(record.runs) != 1 or record.runs[0].bvp_diagnostics is None:
        raise ScientificValidationError(
            "source BVP record must contain exactly one run with BVP diagnostics"
        )
    source_run = record.runs[0]
    diagnostics = source_run.bvp_diagnostics
    if source_run.status != ExperimentRunStatus.SUCCEEDED or diagnostics is None:
        raise ScientificValidationError("source BVP run did not succeed")
    if record.specification.bvp_problem.problem_id != "steady_linear_heat_conduction":
        raise ScientificValidationError("source BVP is not the built-in heat model")

    problem = record.specification.bvp_problem
    source_parameters = {
        parameter.name: parameter.value for parameter in problem.parameters
    }
    required_parameters = {"left_temperature_k", "right_temperature_k"}
    if not required_parameters.issubset(source_parameters):
        raise ScientificValidationError(
            "source heat BVP is missing endpoint temperature parameters"
        )
    if (
        not np.isclose(problem.domain[0], 0.0)
        or not np.isclose(problem.domain[1], model.length_m)
        or not np.isclose(
            source_parameters["left_temperature_k"], model.left_temperature_k
        )
        or not np.isclose(
            source_parameters["right_temperature_k"], model.right_temperature_k
        )
        or "No internal heat generation." not in problem.assumptions
        or not any(
            "constant thermal conductivity" in item for item in problem.assumptions
        )
    ):
        raise ScientificValidationError(
            "source BVP conditions must match the case-study steady baseline"
        )

    source_x = np.asarray(diagnostics.mesh, dtype=np.float64)
    expected = (
        model.left_temperature_k
        + (model.right_temperature_k - model.left_temperature_k)
        * source_x
        / model.length_m
    )
    profile_error = float(
        np.max(np.abs(np.asarray(diagnostics.solution[0]) - expected))
    )
    positions = tuple(float(x) for x in np.linspace(0.0, model.length_m, sample_count))
    evaluation = evaluate_heat_diffusion(model, positions, times)
    signal = _build_signal(model, positions, times)
    discovery = run_discovery((record,))
    experiment_report = build_experiment_report(
        (record,), title="Phase 19 steady heat-conduction source experiment"
    ).body
    limitations = (
        "The built-in equation registry has no heat-conduction equation record; "
        "the source BVP declares no equation_record_ids, so discovery must not "
        "manufacture an equation identity or relationship.",
        "The transient model is an exact finite sine-series solution, not a "
        "general numerical PDE solver or a validated heat-transfer simulation.",
        "The signal is synthetic and uses arbitrary units; reduced high-frequency "
        "content is not evidence of denoising, improved prediction, or utility.",
        "Finite spatial samples display the analytic solution only at the listed "
        "grid positions; modal attenuation and RMS gradient are analytic quantities.",
    )
    falsification_tests = (
        "Compare the modal expression against an independent transient PDE solver "
        "under the same homogeneous Dirichlet conditions and refine its "
        "mesh/time step.",
        "Add nonconstant diffusivity, a source term, or changing boundaries and "
        "test whether the stated modal correspondence still applies.",
        "For an application claim, predefine an out-of-sample task and compare "
        "smoothed and unsmoothed synthetic signals against known truth.",
    )
    relationship = (
        "For constant diffusivity and fixed endpoint temperatures, subtracting "
        "the linear steady profile T_base(x) reduces the transient equation to "
        "u_t = alpha*u_xx with homogeneous Dirichlet boundaries. Each sine mode "
        "is an eigenfunction and decays by exp(-alpha*(n*pi/L)^2*t). Thus the "
        "time-independent profile solves the steady equation, while the transient "
        "solution approaches it for this finite-mode perturbation."
    )
    report = _render_report(
        record=record,
        diagnostics_accepted=diagnostics.numerical_accepted,
        profile_error=profile_error,
        model=model,
        evaluation=evaluation,
        signal=signal,
        relationship=relationship,
        source_experiment_report=experiment_report,
        limitations=limitations,
        falsification_tests=falsification_tests,
    )
    return HeatDiffusionSmoothingCaseStudy(
        source_experiment=record,
        source_bvp_run_index=source_run.run_index,
        source_bvp_numerically_accepted=diagnostics.numerical_accepted,
        source_bvp_max_profile_error_k=profile_error,
        source_equation_record_ids=record.specification.equation_record_ids,
        source_experiment_report=experiment_report,
        discovery=discovery,
        diffusion_model=model,
        evaluation=evaluation,
        synthetic_signal=signal,
        mathematical_relationship=relationship,
        limitations=limitations,
        falsification_tests=falsification_tests,
        report_markdown=report,
    )


def _build_signal(
    model: HeatDiffusionModel,
    positions: tuple[float, ...],
    times: tuple[float, ...],
) -> SyntheticSignalComparison:
    amplitudes_initial = (1.0, 0.4)
    mode_numbers = (1, 4)
    x = np.asarray(positions, dtype=np.float64)
    rows: list[tuple[float, ...]] = []
    amplitude_rows: list[tuple[float, ...]] = []
    for time_s in times:
        factors = tuple(
            exp(
                -model.diffusivity_m2_per_s
                * (number * pi / model.length_m) ** 2
                * time_s
            )
            for number in mode_numbers
        )
        amplitudes = tuple(
            initial * factor
            for initial, factor in zip(amplitudes_initial, factors, strict=True)
        )
        values = np.zeros_like(x)
        for number, amplitude in zip(mode_numbers, amplitudes, strict=True):
            values += amplitude * np.sin(number * pi * x / model.length_m)
        rows.append(tuple(float(value) for value in values))
        amplitude_rows.append(amplitudes)
    return SyntheticSignalComparison(
        initial_mode_amplitudes=amplitudes_initial,
        mode_numbers=mode_numbers,
        times_s=times,
        modal_amplitudes=tuple(amplitude_rows),
        sampled_values=tuple(rows),
    )


def _finite_sequence(
    values: Sequence[float], name: str, *, minimum_length: int
) -> tuple[float, ...]:
    if any(isinstance(value, bool) for value in values):
        raise ScientificValidationError(f"{name} cannot contain booleans")
    try:
        result = tuple(float(value) for value in values)
    except (TypeError, ValueError) as error:
        raise ScientificValidationError(
            f"{name} must contain finite numbers"
        ) from error
    if len(result) < minimum_length or not np.isfinite(result).all():
        raise ScientificValidationError(
            f"{name} must contain at least {minimum_length} finite values"
        )
    return result


def _render_report(
    *,
    record: ExperimentRecord,
    diagnostics_accepted: bool,
    profile_error: float,
    model: HeatDiffusionModel,
    evaluation: HeatDiffusionEvaluation,
    signal: SyntheticSignalComparison,
    relationship: str,
    source_experiment_report: str,
    limitations: tuple[str, ...],
    falsification_tests: tuple[str, ...],
) -> str:
    diagnostics = record.runs[0].bvp_diagnostics
    assert diagnostics is not None
    lines = [
        "# Phase 19 — Heat conduction, diffusion, and smoothing",
        "",
        "## Research question",
        "",
        "Under what assumptions does the existing steady heat BVP connect to "
        "transient diffusion, and does its modal attenuation reduce variation "
        "in a controlled synthetic signal?",
        "",
        "## Existing steady heat-conduction model",
        "",
        "The source is Newton Lab's existing `steady_linear_heat_conduction` "
        "BVP. For constant conductivity k and no internal heat generation, its "
        "equation is d/dx(k dT/dx)=0, or T_xx=0. It uses domain "
        f"[0, {model.length_m:g}] m, temperature in K, and fixed endpoints "
        f"{model.left_temperature_k:g} K and {model.right_temperature_k:g} K. "
        "The model assumes a homogeneous material and constant conductivity; "
        "the source term is absent. Conductivity and thermal capacity are not "
        "estimated by this BVP. The solver returned a spatial mesh of "
        f"{len(diagnostics.mesh)} points, solver_success={diagnostics.solver_success}, "
        f"numerical_accepted={diagnostics_accepted}, maximum boundary residual "
        f"{max(abs(value) for value in diagnostics.boundary_residuals):.3g}, "
        f"and maximum differential residual "
        f"{max(diagnostics.differential_residuals):.3g}. Its maximum difference "
        f"from the analytical steady line on that mesh is {profile_error:.3g} K. "
        "This result is spatial and steady; it is not a time simulation.",
        "",
        "## Transient diffusion model and relationship",
        "",
        "The transient model is T_t = alpha T_xx, with alpha > 0 in m^2/s, "
        "constant material properties, no internal source, and fixed endpoint "
        "temperatures. Its steady equation T_xx=0 has the linear profile "
        "T_base(x)=T_left+(T_right-T_left)x/L. Writing u=T-T_base gives "
        "u_t=alpha u_xx and u(0,t)=u(L,t)=0. For the declared finite sine-series "
        "initial perturbation, the exact solution is "
        "T=T_base+sum(a_n sin(n*pi*x/L) exp(-alpha*(n*pi/L)^2*t)). "
        "This establishes the equilibrium relation and convergence of these "
        "finite modes under the stated conditions; it does not cover arbitrary "
        "media, boundary histories, or forcing. This is a mathematical "
        "correspondence under matched assumptions, not evidence that the steady "
        "BVP alone validates transient material properties or physical time "
        "scales. " + relationship,
        "",
        "## Modal attenuation experiment",
        "",
        f"Length={model.length_m:g} m; alpha={model.diffusivity_m2_per_s:g} m^2/s; "
        f"sampled positions={len(evaluation.positions_m)}; times (s)="
        f"{', '.join(f'{item:g}' for item in evaluation.times_s)}. The listed "
        "temperature values are analytical evaluations on a finite spatial grid. "
        "The attenuation factor is exact for each mode.",
        "",
        "| Mode n | Initial amplitude (K) | "
        + " | ".join(f"factor at t={time:g} s" for time in evaluation.times_s)
        + " |",
        "|---:|---:|" + "---:|" * len(evaluation.times_s),
    ]
    for index, mode in enumerate(model.modes):
        lines.append(
            f"| {mode.mode_number} | {mode.amplitude:.6g} | "
            + " | ".join(f"{row[index]:.6g}" for row in evaluation.attenuation_factors)
            + " |"
        )
    lines.extend(
        [
            "",
            "The exact perturbation spatial RMS-gradient metric is "
            "sqrt((1/2) sum((a_n n pi/L)^2 factor_n^2)) in K/m. This uses "
            "orthogonality over the continuous interval, not a finite-difference "
            "grid estimate:",
            "",
            "| Time (s) | Perturbation RMS gradient (K/m) |",
            "|---:|---:|",
        ]
    )
    lines.extend(
        f"| {time:g} | {variation:.8g} |"
        for time, variation in zip(
            evaluation.times_s,
            evaluation.perturbation_rms_gradient_k_per_m,
            strict=True,
        )
    )
    lines.extend(
        [
            "",
            "For each mode, the decay rate is proportional to n squared, so the "
            "higher-frequency mode decays faster under this model. The metric "
            "tracks only the transient perturbation; the steady baseline gradient "
            "is not counted.",
            "",
            "## Synthetic signal comparison",
            "",
            "The target is a deterministic synthetic sum of the same n=1 and n=4 "
            "sine modes, sampled at the same positions and transformed with the "
            "same exponential factors. Sampling evaluates the continuous "
            "analytical solution on a finite grid; it is not a discrete "
            "finite-difference filter. Its units are arbitrary. This is a "
            "mathematical smoothing analogy, not market data, a trading strategy, "
            "or a forecasting test.",
            "",
            "| Mode n | Initial signal amplitude | "
            + " | ".join(f"amplitude at t={time:g} s" for time in signal.times_s)
            + " |",
            "|---:|---:|" + "---:|" * len(signal.times_s),
        ]
    )
    for index, number in enumerate(signal.mode_numbers):
        lines.append(
            f"| {number} | {signal.initial_mode_amplitudes[index]:.6g} | "
            + " | ".join(f"{row[index]:.6g}" for row in signal.modal_amplitudes)
            + " |"
        )
    lines.extend(
        [
            "",
            "The shared mathematics is the mode-wise linear attenuation operator. "
            "Smoothing is not denoising: high-frequency components can carry "
            "meaningful signal, and reducing variation does not establish better "
            "recovery of an underlying truth.",
            "",
            "## Newton Lab discovery and evidence",
            "",
            f"Source experiment: `{record.specification.experiment_id}`; model: "
            f"`{record.specification.model_id}`; source equation IDs declared: "
            f"{record.specification.equation_record_ids or 'none'}. The existing "
            "discovery pipeline was run on the unchanged source experiment. "
            "Because the registry currently contains no heat-conduction equation "
            "record and the source has no equation ID, no registry structure match "
            "or cross-domain relationship is claimed. The analytical derivation "
            "above is classified as an analytical model result; the signal result "
            "is a synthetic structural analogy only. No real-world or financial "
            "application is validated.",
            "",
            "### Source experiment report",
            "",
            source_experiment_report,
            "",
            "## Limitations",
            "",
        ]
    )
    lines.extend(f"- {item}" for item in limitations)
    lines.extend(["", "## Tests that could challenge the connection", ""])
    lines.extend(f"- {item}" for item in falsification_tests)
    lines.extend(
        [
            "",
            "## Reproduce",
            "",
            "From the Newton Lab workspace, run:",
            "",
            "```python",
            "from newton_lab.heat_diffusion_case_study import (",
            "    run_heat_diffusion_smoothing_case_study,",
            ")",
            "",
            "case = run_heat_diffusion_smoothing_case_study()",
            "print(case.report_markdown)",
            "```",
            "",
        ]
    )
    return "\n".join(lines)

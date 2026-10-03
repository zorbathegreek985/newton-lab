"""Focused, reproducible comparison of Newton Lab case studies.

This module keeps physical models and inference procedures distinct.
"""

from dataclasses import dataclass
from enum import StrEnum
from math import isfinite, pi, sqrt

import numpy as np

from newton_lab.dynamics import simulate_damped_oscillator
from newton_lab.exceptions import ScientificValidationError
from newton_lab.heat_diffusion_case_study import (
    HeatDiffusionModel,
    SpatialMode,
    evaluate_heat_diffusion,
    steady_temperature_profile,
)
from newton_lab.modal_research import heat_mode_decay_rate
from newton_lab.pendulum import simulate_damped_pendulum


class EvidenceStatus(StrEnum):
    """Evidence status for one narrowly scoped relationship."""

    ESTABLISHED = "established mathematical relationship"
    SUPPORTED = "numerically supported relationship"
    ANALOGY = "conditional structural analogy"
    CANDIDATE = "candidate requiring further testing"
    CONTRADICTED = "rejected or contradicted candidate"
    INCONCLUSIVE = "inconclusive"


class StudyRole(StrEnum):
    """Physical equations and inference methods have different roles."""

    PHYSICAL_MODEL = "physical model"
    INFERENCE = "inference procedure"


@dataclass(frozen=True, slots=True)
class ComparisonEntry:
    study_id: str
    phase: int
    name: str
    role: StudyRole
    equation: str
    assumptions: tuple[str, ...]
    states_units: tuple[str, ...]
    parameters_units: tuple[str, ...]
    timescale: str
    dimensionless_groups: str
    observables: tuple[str, ...]
    analytical_reference: str
    numerical_evidence: str
    approximation_limits: str
    evidence_status: EvidenceStatus
    sources: tuple[str, ...]
    limitations: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class CrossSystemFinding:
    finding_id: str
    hypothesis: str
    status: EvidenceStatus
    support: str
    counterexample: str
    assumptions: tuple[str, ...]
    tested_range: str
    evidence: tuple[str, ...]
    untested: tuple[str, ...]
    sources: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class ObservableCompatibility:
    compatible: bool
    left_unit: str
    right_unit: str
    left_meaning: str
    right_meaning: str
    status: str


@dataclass(frozen=True, slots=True)
class TrajectoryComparison:
    time_s: tuple[float, ...]
    oscillator_x_normalized: tuple[float, ...]
    linear_pendulum_angle: tuple[float, ...]
    nonlinear_pendulum_angle: tuple[float, ...]
    linear_max_error: float
    linear_normalized_error: float
    nonlinear_max_error: float
    nonlinear_normalized_error: float
    normalization_scale: float
    initial_angle_rad: float
    initial_angular_velocity_rad_per_s: float
    duration_tau: float
    sample_count: int
    relative_tolerance: float
    absolute_tolerance: float
    mapping_assumptions: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class CrossSystemReport:
    matrix: tuple[ComparisonEntry, ...]
    findings: tuple[CrossSystemFinding, ...]
    trajectory: TrajectoryComparison
    markdown: str


def build_comparison_matrix() -> tuple[ComparisonEntry, ...]:
    """Build a source-linked matrix; Phase 21 JSON artifacts are not present."""
    return (
        ComparisonEntry(
            "phase18_damping_transfer",
            18,
            "Damped linear oscillator",
            StudyRole.PHYSICAL_MODEL,
            "m x¨ + c x˙ + k x = 0",
            (
                "linear spring",
                "linear viscous damping",
                "constant parameters",
                "no force",
            ),
            ("x: m", "x˙: m/s"),
            ("m: kg", "c: kg/s", "k: N/m"),
            "ωn=√(k/m); underdamped envelope time 1/γ, γ=c/(2m)",
            "ζ=c/(2√(mk)); τ=ωn t",
            ("displacement", "velocity", "sampled envelope"),
            "Characteristic roots and oscillator closed form.",
            "Source report: m=1 kg, k=4 N/m, x0=0.1 m, v0=0; c=0,0.1,0.2 kg/s.",
            "Envelope rate applies only when underdamped; zero damping does not decay.",
            EvidenceStatus.ANALOGY,
            (
                "docs/cross_domain_case_study.md",
                "src/newton_lab/dynamics.py",
                "src/newton_lab/cross_domain.py",
            ),
            (
                "Phase 21 source JSON is absent; values are report-derived.",
                "No empirical application data.",
            ),
        ),
        ComparisonEntry(
            "phase26_nonlinear_pendulum",
            26,
            "Nonlinear pendulum",
            StudyRole.PHYSICAL_MODEL,
            "θ¨ + (g/L) sin θ = 0",
            (
                "ideal point mass",
                "fixed rigid length",
                "uniform gravity",
                "undamped, unforced",
            ),
            ("θ: rad", "θ˙: rad/s"),
            ("m: kg", "L: m", "g: m/s²"),
            "√(L/g); libration period also depends on energy",
            "τ=√(g/L)t; energy h=θ˙²L/(2g)+1−cosθ",
            ("turning amplitude", "period", "angle error", "state-coordinate phase"),
            "Elliptic-integral period for librations; full-sine numerical solver.",
            "Phase 26 joint initial-angle/velocity grid; report covers windows "
            "to 10 linear periods; rotations and separatrix excluded from "
            "libration-period formula.",
            "Thresholds depend on state, observable, tolerance and window.",
            EvidenceStatus.SUPPORTED,
            (
                "docs/research/phase_26_pendulum_initial_velocity.md",
                "src/newton_lab/pendulum_initial_velocity.py",
                "src/newton_lab/pendulum.py",
            ),
            ("Ideal equations only; no measurements. Phase 21 JSON record absent.",),
        ),
        ComparisonEntry(
            "phase25_small_angle",
            25,
            "Small-angle pendulum",
            StudyRole.PHYSICAL_MODEL,
            "θ¨ + (g/L) θ = 0",
            ("sin θ ≈ θ", "fixed length", "uniform gravity", "undamped, unforced"),
            ("θ: rad", "θ˙: rad/s"),
            ("L: m", "g: m/s²"),
            "ω0=√(g/L); T0=2π√(L/g)",
            "τ=√(g/L)t",
            ("trajectory", "period", "phase"),
            "Exact harmonic solution of the linearized equation.",
            "Phase 25/26 compare DOP853 finite-grid trajectories with analytical "
            "references.",
            "Linearization error accumulates and depends on metric and horizon.",
            EvidenceStatus.ESTABLISHED,
            (
                "docs/research/phase_25_pendulum_approximation.md",
                "docs/research/phase_26_pendulum_initial_velocity.md",
            ),
            ("Linear model is not the full nonlinear pendulum.",),
        ),
        ComparisonEntry(
            "phase22_modal_dynamics",
            22,
            "Oscillator and heat modal dynamics",
            StudyRole.PHYSICAL_MODEL,
            "Heat: a˙n=−α(nπ/L)²an; oscillator r=−γ±iωd",
            (
                "linear time-invariant equations",
                "fixed-end heat perturbation",
                "underdamped oscillator",
            ),
            ("heat coefficient: K", "oscillator displacement: m"),
            ("α: m²/s", "L: m", "m: kg", "c: kg/s", "k: N/m"),
            "heat τn=1/λn; oscillator envelope τ=1/γ and period=2π/ωd",
            "heat λn t; oscillator ζ",
            ("mode coefficient", "displacement envelope", "sampled peaks"),
            "Heat mode is exponential; oscillator envelope has an oscillatory carrier.",
            "Phase 22 uses Phase 18/19 APIs and reports analytical rates and "
            "sampled estimates.",
            "Shared exponential factor does not imply same spectrum, observable "
            "or mechanism.",
            EvidenceStatus.ANALOGY,
            (
                "src/newton_lab/modal_research.py",
                "docs/research/phase_22_modal_dynamics.md",
            ),
            ("Generated report is the available result record.",),
        ),
        ComparisonEntry(
            "phase23_modal_estimation",
            23,
            "Single-rate estimation",
            StudyRole.INFERENCE,
            "Fit y(t)=A exp(−λt); estimator assumption, not physical law",
            ("synthetic truth", "finite cadence/window", "declared noise and method"),
            ("signal: declared units", "time: s"),
            ("A: signal units", "λ: s⁻¹", "σ: signal units"),
            "window and cadence relative to 1/λ",
            "λT and λΔt",
            ("rate estimate", "fit error", "status"),
            "Synthetic truth is kept separate from fitted observations.",
            "Phase 23 reports cadence, duration, noise, mixed heat modes and "
            "oscillator envelope.",
            "Fits can fail at noise floors, sparse cadence or short windows.",
            EvidenceStatus.INCONCLUSIVE,
            (
                "src/newton_lab/modal_estimation.py",
                "docs/research/phase_23_modal_estimation.md",
            ),
            (
                "Report-derived synthetic results do not guarantee behavior on "
                "measured data.",
            ),
        ),
        ComparisonEntry(
            "phase24_mode_identification",
            24,
            "Single/two-mode identification",
            StudyRole.INFERENCE,
            "Fit ΣAi exp(−λi t); compare models under declared criteria",
            (
                "synthetic positive modes",
                "known-noise BIC where declared",
                "multiple starts",
            ),
            ("observation: signal units", "time: s"),
            ("Ai: signal units", "λi: s⁻¹", "σ: signal units"),
            "window and rate separation jointly affect recovery",
            "Δ=|λ2−λ1|/max(λ1,λ2)",
            ("model selection", "rate recovery", "residual", "conditioning"),
            "Known synthetic truth is used only after fitting.",
            "Phase 24 tests noise repetitions, irregular sampling, weak and "
            "close modes.",
            "Extra parameters lower residual; low residual alone is not rate recovery.",
            EvidenceStatus.SUPPORTED,
            (
                "src/newton_lab/modal_identification.py",
                "docs/research/phase_24_mode_identification.md",
            ),
            ("Finite synthetic design does not generalize to untested observations.",),
        ),
        ComparisonEntry(
            "phase19_steady_heat",
            19,
            "Steady heat conduction",
            StudyRole.PHYSICAL_MODEL,
            "d/dx(k dT/dx)=0; constant k gives Tₓₓ=0",
            (
                "steady",
                "homogeneous",
                "constant conductivity",
                "no generation",
                "fixed endpoints",
            ),
            ("T: K", "x: m"),
            ("k: W/(m K)", "endpoint temperatures: K"),
            "No time variable in this BVP.",
            "not applicable",
            ("spatial profile", "boundary/differential residual"),
            "Linear profile between fixed endpoint temperatures.",
            "Report: 21 mesh points, accepted solve, zero boundary residual and "
            "1.35e-15 differential residual.",
            "BVP acceptance says nothing alone about transient rates.",
            EvidenceStatus.ESTABLISHED,
            (
                "docs/heat_diffusion_smoothing_case_study.md",
                "src/newton_lab/heat_diffusion_case_study.py",
            ),
            ("No heat equation registry identity; no Phase 21 JSON artifact.",),
        ),
        ComparisonEntry(
            "phase19_transient_heat",
            19,
            "Transient heat diffusion",
            StudyRole.PHYSICAL_MODEL,
            "Tt=αTxx; u=T−Tbase obeys ut=αuxx",
            (
                "1D homogeneous medium",
                "constant α",
                "no source",
                "fixed endpoints",
                "finite sine perturbation",
            ),
            ("T,u: K", "x: m", "t: s"),
            ("α: m²/s", "L: m"),
            "τn=1/[α(nπ/L)²]",
            "λn=α(nπ/L)²; λnt",
            ("modal amplitudes", "grid temperatures", "analytic RMS gradient"),
            "an(t)=an(0)e^(−λnt); Tbase is the steady equilibrium.",
            "Phase 19: L=1 m, α=1e-4 m²/s, modes 1/4, four times and 101 positions.",
            "Exact finite-series case, not a general PDE solver or denoising proof.",
            EvidenceStatus.SUPPORTED,
            (
                "docs/heat_diffusion_smoothing_case_study.md",
                "src/newton_lab/heat_diffusion_case_study.py",
            ),
            ("No independent transient PDE solver cross-check.",),
        ),
    )


def oscillator_dimensionless_parameters(
    mass_kg: float, damping_coefficient_kg_per_s: float, stiffness_n_per_m: float
) -> tuple[float, float]:
    """Return ωn and ζ for the dimensionless oscillator equation."""
    values = (mass_kg, damping_coefficient_kg_per_s, stiffness_n_per_m)
    if not all(isfinite(value) for value in values):
        raise ScientificValidationError("oscillator parameters must be finite")
    if mass_kg <= 0.0 or stiffness_n_per_m <= 0.0:
        raise ScientificValidationError("mass and stiffness must be positive")
    if damping_coefficient_kg_per_s < 0.0:
        raise ScientificValidationError("damping must be non-negative")
    omega = sqrt(stiffness_n_per_m / mass_kg)
    zeta = damping_coefficient_kg_per_s / (2.0 * sqrt(mass_kg * stiffness_n_per_m))
    return omega, zeta


def assess_observable_compatibility(
    *, left_unit: str, right_unit: str, left_meaning: str, right_meaning: str
) -> ObservableCompatibility:
    "Compare declared units and meanings; this caller declaration is"
    "unverified."
    values = (
        left_unit.strip(),
        right_unit.strip(),
        left_meaning.strip(),
        right_meaning.strip(),
    )
    if not all(values):
        raise ScientificValidationError("units and meanings must be non-empty")
    compatible = values[0] == values[1] and values[2] == values[3]
    return ObservableCompatibility(
        compatible,
        *values,
        "caller-declared compatibility assumption; not independently validated",
    )


def normalized_discrepancy(
    reference: tuple[float, ...],
    candidate: tuple[float, ...],
    *,
    normalization_scale: float,
) -> float:
    """Compute max absolute discrepancy using an explicit safe scale."""
    if not reference or len(reference) != len(candidate):
        raise ScientificValidationError("sequences must be non-empty and aligned")
    if not isfinite(normalization_scale) or normalization_scale <= np.finfo(float).tiny:
        raise ScientificValidationError("normalization scale must be safely positive")
    combined = np.asarray((*reference, *candidate), dtype=np.float64)
    if not np.isfinite(combined).all():
        raise ScientificValidationError("sequences must be finite")
    return float(
        np.max(np.abs(np.subtract(reference, candidate))) / normalization_scale
    )


def record_finding(
    *,
    finding_id: str,
    hypothesis: str,
    status: EvidenceStatus,
    support: str,
    counterexample: str,
    assumptions: tuple[str, ...],
    tested_range: str,
    evidence: tuple[str, ...],
    untested: tuple[str, ...],
    sources: tuple[str, ...],
    mapping_declaration: str | None = None,
) -> CrossSystemFinding:
    "Record falsifiable scope and retain user mappings as unvalidated"
    "assumptions."
    required = (finding_id, hypothesis, support, counterexample, tested_range)
    if not all(text.strip() for text in required) or not assumptions or not sources:
        raise ScientificValidationError(
            "finding requires text, assumptions and sources"
        )
    if mapping_declaration is not None:
        if not mapping_declaration.strip():
            raise ScientificValidationError("mapping declaration cannot be blank")
        assumptions += (
            "Caller mapping hypothesis (not independently validated): "
            + mapping_declaration.strip(),
        )
    return CrossSystemFinding(
        finding_id,
        hypothesis,
        status,
        support,
        counterexample,
        assumptions,
        tested_range,
        evidence,
        untested,
        sources,
    )


def compare_linear_oscillator_and_pendulum(
    *,
    initial_angle_rad: float = 0.2,
    initial_angular_velocity_rad_per_s: float = 0.3,
    duration_tau: float = 4.0 * pi,
    sample_count: int = 2001,
    relative_tolerance: float = 1e-10,
    absolute_tolerance: float = 1e-12,
) -> TrajectoryComparison:
    """Compare numerical oscillator, analytic linear pendulum and full pendulum.

    Use m=1 kg, k=1 N/m, L=1 m, g=1 m/s² and coordinate scale 1 m/rad.
    Both natural frequencies are 1 s⁻¹, so τ=t. Full-sine pendulum error is
    kept distinct from numerical agreement with its linearized exact solution.
    """
    scalars = (
        initial_angle_rad,
        initial_angular_velocity_rad_per_s,
        duration_tau,
        relative_tolerance,
        absolute_tolerance,
    )
    if not all(isfinite(value) for value in scalars):
        raise ScientificValidationError("comparison inputs must be finite")
    if duration_tau <= 0 or relative_tolerance <= 0 or absolute_tolerance <= 0:
        raise ScientificValidationError("duration and tolerances must be positive")
    if (
        isinstance(sample_count, bool)
        or not isinstance(sample_count, int)
        or sample_count < 2
    ):
        raise ScientificValidationError("sample_count must be an integer of at least 2")
    oscillator = simulate_damped_oscillator(
        1.0,
        0.0,
        1.0,
        initial_angle_rad,
        initial_angular_velocity_rad_per_s,
        duration_tau,
        num_points=sample_count,
        method="DOP853",
        relative_tolerance=relative_tolerance,
        absolute_tolerance=absolute_tolerance,
    )
    pendulum = simulate_damped_pendulum(
        1.0,
        1.0,
        0.0,
        1.0,
        initial_angle_rad,
        initial_angular_velocity_rad_per_s,
        duration_tau,
        num_points=sample_count,
        method="DOP853",
        relative_tolerance=relative_tolerance,
        absolute_tolerance=absolute_tolerance,
    )
    times = oscillator.time_s
    x_normalized = oscillator.displacement_m
    linear = initial_angle_rad * np.cos(
        times
    ) + initial_angular_velocity_rad_per_s * np.sin(times)
    linear_delta = np.abs(x_normalized - linear)
    nonlinear_delta = np.abs(x_normalized - pendulum.angle_rad)
    linear_scale = max(
        float(np.max(np.abs(x_normalized))), float(np.max(np.abs(linear)))
    )
    scale = max(
        float(np.max(np.abs(x_normalized))), float(np.max(np.abs(pendulum.angle_rad)))
    )
    lmax = float(np.max(linear_delta))
    nmax = float(np.max(nonlinear_delta))
    return TrajectoryComparison(
        tuple(float(x) for x in times),
        tuple(float(x) for x in x_normalized),
        tuple(float(x) for x in linear),
        tuple(float(x) for x in pendulum.angle_rad),
        lmax,
        lmax / linear_scale if linear_scale > np.finfo(float).tiny else lmax,
        nmax,
        nmax / scale if scale > np.finfo(float).tiny else nmax,
        scale,
        initial_angle_rad,
        initial_angular_velocity_rad_per_s,
        duration_tau,
        sample_count,
        relative_tolerance,
        absolute_tolerance,
        (
            "Assumed x/(1 m)=θ/(1 rad): numeric coordinates map, physical "
            "meanings do not.",
            "Oscillator m=1 kg,k=1 N/m,c=0; pendulum L=1 m,g=1 m/s², undamped.",
            "Full pendulum retains sin(θ); its finite-angle discrepancy is not "
            "an integration error.",
        ),
    )


def build_findings(result: TrajectoryComparison) -> tuple[CrossSystemFinding, ...]:
    "Return bounded H1–H5 evidence with support, refutation and untested"
    "scope."
    rate1 = heat_mode_decay_rate(1e-4, 1, 1.0)
    rate4 = heat_mode_decay_rate(1e-4, 4, 1.0)
    mixture_times, mixture_values = heat_signed_aggregate_counterexample()
    selected_indices = (0, 2, 3, 4)
    mixture_summary = ", ".join(
        f"{mixture_values[index]:+.5f} K at {mixture_times[index]:g} s"
        for index in selected_indices
    )
    return (
        CrossSystemFinding(
            "H1_dimensionless_structure",
            "Undamped linear oscillator and small-angle pendulum share X″+X=0 "
            "after scaling.",
            EvidenceStatus.ESTABLISHED,
            "Identical dimensionless equations plus matched normalized initial "
            "states give identical exact trajectories.",
            "Damping, unmatched initial data or the retained finite-angle sine "
            "term breaks that correspondence.",
            ("c=0", "τ=ωn t=√(g/L)t", "matched coordinate and derivative scales"),
            "Exact for all mapped states; numerical check uses θ0=0.2 rad, "
            "θ˙0=0.3 rad/s, τ=0…4π.",
            (
                "Numerical oscillator vs exact linear pendulum normalized max "
                f"error: {result.linear_normalized_error:.3g}.",
                "Full-sine model deviation on same finite grid: "
                f"{result.nonlinear_normalized_error:.3g}.",
            ),
            ("Damped oscillator mapping and other solver implementations.",),
            (
                "src/newton_lab/dynamics.py",
                "src/newton_lab/pendulum.py",
                "docs/research/phase_25_pendulum_approximation.md",
            ),
        ),
        CrossSystemFinding(
            "H2_modal_decay",
            "Heat modes and underdamped oscillator envelopes have exponential "
            "factors but different rates and observables.",
            EvidenceStatus.ANALOGY,
            "For fixed-end heat, λn=α(nπ/L)²; for the oscillator envelope, γ=c/(2m).",
            "The Phase 22 signed two-mode sum crosses zero although each "
            "component decays, so aggregate magnitude need not decay "
            "monotonically.",
            (
                "linear time-invariant equations",
                "modal amplitudes and observable specified",
            ),
            "Heat α=1e-4 m²/s,L=1 m,n=1,4; oscillator m=1 kg,c=0.05–0.2 kg/s,k=4 N/m.",
            (
                f"lambda1={rate1:.8g}/s, lambda4={rate4:.8g}/s; "
                "gamma=0.05/s at c=0.1 kg/s,m=1 kg.",
                "Phase 22 Phase19-derived signed n=1,2 aggregate at x=L/4: "
                + mixture_summary
                + ".",
                "Phase 23/24 report findings show short/noisy/sparse windows and "
                "weak or close modes limit estimation.",
            ),
            ("Untested real observations and arbitrary sampling/noise regimes.",),
            (
                "src/newton_lab/modal_research.py",
                "src/newton_lab/modal_estimation.py",
                "src/newton_lab/modal_identification.py",
                "docs/research/phase_22_modal_dynamics.md",
                "docs/research/phase_23_modal_estimation.md",
                "docs/research/phase_24_mode_identification.md",
            ),
        ),
        CrossSystemFinding(
            "H3_linearization_limits",
            "Small-angle validity depends on initial state, observable, "
            "tolerance and horizon.",
            EvidenceStatus.SUPPORTED,
            "Phase25 grid limits differ by metric; Phase26 results vary across "
            "initial velocity slices.",
            "At 0.35 rad the period and 3-period angle criteria pass while "
            "the 5-period phase criterion fails, refuting one shared cutoff.",
            ("undamped ideal pendulum", "declared windows and tolerances", "libration"),
            "Phase25 rest grid 0.01–2 rad; Phase26 two-dimensional state grid "
            "and windows to 10 periods.",
            (
                "Phase25 reports 0.35 rad period/3-period-angle limits and 0.2 rad "
                "phase limit.",
                "At 0.5 rad, period relative error 1.56%, angle error 0.1359 rad at "
                "3 periods, phase lag 0.4804 rad at 5 periods.",
                "Phase26 rest controls report same 0.35/0.2 rad limits; (0.5,+3 "
                "rad/s) phase difference is 5.11159 rad at 10 periods.",
            ),
            ("No universal threshold; untested states and real apparatus.",),
            (
                "docs/research/phase_25_pendulum_approximation.md",
                "docs/research/phase_26_pendulum_initial_velocity.md",
            ),
        ),
        CrossSystemFinding(
            "H3b_universal_threshold",
            "One initial-angle threshold predicts validity across observables "
            "and initial angular velocities.",
            EvidenceStatus.CONTRADICTED,
            "A shared threshold would require period, trajectory and phase "
            "criteria to change consistently across velocity slices.",
            "Phase25 and Phase26 report distinct phase, period and trajectory "
            "limits; nonzero signed velocity changes the conditional limits.",
            (
                "same ideal undamped model",
                "same declared accuracy criteria",
                "threshold applies across observables and initial states",
            ),
            "Phase25 release from rest over 0.01–2 rad; Phase26 includes "
            "velocity slices at 0 and ±1 rad/s.",
            (
                "At rest, period/angle largest passing tested angle is 0.35 rad; "
                "phase limit is 0.2 rad.",
                "At ±1 rad/s, period limit is 0.2 rad, angle limit is 0.1 rad, "
                "and no tested angle passes the 5-period phase criterion.",
            ),
            (
                "This rejects a universal threshold only for the tested grids "
                "and criteria; denser grids and other observables remain untested.",
            ),
            (
                "docs/research/phase_25_pendulum_approximation.md",
                "docs/research/phase_26_pendulum_initial_velocity.md",
            ),
        ),
        CrossSystemFinding(
            "H4_steady_transient_heat",
            "The fixed-end steady heat line is the equilibrium baseline of the "
            "documented transient equation.",
            EvidenceStatus.ESTABLISHED,
            "Subtracting Tbase leaves homogeneous diffusion; each finite sine "
            "mode decays as exp(-λn t).",
            "Changing boundaries, adding a source, or varying material "
            "properties requires a new derivation.",
            (
                "same endpoints",
                "constant properties",
                "no source",
                "homogeneous Dirichlet perturbation",
            ),
            "Phase19 finite modes, L=1 m, α=1e-4 m²/s, t=0…1200 s, n=1,4.",
            (
                "Exact modal factors and analytic RMS gradient; temperatures sampled "
                "at 101 grid positions.",
            ),
            (
                "Independent transient PDE solver and arbitrary forcing "
                "remain untested.",
            ),
            (
                "docs/heat_diffusion_smoothing_case_study.md",
                "src/newton_lab/heat_diffusion_case_study.py",
            ),
        ),
        CrossSystemFinding(
            "H5_inference_reliability",
            "A low residual is insufficient evidence that a fit recovered "
            "physical modes.",
            EvidenceStatus.INCONCLUSIVE,
            "Phase23/24 compare fits with held-out synthetic truth under varied "
            "designs.",
            "Wrong model selection, ill-conditioning or rates outside the truth "
            "tolerance refute residual-only validation.",
            ("truth separate from fit", "noise model declared for BIC"),
            "Finite synthetic designs, seeds and noise assumptions in Phase23/24 "
            "reports.",
            (
                "Inference code is not a physical model; weak/close modes and sparse "
                "data affect identifiability.",
            ),
            (
                "Structured run artifacts are absent; this synthesis uses report "
                "records and inspected source.",
            ),
            (
                "src/newton_lab/modal_estimation.py",
                "src/newton_lab/modal_identification.py",
                "docs/research/phase_23_modal_estimation.md",
                "docs/research/phase_24_mode_identification.md",
            ),
        ),
    )


def heat_signed_aggregate_counterexample() -> tuple[
    tuple[float, ...], tuple[float, ...]
]:
    """Evaluate Phase 22's signed heat mixture with Phase 19's analytic API.

    This does not rerun the steady BVP or reconstruct values from report text.
    """
    model = HeatDiffusionModel(
        length_m=1.0,
        left_temperature_k=300.0,
        right_temperature_k=400.0,
        diffusivity_m2_per_s=1e-4,
        modes=(
            SpatialMode(mode_number=1, amplitude=0.5),
            SpatialMode(mode_number=2, amplitude=-1.0),
        ),
    )
    times = (0.0, 100.0, 300.0, 400.0, 600.0)
    position = 0.25
    baseline = steady_temperature_profile(model, (position,))[0]
    evaluation = evaluate_heat_diffusion(model, (position, 0.75), times)
    values = tuple(row[0] - baseline for row in evaluation.temperature_k)
    return times, values


def render_cross_system_report(
    matrix: tuple[ComparisonEntry, ...],
    findings: tuple[CrossSystemFinding, ...],
    trajectory: TrajectoryComparison | None = None,
) -> str:
    """Render stable Markdown from existing records; no simulations are run."""
    lines = [
        "# Phase 27 — Cross-system structural discovery and falsification",
        "",
        "## Executive summary",
        "",
        "The undamped linear oscillator and small-angle pendulum have "
        "identical dimensionless equations under matched scaling and initial "
        "states. That is mathematical, not physical, equivalence. Heat modes "
        "and oscillator envelopes share conditional exponential forms but "
        "not spectra or observables. Phase25–26 reject a universal "
        "small-angle threshold. The steady heat BVP is the transient model's "
        "equilibrium only under matched assumptions. Estimation results are "
        "methods evidence, not new physical laws. No application or "
        "financial performance is tested.",
        "",
        "## Questions",
        "",
        "H1: dimensionless oscillator/pendulum correspondence. H2: modal "
        "rates and inference. H3: nonlinear approximation limits. H4: steady "
        "versus transient heat. H5: estimation reliability.",
        "",
        "## Source inventory",
        "",
        "Phases 18–26 source modules, tests, and relevant checked-in reports "
        "were inspected. Phase18/19 models are available from source APIs; "
        "no Phase21 serialized results were found. The expected "
        "the expected Phase21 JSON artifact directory "
        "`.newton_lab/artifacts` directory is absent. Phase20 integration, "
        "registry, workflow, and Phase21 schemas (`knowledge/integration.py`, "
        "`knowledge/registry.py`, `knowledge/workflow.py`, "
        "`knowledge/artifacts.py`) were inspected; this module "
        "does not mutate them. Phase22–26 computed values below come from "
        "their checked-in generated reports. The focused H1 solver check is "
        "calculated here from existing solvers.",
        "",
        "## Comparison matrix",
        "",
        "| Phase / role | Model | Assumptions, states and units | Parameters "
        "/ scales | Observables / evidence | Limits |",
        "|---|---|---|---|---|---|",
    ]
    for row in matrix:
        cells = (
            f"Phase {row.phase}: {row.name} ({row.role.value}; `{row.study_id}`)",
            row.equation,
            "; ".join((*row.assumptions, *row.states_units)),
            "; ".join((*row.parameters_units, row.timescale, row.dimensionless_groups)),
            f"{'; '.join(row.observables)}. {row.evidence_status.value}: "
            f"{row.numerical_evidence}",
            row.approximation_limits + " " + " ".join(row.limitations),
        )
        lines.append("| " + " | ".join(cells) + " |")
    lines.extend(("", "## Findings and falsification", ""))
    for finding in findings:
        lines.extend(
            (
                f"### {finding.finding_id}: {finding.status.value}",
                "",
                finding.hypothesis,
                "",
                f"**Support:** {finding.support}",
                "",
                f"**Counterexample:** {finding.counterexample}",
                "",
                f"**Assumptions:** {'; '.join(finding.assumptions)}",
                "",
                f"**Tested range:** {finding.tested_range}",
                "",
                "**Evidence:**",
                "",
            )
        )
        lines.extend(f"- {fact}" for fact in finding.evidence)
        lines.extend(("", "**Untested:**", ""))
        lines.extend(f"- {fact}" for fact in finding.untested)
        lines.extend(
            ("", "**Sources:** " + ", ".join(f"`{src}`" for src in finding.sources), "")
        )
    if trajectory is not None:
        lines.extend(
            (
                "## Dimensionless trajectory experiment",
                "",
                "Mapping: m=1 kg, k=1 N/m, c=0; L=1 m, g=1 m/s²; x/(1 m) maps "
                "numerically to θ/(1 rad). Thus both natural frequencies are 1 s⁻¹ "
                "and τ=t. This is an assumed coordinate mapping, not shared physical "
                "meaning.",
                f"Initial state: theta0={trajectory.initial_angle_rad:g} rad, "
                f"theta_dot0={trajectory.initial_angular_velocity_rad_per_s:g} rad/s; "
                f"tau=0 to {trajectory.duration_tau:g}; "
                f"{trajectory.sample_count} samples; "
                f"DOP853 rtol={trajectory.relative_tolerance:g}, "
                f"atol={trajectory.absolute_tolerance:g}.",
                "Oscillator vs exact linearized pendulum normalized max error: "
                f"{trajectory.linear_normalized_error:.8g}. Solver error only.",
                "Oscillator vs full-sine pendulum normalized max deviation: "
                f"{trajectory.nonlinear_normalized_error:.8g}; scale="
                f"{trajectory.normalization_scale:.8g}. Finite-angle difference, "
                "not physical validation.",
                "",
            )
        )
    lines.extend(
        (
            "## Reproducibility and limitations",
            "",
            "Run `build_cross_system_report()` to compute the focused ODE "
            "comparison and render the report. `render_cross_system_report()` "
            "itself performs no simulations. Report-only historical values are "
            "not presented as structured artifacts. Phase19 uses an analytic "
            "finite-mode solution without an independent transient PDE solver; "
            "Phase23/24 run-level output is not persisted. Finite grids do not "
            "establish universal behavior. Smoothing is not denoising or "
            "prediction.",
            "",
            "## Next research questions",
            "",
            "1. Independently integrate the linearized pendulum over multiple "
            "matched states and solver settings.",
            "2. Cross-check Phase19 finite-mode diffusion against a refined "
            "transient PDE discretization.",
            "3. Preserve structured Phase23/24 result records to query recovery "
            "and failure cases directly.",
            "",
            "No real-world control, predictive, financial, or trading claim is "
            "established.",
            "",
        )
    )
    return "\n".join(lines)


def build_cross_system_report() -> CrossSystemReport:
    """Compute the bounded comparison and return deterministic report text."""
    matrix = build_comparison_matrix()
    trajectory = compare_linear_oscillator_and_pendulum()
    findings = build_findings(trajectory)
    return CrossSystemReport(
        matrix,
        findings,
        trajectory,
        render_cross_system_report(matrix, findings, trajectory),
    )

"""Scientific meaning and provenance tests for Phase 27."""

from __future__ import annotations

import pytest

from newton_lab.cross_system_discovery import (
    EvidenceStatus,
    StudyRole,
    assess_observable_compatibility,
    build_comparison_matrix,
    build_cross_system_report,
    build_findings,
    compare_linear_oscillator_and_pendulum,
    heat_signed_aggregate_counterexample,
    normalized_discrepancy,
    oscillator_dimensionless_parameters,
    record_finding,
)
from newton_lab.exceptions import ScientificValidationError


def test_matrix_covers_required_models_and_inference_separately() -> None:
    rows = build_comparison_matrix()
    assert len(rows) == 8
    assert {row.phase for row in rows} >= {18, 19, 22, 23, 24, 25, 26}
    assert sum(row.role is StudyRole.INFERENCE for row in rows) == 2
    assert sum(row.role is StudyRole.PHYSICAL_MODEL for row in rows) == 6


def test_matrix_retains_provenance_and_missing_artifact_limitation() -> None:
    rows = {row.study_id: row for row in build_comparison_matrix()}
    assert "src/newton_lab/dynamics.py" in rows["phase18_damping_transfer"].sources
    assert (
        "docs/research/phase_26_pendulum_initial_velocity.md"
        in rows["phase26_nonlinear_pendulum"].sources
    )
    assert any(
        "Phase 21" in item for item in rows["phase18_damping_transfer"].limitations
    )


def test_modal_fit_is_not_mislabeled_as_physical_model() -> None:
    row = next(item for item in build_comparison_matrix() if item.phase == 24)
    assert row.role is StudyRole.INFERENCE
    assert "fit" in row.equation.lower()


def test_dimensionless_oscillator_coefficients() -> None:
    omega, zeta = oscillator_dimensionless_parameters(2.0, 6.0, 18.0)
    assert omega == pytest.approx(3.0)
    assert zeta == pytest.approx(0.5)
    assert 6.0 / (2.0 * omega) == pytest.approx(2.0 * zeta)
    assert 18.0 / (2.0 * omega**2) == pytest.approx(1.0)


@pytest.mark.parametrize("values", [(0.0, 0.0, 1.0), (1.0, -1.0, 1.0)])
def test_dimensionless_transform_rejects_invalid_parameters(
    values: tuple[float, float, float],
) -> None:
    with pytest.raises(ScientificValidationError):
        oscillator_dimensionless_parameters(*values)


def test_linearized_trajectory_agrees_and_nonlinear_limit_is_separate() -> None:
    result = compare_linear_oscillator_and_pendulum(
        initial_angle_rad=1e-4,
        initial_angular_velocity_rad_per_s=2e-4,
        duration_tau=2.0,
        sample_count=201,
    )
    assert result.linear_normalized_error < 1e-8
    assert result.nonlinear_normalized_error < 1e-8
    assert result.mapping_assumptions


def test_finite_amplitude_full_sine_pendulum_departs_from_harmonic_model() -> None:
    result = compare_linear_oscillator_and_pendulum(
        initial_angle_rad=0.8, duration_tau=4.0 * 3.141592653589793
    )
    assert result.linear_normalized_error < 1e-8
    assert result.nonlinear_normalized_error > 1e-3


def test_identical_units_do_not_imply_identical_observable() -> None:
    result = assess_observable_compatibility(
        left_unit="K",
        right_unit="K",
        left_meaning="temperature",
        right_meaning="modal coefficient",
    )
    assert not result.compatible
    assert "not independently validated" in result.status


def test_observable_assessment_requires_units_and_meanings_to_match() -> None:
    result = assess_observable_compatibility(
        left_unit="m",
        right_unit="rad",
        left_meaning="coordinate",
        right_meaning="coordinate",
    )
    assert not result.compatible


def test_math_equivalence_does_not_claim_physical_equivalence() -> None:
    h1 = build_findings(compare_linear_oscillator_and_pendulum())[0]
    assert h1.status is EvidenceStatus.ESTABLISHED
    assert "physical meanings do not" in " ".join(
        compare_linear_oscillator_and_pendulum().mapping_assumptions
    )


def test_modal_relationship_is_conditional_and_has_counterexample() -> None:
    h2 = build_findings(compare_linear_oscillator_and_pendulum())[1]
    assert h2.status is EvidenceStatus.ANALOGY
    assert "crosses zero" in h2.counterexample
    assert any("lambda1=" in item for item in h2.evidence)


def test_signed_heat_mode_mixture_can_cross_and_reverse_aggregate_magnitude() -> None:
    times, values = heat_signed_aggregate_counterexample()
    assert len(times) == len(values) == 5
    assert values[0] < 0.0 < values[3]
    assert abs(values[4]) > abs(values[3])


def test_candidate_contradicted_and_inconclusive_classifications_are_recordable() -> (
    None
):
    for status in (
        EvidenceStatus.CANDIDATE,
        EvidenceStatus.CONTRADICTED,
        EvidenceStatus.INCONCLUSIVE,
    ):
        finding = record_finding(
            finding_id="test",
            hypothesis="bounded hypothesis",
            status=status,
            support="support criterion",
            counterexample="refutation criterion",
            assumptions=("stated assumption",),
            tested_range="declared finite range",
            evidence=(),
            untested=("outside range",),
            sources=("source.py",),
        )
        assert finding.status is status
        assert finding.tested_range == "declared finite range"


def test_findings_declare_assumptions_range_and_provenance() -> None:
    result = compare_linear_oscillator_and_pendulum()
    finding = build_findings(result)[2]
    assert "0.35 rad" in finding.counterexample
    assert "Phase25" in finding.tested_range
    assert finding.assumptions and finding.sources


def test_counterexample_provenance_is_preserved() -> None:
    finding = build_findings(compare_linear_oscillator_and_pendulum())[2]
    assert "docs/research/phase_25_pendulum_approximation.md" in finding.sources
    assert "docs/research/phase_26_pendulum_initial_velocity.md" in finding.sources


def test_tested_universal_angle_threshold_candidate_is_contradicted() -> None:
    finding = next(
        item
        for item in build_findings(compare_linear_oscillator_and_pendulum())
        if item.finding_id == "H3b_universal_threshold"
    )
    assert finding.status is EvidenceStatus.CONTRADICTED
    assert "0.2 rad" in " ".join(finding.evidence)


def test_normalized_error_uses_declared_denominator() -> None:
    assert normalized_discrepancy(
        (2.0, 0.0), (1.0, 1.0), normalization_scale=2.0
    ) == pytest.approx(0.5)


@pytest.mark.parametrize("scale", [0.0, 1e-320, -1.0])
def test_normalized_error_rejects_zero_or_unsafe_denominator(scale: float) -> None:
    with pytest.raises(ScientificValidationError):
        normalized_discrepancy((0.0,), (0.0,), normalization_scale=scale)


def test_normalized_error_rejects_missing_or_misaligned_data() -> None:
    with pytest.raises(ScientificValidationError):
        normalized_discrepancy((), (), normalization_scale=1.0)
    with pytest.raises(ScientificValidationError):
        normalized_discrepancy((1.0,), (1.0, 2.0), normalization_scale=1.0)


def test_caller_mapping_is_kept_as_unverified_assumption() -> None:
    finding = record_finding(
        finding_id="mapping",
        hypothesis="mapping test",
        status=EvidenceStatus.CANDIDATE,
        support="if",
        counterexample="unless",
        assumptions=("source assumptions",),
        tested_range="none",
        evidence=(),
        untested=("all",),
        sources=("source.py",),
        mapping_declaration="map q to price deviation",
    )
    assert any("not independently validated" in item for item in finding.assumptions)


def test_failed_or_unreliable_inference_scope_is_not_hidden() -> None:
    row = next(item for item in build_comparison_matrix() if item.phase == 24)
    assert "low residual alone is not rate recovery" in row.approximation_limits
    assert "Finite synthetic design" in row.limitations[0]


def test_empty_compatibility_or_finding_inputs_are_rejected() -> None:
    with pytest.raises(ScientificValidationError):
        assess_observable_compatibility(
            left_unit="", right_unit="m", left_meaning="x", right_meaning="x"
        )
    with pytest.raises(ScientificValidationError):
        record_finding(
            finding_id="",
            hypothesis="",
            status=EvidenceStatus.CANDIDATE,
            support="",
            counterexample="",
            assumptions=(),
            tested_range="",
            evidence=(),
            untested=(),
            sources=(),
        )


def test_report_is_deterministic_and_preserves_scope_disclaimers() -> None:
    first = build_cross_system_report()
    second = build_cross_system_report()
    assert first.markdown == second.markdown
    assert first.trajectory == second.trajectory
    assert "No application or financial performance is tested" in first.markdown


def test_phase18_to_26_source_api_matrix_is_read_only() -> None:
    rows = build_comparison_matrix()
    assert {18, 19, 22, 23, 24, 25, 26}.issubset({row.phase for row in rows})
    assert all(row.equation and row.analytical_reference for row in rows)

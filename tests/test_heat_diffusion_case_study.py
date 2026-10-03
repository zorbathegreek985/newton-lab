"""Deterministic checks for the bounded Phase 19 analytical case study."""

from __future__ import annotations

import math

import pytest
from pydantic import ValidationError

from newton_lab.exceptions import ScientificValidationError
from newton_lab.heat_diffusion_case_study import (
    HeatDiffusionModel,
    SpatialMode,
    evaluate_heat_diffusion,
    run_heat_diffusion_smoothing_case_study,
    steady_temperature_profile,
)


def _model() -> HeatDiffusionModel:
    return HeatDiffusionModel(
        length_m=2.0,
        left_temperature_k=10.0,
        right_temperature_k=20.0,
        diffusivity_m2_per_s=0.25,
        modes=(
            SpatialMode(mode_number=1, amplitude=2.0),
            SpatialMode(mode_number=3, amplitude=-0.5),
        ),
    )


def test_steady_profile_is_the_linear_fixed_boundary_solution() -> None:
    model = _model()

    assert steady_temperature_profile(model, (0.0, 0.5, 1.0, 2.0)) == pytest.approx(
        (10.0, 12.5, 15.0, 20.0)
    )


def test_transient_modal_solution_matches_closed_form_and_boundaries() -> None:
    model = _model()
    positions = (0.0, 0.5, 1.0, 2.0)
    result = evaluate_heat_diffusion(model, positions, (0.0, 0.4))
    base = steady_temperature_profile(model, positions)
    expected_initial = tuple(
        base[index]
        + 2.0 * math.sin(math.pi * x / 2.0)
        - 0.5 * math.sin(3.0 * math.pi * x / 2.0)
        for index, x in enumerate(positions)
    )
    assert result.temperature_k[0] == pytest.approx(expected_initial)
    assert result.temperature_k[1][0] == pytest.approx(10.0)
    assert result.temperature_k[1][-1] == pytest.approx(20.0)
    assert result.modal_amplitudes[0] == pytest.approx((2.0, -0.5))
    assert result.modal_amplitudes[1][0] == pytest.approx(
        2.0 * math.exp(-0.25 * (math.pi / 2.0) ** 2 * 0.4)
    )
    assert result.modal_amplitudes[1][1] == pytest.approx(
        -0.5 * math.exp(-0.25 * (3.0 * math.pi / 2.0) ** 2 * 0.4)
    )


def test_higher_frequency_attenuates_faster_and_variation_falls() -> None:
    model = _model()
    result = evaluate_heat_diffusion(model, (0.0, 1.0, 2.0), (0.0, 0.5, 1.0))

    assert result.attenuation_factors[1][1] < result.attenuation_factors[1][0]
    assert (
        result.perturbation_rms_gradient_k_per_m[2]
        < result.perturbation_rms_gradient_k_per_m[1]
        < result.perturbation_rms_gradient_k_per_m[0]
    )
    expected_initial = math.sqrt(
        ((2.0 * math.pi / 2.0) ** 2 + (-0.5 * 3.0 * math.pi / 2.0) ** 2) / 2.0
    )
    assert result.perturbation_rms_gradient_k_per_m[0] == pytest.approx(
        expected_initial
    )


@pytest.mark.parametrize(
    "updates",
    [
        {"length_m": 0.0},
        {"length_m": math.inf},
        {"diffusivity_m2_per_s": -1.0},
        {"diffusivity_m2_per_s": math.nan},
        {"left_temperature_k": math.inf},
        {"modes": ()},
        {"modes": (SpatialMode(mode_number=1, amplitude=1.0),) * 2},
    ],
)
def test_invalid_model_parameters_are_rejected(updates: dict[str, object]) -> None:
    values: dict[str, object] = {
        "length_m": 1.0,
        "left_temperature_k": 0.0,
        "right_temperature_k": 1.0,
        "diffusivity_m2_per_s": 1.0,
        "modes": (SpatialMode(mode_number=1, amplitude=1.0),),
    }
    values.update(updates)
    with pytest.raises(ValidationError):
        HeatDiffusionModel.model_validate(values)


@pytest.mark.parametrize(
    ("positions", "times"),
    [
        ((0.0,), (0.0,)),
        ((-0.1, 0.5), (0.0,)),
        ((0.0, 2.1), (0.0,)),
        ((0.0, 1.0), (-0.1,)),
        ((0.0, math.nan), (0.0,)),
        ((0.0, 1.0), (math.inf,)),
    ],
)
def test_invalid_evaluation_grids_are_rejected(
    positions: tuple[float, ...], times: tuple[float, ...]
) -> None:
    with pytest.raises(ScientificValidationError):
        evaluate_heat_diffusion(_model(), positions, times)


def test_synthetic_signal_uses_same_attenuation_and_is_reproducible() -> None:
    first = run_heat_diffusion_smoothing_case_study(sample_count=21)
    second = run_heat_diffusion_smoothing_case_study(sample_count=21)

    assert first.report_markdown == second.report_markdown
    assert first.evaluation == second.evaluation
    assert first.synthetic_signal == second.synthetic_signal
    assert first.synthetic_signal.mode_numbers == (1, 4)
    assert first.synthetic_signal.modal_amplitudes[0] == pytest.approx((1.0, 0.4))
    assert first.synthetic_signal.modal_amplitudes[-1][1] < (
        first.synthetic_signal.modal_amplitudes[-1][0] * 0.4
    )
    assert len(first.synthetic_signal.sampled_values[0]) == 21


def test_source_bvp_provenance_is_retained_without_invented_equation_match() -> None:
    case = run_heat_diffusion_smoothing_case_study(sample_count=21)
    source = case.source_experiment

    assert case.source_bvp_numerically_accepted
    assert case.source_bvp_run_index == 0
    assert case.source_bvp_max_profile_error_k < 1e-9
    assert case.source_equation_record_ids == ()
    assert source.specification.bvp_problem is not None
    assert (
        source.specification.bvp_problem.problem_id == "steady_linear_heat_conduction"
    )
    assumptions = source.specification.bvp_problem.assumptions
    assert "Steady state." in assumptions
    assert "No internal heat generation." in assumptions
    assert any("constant thermal conductivity" in item for item in assumptions)
    assert not case.discovery.structure_matches
    assert not case.discovery.hypotheses
    assert case.discovery.observations
    assert all(not item.equation_record_ids for item in case.discovery.observations)


def test_evidence_report_names_scope_and_does_not_claim_application_validation() -> (
    None
):
    case = run_heat_diffusion_smoothing_case_study(sample_count=21)
    report = case.report_markdown.lower()

    assert case.evidence_classification == "analytical_model_result"
    assert case.cross_domain_status == "synthetic_structural_analogy_only"
    assert case.real_world_or_financial_validation is False
    for heading in (
        "research question",
        "existing steady heat-conduction model",
        "transient diffusion model",
        "modal attenuation experiment",
        "synthetic signal comparison",
        "newton lab discovery and evidence",
        "limitations",
        "tests that could challenge the connection",
    ):
        assert heading in report
    assert "does not establish better" in report
    assert "no equation id" in report
    assert "numerical_accepted=true" in report


def test_source_mismatch_fails_without_changing_the_source_record() -> None:
    source = run_heat_diffusion_smoothing_case_study(sample_count=11).source_experiment
    before = source.model_copy(deep=True)

    with pytest.raises(ScientificValidationError, match="conditions must match"):
        run_heat_diffusion_smoothing_case_study(
            source_record=source, length_m=2.0, sample_count=11
        )

    assert source == before

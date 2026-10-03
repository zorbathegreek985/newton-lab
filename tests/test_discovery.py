"""End-to-end evidence-aware experiment discovery tests."""

import pytest

from newton_lab.algebraic import AlgebraicSolverConfiguration
from newton_lab.discovery import (
    DiscoveryEvidenceStatus,
    ObservationAvailability,
    ObservationBasis,
    run_discovery,
)
from newton_lab.experiments import (
    ExperimentRecord,
    make_heat_conduction_experiment,
    make_nonlinear_spring_experiment,
    make_oscillator_experiment,
    make_pendulum_experiment,
    run_experiment,
)
from newton_lab.knowledge import build_example_knowledge_base
from newton_lab.simulation import OutputSampling


def _damping_record() -> ExperimentRecord:
    return run_experiment(
        make_oscillator_experiment(
            experiment_id="discovery_damping_study",
            sweep_values=(0.0, 0.1, 0.2),
            metrics=("maximum_absolute_velocity_m_per_s",),
        )
    )


def test_end_to_end_discovery_preserves_evidence_and_is_deterministic() -> None:
    record = _damping_record()
    before = record.model_dump(mode="python")

    result = run_discovery((record,))
    repeated = run_discovery((record,))

    run_observations = [
        item for item in result.observations if item.descriptor == "recorded_metric"
    ]
    change = next(
        item
        for item in result.observations
        if item.descriptor == "sweep_endpoint_change"
    )
    assert len(run_observations) == 3
    assert all(item.basis == ObservationBasis.COMPUTED for item in run_observations)
    assert all(
        item.equation_record_ids == ("damped_harmonic_oscillator",)
        for item in run_observations
    )
    assert all(item.source_record_index == 0 for item in run_observations)
    assert all(item.solver_success is True for item in run_observations)
    assert all(
        "returned samples" in " ".join(item.limitations) for item in run_observations
    )
    assert change.availability == ObservationAvailability.AVAILABLE
    assert change.basis == ObservationBasis.DERIVED
    assert change.source_run_ids == tuple(
        f"0:discovery_damping_study#{index}" for index in range(3)
    )
    assert change.swept_parameter_values == (0.0, 0.1, 0.2)
    assert change.swept_parameter_unit == "kg/s"
    first_value = run_observations[0].value
    last_value = run_observations[-1].value
    assert first_value is not None and last_value is not None
    assert change.value == pytest.approx(last_value - first_value)
    assert change.value is not None and change.value < 0.0

    established = next(
        item
        for item in result.structure_matches
        if item.target_equation_id == "newton_second_law"
    )
    assert established.evidence_status == (
        DiscoveryEvidenceStatus.ESTABLISHED_RELATIONSHIP
    )
    assert established.relationship_type.value == "derived_from"
    assert established.relationship_source_equation_id == ("damped_harmonic_oscillator")
    assert established.relationship_target_equation_id == "newton_second_law"

    match = next(
        item
        for item in result.structure_matches
        if item.target_equation_id == "damped_nonlinear_pendulum"
    )
    assert match.evidence_status == (
        DiscoveryEvidenceStatus.CANDIDATE_STRUCTURAL_CONNECTION
    )
    assert match.relationship_id == "small_angle_oscillator_approximation"
    assert match.relationship_review_status is not None
    assert match.relationship_review_status.value == "candidate_requires_review"
    assert any(item.startswith("ode_order=") for item in match.shared_features)
    assert any(
        item.startswith("derivative_operators=") for item in match.shared_features
    )
    assert match.conditions
    assert match.supporting_observation_ids

    assert len(result.hypotheses) == 1
    hypothesis = result.hypotheses[0]
    assert hypothesis.evidence_status == (
        DiscoveryEvidenceStatus.PROPOSED_APPLICATION_HYPOTHESIS
    )
    assert hypothesis.validated_application is False
    assert hypothesis.target_domain == "quantitative_finance_risk_control_research"
    assert hypothesis.test_plan
    assert any("transaction costs" in item for item in hypothesis.test_plan)
    assert any("nonstationary" in item for item in hypothesis.differences)
    assert "not validated" in result.report_markdown
    assert "candidate_requires_review" in result.report_markdown
    assert result.report_markdown == repeated.report_markdown
    assert record.model_dump(mode="python") == before


def test_failed_runs_remain_unavailable_and_are_not_zero_filled() -> None:
    record = run_experiment(
        make_oscillator_experiment(
            experiment_id="discovery_failed_sweep",
            sweep_values=(0.0, -0.1, 0.2),
            metrics=("maximum_absolute_velocity_m_per_s",),
        )
    )

    result = run_discovery((record,))

    failed = next(
        item
        for item in result.observations
        if item.descriptor == "recorded_metric" and "#1" in item.source_run_ids[0]
    )
    endpoint_change = next(
        item
        for item in result.observations
        if item.descriptor == "sweep_endpoint_change"
    )
    assert failed.availability == ObservationAvailability.FAILED_RUN
    assert failed.value is None
    assert failed.basis == ObservationBasis.UNAVAILABLE
    assert endpoint_change.value is not None
    assert "One or more sweep runs failed" in " ".join(endpoint_change.limitations)
    assert not result.hypotheses
    assert "| failed |" in result.report_markdown


def test_unreliable_sweep_metadata_yields_unavailable_derived_descriptor() -> None:
    record = run_experiment(
        make_oscillator_experiment(
            experiment_id="discovery_unreliable_sweep",
            sweep_values=(0.0, 0.1, 0.2),
            metrics=("maximum_absolute_velocity_m_per_s",),
        )
    )
    changed_run = record.runs[1].model_copy(update={"swept_parameter_value": 0.15})
    unreliable = record.model_copy(
        update={"runs": (record.runs[0], changed_run, record.runs[2])}
    )

    result = run_discovery((unreliable,))
    change = next(
        item
        for item in result.observations
        if item.descriptor == "sweep_endpoint_change"
    )

    assert change.availability == ObservationAvailability.UNRELIABLE_SWEEP
    assert change.value is None
    assert not result.hypotheses
    assert "unreliable" in " ".join(change.limitations).lower()


def test_algebraic_and_bvp_metrics_are_described_without_guessed_structures() -> None:
    algebraic = run_experiment(
        make_nonlinear_spring_experiment(
            experiment_id="discovery_algebraic",
            sweep_values=None,
            metrics=("residual_norm",),
            solver_configuration=AlgebraicSolverConfiguration(),
        )
    )
    bvp = run_experiment(
        make_heat_conduction_experiment(
            experiment_id="discovery_bvp",
            sweep_values=None,
            metrics=("right_endpoint:temperature",),
        )
    )

    result = run_discovery((algebraic, bvp))

    assert result.collection_summary.counts_by_family
    assert {item.model_id for item in result.observations} == {
        "nonlinear_algebraic_system",
        "boundary_value_problem",
    }
    assert all(not item.equation_record_ids for item in result.observations)
    assert all(item.value is not None for item in result.observations)
    assert algebraic.runs[0].algebraic_diagnostics is not None
    assert algebraic.runs[0].ode_trajectory is None
    assert bvp.runs[0].bvp_diagnostics is not None
    assert bvp.runs[0].ode_trajectory is None
    assert not result.structure_matches
    assert not result.hypotheses
    assert "did not guess from names" in " ".join(result.limitations)
    assert "No application is validated" in result.report_markdown


def test_empty_and_unmatched_metric_inputs_report_insufficient_evidence() -> None:
    empty = run_discovery(())
    record = run_experiment(
        make_oscillator_experiment(
            experiment_id="discovery_selected_metric",
            sweep_values=None,
            metrics=("final_displacement_m",),
        )
    )
    unmatched = run_discovery((record,), selected_metrics=("unknown_metric",))

    assert empty.observations == ()
    assert empty.structure_matches == ()
    assert empty.hypotheses == ()
    assert "No experiment records" in " ".join(empty.limitations)
    assert "No metric observations" in empty.report_markdown
    assert not unmatched.observations
    assert not unmatched.hypotheses
    assert "not requested" in " ".join(unmatched.limitations)


def test_explicit_equation_provenance_is_used_and_unknown_ids_are_not_guessed() -> None:
    original = _damping_record()
    spec = original.specification.model_copy(
        update={"equation_record_ids": ("missing_equation_record",)}
    )
    declared_unknown = original.model_copy(update={"specification": spec})

    result = run_discovery(
        (declared_unknown,), knowledge_base=build_example_knowledge_base()
    )

    assert all(not item.equation_record_ids for item in result.observations)
    assert not result.structure_matches
    assert not result.hypotheses
    assert "Unrecognized equation record IDs" in " ".join(result.limitations)


def test_retained_trajectory_yields_qualified_period_and_amplitude_descriptors() -> (
    None
):
    record = run_experiment(
        make_oscillator_experiment(
            experiment_id="temporal_descriptors",
            sweep_values=None,
            duration_s=10.0,
            output_sampling=OutputSampling(point_count=2001),
            retain_ode_trajectory=True,
        )
    )
    result = run_discovery((record,))
    period = next(
        item for item in result.observations if item.descriptor == "observed_period"
    )
    amplitude = next(
        item for item in result.observations if item.descriptor == "amplitude_change"
    )
    return_indicator = next(
        item
        for item in result.observations
        if item.descriptor == "finite_window_return_indicator"
    )

    assert period.availability == ObservationAvailability.AVAILABLE
    assert period.unit == "s"
    assert period.value == pytest.approx(3.143, abs=0.02)
    assert period.sample_count == 2001
    assert period.time_interval_s == (0.0, 10.0)
    assert period.state_name == "displacement"
    assert amplitude.availability == ObservationAvailability.AVAILABLE
    assert amplitude.unit == "m"
    assert amplitude.value is not None and amplitude.value < 0.0
    assert "not a fitted decay law" in amplitude.method
    assert return_indicator.availability == ObservationAvailability.AVAILABLE
    assert return_indicator.unit == "1"
    assert return_indicator.value in (0.0, 1.0)


def test_temporal_descriptors_preserve_sweep_order_and_failed_run_provenance() -> None:
    record = run_experiment(
        make_oscillator_experiment(
            experiment_id="temporal_sweep",
            sweep_values=(0.0, -0.1, 0.2, 0.2),
            metrics=("final_displacement_m",),
            retain_ode_trajectory=True,
        )
    )
    result = run_discovery((record,))
    periods = [
        item for item in result.observations if item.descriptor == "observed_period"
    ]
    assert [item.swept_parameter_values for item in periods] == [
        (0.0,),
        (-0.1,),
        (0.2,),
        (0.2,),
    ]
    assert periods[1].availability == ObservationAvailability.FAILED_RUN
    assert periods[1].value is None
    assert periods[0].availability == ObservationAvailability.AVAILABLE


def test_ode_temporal_descriptors_are_explicitly_unavailable_without_retention() -> (
    None
):
    result = run_discovery((_damping_record(),))
    periods = [
        item for item in result.observations if item.descriptor == "observed_period"
    ]
    assert len(periods) == 3
    assert all(
        item.availability == ObservationAvailability.INSUFFICIENT_SAMPLES
        for item in periods
    )
    assert all(item.value is None for item in periods)


def test_nonlinear_pendulum_trajectory_uses_angle_units() -> None:
    record = run_experiment(
        make_pendulum_experiment(
            experiment_id="pendulum_temporal",
            sweep_values=None,
            duration_s=8.0,
            output_sampling=OutputSampling(point_count=1601),
            retain_ode_trajectory=True,
        )
    )
    result = run_discovery((record,))
    amplitude = next(
        item for item in result.observations if item.descriptor == "amplitude_change"
    )
    assert amplitude.state_name == "angle"
    assert amplitude.unit == "rad"
    assert amplitude.availability == ObservationAvailability.AVAILABLE

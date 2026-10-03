"""Tests for read-only experiment summaries and cross-family analysis."""

import pytest

from newton_lab.algebraic import AlgebraicSolverConfiguration
from newton_lab.experiment_analysis import (
    MetricCompatibilityDeclaration,
    MetricSourceDeclaration,
    analyze_parameter_sweep,
    compare_metric_observations,
    extract_metric_observations,
    project_experiment_runs,
    summarize_experiments,
    summarize_metric,
)
from newton_lab.experiments import (
    ExperimentRunStatus,
    MetricValue,
    make_heat_conduction_experiment,
    make_nonlinear_spring_experiment,
    make_oscillator_experiment,
    run_experiment,
)


def test_empty_collection_summary_is_well_defined() -> None:
    summary = summarize_experiments(())

    assert summary.record_count == 0
    assert summary.run_count == 0
    assert summary.successful_run_count == 0
    assert summary.failed_run_count == 0
    assert summary.counts_by_family == ()
    assert summary.duplicate_experiment_ids == ()


def test_collection_counts_status_acceptance_duplicates_and_provenance() -> None:
    specification = make_oscillator_experiment(
        experiment_id="analysis_sweep",
        sweep_values=(0.0, -0.1, 0.0),
        metrics=("final_displacement_m",),
    )
    record = run_experiment(specification)
    before = record.model_dump(mode="python")

    summary = summarize_experiments((record, record))
    projection = extract_metric_observations((record,), "final_displacement_m")

    assert summary.record_count == 2
    assert summary.run_count == 6
    assert summary.successful_run_count == 4
    assert summary.failed_run_count == 2
    assert summary.duplicate_experiment_ids[0].experiment_id == "analysis_sweep"
    assert summary.duplicate_experiment_ids[0].record_indices == (0, 1)
    assert [(item.category, item.count) for item in summary.counts_by_family] == [
        ("ode", 6)
    ]
    assert summary.solver_success_count == 4
    assert summary.solver_failure_count == 0
    assert summary.solver_success_unavailable_count == 2
    assert summary.numerical_acceptance_unavailable_count == 6
    assert len(summary.parameter_configurations) == 2
    assert [config.run_count for config in summary.parameter_configurations] == [4, 2]
    assert [item.record_index for item in projection] == [0, 0, 0]
    assert projection[1].execution_status == ExperimentRunStatus.FAILED
    assert projection[0].reproducibility is not None
    assert projection[0].solver_id == "scipy_solve_ivp"
    assert projection[0].swept_parameter_name == "damping_coefficient_kg_per_s"
    assert record.model_dump(mode="python") == before


def test_missing_zero_and_nonfinite_metric_observations_remain_distinct() -> None:
    specification = make_oscillator_experiment(
        experiment_id="metric_integrity",
        sweep_values=None,
        metrics=("final_displacement_m",),
    )
    record = run_experiment(specification)
    run = record.runs[0]
    measured_zero = run.model_copy(
        update={
            "metrics": (
                MetricValue(metric_id="final_displacement_m", value=0.0, unit="m"),
            )
        }
    )
    zero_record = record.model_copy(update={"runs": (measured_zero,)})
    missing_record = record.model_copy(
        update={"runs": (run.model_copy(update={"metrics": ()}),)}
    )
    nonfinite_metric = MetricValue.model_construct(
        metric_id="final_displacement_m", value=float("nan"), unit="m"
    )
    nonfinite_record = record.model_copy(
        update={"runs": (run.model_copy(update={"metrics": (nonfinite_metric,)}),)}
    )

    zero_observation, missing_observation, nonfinite_observation = (
        extract_metric_observations(
            (zero_record, missing_record, nonfinite_record), "final_displacement_m"
        )
    )

    assert zero_observation.value == 0.0
    assert zero_observation.observation_status == "valid"
    assert missing_observation.value is None
    assert missing_observation.metric_status == "missing"
    assert nonfinite_observation.value is None
    assert nonfinite_observation.metric_status == "non_finite"
    comparison = compare_metric_observations(
        (zero_observation, missing_observation, nonfinite_observation)
    )
    stats = summarize_metric(comparison)
    assert stats.valid_count == 1
    assert stats.mean == 0.0
    assert stats.missing_metric_count == 1
    assert stats.non_finite_metric_count == 1

    duplicated_run = run.model_copy(update={"metrics": run.metrics * 2})
    invalid_record = record.model_copy(update={"runs": (duplicated_run,)})
    invalid_observation = extract_metric_observations(
        (invalid_record,), "final_displacement_m"
    )[0]
    assert invalid_observation.metric_status == "invalid"
    assert (
        compare_metric_observations((zero_observation, invalid_observation)).status
        == "not_comparable"
    )


def test_comparison_requires_declaration_across_different_physical_models() -> None:
    ode_record = run_experiment(
        make_oscillator_experiment(
            experiment_id="oscillator_observation",
            sweep_values=None,
            metrics=("final_displacement_m",),
        )
    )
    algebraic_record = run_experiment(
        make_nonlinear_spring_experiment(
            experiment_id="spring_observation",
            sweep_values=None,
            metrics=("solution:displacement",),
            solver_configuration=AlgebraicSolverConfiguration(),
        )
    )
    ode_observation = extract_metric_observations(
        (ode_record,), "final_displacement_m"
    )[0]
    algebraic_observation = extract_metric_observations(
        (algebraic_record,), "solution:displacement"
    )[0]
    observations = (ode_observation, algebraic_observation)

    automatic = compare_metric_observations(observations)
    explicit = compare_metric_observations(
        observations,
        compatibility=MetricCompatibilityDeclaration(
            shared_definition=(
                "Displacement in metres at the declared equilibrium/state."
            ),
            unit="m",
            sources=(
                MetricSourceDeclaration(
                    model_id=ode_observation.model_id,
                    problem_id=ode_observation.problem_id,
                    metric_id=ode_observation.metric_id,
                    unit="m",
                    source_definition=ode_observation.definition
                    or "ODE final displacement",
                ),
                MetricSourceDeclaration(
                    model_id=algebraic_observation.model_id,
                    problem_id=algebraic_observation.problem_id,
                    metric_id=algebraic_observation.metric_id,
                    unit="m",
                    source_definition=(
                        algebraic_observation.definition or "Algebraic displacement"
                    ),
                ),
            ),
            assumptions=(
                "The compared coordinates represent the same physical displacement.",
            ),
        ),
    )

    assert automatic.status == "not_comparable"
    assert explicit.status == "conditional"
    assert explicit.unit == "m"
    assert explicit.assumptions
    bad_units = MetricCompatibilityDeclaration.model_construct(
        shared_definition="incorrect mixed-unit mapping",
        unit="m",
        sources=(
            MetricSourceDeclaration(
                model_id=ode_observation.model_id,
                problem_id=ode_observation.problem_id,
                metric_id=ode_observation.metric_id,
                unit="m",
                source_definition=ode_observation.definition or "ODE displacement",
            ),
            MetricSourceDeclaration(
                model_id=algebraic_observation.model_id,
                problem_id=algebraic_observation.problem_id,
                metric_id=algebraic_observation.metric_id,
                unit="cm",
                source_definition=algebraic_observation.definition
                or "Spring displacement",
            ),
        ),
    )
    assert (
        compare_metric_observations(observations, compatibility=bad_units).status
        == "not_comparable"
    )


def test_parameter_sweep_analysis_preserves_order_duplicates_and_failures() -> None:
    record = run_experiment(
        make_oscillator_experiment(
            experiment_id="ordered_analysis_sweep",
            sweep_values=(0.0, -0.1, 0.0),
            metrics=("final_displacement_m",),
        )
    )

    analysis = analyze_parameter_sweep(
        record,
        parameter_name="damping_coefficient_kg_per_s",
        metric_id="final_displacement_m",
    )

    assert analysis.ordering_reliable
    assert [point.parameter_value for point in analysis.points] == [0.0, -0.1, 0.0]
    assert [point.repeated_parameter_value for point in analysis.points] == [
        False,
        False,
        True,
    ]
    assert [point.observation.execution_status for point in analysis.points] == [
        ExperimentRunStatus.SUCCEEDED,
        ExperimentRunStatus.FAILED,
        ExperimentRunStatus.SUCCEEDED,
    ]
    unavailable = analyze_parameter_sweep(
        record, parameter_name="stiffness_n_per_m", metric_id="final_displacement_m"
    )
    assert not unavailable.ordering_reliable
    assert (
        unavailable.points[1].observation.execution_status == ExperimentRunStatus.FAILED
    )


def test_statistics_use_population_standard_deviation_and_handle_small_samples() -> (
    None
):
    record = run_experiment(
        make_oscillator_experiment(
            experiment_id="stats_sweep",
            sweep_values=(0.0, 0.1),
            metrics=("final_displacement_m",),
        )
    )
    first, second = extract_metric_observations((record,), "final_displacement_m")
    first = first.model_copy(update={"value": 1.0})
    second = second.model_copy(update={"value": 3.0})
    comparison = compare_metric_observations((first, second))

    stats = summarize_metric(comparison)
    one = summarize_metric(comparison.model_copy(update={"observations": (first,)}))
    no_valid = summarize_metric(comparison.model_copy(update={"observations": ()}))

    assert stats.mean == pytest.approx(2.0)
    assert stats.population_standard_deviation == pytest.approx(1.0)
    assert stats.minimum == 1.0 and stats.maximum == 3.0
    assert one.mean == 1.0
    assert one.population_standard_deviation == 0.0
    assert one.status == "computed"
    assert no_valid.mean is None
    assert no_valid.status == "no_valid_observations"


def test_bvp_metric_provenance_and_family_specific_diagnostics_are_retained() -> None:
    record = run_experiment(
        make_heat_conduction_experiment(
            experiment_id="analysis_heat_case",
            sweep_values=None,
            metrics=("sampled_minimum:temperature",),
        )
    )
    observation = extract_metric_observations((record,), "sampled_minimum:temperature")[
        0
    ]
    summary = summarize_experiments((record,))

    assert observation.problem_family == "boundary_value"
    assert observation.unit == "K"
    assert observation.definition is not None
    assert observation.reproducibility is not None
    projection = project_experiment_runs((record,))[0]
    assert projection.bvp_diagnostics is not None
    assert observation.bvp_diagnostics is not None
    assert projection.algebraic_diagnostics is None
    assert summary.numerical_acceptance_count == 1
    assert summary.numerical_rejection_count == 0
    assert summary.numerical_acceptance_unavailable_count == 0

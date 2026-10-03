"""Tests for reusable simulation contracts and existing model adapters."""

from dataclasses import replace

import numpy as np
import pytest
from pydantic import ValidationError

from newton_lab.dynamics import simulate_damped_oscillator
from newton_lab.exceptions import IntegrationError, ScientificValidationError
from newton_lab.knowledge.examples import build_example_registry
from newton_lab.pendulum import simulate_damped_pendulum
from newton_lab.simulation import (
    KnowledgeMappingKind,
    KnowledgeMappingStatus,
    ODESolverConfiguration,
    OutputSampling,
    SimulationFailure,
    SimulationProblem,
    SimulationResult,
    SimulationStatus,
    validate_knowledge_mappings,
)
from newton_lab.simulation_adapters import (
    DampedOscillatorAdapter,
    DampedPendulumAdapter,
    make_damped_oscillator_problem,
    make_damped_pendulum_problem,
)


def test_problem_contract_validates_intervals_sampling_and_named_values() -> None:
    spatial_problem = SimulationProblem(
        model_id="spatial_field_example",
        name="Spatial field problem",
        family="boundary_value_pde",
        model_version="0.1",
        boundary_conditions=("field vanishes at boundary",),
    )
    assert spatial_problem.interval is None
    assert spatial_problem.independent_variable is None
    assert spatial_problem.states == ()
    with pytest.raises(ValidationError, match="finite increasing endpoints"):
        make_damped_oscillator_problem(
            mass_kg=1.0,
            damping_coefficient_kg_per_s=0.1,
            stiffness_n_per_m=1.0,
            initial_displacement_m=0.1,
            initial_velocity_m_per_s=0.0,
            duration_s=float("inf"),
        )
    with pytest.raises(ValidationError):
        OutputSampling(point_count=1)
    with pytest.raises(ValidationError):
        OutputSampling(point_count=True)
    with pytest.raises(ValidationError):
        ODESolverConfiguration(maximum_step=float("inf"))


def test_solver_configuration_and_capabilities_are_explicit() -> None:
    configuration = ODESolverConfiguration(
        method="DOP853", relative_tolerance=1e-9, maximum_step=0.1
    )
    assert configuration.maximum_step == 0.1
    with pytest.raises(ValidationError):
        ODESolverConfiguration.model_validate({"method": "not-a-solver"})
    assert "initial_value_ode" in DampedOscillatorAdapter.capabilities.problem_families
    assert DampedOscillatorAdapter.capabilities.adaptive_internal_steps
    assert "PDE" not in DampedOscillatorAdapter.capabilities.problem_families


def test_oscillator_adapter_preserves_legacy_values_and_adds_common_metadata() -> None:
    problem = make_damped_oscillator_problem(
        mass_kg=1.0,
        damping_coefficient_kg_per_s=0.2,
        stiffness_n_per_m=4.0,
        initial_displacement_m=0.1,
        initial_velocity_m_per_s=0.0,
        duration_s=2.0,
        point_count=101,
    )
    legacy = simulate_damped_oscillator(1.0, 0.2, 4.0, 0.1, 0.0, 2.0, num_points=101)
    result = DampedOscillatorAdapter().run(problem, ODESolverConfiguration())

    assert isinstance(result, SimulationResult)
    assert result.status == SimulationStatus.SUCCESS
    assert result.independent_variable_name == "time"
    assert result.states[0].name == "displacement"
    assert result.states[0].unit == "m"
    np.testing.assert_array_equal(result.independent_variable, legacy.time_s)
    np.testing.assert_array_equal(result.states[0].values, legacy.displacement_m)
    assert result.reproducibility.model_parameters == problem.parameters
    assert result.reproducibility.numpy_version
    assert result.diagnostics.function_evaluations == legacy.function_evaluations
    assert result.diagnostics.accepted_internal_steps is None


def test_pendulum_adapter_preserves_full_sine_solution_and_units() -> None:
    problem = make_damped_pendulum_problem(
        mass_kg=0.8,
        length_m=1.2,
        damping_coefficient_kg_m2_per_s=0.1,
        gravity_m_per_s2=9.81,
        initial_angle_rad=1.4,
        initial_angular_velocity_rad_per_s=-0.1,
        duration_s=1.5,
        point_count=101,
    )
    legacy = simulate_damped_pendulum(
        0.8, 1.2, 0.1, 9.81, 1.4, -0.1, 1.5, num_points=101
    )
    result = DampedPendulumAdapter().run(problem, ODESolverConfiguration())

    assert isinstance(result, SimulationResult)
    assert result.states[0].unit == "rad"
    np.testing.assert_array_equal(result.states[0].values, legacy.angle_rad)
    np.testing.assert_array_equal(
        result.states[1].values, legacy.angular_velocity_rad_per_s
    )
    assert "full nonlinear sine" in result.assumptions[-1]


def test_solver_settings_maximum_step_are_shared_with_legacy_functions() -> None:
    problem = make_damped_oscillator_problem(
        mass_kg=1.0,
        damping_coefficient_kg_per_s=0.2,
        stiffness_n_per_m=4.0,
        initial_displacement_m=0.1,
        initial_velocity_m_per_s=0.0,
        duration_s=2.0,
    )
    config = ODESolverConfiguration(maximum_step=0.05)
    common = DampedOscillatorAdapter().run(problem, config)
    legacy = simulate_damped_oscillator(
        1.0,
        0.2,
        4.0,
        0.1,
        0.0,
        2.0,
        maximum_step_s=0.05,
    )

    assert isinstance(common, SimulationResult)
    assert common.reproducibility.solver_maximum_step == 0.05
    np.testing.assert_array_equal(common.states[0].values, legacy.displacement_m)


def test_integration_failure_is_explicit_in_adapter_outcome(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    import newton_lab.simulation_adapters as adapters

    def fail(*_args: object, **_kwargs: object) -> None:
        raise IntegrationError("forced failure")

    monkeypatch.setattr(adapters, "simulate_damped_oscillator", fail)
    problem = make_damped_oscillator_problem(
        mass_kg=1.0,
        damping_coefficient_kg_per_s=0.2,
        stiffness_n_per_m=4.0,
        initial_displacement_m=0.1,
        initial_velocity_m_per_s=0.0,
        duration_s=1.0,
    )
    result = DampedOscillatorAdapter().run(problem, ODESolverConfiguration())

    assert isinstance(result, SimulationFailure)
    assert result.status == SimulationStatus.FAILED
    assert result.error_type == "IntegrationError"
    assert result.message == "forced failure"


def test_adapters_reject_wrong_problem_family_without_guessing() -> None:
    problem = make_damped_oscillator_problem(
        mass_kg=1.0,
        damping_coefficient_kg_per_s=0.2,
        stiffness_n_per_m=4.0,
        initial_displacement_m=0.1,
        initial_velocity_m_per_s=0.0,
        duration_s=1.0,
    )
    pde_problem = problem.model_copy(update={"family": "pde"})
    with pytest.raises(ScientificValidationError, match="only initial_value_ode"):
        DampedOscillatorAdapter().run(pde_problem, ODESolverConfiguration())


def test_knowledge_mappings_are_explicit_and_validated_separately() -> None:
    problem = make_damped_oscillator_problem(
        mass_kg=1.0,
        damping_coefficient_kg_per_s=0.2,
        stiffness_n_per_m=4.0,
        initial_displacement_m=0.1,
        initial_velocity_m_per_s=0.0,
        duration_s=1.0,
    )
    mapping = problem.knowledge_mappings[0]
    record = build_example_registry().get("damped_harmonic_oscillator")
    validate_knowledge_mappings(problem, (record,))

    assert mapping.kind == KnowledgeMappingKind.NUMERICAL_REALIZATION
    assert mapping.status == KnowledgeMappingStatus.CANDIDATE
    assert mapping.record_review(reviewer="reviewer", notes="Checked").status == (
        KnowledgeMappingStatus.REVIEWED
    )
    bad_problem = problem.model_copy(
        update={
            "knowledge_mappings": (
                mapping.model_copy(update={"equation_record_id": "unknown_record"}),
            )
        }
    )
    with pytest.raises(ScientificValidationError, match="unknown equation"):
        validate_knowledge_mappings(bad_problem, (record,))
    source_problem = problem.model_copy(
        update={
            "knowledge_mappings": (
                mapping.model_copy(
                    update={"source_reference_ids": ("unlisted_source",)}
                ),
            )
        }
    )
    with pytest.raises(ScientificValidationError, match="unknown source IDs"):
        validate_knowledge_mappings(source_problem, (record,))


def test_result_contract_rejects_inconsistent_sampled_state() -> None:
    problem = make_damped_oscillator_problem(
        mass_kg=1.0,
        damping_coefficient_kg_per_s=0.2,
        stiffness_n_per_m=4.0,
        initial_displacement_m=0.1,
        initial_velocity_m_per_s=0.0,
        duration_s=1.0,
        point_count=4,
    )
    result = DampedOscillatorAdapter().run(problem, ODESolverConfiguration())
    assert isinstance(result, SimulationResult)
    with pytest.raises(ScientificValidationError, match="state names and units"):
        replace(result, states=(replace(result.states[0], unit="s"), result.states[1]))


def test_model_mapping_does_not_change_legacy_error_semantics() -> None:
    with pytest.raises(ScientificValidationError):
        simulate_damped_oscillator(0.0, 0.1, 1.0, 0.1, 0.0, 1.0)
    with pytest.raises(ScientificValidationError):
        simulate_damped_pendulum(1.0, 1.0, 0.1, 0.0, 0.1, 0.0, 1.0)

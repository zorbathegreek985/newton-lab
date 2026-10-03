"""Adapters from reusable simulation contracts to existing physical models."""

from __future__ import annotations

from importlib.metadata import PackageNotFoundError, version

import numpy as np

from newton_lab.dynamics import simulate_damped_oscillator
from newton_lab.exceptions import IntegrationError, ScientificValidationError
from newton_lab.pendulum import simulate_damped_pendulum
from newton_lab.simulation import (
    IntegrationDiagnostics,
    KnowledgeMappingKind,
    KnowledgeModelMapping,
    NamedValue,
    ODESolverConfiguration,
    OutputSampling,
    ParameterCorrespondence,
    SimulationFailure,
    SimulationOutcome,
    SimulationProblem,
    SimulationResult,
    SimulationStatus,
    SolverCapabilities,
    StateSeries,
    StateVariable,
    VariableCorrespondence,
    make_reproducibility_record,
)

_SCIPY_METHODS = ("RK23", "RK45", "DOP853", "Radau", "BDF", "LSODA")
_DEFAULT_CONFIGURATION = ODESolverConfiguration()


def _scipy_version() -> str:
    try:
        return version("scipy")
    except PackageNotFoundError:
        return "unknown"


def _validate_adapter_input(
    problem: SimulationProblem,
    configuration: ODESolverConfiguration,
    *,
    expected_model_id: str,
    expected_parameters: tuple[str, ...],
    expected_states: tuple[str, ...],
    expected_initial_conditions: tuple[str, ...],
) -> None:
    if problem.model_id != expected_model_id:
        raise ScientificValidationError(
            f"adapter for {expected_model_id!r} cannot run {problem.model_id!r}"
        )
    if problem.family != "initial_value_ode":
        raise ScientificValidationError(
            "the SciPy model adapters support only initial_value_ode problems"
        )
    if (
        problem.interval is None
        or problem.sampling is None
        or problem.independent_variable is None
        or not problem.states
    ):
        raise ScientificValidationError(
            "initial_value_ode adapters require an interval, time variable, "
            "states, and output sampling"
        )
    if problem.interval[0] != 0.0:
        raise ScientificValidationError(
            "these model adapters require a simulation interval starting at zero"
        )
    if configuration.method not in _SCIPY_METHODS:
        raise ScientificValidationError(
            f"unsupported SciPy solve_ivp method: {configuration.method!r}"
        )
    actual_parameters = {item.name for item in problem.parameters}
    if actual_parameters != set(expected_parameters):
        raise ScientificValidationError(
            f"problem parameters must be exactly {expected_parameters}"
        )
    actual_states = tuple(item.name for item in problem.states)
    if actual_states != expected_states:
        raise ScientificValidationError(
            f"problem states must be ordered as {expected_states}"
        )
    actual_initial = {item.name for item in problem.initial_conditions}
    if actual_initial != set(expected_initial_conditions):
        raise ScientificValidationError(
            f"initial conditions must be exactly {expected_initial_conditions}"
        )


def _values_by_name(values: tuple[NamedValue, ...]) -> dict[str, float]:
    return {item.name: item.value for item in values}


def _result_or_failure(
    *,
    problem: SimulationProblem,
    configuration: ODESolverConfiguration,
    capabilities: SolverCapabilities,
    independent_variable: np.ndarray | None,
    state_arrays: tuple[np.ndarray, ...] = (),
    function_evaluations: int = 0,
    error: IntegrationError | None = None,
) -> SimulationOutcome:
    reproducibility = make_reproducibility_record(problem, configuration, capabilities)
    if error is not None:
        return SimulationFailure(
            problem=problem,
            solver_configuration=configuration,
            error_type=type(error).__name__,
            message=str(error),
            reproducibility=reproducibility,
        )
    if (
        problem.independent_variable is None
        or problem.sampling is None
        or problem.interval is None
    ):
        raise ScientificValidationError(
            "successful ODE results require interval, sampling, "
            "and independent variable"
        )
    if independent_variable is None:
        raise RuntimeError("successful simulation requires output samples")
    states = tuple(
        StateSeries(
            name=declaration.name,
            unit=declaration.unit,
            values=array.copy(),
        )
        for declaration, array in zip(problem.states, state_arrays, strict=True)
    )
    return SimulationResult(
        independent_variable=independent_variable.copy(),
        independent_variable_name=problem.independent_variable.name,
        independent_variable_unit=problem.independent_variable.unit,
        states=states,
        parameters=problem.parameters,
        problem=problem,
        solver_configuration=configuration,
        solver_capabilities=capabilities,
        status=SimulationStatus.SUCCESS,
        diagnostics=IntegrationDiagnostics(
            function_evaluations=function_evaluations,
            # solve_ivp's t_eval output does not expose internal accepted steps.
            accepted_internal_steps=None,
            message="Integration completed successfully.",
        ),
        reproducibility=reproducibility,
        assumptions=problem.assumptions,
        validity_limits=problem.validity_limits,
    )


def _scipy_capabilities() -> SolverCapabilities:
    return SolverCapabilities(
        solver_id="scipy_solve_ivp",
        solver_version=_scipy_version(),
        problem_families=("initial_value_ode",),
        supported_methods=_SCIPY_METHODS,
        adaptive_internal_steps=True,
        supports_requested_output_grid=True,
        supports_events=False,
        stiffness_support=(
            "Radau, BDF, and LSODA can address some stiff IVPs; method selection "
            "requires problem-specific judgment."
        ),
        limitations=(
            "This adapter handles initial-value ODEs only.",
            "Completion does not establish model validity or global error bounds.",
            "Requested output samples are interpolated from adaptive internal steps.",
        ),
    )


class DampedOscillatorAdapter:
    """Run the existing oscillator function through the common result contract."""

    capabilities = _scipy_capabilities()

    def run(
        self,
        problem: SimulationProblem,
        configuration: ODESolverConfiguration = _DEFAULT_CONFIGURATION,
    ) -> SimulationOutcome:
        """Adapt a validated oscillator problem; integration errors become outcomes."""
        _validate_adapter_input(
            problem,
            configuration,
            expected_model_id="damped_harmonic_oscillator",
            expected_parameters=(
                "mass_kg",
                "damping_coefficient_kg_per_s",
                "stiffness_n_per_m",
            ),
            expected_states=("displacement", "velocity"),
            expected_initial_conditions=("displacement", "velocity"),
        )
        assert problem.interval is not None and problem.sampling is not None
        parameters = _values_by_name(problem.parameters)
        initial = _values_by_name(problem.initial_conditions)
        try:
            result = simulate_damped_oscillator(
                mass_kg=parameters["mass_kg"],
                damping_coefficient_kg_per_s=parameters["damping_coefficient_kg_per_s"],
                stiffness_n_per_m=parameters["stiffness_n_per_m"],
                initial_displacement_m=initial["displacement"],
                initial_velocity_m_per_s=initial["velocity"],
                duration_s=problem.interval[1] - problem.interval[0],
                num_points=problem.sampling.point_count,
                method=configuration.method,
                relative_tolerance=configuration.relative_tolerance,
                absolute_tolerance=configuration.absolute_tolerance,
                maximum_step_s=configuration.maximum_step,
            )
        except IntegrationError as exc:
            return _result_or_failure(
                problem=problem,
                configuration=configuration,
                capabilities=self.capabilities,
                independent_variable=None,
                error=exc,
            )
        return _result_or_failure(
            problem=problem,
            configuration=configuration,
            capabilities=self.capabilities,
            independent_variable=result.time_s,
            state_arrays=(result.displacement_m, result.velocity_m_per_s),
            function_evaluations=result.function_evaluations,
        )


class DampedPendulumAdapter:
    """Run the existing full-sine pendulum through the common result contract."""

    capabilities = _scipy_capabilities()

    def run(
        self,
        problem: SimulationProblem,
        configuration: ODESolverConfiguration = _DEFAULT_CONFIGURATION,
    ) -> SimulationOutcome:
        """Adapt a pendulum problem; retain its nonlinear sine restoring term."""
        _validate_adapter_input(
            problem,
            configuration,
            expected_model_id="damped_nonlinear_pendulum",
            expected_parameters=(
                "mass_kg",
                "length_m",
                "damping_coefficient_kg_m2_per_s",
                "gravity_m_per_s2",
            ),
            expected_states=("angle", "angular_velocity"),
            expected_initial_conditions=("angle", "angular_velocity"),
        )
        assert problem.interval is not None and problem.sampling is not None
        parameters = _values_by_name(problem.parameters)
        initial = _values_by_name(problem.initial_conditions)
        try:
            result = simulate_damped_pendulum(
                mass_kg=parameters["mass_kg"],
                length_m=parameters["length_m"],
                damping_coefficient_kg_m2_per_s=parameters[
                    "damping_coefficient_kg_m2_per_s"
                ],
                gravity_m_per_s2=parameters["gravity_m_per_s2"],
                initial_angle_rad=initial["angle"],
                initial_angular_velocity_rad_per_s=initial["angular_velocity"],
                duration_s=problem.interval[1] - problem.interval[0],
                num_points=problem.sampling.point_count,
                method=configuration.method,
                relative_tolerance=configuration.relative_tolerance,
                absolute_tolerance=configuration.absolute_tolerance,
                maximum_step_s=configuration.maximum_step,
            )
        except IntegrationError as exc:
            return _result_or_failure(
                problem=problem,
                configuration=configuration,
                capabilities=self.capabilities,
                independent_variable=None,
                error=exc,
            )
        return _result_or_failure(
            problem=problem,
            configuration=configuration,
            capabilities=self.capabilities,
            independent_variable=result.time_s,
            state_arrays=(result.angle_rad, result.angular_velocity_rad_per_s),
            function_evaluations=result.function_evaluations,
        )


def make_damped_oscillator_problem(
    *,
    mass_kg: float,
    damping_coefficient_kg_per_s: float,
    stiffness_n_per_m: float,
    initial_displacement_m: float,
    initial_velocity_m_per_s: float,
    duration_s: float,
    point_count: int = 1001,
) -> SimulationProblem:
    """Describe an oscillator run and its explicit mapping to Phase 4 metadata."""
    mapping = KnowledgeModelMapping(
        equation_record_id="damped_harmonic_oscillator",
        kind=KnowledgeMappingKind.NUMERICAL_REALIZATION,
        variable_correspondences=(
            VariableCorrespondence(
                equation_symbol="x",
                model_state="displacement",
                equation_unit="m",
                model_unit="m",
            ),
            VariableCorrespondence(
                equation_symbol="x_dot",
                model_state="velocity",
                equation_unit="m/s",
                model_unit="m/s",
            ),
        ),
        parameter_correspondences=(
            ParameterCorrespondence(
                equation_symbol="m",
                model_parameter="mass_kg",
                equation_unit="kg",
                model_unit="kg",
            ),
            ParameterCorrespondence(
                equation_symbol="c",
                model_parameter="damping_coefficient_kg_per_s",
                equation_unit="kg/s",
                model_unit="kg/s",
            ),
            ParameterCorrespondence(
                equation_symbol="k",
                model_parameter="stiffness_n_per_m",
                equation_unit="N/m",
                model_unit="N/m",
            ),
        ),
        assumptions=("Free oscillator with linear spring and viscous damping.",),
        validity_limits=("No external drive or nonlinear restoring force is modeled.",),
    )
    return SimulationProblem(
        model_id="damped_harmonic_oscillator",
        name="Damped harmonic oscillator",
        family="initial_value_ode",
        model_version="1",
        independent_variable=StateVariable(name="time", unit="s"),
        interval=(0.0, duration_s),
        states=(
            StateVariable(name="displacement", unit="m"),
            StateVariable(name="velocity", unit="m/s"),
        ),
        parameters=(
            NamedValue(name="mass_kg", value=mass_kg, unit="kg"),
            NamedValue(
                name="damping_coefficient_kg_per_s",
                value=damping_coefficient_kg_per_s,
                unit="kg/s",
            ),
            NamedValue(name="stiffness_n_per_m", value=stiffness_n_per_m, unit="N/m"),
        ),
        initial_conditions=(
            NamedValue(name="displacement", value=initial_displacement_m, unit="m"),
            NamedValue(name="velocity", value=initial_velocity_m_per_s, unit="m/s"),
        ),
        sampling=OutputSampling(point_count=point_count),
        assumptions=(
            "Linear spring, linear viscous damping, constant parameters, no drive.",
        ),
        validity_limits=(
            "Real oscillators may have nonlinear or non-viscous effects.",
        ),
        equation_record_ids=("damped_harmonic_oscillator",),
        knowledge_mappings=(mapping,),
    )


def make_damped_pendulum_problem(
    *,
    mass_kg: float,
    length_m: float,
    damping_coefficient_kg_m2_per_s: float,
    gravity_m_per_s2: float,
    initial_angle_rad: float,
    initial_angular_velocity_rad_per_s: float,
    duration_s: float,
    point_count: int = 1001,
) -> SimulationProblem:
    """Describe a full nonlinear pendulum run and its equation correspondence."""
    mapping = KnowledgeModelMapping(
        equation_record_id="damped_nonlinear_pendulum",
        kind=KnowledgeMappingKind.NUMERICAL_REALIZATION,
        variable_correspondences=(
            VariableCorrespondence(
                equation_symbol="theta",
                model_state="angle",
                equation_unit="rad",
                model_unit="rad",
            ),
            VariableCorrespondence(
                equation_symbol="theta_dot",
                model_state="angular_velocity",
                equation_unit="rad/s",
                model_unit="rad/s",
            ),
        ),
        parameter_correspondences=(
            ParameterCorrespondence(
                equation_symbol="m",
                model_parameter="mass_kg",
                equation_unit="kg",
                model_unit="kg",
            ),
            ParameterCorrespondence(
                equation_symbol="L",
                model_parameter="length_m",
                equation_unit="m",
                model_unit="m",
            ),
            ParameterCorrespondence(
                equation_symbol="b",
                model_parameter="damping_coefficient_kg_m2_per_s",
                equation_unit="kg·m²/s",
                model_unit="kg·m²/s",
            ),
            ParameterCorrespondence(
                equation_symbol="g",
                model_parameter="gravity_m_per_s2",
                equation_unit="m/s²",
                model_unit="m/s²",
            ),
        ),
        assumptions=(
            "Point-mass bob, massless rigid rod, fixed pivot, uniform gravity.",
            "Linear viscous torque and no external torque.",
            "The full nonlinear sine restoring term is retained.",
        ),
        validity_limits=("Rod flexibility and non-viscous damping are not modeled.",),
    )
    return SimulationProblem(
        model_id="damped_nonlinear_pendulum",
        name="Damped nonlinear pendulum",
        family="initial_value_ode",
        model_version="1",
        independent_variable=StateVariable(name="time", unit="s"),
        interval=(0.0, duration_s),
        states=(
            StateVariable(name="angle", unit="rad"),
            StateVariable(name="angular_velocity", unit="rad/s"),
        ),
        parameters=(
            NamedValue(name="mass_kg", value=mass_kg, unit="kg"),
            NamedValue(name="length_m", value=length_m, unit="m"),
            NamedValue(
                name="damping_coefficient_kg_m2_per_s",
                value=damping_coefficient_kg_m2_per_s,
                unit="kg·m²/s",
            ),
            NamedValue(name="gravity_m_per_s2", value=gravity_m_per_s2, unit="m/s²"),
        ),
        initial_conditions=(
            NamedValue(name="angle", value=initial_angle_rad, unit="rad"),
            NamedValue(
                name="angular_velocity",
                value=initial_angular_velocity_rad_per_s,
                unit="rad/s",
            ),
        ),
        sampling=OutputSampling(point_count=point_count),
        assumptions=(
            "Point-mass bob, massless rigid rod, fixed pivot, uniform gravity.",
            "Linear viscous torque and no external torque.",
            "The full nonlinear sine restoring term is retained.",
        ),
        validity_limits=("Rod flexibility and non-viscous damping are not modeled.",),
        equation_record_ids=("damped_nonlinear_pendulum",),
        knowledge_mappings=(mapping,),
    )

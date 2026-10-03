# Reusable simulation architecture

Newton Lab now separates four concerns: descriptive physics knowledge,
model-specific numerical functions, a family-aware simulation request/result
contract, and adapters that connect a supported model to a supported solver.
The current adapters reuse the oscillator and pendulum functions; they do not
interpret registry expression strings.

## Problem and result contracts

`SimulationProblem` describes a model ID and version, an open solver-family
label, independent variable and interval, named states and units, parameters,
initial and boundary conditions, output sampling, assumptions, validity limits,
equation IDs, and source IDs. Not every problem family needs every field. The
shared container does not assert that every solver can handle every family.
Family-specific adapters validate their own requirements.

`SimulationResult` carries requested independent-variable samples, named state
arrays and units, parameters, problem and solver metadata, diagnostics,
assumptions, validity notes, and a reproducibility record. Sampled arrays
contain the requested output grid; adaptive internal integration steps remain
controlled by the solver. SciPy's `t_eval` results do not expose the exact
number of accepted internal steps through the current model functions, so this
diagnostic is left unknown rather than inferred from output sample count.
`SimulationFailure` carries a clear error type and message for adapter calls.
The existing direct functions continue to raise `IntegrationError` and
`ScientificValidationError` as before.

## Solver family and capabilities

`SimulationAdapter` is a small protocol with an explicit `SolverCapabilities`
declaration. The current SciPy `solve_ivp` adapters support only
`initial_value_ode`, six declared integration methods, adaptive internal
steps, and requested output grids. The adapter does not expose events. `RK23`,
`RK45`, and `DOP853` are explicit Runge–Kutta methods; `Radau`, `BDF`, and
`LSODA` can address some stiff problems. A method name is not an automatic
suitability guarantee; users must judge stiffness, scale, accuracy, and model
domain. No automatic solver selection is implemented.

`ODESolverConfiguration` records method, relative tolerance, absolute
tolerance, and optional maximum step. Tolerances control local error estimates
and do not guarantee a global error bound. Output point count controls only
requested samples. Units are recorded and never converted silently.

Future PDE, boundary-value, algebraic, optimization, stochastic, or
event-driven solvers can define family-specific problem/result data and
implement the adapter protocol or a distinct protocol when their data do not
fit these time-series arrays. They do not need to change the current ODE
adapters.

## Existing model adapters

`make_damped_oscillator_problem()` and
`make_damped_pendulum_problem()` construct explicit problem records. Pass those
to `DampedOscillatorAdapter` or `DampedPendulumAdapter`, optionally with an
`ODESolverConfiguration`.

```python
from newton_lab.simulation import ODESolverConfiguration, SimulationResult
from newton_lab.simulation_adapters import (
    DampedPendulumAdapter,
    make_damped_pendulum_problem,
)

problem = make_damped_pendulum_problem(
    mass_kg=0.2,
    length_m=0.75,
    damping_coefficient_kg_m2_per_s=0.01,
    gravity_m_per_s2=9.81,
    initial_angle_rad=0.8,
    initial_angular_velocity_rad_per_s=0.0,
    duration_s=8.0,
)
outcome = DampedPendulumAdapter().run(
    problem, ODESolverConfiguration(method="DOP853", maximum_step=0.02)
)
if isinstance(outcome, SimulationResult):
    print(outcome.states[0].values[-1])
else:
    print(outcome.error_type, outcome.message)
```

The pendulum adapter delegates to the existing full-sine implementation; it
does not substitute the small-angle equation. The existing
`simulate_damped_oscillator` and `simulate_damped_pendulum` signatures and
result data remain available. Both now also accept optional `maximum_step_s`;
the default leaves SciPy unrestricted as before.

## Knowledge mapping

Each adapter problem includes an explicit candidate
`KnowledgeModelMapping`: equation-record ID, mapping kind, typed variable and
parameter correspondences with units, assumptions, and validity limits.
`validate_knowledge_mappings(problem, records)` verifies equation IDs and
referenced source IDs against the knowledge records. Mapping review is separate
from equation scientific verification and software implementation status.
A candidate equation relationship does not create or approve an
implementation mapping. Review a mapping explicitly by recording reviewer and
notes. No mapping runs the equation text.

The current model mappings are local metadata on simulation problems and
results; no new persistent store record was introduced. Environment metadata
records Python/platform, NumPy and SciPy versions, model and solver versions,
parameters, initial/boundary conditions, tolerances, method, maximum step,
sample grid configuration, and equation/source IDs. This aids reproduction but
does not promise bit-for-bit results across platforms or library versions.
Deterministic models have no random seed; stochastic models can populate that
field in a future adapter.

## Adding a model or solver

For a new model, retain a focused scientific function with its own meaningful
validation and exception behavior. Add a `SimulationProblem` constructor and
an adapter that checks family, named inputs, units, and supported solver
configuration, delegates to the model, and builds the common result or
failure. Document assumptions and validity limits, then test numerical
compatibility against the model function.

For a new solver family, first define that family's capabilities, configuration,
and result shape. Reuse the shared envelope only where it preserves the
scientific meaning of the outputs. Do not present solver completion as proof of
physical validity, and do not generalize the ODE methods to unrelated problem
families.

# Reproducible experiments

An experiment is a versioned, explicit record of model parameters, requested
outputs, solver settings, and one or more ordered executions. The framework
currently demonstrates this with the damped harmonic oscillator, nonlinear
pendulum, nonlinear algebraic root problem, and spatial boundary-value problems.
These are initial adapters, not the intended limit of Newton Lab's
Physics-to-Systems research scope.

## Parameter sweeps and metrics

`ExperimentSpecification` records named parameters with units and an optional
`ParameterSweep`. A sweep substitutes values for one declared parameter in the
order supplied by the caller. Each point receives a fresh parameter tuple. The
runner is sequential; expected input validation, solver, and metric extraction
failures are stored on that run, and the sweep continues to later points.
Unexpected programming exceptions remain visible.

The oscillator adapter supports four explicitly defined metrics:

| Metric | Definition | Unit |
| --- | --- | --- |
| `maximum_absolute_displacement_m` | Maximum `abs(displacement_m)` over the requested output samples | m |
| `final_displacement_m` | Displacement at the final requested sample | m |
| `final_velocity_m_per_s` | Velocity at the final requested sample | m/s |
| `maximum_absolute_velocity_m_per_s` | Maximum `abs(velocity_m_per_s)` over the requested output samples | m/s |

The maximum metrics are over the requested output grid, not a guarantee of the
continuous-time maximum between samples. Metrics are finite scalar records,
not opaque dictionaries. Other solver families do not get coerced into these
time-series definitions: algebraic results have their own scalar metrics and
diagnostics, and boundary-value results retain their spatial mesh and profile.

The pendulum adapter separately defines its metric IDs and units:

| Metric | Definition | Unit |
| --- | --- | --- |
| `final_angle_rad` | Angle at the final requested sample | rad |
| `final_angular_velocity_rad_per_s` | Angular velocity at the final requested sample | rad/s |
| `maximum_absolute_angle_rad` | Maximum `abs(angle_rad)` over requested samples, including both endpoints | rad |
| `maximum_absolute_angular_velocity_rad_per_s` | Maximum `abs(angular_velocity_rad_per_s)` over requested samples | rad/s |

The pendulum parameters use the existing model's API: mass in kg, length in m,
damping torque coefficient in kg*m^2/s, gravity in m/s^2, initial angle in
radians, initial angular velocity in rad/s, and duration in seconds. The model
is a point mass on a massless rigid rod with fixed pivot, uniform gravity,
linear viscous damping torque, constant parameters, and no external drive. Its
nonlinear `sin(theta)` restoring term is retained by the adapter.

## Nonlinear algebraic experiments

`make_nonlinear_spring_experiment` uses the existing
`make_nonlinear_spring_equilibrium_problem` factory and
`ScipyHybridRootSolver`; the experiment layer does not implement root finding.
The example residual is `F(x) = k*x + a*x**3 - f`, where displacement is in m,
`k` is in N/m, `a` is in N/m^3, and applied force is in N. Its default sweep
varies force while using the same configured initial guess independently at
every point. It never seeds one run from the previous run's candidate.

For an independently checkable reference, `k=100 N/m`, `a=1000 N/m^3`, and
`f=1.001 N` give `x=0.01 m` because `100*0.01 + 1000*(0.01)^3 = 1.001`.
The sometimes-cited combination `f=1.1 N`, `x=0.01 m` does **not** satisfy
this equation: its force-balance residual is `-0.099 N`.

Algebraic metrics are separate from the ODE time-series metrics:

| Metric ID | Definition | Unit |
| --- | --- | --- |
| `solution:<unknown-name>` | Corresponding component of the solver's candidate vector | Declared variable unit, or `unknown` |
| `residual:<residual-name>` | Signed raw equation residual component | Declared residual unit, or `unknown` |
| `residual_norm` | Euclidean L2 norm of raw residual entries, without scaling | Common declared residual unit, or `unknown` |

Discrete evidence stays in `algebraic_diagnostics`: solver convergence,
residual acceptance, initial guess, termination message, iteration count when
available, and function evaluations. `solver_success` means the algebraic
solver reported convergence; `numerical_acceptance` reflects the solver's
separate residual-tolerance check. A nonconverged candidate or a converged
candidate with an unacceptable residual is a failed experiment run, but its
candidate, metrics, residuals, and diagnostics are retained when available.
The solver's raw L2 residual tolerance is not converted between units or
scaled. The adapter does not mark a candidate physically valid.

The experiment checks a caller-supplied reference vector against the current
problem's residual tolerance before comparing a selected solution component.
The comparison records the candidate value, reference value, absolute error,
and a separate metric comparison tolerance. It neither changes residual
acceptance nor proves the model correct.

Initial guess can affect a local solve and the root it returns. Diagnostics
record it on every run. Multiple roots may exist; neither success nor two
different candidates establishes uniqueness or root-set completeness. The
cubic spring's stated calibration-range constraint is not enforced by the
solver; linear and cubic stiffness are constrained to positive values by the
spring experiment adapter. Physical validity remains unassessed.

Central sensitivity is available for a scalar algebraic metric, but returned
candidate roots may not follow a verified common branch. Algebraic sensitivity
records are therefore labeled `candidate_difference` with branch consistency
`not_established`. Interpret the central difference as a comparison of returned
candidates, not a validated local derivative. Solver tolerances and perturbation
size can materially affect it; no continuation strategy is implied.

## Use a parameter sweep

```python
from newton_lab.experiments import make_oscillator_experiment, run_experiment
from newton_lab.simulation import ODESolverConfiguration, OutputSampling

specification = make_oscillator_experiment(
    sweep_values=(0.0, 0.1, 0.2),
    metrics=("final_displacement_m", "maximum_absolute_displacement_m"),
    solver_configuration=ODESolverConfiguration(
        relative_tolerance=1e-9,
        absolute_tolerance=1e-11,
    ),
    output_sampling=OutputSampling(point_count=501),
)
record = run_experiment(specification)
for run in record.runs:
    print(run.swept_parameter_value, run.status, run.metrics, run.message)
```

The same runner dispatches the pendulum adapter from its explicit model ID:

```python
from newton_lab.experiments import make_pendulum_experiment, run_experiment

pendulum = make_pendulum_experiment(
    sweep_parameter_name="initial_angle_rad",
    sweep_values=(0.05, 0.2, 0.5),
    metrics=("final_angle_rad", "maximum_absolute_angle_rad"),
)
record = run_experiment(pendulum)
```

Each run calls the existing `simulate_damped_pendulum` function with a fresh
parameter tuple and the requested ODE configuration/output grid. Parameter
order and sweep order are preserved, and expected failures follow the same
record-and-continue policy as oscillator runs.

An invalid damping value, for example, produces a model-validation failure
record for that point without clamping it or preventing subsequent points from
running. Solver failures retain their exception type and message. Metric
extraction failures are distinct from input and integration failures.

## Finite-difference sensitivity

`run_sensitivity_analysis` evaluates the base point, `p-h`, then `p+h`, and
estimates the absolute derivative with the central difference
`[f(p+h)-f(p-h)]/(2h)`. The result records the base value, perturbation,
perturbed parameter values, all three metric values, derivative units, and the
three run records. If any run fails, a `SensitivityFailure` retains all
outcomes; no one-sided substitute is used. The utility
`estimate_central_difference` also supports a caller-supplied scalar function
and is tested against a function with a known derivative.

For example, select `final_angle_rad` and perturb `initial_angle_rad` on a
pendulum specification with `sweep_values=None`. The derivative has units
`rad/rad`; it remains an absolute sensitivity estimate and is not normalized.

This estimate is not an exact derivative. It can be sensitive to perturbation
size, round-off, discretization, solver tolerances, and model validity limits.
The perturbation must be positive and finite, and every parameter value must
remain valid for the model.

## Analytical reference comparisons

An experiment may request an absolute-tolerance comparison for a scalar metric.
The oscillator adapter supplies an exact analytical final-displacement
reference only for zero damping, using the undamped harmonic-oscillator
solution. Other model points or metrics record the reference as unavailable.
The comparison records the computed and reference values, common unit, absolute
error, tolerance, and pass/fail status. Agreement for one metric verifies only
that comparison under its assumptions; it does not establish physical validity.

The pendulum adapter does not reuse the oscillator reference. For an undamped
pendulum it can compare the nonlinear solver's final angle to the analytical
solution of the **small-angle linearized equation**, whose angular frequency is
`sqrt(g/L)`. The reported absolute error is the discrepancy in final angle at
the same requested end time and initial conditions. This compares the full
`sin(theta)` model with its `sin(theta) ~= theta` approximation; the reference
is not an exact solution of the nonlinear pendulum. Approximation accuracy
depends on angular regime and measured quantity. This implies no universal
angle threshold or general model equivalence. With damping, the comparison
reference is unavailable in the current adapter.

## Boundary-value experiments

`make_heat_conduction_experiment` adapts the existing `solve_bvp` heat example
to this experiment framework. It exposes left and right endpoint temperatures
as ordered parameters. Each sweep point substitutes one value into a fresh BVP
problem copy and retains its own parameter record. This example constrains
temperatures to non-negative kelvins. Expected parameter and solver failures
become failed run records, and later sweep points continue.

The supported spatial metrics are:

| Metric | Definition | Unit |
| --- | --- | --- |
| `domain_start`, `domain_end` | First and last coordinates of the returned mesh | Problem domain unit |
| `mesh_point_count` | Number of returned spatial nodes | count |
| `sampled_minimum:<state>`, `sampled_maximum:<state>` | Extrema over returned nodal samples | State unit |
| `sampled_minimum_location:<state>`, `sampled_maximum_location:<state>` | Returned mesh coordinate of the first sampled extremum | Problem domain unit |
| `left_endpoint:<state>`, `right_endpoint:<state>` | State values at returned first and last nodes | State unit |
| `boundary_residual_max_abs`, `boundary_residual_l2_norm` | Raw boundary callback residual norm | Common declared residual unit, otherwise `unknown` |
| `boundary_residual_component_0`, `boundary_residual_component_1` | Signed raw endpoint residual components | Corresponding declared residual unit, otherwise `unknown` |

Extrema are sampled mesh extrema; they do not assert a continuous-domain
extremum. The run's `bvp_diagnostics` stores the actual mesh, both aligned
solution components, differential residual estimates, boundary residuals,
solver status/message, iteration count, tolerances, and the solver's
`numerical_accepted` result. The heat example declares both endpoint residuals
in kelvins. A generic BVP that omits residual units records them as `unknown`;
the raw L2 norm is not scaled. The solver's `solver_success`, its boundary
residual check, and numerical acceptance remain distinct fields. Neither
acceptance nor a small boundary residual establishes interior accuracy,
physical validity, uniqueness, or completeness.

The heat example offers an explicit `BVPReferenceRequest`. For constant
conductivity, no internal generation, and fixed endpoint temperatures, the
exact steady profile is linear. The adapter evaluates that reference directly
at every coordinate in the returned mesh, so it needs no interpolation and
does not compare arrays by index across different grids. It reports maximum
absolute and RMS sampled profile errors in `bvp_reference_comparison`.
Unsupported models report the comparison as unavailable. Agreement checks one
known case; it does not validate the solver generally. A reference comparison
is separate evidence and does not change the run's solver or acceptance status.

```python
from newton_lab.experiments import (
    BVPReferenceRequest,
    make_heat_conduction_experiment,
    run_experiment,
)

specification = make_heat_conduction_experiment(
    left_temperature_k=300.0,
    sweep_parameter_name="right_temperature_k",
    sweep_values=(350.0, 400.0, 450.0),
    bvp_reference=BVPReferenceRequest(maximum_absolute_tolerance=1e-5),
)
record = run_experiment(specification)
for run in record.runs:
    print(run.parameters, run.status, run.bvp_reference_comparison)
```

The existing central-difference utility also accepts BVP scalar metrics. Such
results are labeled `mesh_dependent_metric_difference`: adaptive mesh changes
can affect sampled quantities. The method and perturbation are recorded, all
base/lower/upper runs are retained, and an invalid parameter point yields a
`SensitivityFailure`; no one-sided estimate is substituted. This numerical
difference alone does not establish a physically meaningful sensitivity.

## Reproducibility and scope

Every run records its actual ordered model parameters, experiment/configuration
version, assumptions, available equation and source IDs, solver identity and
method, tolerances, maximum step, requested output count, Python/NumPy/SciPy
versions, and random seed when applicable. This supports configuration
reconstruction and investigation. It does not promise bitwise-identical output
across machines, platforms, or dependency versions. No experiment persistence
format is added in this phase; records remain in memory.

The nonlinear spring experiment can be run directly:

```python
from newton_lab.algebraic import AlgebraicSolverConfiguration
from newton_lab.experiments import (
    AnalyticalReferenceRequest,
    make_nonlinear_spring_experiment,
    run_experiment,
)

specification = make_nonlinear_spring_experiment(
    applied_force_n=1.001,
    sweep_values=(0.8, 1.001, 2.008),
    initial_displacement_m=0.005,
    metrics=("solution:displacement", "residual:force_balance", "residual_norm"),
    solver_configuration=AlgebraicSolverConfiguration(residual_tolerance=1e-10),
    analytical_reference=AnalyticalReferenceRequest(
        metric_id="solution:displacement", absolute_tolerance=1e-8
    ),
    reference_solution=(0.01,),
)
record = run_experiment(specification)
for run in record.runs:
    print(run.swept_parameter_value, run.status, run.metrics)
```

The supplied reference is a root only at the force value `1.001 N`; other
sweep points report that reference comparison as unavailable because it does
not satisfy their residual tolerance.

The framework supports oscillator and nonlinear pendulum ODEs, the nonlinear
algebraic root family, and the boundary-value heat example. The BVP adapter
keeps spatial grids and solution profiles separate from ODE time-series and
algebraic results. No registry equation text is executed.

Read-only collection summaries, metric compatibility checks, sweep projections,
and descriptive statistics are documented in the
[experiment analysis guide](experiment_analysis.md).

Headless scientific plots and Markdown reports over these stored experiment
records are covered by the [visualization and reporting guide](visualization.md).

# Architecture

## Current foundation

Newton Lab starts as a small Python package using a `src` layout. Public
configuration lives in `config.py`; project-wide exception types live in
`exceptions.py`; reusable scientific input checks live under `utils/`.
Dependencies flow inward: utility and domain modules may depend on shared
exceptions, while configuration does not depend on application features.
Tests mirror the package's public behavior under `tests/`.

NumPy, SciPy, pandas, and Pydantic are runtime dependencies selected for
scientific data and validated configuration needs. Feature modules should add
dependencies only when a concrete, tested use case requires them.

## Adding research areas

Add a focused module when there is a well-defined concept or workflow, with
documented inputs, outputs, assumptions, and tests. Keep domain calculations
independent of user interfaces and external services. Models should expose
assumptions and units where applicable; simulations should make random seeds
controllable. Cross-cutting abstractions should be introduced only after more
than one real use case demonstrates the need.

The first dynamical-systems module, `dynamics.py`, contains a damped harmonic
oscillator and its result type. Further models should remain focused and expose
their assumptions, parameter units, solver settings, and failure behavior.
Stochastic processes or network models can be added as separate modules when
there is a concrete, tested research need. Trading, interfaces, persistence, and
integrations remain outside this foundation and require separate scope.

The `knowledge` package adds an extensible physics taxonomy and validated
equation registry. Registry records hold scientific and source metadata;
solver modules remain separate. See [the registry documentation](physics_registry.md)
for schema fields, evidence categories, example records, and extension guidance.

Candidate ingestion and review live alongside, but remain separate from, the
accepted equation registry. Human acceptance does not promote scientific
verification. `JsonKnowledgeStore` persists the complete workspace in a
versioned document and protects existing files by default. See
[the knowledge workflow guide](knowledge_workflow.md) for state transitions,
supported import paths, and persistence behavior.

Equation structure metadata and graph relationships extend the same registry
and workspace. Comparisons consume declared metadata only; candidate discovery
never accepts a relationship. See
[the mathematical structure guide](mathematical_structure.md) for comparison,
directionality, compatibility, and scientific limits.

Reusable simulation contracts describe model-independent requests and common
results, while the current SciPy adapters support only initial-value ODEs.
They delegate to the existing oscillator and full-sine pendulum implementations
and retain their public APIs. Future solver families can introduce suitable
family-specific data without inheriting an ODE-only abstraction. See
[the simulation architecture guide](simulation_architecture.md) for solver
capabilities, mappings, reproducibility metadata, and limitations.

The nonlinear algebraic root family is separate from the ODE time-series
contracts. It uses a trusted residual callable, SciPy's local hybrid root
method, and an algebraic-specific result that distinguishes solver convergence
from residual acceptance. See
[the nonlinear algebraic solver guide](nonlinear_algebraic_solver.md) for its
supported square systems, equilibrium example, and limitations.

The boundary-value family represents endpoint-constrained spatial problems
with its own mesh, collocation settings, candidate result, and residual checks.
It uses SciPy `solve_bvp`, and does not inherit time-series output contracts or
algebraic root semantics. See [the boundary-value solver guide](boundary_value_solver.md)
for the heat-conduction example and interpretation of its numerical evidence.

The experiment layer records explicit parameters, solver settings, ordered
sweep outcomes, and model-specific metrics through adapters. Current adapters
call the existing oscillator and nonlinear pendulum functions and the
algebraic and boundary-value solvers. Algebraic records retain candidates and
residual acceptance separately from ODE time-series data. BVP records retain
the returned spatial mesh, profile, and boundary/differential residual
diagnostics separately from both ODE and algebraic results. Sensitivity
estimates and analytical comparisons are experiment utilities, not new solver
families. See
[the experiment guide](experiments.md) for failure policy, reproducibility,
root limitations, and extension points. The separate read-only analysis layer
projects these records into collection summaries, qualified metric comparisons,
sweep views, and descriptive statistics without combining unrelated metrics
into a score. See [the experiment analysis guide](experiment_analysis.md).

The optional visualization and reporting module consumes these stored results
and analysis projections without rerunning solvers. It preserves ODE time,
algebraic candidate, and BVP spatial structures, and it is not imported by the
core solver modules. See [the visualization guide](visualization.md).

The discovery workflow reuses Phase 14 metric and sweep analysis, the declared
equation structures and review-required relationship graph, and the Phase 15
Markdown experiment report. Optional sampled ODE trajectories support
qualified temporal descriptors while algebraic and BVP result types remain
distinct. It keeps computed observations, structural candidates, and
application hypotheses at separate evidence levels. It does not modify the
knowledge base or run solvers. See
[the discovery guide](discovery.md).

The Phase 18 oscillator-to-adjustment case study reuses an existing discovery
hypothesis and oscillator solver while keeping source and target measurements
distinct. Its restricted equation mapping and evidence boundary are documented
in [the cross-domain case-study guide](cross_domain_case_study.md).

The Phase 19 case study reuses the steady heat-conduction BVP and discovery
reporting APIs, then evaluates a separate exact finite sine-mode solution of
the transient diffusion equation. At the time it was created, the registry had
no heat equation record; its experiment and observations therefore retain their
original missing equation identity. Phase 20 adds a separate retrospective
association instead of rewriting those records. See
[the heat diffusion and smoothing case study](heat_diffusion_smoothing_case_study.md).

The Phase 20 integration adds source-backed steady and transient heat equation
records, expresses their conditional reduction with the existing relationship
types, and connects the Phase 18/19 experiments through immutable retrospective
association records. It preserves original run snapshots and discovery identity
data, while oscillator/diffusion structure sharing stays a review candidate.
See [the knowledge integration guide](knowledge_integration.md).

Phase 21 adds versioned local JSON persistence for completed Phase 18 and
Phase 19 case studies. Loaded artifacts rebuild the Phase 20 index and report
without rerunning experiments. See [persistent discovery artifacts](discovery_artifacts.md).

Phase 22 adds a bounded hypothesis-driven analysis of modal decay in those
existing oscillator and heat models. It reuses their solver and analytical APIs
and reports parameter dependence and observation limits without claiming
physical equivalence. See [the Phase 22 research report](research/phase_22_modal_dynamics.md).

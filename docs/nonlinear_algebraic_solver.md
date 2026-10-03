# Nonlinear algebraic root solver

Newton Lab's algebraic solver handles one narrow family: finite-dimensional,
real, square systems of nonlinear residual equations

\[
F(x; p)=0,
\]

where the unknown vector `x` has one component per residual and `p` contains
declared scalar parameters. A problem supplies a trusted Python callable for
`F`; Newton Lab does not parse or execute equation-expression strings from the
knowledge registry. Units and residual labels are optional metadata and remain
unknown when not supplied. Variable and parameter units may also be unknown.
Omitted assumptions, constraints, and validity limits remain unknown (`None`)
rather than asserting that none apply.

## Method and configuration

`ScipyHybridRootSolver` uses `scipy.optimize.root(method="hybr")`, the Powell
hybrid method for local root finding. The configuration exposes:

- `step_tolerance`, passed to SciPy as `xtol` for relative iterate-change
  convergence;
- `residual_tolerance`, independently used to accept the final residual's
  Euclidean L2 norm;
- `maximum_function_evaluations`, passed as `maxfev`.

There is no automatic variable/residual scaling or unit conversion. If residual
units are declared, they must all be the same; alternatively leave all
residual units unknown. The reported L2 norm uses the raw residual entries.
Residual magnitudes therefore depend on the model's units and scaling, and users
must choose an appropriate residual tolerance. Variables with very different
scales can also make local convergence more difficult.

The solver returns an `AlgebraicResult` when SciPy returns a candidate, even if
SciPy reports non-convergence. Its `solver_converged` and `residual_accepted`
fields are separate; `success` requires both. Iteration count is unknown for
the current SciPy method and is reported as `None`; function evaluations are
counted. A malformed/non-finite residual or solver exception without a usable
candidate returns `AlgebraicFailure` with diagnostic text.

## Equilibrium example

```python
from newton_lab.algebraic import (
    AlgebraicResult,
    ScipyHybridRootSolver,
    make_nonlinear_spring_equilibrium_problem,
)

problem = make_nonlinear_spring_equilibrium_problem()
result = ScipyHybridRootSolver().solve(problem)
if isinstance(result, AlgebraicResult) and result.success:
    print(result.solution_candidate[0], result.residual_norm)
```

The example residual is `k*x + a*x**3 - f`, with displacement `x` in metres,
linear stiffness `k` in N/m, cubic stiffness `a` in N/m³, and applied force `f`
in newtons. Its default parameters make `x=0.1 m` directly checkable because
`100 N/m * 0.1 m + 1000 N/m³ * (0.1 m)^3 = 11 N`. The model assumes one-degree-
of-freedom quasi-static force balance and a calibrated cubic restoring law.
The included validity note says that the illustrative constitutive model needs
calibration before use for a real spring. The example does not create a source
or claim experimental validation.

## Interpretation and limits

This is local root finding. The initial guess can affect convergence and which
root is found. Multiple roots may exist; one returned candidate does not prove
uniqueness or completeness. A small residual indicates approximate satisfaction
of the supplied equations, not correctness of those equations, physical
validity of the candidate, or satisfaction of the problem's descriptive
constraints. The current solver does not enforce bounds or constraints, handle
complex unknowns, solve over/underdetermined least-squares systems, or promise
convergence for every square system. Singular or poorly conditioned systems may
fail, stagnate, or be sensitive to numerical perturbations. `hybr` does not
provide an iteration count or a condition estimate through this interface.

## Extending the algebraic family

Add another solver as an explicit implementation of the `AlgebraicSolver`
protocol, with its own typed configuration and `AlgebraicSolverCapabilities`.
Keep its supported system shapes, bounds, derivatives, scaling, and termination
criteria explicit. Return the same result only when its convergence semantics
can be represented honestly; otherwise define a family-specific result.
Optimization, boundary-value, PDE, and stochastic problems are not treated as
nonlinear algebraic roots here and should receive separate families and
scientifically appropriate results.

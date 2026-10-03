# Boundary-value solver family

Newton Lab's boundary-value API represents spatial endpoint-constrained
problems. It is separate from initial-value ODE time series and from
finite-dimensional algebraic root solving. The first example is steady,
one-dimensional heat conduction in a homogeneous rod:

\[
\frac{d^2 T}{dx^2}=0,\qquad T(0)=T_0,\quad T(L)=T_L.
\]

Here `x` and `L` are metres, `T` and endpoint values are kelvins, and the
temperature gradient is K/m. This equation follows from steady Fourier
conduction for constant thermal conductivity and no internal heat generation.
It assumes a continuum, one-dimensional geometry and fixed endpoint
temperatures. Real materials, sources, contact resistance, and multidimensional
effects can require different models; this example does not represent every
heat-transfer system.

## Numerical method and evidence

The solver rewrites a scalar second-order equation as two first-order equations
and passes them to SciPy `solve_bvp`, which uses adaptive fourth-order
collocation. The returned nodal values are a numerical approximation. For this
special equation the exact analytical profile is linear,
`T(x) = T0 + (TL-T0)x/L`; that expression is useful for verification but is not
how the solver computes its output.

`solver_success` reports SciPy's termination status. `boundary_conditions_satisfied`
reports the independently evaluated endpoint residual check, and
`numerical_accepted` requires solver success plus endpoint and differential
residuals within the configured acceptance tolerance. These flags do not
establish physical validity, uniqueness, or completeness of a solution set.
Collocation error estimates and tolerance choices are numerical controls, not
guarantees about model error. Refine tolerances and compare with known solutions
when accuracy matters.

## Usage

```python
from newton_lab.boundary_value import (
    ScipyBoundaryValueSolver,
    make_linear_heat_conduction_problem,
)

problem = make_linear_heat_conduction_problem(
    length_m=0.5,
    left_temperature_k=300.0,
    right_temperature_k=350.0,
)
result = ScipyBoundaryValueSolver().solve(problem)
if result.numerical_accepted:
    print(result.mesh)  # position in metres
    print(result.solution[0])  # temperature in kelvins
    print(result.solution[1])  # temperature gradient in K/m
else:
    print(result.message)
```

Callback-based problems are trusted Python functions, not parsed equation
strings. This initial family supports two-state first-order systems and two
endpoint residuals without unknown solver parameters. Strong nonlinearities,
poor initial guesses, singular systems, and multiple solutions can prevent
convergence or make a candidate difficult to interpret.

# Phase 19 — Heat conduction, diffusion, and smoothing

## Research question

Under what assumptions does the existing steady heat BVP connect to transient diffusion, and does its modal attenuation reduce variation in a controlled synthetic signal?

## Existing steady heat-conduction model

The source is Newton Lab's existing `steady_linear_heat_conduction` BVP. For constant conductivity k and no internal heat generation, its equation is d/dx(k dT/dx)=0, or T_xx=0. It uses domain [0, 1] m, temperature in K, and fixed endpoints 300 K and 400 K. The model assumes a homogeneous material and constant conductivity; the source term is absent. Conductivity and thermal capacity are not estimated by this BVP. The solver returned a spatial mesh of 21 points, solver_success=True, numerical_accepted=True, maximum boundary residual 0, and maximum differential residual 1.35e-15. Its maximum difference from the analytical steady line on that mesh is 0 K. This result is spatial and steady; it is not a time simulation.

## Transient diffusion model and relationship

The transient model is T_t = alpha T_xx, with alpha > 0 in m^2/s, constant material properties, no internal source, and fixed endpoint temperatures. Its steady equation T_xx=0 has the linear profile T_base(x)=T_left+(T_right-T_left)x/L. Writing u=T-T_base gives u_t=alpha u_xx and u(0,t)=u(L,t)=0. For the declared finite sine-series initial perturbation, the exact solution is T=T_base+sum(a_n sin(n*pi*x/L) exp(-alpha*(n*pi/L)^2*t)). This establishes the equilibrium relation and convergence of these finite modes under the stated conditions; it does not cover arbitrary media, boundary histories, or forcing. This is a mathematical correspondence under matched assumptions, not evidence that the steady BVP alone validates transient material properties or physical time scales. For constant diffusivity and fixed endpoint temperatures, subtracting the linear steady profile T_base(x) reduces the transient equation to u_t = alpha*u_xx with homogeneous Dirichlet boundaries. Each sine mode is an eigenfunction and decays by exp(-alpha*(n*pi/L)^2*t). Thus the time-independent profile solves the steady equation, while the transient solution approaches it for this finite-mode perturbation.

## Modal attenuation experiment

Length=1 m; alpha=0.0001 m^2/s; sampled positions=101; times (s)=0, 60, 300, 1200. The listed temperature values are analytical evaluations on a finite spatial grid. The attenuation factor is exact for each mode.

| Mode n | Initial amplitude (K) | factor at t=0 s | factor at t=60 s | factor at t=300 s | factor at t=1200 s |
|---:|---:|---:|---:|---:|---:|
| 1 | 1 | 1 | 0.942502 | 0.743722 | 0.305944 |
| 4 | 0.5 | 1 | 0.387716 | 0.00876131 | 5.89218e-09 |

The exact perturbation spatial RMS-gradient metric is sqrt((1/2) sum((a_n n pi/L)^2 factor_n^2)) in K/m. This uses orthogonality over the continuous interval, not a finite-difference grid estimate:

| Time (s) | Perturbation RMS gradient (K/m) |
|---:|---:|
| 0 | 4.9672941 |
| 60 | 2.7112548 |
| 300 | 1.6525931 |
| 1200 | 0.67963715 |

For each mode, the decay rate is proportional to n squared, so the higher-frequency mode decays faster under this model. The metric tracks only the transient perturbation; the steady baseline gradient is not counted.

## Synthetic signal comparison

The target is a deterministic synthetic sum of the same n=1 and n=4 sine modes, sampled at the same positions and transformed with the same exponential factors. Sampling evaluates the continuous analytical solution on a finite grid; it is not a discrete finite-difference filter. Its units are arbitrary. This is a mathematical smoothing analogy, not market data, a trading strategy, or a forecasting test.

| Mode n | Initial signal amplitude | amplitude at t=0 s | amplitude at t=60 s | amplitude at t=300 s | amplitude at t=1200 s |
|---:|---:|---:|---:|---:|---:|
| 1 | 1 | 1 | 0.942502 | 0.743722 | 0.305944 |
| 4 | 0.4 | 0.4 | 0.155086 | 0.00350452 | 2.35687e-09 |

The shared mathematics is the mode-wise linear attenuation operator. Smoothing is not denoising: high-frequency components can carry meaningful signal, and reducing variation does not establish better recovery of an underlying truth.

## Newton Lab discovery and evidence

Source experiment: `phase19_steady_heat_source`; model: `boundary_value_problem`; source equation IDs declared: none. The existing discovery pipeline was run on the unchanged source experiment. Because the registry currently contains no heat-conduction equation record and the source has no equation ID, no registry structure match or cross-domain relationship is claimed. The analytical derivation above is classified as an analytical model result; the signal result is a synthetic structural analogy only. No real-world or financial application is validated.

### Source experiment report

# Phase 19 steady heat-conduction source experiment

Generated from existing experiment records; no simulations were run.

## Collection summary

- Records: 1
- Runs: 1
- Succeeded: 1
- Failed: 0
- Solver success / failure / unavailable: 1 / 0 / 0
- Numerical acceptance / rejection / unavailable: 1 / 0 / 0

## Source records

### Record 0: `phase19_steady_heat_source`

- Model: `boundary_value_problem`
- Record status: `succeeded`
- Problem family: `boundary_value`
- Parameters and units are recorded for each run below.

Assumptions:
- Steady state.
- One-dimensional homogeneous material with constant thermal conductivity.
- No internal heat generation.

| Run | Outcome | Solver | Solver success | Numerical acceptance | Parameters | Failure |
| ---: | --- | --- | --- | --- | --- | --- |
| 0 | succeeded | scipy_solve_bvp / fourth_order_collocation | True | True | left_temperature_k=300 K; right_temperature_k=400 K |  |

Run source `phase19_steady_heat_source#0`; solver metadata: `configuration=1; Python=3.12.10; NumPy=2.5.3; SciPy=1.18.1; seed=None; rtol=None; atol=None; residual tolerance=None; BVP tolerance=1e-06`.

Run assumptions:
- Steady state.
- One-dimensional homogeneous material with constant thermal conductivity.
- No internal heat generation.

## Provenance and limitations

- Runs without reproducibility metadata: 0 of 1.
- Run parameters, solver status, and family-specific diagnostics remain in the source records.
- Numerical acceptance does not establish physical validity, uniqueness, or model accuracy.
- Descriptive summaries and visual appearance do not establish causality or predictive power.
- A smaller metric value is not inherently better; interpret each quantity using its definition and assumptions.


## Limitations

- The built-in equation registry has no heat-conduction equation record; the source BVP declares no equation_record_ids, so discovery must not manufacture an equation identity or relationship.
- The transient model is an exact finite sine-series solution, not a general numerical PDE solver or a validated heat-transfer simulation.
- The signal is synthetic and uses arbitrary units; reduced high-frequency content is not evidence of denoising, improved prediction, or utility.
- Finite spatial samples display the analytic solution only at the listed grid positions; modal attenuation and RMS gradient are analytic quantities.

## Tests that could challenge the connection

- Compare the modal expression against an independent transient PDE solver under the same homogeneous Dirichlet conditions and refine its mesh/time step.
- Add nonconstant diffusivity, a source term, or changing boundaries and test whether the stated modal correspondence still applies.
- For an application claim, predefine an out-of-sample task and compare smoothed and unsmoothed synthetic signals against known truth.

## Reproduce

From the Newton Lab workspace, run:

```python
from newton_lab.heat_diffusion_case_study import (
    run_heat_diffusion_smoothing_case_study,
)

case = run_heat_diffusion_smoothing_case_study()
print(case.report_markdown)
```
# Phase 18: oscillator to second-order adjustment

This case study tests the Phase 16 damping risk-control hypothesis against one
explicit toy model. It is a mathematical comparison, not a market model,
trading strategy, financial recommendation, or profitability claim.

## Source hypothesis and experiment

The case links the unchanged Phase 16 `damping_risk_control` hypothesis and
source experiment. That hypothesis is proposed, not validated. Its source
metric is `maximum_absolute_velocity_m_per_s`: the largest absolute returned
oscillator velocity sample, in m/s. The Phase 16 sweep uses damping values
`(0.0, 0.1, 0.2) kg/s`, in that order. The resulting source values are retained
in the case report and result object. Its endpoint change is the last sampled
metric minus the first; it is one endpoint comparison, not a trend proof.

The case accepts the source mapping only when the record has the damped
oscillator model ID; mass `1 kg`; stiffness `4 N/m`; initial displacement
`0.1 m`; initial velocity `0 m/s`; a damping-coefficient sweep in `kg/s`; the
requested velocity metric in `m/s`; the no-driving-force assumption; matching
ordered run values and model provenance; and successful finite source metrics.
This is the supported source configuration for the built-in hypothesis. A
source with failures cannot yield the existing complete Phase 16 hypothesis,
so the case rejects it clearly and does not substitute a generated source.
Successful duplicate sweep values remain duplicated and ordered. The embedded
source `ExperimentRecord` preserves its original run and equation provenance.

## Equations and exact transformation

The physical equation is

\[
m\ddot{x}+c\dot{x}+kx=0,
\]

with `m` in kg, `c` in kg/s, `k` in N/m, and `x` in m. It assumes constant
parameters, a linear spring, linear viscous damping, and no external force.

The target is a post-step residual error `q(t)` from a fixed setpoint:

\[
\ddot{q}+\frac{2\zeta}{\tau}\dot{q}+\frac{1}{\tau^2}q=0,
\qquad q(0)=1,\quad \dot{q}(0)=0.
\]

`q` is dimensionless; `q'` has units s^-1 by convention; `tau` is seconds;
and `zeta` is dimensionless. Equilibrium is `(q,q')=(0,0)`. A unit setpoint
step establishes the initial error, and then the reference is held constant.
There are no later commands, disturbances, noise, delays, limits, or
saturation.

Set `x=Aq`, where `A=0.1 m` per unit of `q`, and divide the physical equation
by `m A`. The result is

\[
\ddot{q}+\frac{c}{m}\dot{q}+\frac{k}{m}q=0.
\]

The coefficients match exactly when

\[
\tau=\sqrt{m/k},\qquad
\zeta=\frac{c}{2\sqrt{mk}}=\frac{c\tau}{2m}.
\]

The units check: `sqrt(m/k)` is seconds; `sqrt(m*k)` is kg/s, so `zeta` is
dimensionless. For the source `m=1 kg`, `k=4 N/m`, and `tau=0.5 s`. Its
source damping sweep maps to target damping ratios `(0, 0.025, 0.05)`. The
initial state maps as `x(0)=0.1 m` to `q(0)=1`; initial rates are both zero
after scaling. This is exact mathematical equivalence under the declared
linear constant-parameter transformation. The interpretation of `q` as a
financial adjustment variable remains only a structural analogy.

The target sweep adds lower, critical, and overdamped settings so the
comparison can expose transient tradeoffs. By default it tests
`(0, 0.025, 0.05, 0.7, 1, 2)` in that order. Target samples are produced by the
existing oscillator integrator with the mapped coefficients, scaled back to
`q=x/A`, and compared with exact underdamped, critically damped, and
overdamped closed forms. Each `maximum_analytical_absolute_error` is computed
on the returned grid; it checks numerical agreement, not model validity.

## Metrics and interpretation

All following metrics are finite-grid measurements, except the separately
identified closed-form response and coefficient transformation:

- Peak adjustment rate: `max(abs(q'))`, s^-1.
- Overshoot beyond equilibrium: `max(0, -min(q))`, dimensionless.
- Zero crossings: count of sign changes outside a small numerical zero band.
- Sampled settling time: first returned time after which every remaining
  sample satisfies `abs(q) <= 0.02`, seconds. `not reached` means the finite
  horizon did not establish settling.
- Maximum analytical error: `max(abs(q_numerical - q_exact))`, dimensionless.

For the documented 1201-point output grid, the target results are:

| `zeta` | Peak `|q'|` (s^-1) | Overshoot | Crossings | Sampled settling (s) |
| ---: | ---: | ---: | ---: | ---: |
| 0 | 2.0000 | 1.0000 | 38 | not reached |
| 0.025 | 1.9221 | 0.9236 | 38 | not reached |
| 0.05 | 1.8529 | 0.8536 | 38 | 38.05 |
| 0.7 | 0.9170 | 0.0460 | 5 | 3.00 |
| 1 | 0.7358 | approximately 0 | 0 | 2.95 |
| 2 | 0.4368 | approximately 0 | 0 | 7.45 |

Maximum absolute numerical-versus-analytical state errors are below `3e-9`
on this grid. The normalized source peak rates `max|x'|/A` agree with the
corresponding target sampled peak rates within `0.02 s^-1`; the residual
reflects finite output sampling and solver settings.

Peak rate and overshoot are nonincreasing on this tested grid. The overdamped
`zeta=2` setting settles later than the critical `zeta=1` setting. The zero
and `zeta=0.025` cases do not enter and remain inside the 2% band within the
finite 60-second horizon. Thus reduced oscillation or overshoot does not imply
faster settling. These are observations
on this grid, not universal monotonicity claims. At `zeta >= 1`, the exact
response is non-oscillatory for these initial conditions; a finite count of
sampled crossings alone is not used to establish that analytical fact.

The source and target peak-rate numbers have different units (`m/s` and
dimensionless state units/s) and are not compared by raw magnitude. Their
parameter coordinates are compared only through the derived damping-ratio
mapping.

## Reproduce

From the project root using the existing Windows environment:

```powershell
.venv\Scripts\python.exe -c "from newton_lab.cross_domain import run_damping_transfer_case_study; case = run_damping_transfer_case_study(); print(case.report_markdown)"
```

For a smaller output grid in a script:

```python
from newton_lab.cross_domain import run_damping_transfer_case_study

case = run_damping_transfer_case_study(output_point_count=1201)
print(case.report_markdown)
```

Pass `source_record=...` to connect another source record. It must meet the
compatibility requirements above. The input is not mutated. Target damping
ratios can be configured, but must include each mapped source value, a declared
baseline, and settings on both sides of that baseline. The case result keeps
source values and target values in separate fields and tables.

## Evidence boundary

The case establishes coefficient equivalence for the stated toy equations and
reports controlled synthetic responses. It supports testing a second-order
damping analogy under those assumptions. It does not show that a financial
process follows this equation, predict prices or returns, validate a trading
rule, or establish profitability. Solver success is not model validation.

A real-world study would require suitable empirical data, prespecified
defensible baselines, transaction costs and market impact where relevant,
appropriate risk measures, chronological out-of-sample evaluation, and
controls against overfitting. That study is outside this phase.

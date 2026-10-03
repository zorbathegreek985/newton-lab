# Phase 29: Diffusion-inspired signal smoothing

## Research question and hypothesis

Can transient diffusion attenuate rapidly varying synthetic components
more than slow components, while making responsiveness loss measurable?
The hypothesis is limited to the declared modal model and synthetic
signals; it makes no claim about price prediction or trading utility.

Hypotheses specified before measurement: H1, evaluated attenuation
matches the closed form; H2, higher modes attenuate more; H3, removing
controlled high-frequency contamination can reduce error to clean truth
for that case but need not improve every signal; H4, pulse smoothing
broadens the transition and creates pre-edge influence.

## Mathematical mapping and assumptions

For fixed Dirichlet endpoints, `du/dt = alpha*d2u/dx2` has modal
solution `a_n sin(n*pi*x/L) exp(-alpha*(n*pi/L)^2*t)`. With
`x=L*xi`, xi maps to normalized signal time in `[0,1]`; this is a
mathematical coordinate mapping, not physical space in a market. With
`s=alpha*t/L^2`, modal retention is `exp(-(n*pi)^2*s)`. The length
is 1 m and alpha is 1 m^2/s. Exposure is dimensionless; its numeric
value equals diffusion time in seconds for these selected units. Signal
values use arbitrary units, not literal kelvin temperatures.

This is a whole-interval sine-series transform. Coefficients depend on
the complete interval, including observations after a given xi; the
transformation is noncausal. The comparator is also noncausal: a centered
21-sample moving average, with clipped windows
renormalized at both endpoints. No equal-bandwidth claim is made.
The heat model assumes constant positive diffusivity, no source, and
fixed zero endpoint perturbations. These are modeling assumptions, not
properties attributed to a financial time series.

## Cases and metrics

Synthetic cases:

RMSE is `sqrt(mean((estimate-reference)^2))` on the shared sampled xi grid.
Input RMSE is measured to the known clean truth. Change-to-input RMSE
measures how much a transformation alters the raw input. Sine-mode
retention is the discrete DST-I coefficient after/before ratio.
Modal formula error is max absolute difference between evaluated and
closed-form attenuation factors. Step delay is the nearest linearly
interpolated 0.5 crossing within +/-0.15 xi of each known edge, minus
that edge. Transition width is the 10-to-90% (rising) or 90-to-10%
(falling) crossing distance. Plateau retention is mean output on
xi=[0.45,0.60].
Pre-transition influence is max absolute output on xi=[0.25,0.35),
divided by unit step height. Additional diffusion influence is the max
absolute change from the zero-exposure modal projection in that window.
These are operational pulse metrics, not universal responsiveness measures.

Dimensionless exposures retain order and repeats: (0.0, 1e-05, 5e-05, 0.0001, 0.0005, 0.001, 0.002).
Sample count=501; modal truncation=80;
moving-average width=21 samples.
The moving average uses the same raw input realization per case and is
independent of diffusion exposure. At edges, its available window is
clipped and renormalized. The unsmoothed input is the third baseline.

## Results

| Case | s | status | projection RMSE | raw RMSE truth | diffusion RMSE truth | MA RMSE truth | diffusion change to input | MA change to input | formula error |
|---|---:|---|---:|---:|---:|---:|---:|---:|---:|
- `two_frequency`: Ground truth is mode 1 plus mode 4 with amplitudes 1 and 0.5.
- `slow_plus_high_frequency`: Mode 1 truth plus a deterministic mode 12 perturbation of amplitude 0.35; not a random-noise realization.
- `abrupt_pulse`: Unit pulse on [0.35, 0.7), represented by 80 sine modes.
The pulse uses a truncated sine representation; its zero-exposure
projection error is reported separately from smoothing error.
For pulse edges a and b, coefficients are
`b_n=2*(cos(n*pi*a)-cos(n*pi*b))/(n*pi)`.

| two_frequency | 0 | succeeded | 0 | 0 | 0 | 0.0095213 | 0 | 0.0095213 | 0 |
| two_frequency | 1e-05 | succeeded | 0 | 0 | 0.00056166 | 0.0095213 | 0.00056166 | 0.0095213 | 0 |
| two_frequency | 5e-05 | succeeded | 0 | 0 | 0.0027996 | 0.0095213 | 0.0027996 | 0.0095213 | 0 |
| two_frequency | 0.0001 | succeeded | 0 | 0 | 0.0055774 | 0.0095213 | 0.0055774 | 0.0095213 | 0 |
| two_frequency | 0.0005 | succeeded | 0 | 0 | 0.02704 | 0.0095213 | 0.02704 | 0.0095213 | 0 |
| two_frequency | 0.001 | succeeded | 0 | 0 | 0.052059 | 0.0095213 | 0.052059 | 0.0095213 | 0 |
| two_frequency | 0.002 | succeeded | 0 | 0 | 0.096643 | 0.0095213 | 0.096643 | 0.0095213 | 0 |
| slow_plus_high_frequency | 0 | succeeded | 0 | 0.24724 | 0.24724 | 0.22365 | 0 | 0.028877 | 0 |
| slow_plus_high_frequency | 1e-05 | succeeded | 0 | 0.24724 | 0.24375 | 0.22365 | 0.0034897 | 0.028877 | 0 |
| slow_plus_high_frequency | 5e-05 | succeeded | 0 | 0.24724 | 0.23028 | 0.22365 | 0.016963 | 0.028877 | 0 |
| slow_plus_high_frequency | 0.0001 | succeeded | 0 | 0.24724 | 0.21449 | 0.22365 | 0.032763 | 0.028877 | 0 |
| slow_plus_high_frequency | 0.0005 | succeeded | 0 | 0.24724 | 0.12153 | 0.22365 | 0.12581 | 0.028877 | 0 |
| slow_plus_high_frequency | 0.001 | succeeded | 0 | 0.24724 | 0.06009 | 0.22365 | 0.18768 | 0.028877 | 0 |
| slow_plus_high_frequency | 0.002 | succeeded | 0 | 0.24724 | 0.019957 | 0.22365 | 0.23324 | 0.028877 | 0 |
| abrupt_pulse | 0 | succeeded | 0.05132 | 0 | 0.05132 | 0.083488 | 0.05132 | 0.083488 | 0 |
| abrupt_pulse | 1e-05 | succeeded | 0.05132 | 0 | 0.053398 | 0.083488 | 0.053398 | 0.083488 | 0 |
| abrupt_pulse | 5e-05 | succeeded | 0.05132 | 0 | 0.068854 | 0.083488 | 0.068854 | 0.083488 | 0 |
| abrupt_pulse | 0.0001 | succeeded | 0.05132 | 0 | 0.081454 | 0.083488 | 0.081454 | 0.083488 | 0 |
| abrupt_pulse | 0.0005 | succeeded | 0.05132 | 0 | 0.12152 | 0.083488 | 0.12152 | 0.083488 | 0 |
| abrupt_pulse | 0.001 | succeeded | 0.05132 | 0 | 0.14447 | 0.083488 | 0.14447 | 0.083488 | 0 |
| abrupt_pulse | 0.002 | succeeded | 0.05132 | 0 | 0.17179 | 0.083488 | 0.17179 | 0.083488 | 0 |

Mode attenuation and retention; analytic factors are shown next to
the sampled coefficient ratios. Moving-average ratios are measured on
the same input and do not share the diffusion transfer function.

| Case | s | analytic slow H | diffusion slow | MA slow | analytic fast H | diffusion fast | MA fast |
|---|---:|---:|---:|---:|---:|---:|---:|
| two_frequency | 0 | 1 | 1 | 0.9993 | 1 | 1 | 0.98887 |
| two_frequency | 1e-05 | 0.9999 | 0.9999 | 0.9993 | 0.99842 | 0.99842 | 0.98887 |
| two_frequency | 5e-05 | 0.99951 | 0.99951 | 0.9993 | 0.99214 | 0.99214 | 0.98887 |
| two_frequency | 0.0001 | 0.99901 | 0.99901 | 0.9993 | 0.98433 | 0.98433 | 0.98887 |
| two_frequency | 0.0005 | 0.99508 | 0.99508 | 0.9993 | 0.92408 | 0.92408 | 0.98887 |
| two_frequency | 0.001 | 0.99018 | 0.99018 | 0.9993 | 0.85392 | 0.85392 | 0.98887 |
| two_frequency | 0.002 | 0.98045 | 0.98045 | 0.9993 | 0.72919 | 0.72919 | 0.98887 |
| slow_plus_high_frequency | 0 | 1 | 1 | 0.9993 | 1 | 1 | 0.90231 |
| slow_plus_high_frequency | 1e-05 | 0.9999 | 0.9999 | 0.9993 | 0.98589 | 0.98589 | 0.90231 |
| slow_plus_high_frequency | 5e-05 | 0.99951 | 0.99951 | 0.9993 | 0.9314 | 0.9314 | 0.90231 |
| slow_plus_high_frequency | 0.0001 | 0.99901 | 0.99901 | 0.9993 | 0.86752 | 0.86752 | 0.90231 |
| slow_plus_high_frequency | 0.0005 | 0.99508 | 0.99508 | 0.9993 | 0.49134 | 0.49134 | 0.90231 |
| slow_plus_high_frequency | 0.001 | 0.99018 | 0.99018 | 0.9993 | 0.24142 | 0.24142 | 0.90231 |
| slow_plus_high_frequency | 0.002 | 0.98045 | 0.98045 | 0.9993 | 0.058283 | 0.058283 | 0.90231 |
| abrupt_pulse | 0 | 1 | 0.99976 | 0.99928 | 1 | 1.0287 | 0.89898 |
| abrupt_pulse | 1e-05 | 0.9999 | 0.99966 | 0.99928 | 0.98589 | 1.0141 | 0.89898 |
| abrupt_pulse | 5e-05 | 0.99951 | 0.99926 | 0.99928 | 0.9314 | 0.9581 | 0.89898 |
| abrupt_pulse | 0.0001 | 0.99901 | 0.99877 | 0.99928 | 0.86752 | 0.89238 | 0.89898 |
| abrupt_pulse | 0.0005 | 0.99508 | 0.99483 | 0.99928 | 0.49134 | 0.50543 | 0.89898 |
| abrupt_pulse | 0.001 | 0.99018 | 0.98994 | 0.99928 | 0.24142 | 0.24834 | 0.89898 |
| abrupt_pulse | 0.002 | 0.98045 | 0.98022 | 0.99928 | 0.058283 | 0.059953 | 0.89898 |

Pulse responsiveness metrics use the declared edge, threshold and
windows. Pre-projection influence at zero exposure is Fourier/Gibbs
ringing; the additional-change column isolates change caused by the
diffusion evolution relative to that truncated initial projection.

| s | diffusion rise/fall delay | diffusion rise/fall width | MA rise/fall delay | MA rise/fall width | diffusion/MA plateau | diffusion pre / added / MA influence |
|---:|---|---|---|---|---|---|
| 0 | 5.5782e-05/-5.7274e-05 | 0.011226/0.011226 | -0.001/-0.001 | 0.0336/0.0336 | 0.99971/1 | 0.33836/0/0.47619 |
| 1e-05 | 3.5943e-05/-3.6906e-05 | 0.014142/0.014142 | -0.001/-0.001 | 0.0336/0.0336 | 0.99983/1 | 0.36764/0.07358/0.47619 |
| 5e-05 | 4.7622e-06/-4.8896e-06 | 0.025679/0.025679 | -0.001/-0.001 | 0.0336/0.0336 | 0.99998/1 | 0.42148/0.22395/0.47619 |
| 0.0001 | 2.7809e-07/-2.8544e-07 | 0.036264/0.036264 | -0.001/-0.001 | 0.0336/0.0336 | 1/1 | 0.44378/0.30342/0.47619 |
| 0.0005 | 5.5511e-17/0 | 0.081084/0.081084 | -0.001/-0.001 | 0.0336/0.0336 | 0.9999/1 | 0.47479/0.43956/0.47619 |
| 0.001 | 3.3307e-16/-2.2204e-16 | 0.11465/0.11465 | -0.001/-0.001 | 0.0336/0.0336 | 0.99724/1 | 0.48216/0.47863/0.47619 |
| 0.002 | 2.4817e-09/-2.4817e-09 | 0.16213/0.16213 | -0.001/-0.001 | 0.0336/0.0336 | 0.97903/1 | 0.48739/0.50916/0.47619 |

## Hypothesis assessment for this sweep

- H1 attenuation formula agreement: supported.
  Successful runs have maximum formula error <= 1e-14.
- H2 higher-mode attenuation ordering: supported for positive exposures.
- H3 reconstruction improvement: supported for the controlled high-frequency case.
  Two-frequency detail loss observed: True.
- H4 pulse responsiveness cost: supported by increased measured rise width.
  Midpoint delay and pre-edge effects depend on pulse and truncation.

Recorded rows: 21; succeeded: 21; failed: 0.
Failures retain their sweep index and exposure; failed metrics remain
undefined rather than being replaced with zero.

## Findings by evidence class

**Established mathematical result.** For this fixed-boundary model,
mode n has attenuation `exp(-(n*pi)^2*s)`, so larger n has a faster
decay rate. Agreement between the evaluator and this expression checks
the modal scaling and mode indexing. This residual is an algebraic
cross-check, not an independent PDE discretization or convergence test.

**Experimentally supported result.** The table reports measured modal
retention, clean-signal error, input change, moving-average comparison,
and pulse metrics only for the finite cases and exposure values shown.
Lower high-frequency energy does not by itself imply lower error to
truth; compare both metrics for each exposure.

**Candidate cross-domain application.** Diffusion-inspired smoothing is
a candidate whole-series preprocessing or feature-construction method.
Smoothing modifies the input and cannot restore the original signal
automatically. The centered transform can let later values influence
earlier outputs, as quantified by the pulse pre-influence metric.

**Unvalidated financial application.** No financial observations were
downloaded or tested. This study establishes no predictive power,
trading signal, profitability, excess returns, risk-adjusted performance,
asset/regime robustness, or executable performance after costs. A
causal redesign and separate out-of-sample validation would be required
before any real-time financial feature could be assessed.

## Newton Lab discovery integration

Phase 19 and the existing knowledge integration already record the
diffusion-to-synthetic-signal relationship as a synthetic analogy.
Those APIs do not consume this study's per-exposure result rows; changing
their artifact schema would be broader than this focused experiment.
This report therefore preserves the evidence classes locally rather
than promoting the measurements into a registry relationship.

## Reproducibility and limitations

The synthetic inputs, coefficients, sample grid, moving-average window,
modal cutoff and ordered exposures are deterministic. The exact modal
evaluator is reused; no finite-difference PDE solver is introduced.
The modal formula residual checks the dimensionless scaling against the
existing analytical evaluator; this is not a solver accuracy bound.
For the pulse, finite truncation causes Gibbs/projection error even at
zero exposure. The discrete DST coefficient and sampled metrics have
finite-grid effects. The moving average and diffusion operator have
different transfer functions and boundary behavior, so this is not an
equal-bandwidth comparison. Neither method is a causal feature here.

## Conclusion and next question

The study can test attenuation ordering and synthetic trade-offs under
a declared mathematical mapping. It cannot establish that smoothing
improves real data or predicts future outcomes. A next scientific step
would be to derive and validate a causal diffusion-inspired filter,
then compare its transfer function and delay with causal baselines
before considering any real-data study.

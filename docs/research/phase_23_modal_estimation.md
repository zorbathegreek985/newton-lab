# Phase 23 - Modal rate estimation under imperfect observation

## Research question

When can a decay rate be estimated from sampled, finite-duration, noisy, or
mixed observations of the Phase 18 oscillator and Phase 19 heat modes? This
study uses synthetic data generated from those declared models. Hidden
parameters are retained in separate truth records; fit functions receive only
the observation record.

The reproducible experiment grid is available through
`newton_lab.modal_estimation.run_phase23_research()`. Its fixed-seed synthetic
noise makes the listed observations and estimates repeatable. No external or
financial data are used.

## Models and correspondence

For a scalar observable with a single exponential mode, the assumed model is

\[
y(t)=A e^{-\lambda t},\qquad A>0,\quad \lambda>0.
\]

Here `y` and `A` have the observable's units, time is in seconds, and the decay
rate `lambda` is in s^-1. For the homogeneous one-dimensional heat equation
with constant diffusivity `alpha`, fixed endpoints, and no internal source,
subtracting the steady linear profile gives homogeneous Dirichlet conditions.
Its sine mode `n` decays with

\[
\lambda_n=\alpha\left(\frac{n\pi}{L}\right)^2,
\qquad a_n(t)=a_n(0)e^{-\lambda_n t}.
\]

`alpha` is in m^2/s, `L` in m, and modal temperature amplitude in K. These
assumptions and the separation between transient diffusion and the steady BVP
are those of the Phase 19 example; the steady BVP is not used as a transient
solver.

For the Phase 18 underdamped free oscillator,

\[
m\ddot{x}+c\dot{x}+kx=0,\qquad
\gamma=\frac{c}{2m},\qquad
\omega_d=\sqrt{\frac{k}{m}-\gamma^2}.
\]

Its signed displacement is an oscillatory term multiplied by
`exp(-gamma*t)`. The fitted envelope rate is `gamma` in s^-1, not the damped
frequency `omega_d` in rad/s. The tested oscillator has `m=1 kg`, `c=0.1
kg/s`, `k=4 N/m`, `x(0)=0.1 m`, and `v(0)=0 m/s`, so its exact envelope rate
is 0.05 s^-1 and its damped period is about 3.14257 s. Assumptions are a linear
spring, linear viscous damping, constant parameters, no forcing, and an
underdamped response. The sampled absolute-peak method is an approximation to
the analytic envelope; it does not treat signed samples as positive amplitudes.

These mathematical correspondences apply to the specified models and
conditions. They do not establish that real oscillators or arbitrary measured
signals obey these equations.

## Estimation and observation design

The scalar experiments vary one factor at a time:

- cadence: 0.05, 0.25, or 0.8 s over 4 s;
- duration: 1, 4, or 10 s at 0.25 s cadence;
- additive-noise standard deviation: 0.001, 0.01, or 0.05 signal units over
  12 s at 0.25 s cadence.

The scalar fits compare weighted log-linear regression with nonlinear least
squares against the original observations. Log-linear regression excludes
nonpositive values and values at or below three times the declared additive
noise standard deviation. Nonlinear least squares uses a positive amplitude
and rate and fits the signed, untransformed observations. Its fit RMSE is over
all observations; log-linear RMSE is over samples retained for the log fit.
For the oscillator, the fit RMSE is over retained absolute peaks.

Oscillator cadence is varied from 0.05 to about 1.667 s over a 10 s record.
The duration experiment uses 1.5, 3, or 10 s at 0.05 s cadence. At least four
usable local absolute peaks are required for the envelope regression. Fewer
peaks are reported as insufficient. The study flags fitted peak sequences
with fewer than 16 samples per estimated period as unreliable. Both are
declared design screens, not universal statistical thresholds.

The heat experiments fit exact single-mode samples for n=1 and n=4, then fit a
single exponential to two-mode mixtures at full, early, and late windows.
Mixture effective amplitudes describe each component at the fixed observation
position. A final strong-fast mixture adds 0.005 K iid Gaussian noise with
seed 2340 to show the fast component and then the total late signal falling
below the three-sigma log-fit floor. Fixed seeds are recorded alongside the
additive Gaussian noise standard deviations in observations.

For the noisy scalar reference, the nonlinear fit initializes its rate from a
log-linear slope when at least two positive above-floor samples exist, and its
amplitude from the first such sample. Otherwise it starts from the reciprocal
observation span and the largest absolute observation. It constrains amplitude
and rate to be positive. The weighted log-linear fit uses delta-method weights
proportional to `(observed value / noise standard deviation)^2`. The nonlinear
fit uses unweighted least squares in observation space. These are explicit
choices for this synthetic design, not universal optimal estimators.

Local rate-noise sensitivity is a first-order estimate based on the fit
Jacobian and the declared noise scale. The three-times-sensitivity rule is a
screen only; it is not a confidence interval. The log residual screen at 0.15
is likewise a chosen model-mismatch diagnostic. The 10% error field is a
descriptive post-fit score when a single target rate exists. None of these
thresholds proves identifiability in the statistical or structural sense.

## Hypotheses and results

| Hypothesis | Outcome in the tested design | Evidence and limit |
|---|---|---|
| H1 - Sampling cadence affects estimation reliability. | Supported for the oscillator peak estimator. | At 0.05 s cadence the estimate was 0.04992 s^-1; at about 0.833 s it was 0.07648 s^-1, and the coarsest grid did not supply four usable peaks. The scalar cadence cases happened to stay within 1.4% error for the fixed seed. No universal sampling threshold follows. |
| H2 - Observation duration affects identifiability. | Supported for the tested oscillator windows. | At 1.5 and 3 s, there were only two and three usable peaks and no rate was returned; the 10 s window had seven peaks. This does not establish a general minimum duration. |
| H3 - Noise limits observable decay. | Supported under the declared additive-noise design. | At scalar sigma=0.05, the log fit discarded 20 samples and had 14.3% rate error. In the strong-fast heat mixture (sigma=0.005 K), the full fit excluded below-floor samples and the late window had only one usable sample. These synthetic noise conditions are not sensor-noise estimates. |
| H4 - Mode mixing can create misleading effective rates. | Supported for the two-mode heat examples. | The strong-fast mixture produced estimates from about 0.001000 to 0.015406 s^-1 across late and early windows, while the modes' analytical rates were about 0.000987 and 0.015791 s^-1. The full-window single-exponential residual screen also failed. |
| H5 - Fit quality does not establish parameter identifiability. | Supported as a limitation of interpreting the tested fits. | The close n=(10,11) mixture's early-window single-exponential fit had log RMSE about 0.000359, yet its estimate 0.10863 s^-1 lay between the two true component rates and was not a component identification. Separately, the low-rate/high-noise case converged but was marked not identifiable by the declared local sensitivity screen. This screen is diagnostic evidence, not proof of mathematical non-uniqueness. |

### Numerical findings

The scalar rate (0.3 s^-1) was recoverable within the study's descriptive 10%
error tolerance in the tested conditions. A 1 s observation window produced
about 7.4-7.9% error, while 4 s conditions were closer in this particular
fixed-seed design. At noise standard deviation 0.05, the weighted log fit was
flagged unreliable: it excluded 20 of 49 samples and had log RMSE about 0.220
and rate error about 14.3%. The nonlinear fit used all observations and had
about 3.5% rate error. This comparison reflects these generators, seeds, and
methods; it does not show that nonlinear fitting is generally superior.

The oscillator's dense cadence estimate was 0.04992 s^-1, close to the
analytical 0.05 s^-1. The approximately 0.4 s cadence estimate was 0.05536
s^-1 (10.7% error), and the approximately 0.833 s estimate was 0.07648 s^-1
(53.0% error); both are flagged unreliable by the sampling screen. At about
1.667 s cadence there were too few usable peaks. The 1.5 s and 3 s durations
also produced fewer than four usable peaks. This demonstrates that optimizer
or regression convergence does not supply missing temporal resolution or
duration.

The two exact single heat modes recovered the analytical rates: approximately
0.0009869604 s^-1 for n=1 and 0.01579137 s^-1 for n=4. For mixtures, the
single-exponential estimate depends on the observation window and component
amplitudes. For example, the strong-fast n=(1,4) mixture had a full-window
estimate of about 0.003732 s^-1 with log RMSE 0.813, while its early-window
estimate was about 0.015406 s^-1 and its late-window estimate about
0.0010001 s^-1. A close-rate mixture could have a small single-exponential
residual while still yielding an estimate between its two physical component
rates. A visually or numerically adequate single-exponential fit therefore
does not identify a unique physical mode.

In the noisy strong-fast mixture (`sigma=0.005 K`, seed 2340), the full-window
log fit retained 18 of 61 observations and was flagged unreliable (log RMSE
about 1.672). The 0-120 s fit retained seven points and estimated about
0.015384 s^-1, while the 600-1200 s window had only one above-floor positive
observation and was insufficient. For the isolated fast component with
amplitude 1 K, its own amplitude falls below the 0.015 K three-sigma level
after about 266 s under the exact diffusion rate; the noisy mixture's
observed values additionally include the slow component and Gaussian noise.

`fit_converged`, residual metrics, exclusions, and status are separate result
fields. Status values distinguish
`available_interpretable_under_tested_assumptions`,
`fit_failed`, `insufficient_observations`,
`not_identifiable_at_declared_observation_resolution`, and
`available_but_unreliable_under_declared_diagnostics`. “Interpretable” means
only that the stated screens passed; it does not assert truth, uniqueness, or
physical validity. Evaluation against known truth is performed after fitting
and is not supplied to the estimator.

## Limitations

The generated observations are analytic samples with simple iid additive
Gaussian noise, not sensor records. There are no missing-at-random studies,
irregular-cadence designs, colored or multiplicative noise, uncertain time
stamps, calibration drift, or measurement-process models. The oscillator fit
uses simple sampled local maxima of absolute displacement and can be affected
by noise, sparse cadence, endpoint peaks, and peak-selection errors. The heat
mixtures contain only two known sine modes and use exact sampled signals; the
estimator itself is intentionally a single-exponential fit, not a modal
decomposition algorithm.

This study does not validate denoising, a financial application, trading
profitability, or real-world parameter recovery. Synthetic observations make
controlled failure cases possible but do not substitute for empirical
validation against independently collected data.

## Recommended next experiment

Repeat a subset of the cadence, duration, and noise designs over a predeclared
grid of independent fixed seeds and report the distribution of rate errors,
exclusion counts, and status changes. Add irregular sampling and a noise model
whose scale is estimated from observation-only information. Keep heat-mixture
component identification as a separate question, comparing profile residuals
or alternative two-mode parameter combinations on held-out observation times.
Those extensions would test robustness that a single deterministic seed and a
single-exponential fit cannot establish.

## Reproduce

From the project root, with the project environment installed:

```powershell
.venv\Scripts\python.exe -c "from newton_lab.modal_estimation import run_phase23_research, render_phase23_results; print(render_phase23_results(run_phase23_research()))"
```

The generated result table and detailed fit diagnostics follow. It records
the specific observations, exclusions, convergence messages, residuals, local
sensitivity values, and conditioning diagnostics for every designed fit.

\n\n\n# Phase 23 computed modal-estimation results

| Scenario | Condition | Method | Status | Estimate (s^-1) | Truth / nearest mode (s^-1) | Abs. error/gap (s^-1) | Rel. error | n used / excluded | RMSE (observation units) |
|---|---|---|---|---:|---:|---:|---:|---:|---:|
| scalar_cadence_1 | cadence dt=0.05 s; duration=4 s | weighted_log_linear | available_interpretable_under_tested_assumptions | 0.3020586 | 0.3 | 0.002059 | 0.006862 | 81/0 | 0.01925 arbitrary signal units |
| scalar_cadence_1 | cadence dt=0.05 s; duration=4 s | nonlinear_exponential_least_squares | available_interpretable_under_tested_assumptions | 0.3028267 | 0.3 | 0.002827 | 0.009422 | 81/0 | 0.01922 arbitrary signal units |
| scalar_cadence_2 | cadence dt=0.25 s; duration=4 s | weighted_log_linear | available_interpretable_under_tested_assumptions | 0.2961001 | 0.3 | 0.0039 | 0.013 | 17/0 | 0.02136 arbitrary signal units |
| scalar_cadence_2 | cadence dt=0.25 s; duration=4 s | nonlinear_exponential_least_squares | available_interpretable_under_tested_assumptions | 0.2974908 | 0.3 | 0.002509 | 0.008364 | 17/0 | 0.02131 arbitrary signal units |
| scalar_cadence_3 | cadence dt=0.8 s; duration=4 s | weighted_log_linear | available_interpretable_under_tested_assumptions | 0.3025412 | 0.3 | 0.002541 | 0.008471 | 6/0 | 0.01662 arbitrary signal units |
| scalar_cadence_3 | cadence dt=0.8 s; duration=4 s | nonlinear_exponential_least_squares | available_interpretable_under_tested_assumptions | 0.3040611 | 0.3 | 0.004061 | 0.01354 | 6/0 | 0.01657 arbitrary signal units |
| scalar_duration_1 | duration=1 s; cadence dt=0.25 s | weighted_log_linear | available_interpretable_under_tested_assumptions | 0.3221058 | 0.3 | 0.02211 | 0.07369 | 5/0 | 0.01877 arbitrary signal units |
| scalar_duration_1 | duration=1 s; cadence dt=0.25 s | nonlinear_exponential_least_squares | available_interpretable_under_tested_assumptions | 0.323611 | 0.3 | 0.02361 | 0.0787 | 5/0 | 0.01875 arbitrary signal units |
| scalar_duration_2 | duration=4 s; cadence dt=0.25 s | weighted_log_linear | available_interpretable_under_tested_assumptions | 0.299171 | 0.3 | 0.000829 | 0.002763 | 17/0 | 0.02344 arbitrary signal units |
| scalar_duration_2 | duration=4 s; cadence dt=0.25 s | nonlinear_exponential_least_squares | available_interpretable_under_tested_assumptions | 0.3016236 | 0.3 | 0.001624 | 0.005412 | 17/0 | 0.02335 arbitrary signal units |
| scalar_duration_3 | duration=10 s; cadence dt=0.25 s | weighted_log_linear | available_interpretable_under_tested_assumptions | 0.297108 | 0.3 | 0.002892 | 0.00964 | 35/6 | 0.02019 arbitrary signal units |
| scalar_duration_3 | duration=10 s; cadence dt=0.25 s | nonlinear_exponential_least_squares | available_interpretable_under_tested_assumptions | 0.3038721 | 0.3 | 0.003872 | 0.01291 | 41/0 | 0.02117 arbitrary signal units |
| scalar_noise_1 | additive Gaussian sigma=0.001; seed=2330 | weighted_log_linear | available_interpretable_under_tested_assumptions | 0.2997748 | 0.3 | 0.0002252 | 0.0007508 | 49/0 | 0.0009924 arbitrary signal units |
| scalar_noise_1 | additive Gaussian sigma=0.001; seed=2330 | nonlinear_exponential_least_squares | available_interpretable_under_tested_assumptions | 0.299791 | 0.3 | 0.000209 | 0.0006966 | 49/0 | 0.0009923 arbitrary signal units |
| scalar_noise_2 | additive Gaussian sigma=0.01; seed=2330 | weighted_log_linear | available_interpretable_under_tested_assumptions | 0.2961955 | 0.3 | 0.003805 | 0.01268 | 47/2 | 0.01002 arbitrary signal units |
| scalar_noise_2 | additive Gaussian sigma=0.01; seed=2330 | nonlinear_exponential_least_squares | available_interpretable_under_tested_assumptions | 0.2979128 | 0.3 | 0.002087 | 0.006957 | 49/0 | 0.009923 arbitrary signal units |
| scalar_noise_3 | additive Gaussian sigma=0.05; seed=2330 | weighted_log_linear | available_but_unreliable_under_declared_diagnostics | 0.2570585 | 0.3 | 0.04294 | 0.1431 | 29/20 | 0.05628 arbitrary signal units |
| scalar_noise_3 | additive Gaussian sigma=0.05; seed=2330 | nonlinear_exponential_least_squares | available_interpretable_under_tested_assumptions | 0.2896217 | 0.3 | 0.01038 | 0.03459 | 49/0 | 0.04961 arbitrary signal units |
| oscillator_cadence_1 | oscillator dt=0.05 s; duration=10 s; period=3.14257 s; samples/period=62.85 | absolute_sample_peak_log_linear | available_interpretable_under_tested_assumptions | 0.04991892 | 0.05 | 8.108e-05 | 0.001622 | 7/194 | 6.59e-05 m |
| oscillator_cadence_2 | oscillator dt=0.4 s; duration=10 s; period=3.14257 s; samples/period=7.856 | absolute_sample_peak_log_linear | available_but_unreliable_under_declared_diagnostics | 0.05536292 | 0.05 | 0.005363 | 0.1073 | 7/19 | 0.0004715 m |
| oscillator_cadence_3 | oscillator dt=0.833333 s; duration=10 s; period=3.14257 s; samples/period=3.771 | absolute_sample_peak_log_linear | available_but_unreliable_under_declared_diagnostics | 0.07647902 | 0.05 | 0.02648 | 0.5296 | 6/7 | 0.003718 m |
| oscillator_cadence_4 | oscillator dt=1.66667 s; duration=10 s; period=3.14257 s; samples/period=1.886 | absolute_sample_peak_log_linear | insufficient_observations | unavailable | 0.05 | n/a | n/a | 1/6 | n/a m |
| oscillator_duration_1 | oscillator duration=1.5 s; cadence dt=0.05 s; period=3.14257 s | absolute_sample_peak_log_linear | insufficient_observations | unavailable | 0.05 | n/a | n/a | 2/29 | n/a m |
| oscillator_duration_2 | oscillator duration=3 s; cadence dt=0.05 s; period=3.14257 s | absolute_sample_peak_log_linear | insufficient_observations | unavailable | 0.05 | n/a | n/a | 3/58 | n/a m |
| oscillator_duration_3 | oscillator duration=10 s; cadence dt=0.05 s; period=3.14257 s | absolute_sample_peak_log_linear | available_interpretable_under_tested_assumptions | 0.04991892 | 0.05 | 8.108e-05 | 0.001622 | 7/194 | 6.59e-05 m |
| heat_single_mode_1 | individual heat mode n=1; L=1 m; alpha=1e-4 m^2/s; dt=60 s; duration=1200 s | log_linear | available_interpretable_under_tested_assumptions | 0.0009869604 | 0.0009869604 | 4.337e-19 | 4.394e-16 | 21/0 | 2.009e-16 K |
| heat_single_mode_4 | individual heat mode n=4; L=1 m; alpha=1e-4 m^2/s; dt=15 s; duration=300 s | log_linear | available_interpretable_under_tested_assumptions | 0.01579137 | 0.01579137 | 0 | 0 | 21/0 | 7.532e-17 K |
| mix_separated_equal | heat mode mixture n=(1, 4), effective amplitudes K=(1.0, 0.5); early/full/late windows; dt=20 s | log_linear | available_interpretable_under_tested_assumptions | 0.001109661 | 0.0009869604 | 0.0001227 | n/a | 61/0 | 0.06979 K |
| mix_separated_equal | heat mode mixture n=(1, 4), effective amplitudes K=(1.0, 0.5); early/full/late windows; dt=20 s | log_linear | available_interpretable_under_tested_assumptions | 0.00365016 | 0.0009869604 | 0.002663 | n/a | 7/0 | 0.02689 K |
| mix_separated_equal | heat mode mixture n=(1, 4), effective amplitudes K=(1.0, 0.5); early/full/late windows; dt=20 s | log_linear | available_interpretable_under_tested_assumptions | 0.0009870265 | 0.0009869604 | 6.605e-08 | n/a | 31/0 | 5.752e-06 K |
| mix_weak_fast | heat mode mixture n=(1, 4), effective amplitudes K=(1.0, 0.01); early/full/late windows; dt=20 s | log_linear | available_interpretable_under_tested_assumptions | 0.0009897478 | 0.0009869604 | 2.787e-06 | n/a | 61/0 | 0.001379 K |
| mix_weak_fast | heat mode mixture n=(1, 4), effective amplitudes K=(1.0, 0.01); early/full/late windows; dt=20 s | log_linear | available_interpretable_under_tested_assumptions | 0.001053949 | 0.0009869604 | 6.699e-05 | n/a | 7/0 | 0.0006434 K |
| mix_weak_fast | heat mode mixture n=(1, 4), effective amplitudes K=(1.0, 0.01); early/full/late windows; dt=20 s | log_linear | available_interpretable_under_tested_assumptions | 0.0009869618 | 0.0009869604 | 1.321e-09 | n/a | 31/0 | 1.15e-07 K |
| mix_strong_fast | heat mode mixture n=(1, 4), effective amplitudes K=(0.01, 1.0); early/full/late windows; dt=20 s | log_linear | available_but_unreliable_under_declared_diagnostics | 0.003732061 | 0.0009869604 | 0.002745 | n/a | 61/0 | 0.1624 K |
| mix_strong_fast | heat mode mixture n=(1, 4), effective amplitudes K=(0.01, 1.0); early/full/late windows; dt=20 s | log_linear | available_interpretable_under_tested_assumptions | 0.01540576 | 0.01579137 | 0.0003856 | n/a | 7/0 | 0.002099 K |
| mix_strong_fast | heat mode mixture n=(1, 4), effective amplitudes K=(0.01, 1.0); early/full/late windows; dt=20 s | log_linear | available_interpretable_under_tested_assumptions | 0.001000111 | 0.0009869604 | 1.315e-05 | n/a | 31/0 | 1.15e-05 K |
| mix_close_rates | heat mode mixture n=(10, 11), effective amplitudes K=(1.0, 1.0); early/full/late windows; dt=2 s | log_linear | available_interpretable_under_tested_assumptions | 0.1051129 | 0.09869604 | 0.006417 | n/a | 41/0 | 0.02011 K |
| mix_close_rates | heat mode mixture n=(10, 11), effective amplitudes K=(1.0, 1.0); early/full/late windows; dt=2 s | log_linear | available_interpretable_under_tested_assumptions | 0.10863 | 0.09869604 | 0.009934 | n/a | 5/0 | 0.0005179 K |
| mix_close_rates | heat mode mixture n=(10, 11), effective amplitudes K=(1.0, 1.0); early/full/late windows; dt=2 s | log_linear | available_interpretable_under_tested_assumptions | 0.1033721 | 0.09869604 | 0.004676 | n/a | 21/0 | 7.12e-05 K |
| mix_fast_component_below_noise_floor | strong-fast n=(1,4) mixture; sigma=0.005 K; seed=2340; fast-mode tail crosses the three-sigma floor; dt=20 s | weighted_log_linear | available_but_unreliable_under_declared_diagnostics | 0.01408468 | 0.01579137 | 0.001707 | n/a | 18/43 | 0.01638 K |
| mix_fast_component_below_noise_floor | strong-fast n=(1,4) mixture; sigma=0.005 K; seed=2340; fast-mode tail crosses the three-sigma floor; dt=20 s | weighted_log_linear | available_interpretable_under_tested_assumptions | 0.01538412 | 0.01579137 | 0.0004073 | n/a | 7/0 | 0.004344 K |
| mix_fast_component_below_noise_floor | strong-fast n=(1,4) mixture; sigma=0.005 K; seed=2340; fast-mode tail crosses the three-sigma floor; dt=20 s | log_linear | insufficient_observations | unavailable | n/a | n/a | n/a | 1/30 | n/a K |

## Estimation diagnostics

- `scalar_cadence_1` / `weighted_log_linear`: converged=True; message=weighted linear least squares converged; rate-noise-sensitivity=0.0035221492886458056; Jacobian condition=3.2061185278026234; log RMSE=0.037216384767036786; exclusions=none. estimate passes the declared observation and residual screens only
- `scalar_cadence_1` / `nonlinear_exponential_least_squares`: converged=True; message=`ftol` termination condition is satisfied.; rate-noise-sensitivity=0.0035383217882892404; Jacobian condition=2.80262966700863; log RMSE=None; exclusions=none. estimate passes the declared observation and residual screens only
- `scalar_cadence_2` / `weighted_log_linear`: converged=True; message=weighted linear least squares converged; rate-noise-sensitivity=0.0072234085691167925; Jacobian condition=3.0034771343799065; log RMSE=0.041996388914453946; exclusions=none. estimate passes the declared observation and residual screens only
- `scalar_cadence_2` / `nonlinear_exponential_least_squares`: converged=True; message=`gtol` termination condition is satisfied.; rate-noise-sensitivity=0.007282058824211441; Jacobian condition=2.631454066939635; log RMSE=None; exclusions=none. estimate passes the declared observation and residual screens only
- `scalar_cadence_3` / `weighted_log_linear`: converged=True; message=weighted linear least squares converged; rate-noise-sensitivity=0.01119448865815448; Jacobian condition=2.5143109759767084; log RMSE=0.04234187339591705; exclusions=none. estimate passes the declared observation and residual screens only
- `scalar_cadence_3` / `nonlinear_exponential_least_squares`: converged=True; message=`gtol` termination condition is satisfied.; rate-noise-sensitivity=0.011182713172085128; Jacobian condition=2.2334147142156535; log RMSE=None; exclusions=none. estimate passes the declared observation and residual screens only
- `scalar_duration_1` / `weighted_log_linear`: converged=True; message=weighted linear least squares converged; rate-noise-sensitivity=0.030095134903928965; Jacobian condition=3.436306095210826; log RMSE=0.02342334405857526; exclusions=none. estimate passes the declared observation and residual screens only
- `scalar_duration_1` / `nonlinear_exponential_least_squares`: converged=True; message=`gtol` termination condition is satisfied.; rate-noise-sensitivity=0.030147840250672055; Jacobian condition=2.776466482786736; log RMSE=None; exclusions=none. estimate passes the declared observation and residual screens only
- `scalar_duration_2` / `weighted_log_linear`: converged=True; message=weighted linear least squares converged; rate-noise-sensitivity=0.007516555514775933; Jacobian condition=3.001676187517628; log RMSE=0.0579981795405911; exclusions=none. estimate passes the declared observation and residual screens only
- `scalar_duration_2` / `nonlinear_exponential_least_squares`: converged=True; message=`gtol` termination condition is satisfied.; rate-noise-sensitivity=0.007499988689155939; Jacobian condition=2.6235819634847544; log RMSE=None; exclusions=none. estimate passes the declared observation and residual screens only
- `scalar_duration_3` / `weighted_log_linear`: converged=True; message=weighted linear least squares converged; rate-noise-sensitivity=0.004905083442760782; Jacobian condition=3.342380334189853; log RMSE=0.1025971556438317; exclusions=6 x value is at or below the declared three-sigma noise floor (indices [31, 35, 36, 37, 38, 40]). estimate passes the declared observation and residual screens only
- `scalar_duration_3` / `nonlinear_exponential_least_squares`: converged=True; message=`ftol` termination condition is satisfied.; rate-noise-sensitivity=0.004818193242303562; Jacobian condition=2.3324811132340386; log RMSE=None; exclusions=none. estimate passes the declared observation and residual screens only
- `scalar_noise_1` / `weighted_log_linear`: converged=True; message=weighted linear least squares converged; rate-noise-sensitivity=0.000228039767132419; Jacobian condition=3.398695014858744; log RMSE=0.011825238942655731; exclusions=none. estimate passes the declared observation and residual screens only
- `scalar_noise_1` / `nonlinear_exponential_least_squares`: converged=True; message=`gtol` termination condition is satisfied.; rate-noise-sensitivity=0.00022808420868657213; Jacobian condition=2.311783088004478; log RMSE=None; exclusions=none. estimate passes the declared observation and residual screens only
- `scalar_noise_2` / `weighted_log_linear`: converged=True; message=weighted linear least squares converged; rate-noise-sensitivity=0.002264071263865949; Jacobian condition=3.418415397010178; log RMSE=0.10392301073898849; exclusions=2 x value is at or below the declared three-sigma noise floor (indices [45, 48]). estimate passes the declared observation and residual screens only
- `scalar_noise_2` / `nonlinear_exponential_least_squares`: converged=True; message=`gtol` termination condition is satisfied.; rate-noise-sensitivity=0.0022684430837351323; Jacobian condition=2.3131774815458535; log RMSE=None; exclusions=none. estimate passes the declared observation and residual screens only
- `scalar_noise_3` / `weighted_log_linear`: converged=True; message=weighted linear least squares converged; rate-noise-sensitivity=0.011472023673602061; Jacobian condition=3.4498739108589347; log RMSE=0.21971042824127288; exclusions=16 x value is at or below the declared three-sigma noise floor (indices [25, 26, 28, 31, 33, 34, 35, 37] ...); 4 x nonpositive value cannot be log-transformed (indices [32, 38, 43, 45]). log residual exceeds the declared 0.15 model-mismatch screen
- `scalar_noise_3` / `nonlinear_exponential_least_squares`: converged=True; message=`ftol` termination condition is satisfied.; rate-noise-sensitivity=0.011068814648240864; Jacobian condition=2.3195903108112454; log RMSE=None; exclusions=none. estimate passes the declared observation and residual screens only
- `oscillator_cadence_1` / `absolute_sample_peak_log_linear`: converged=True; message=weighted linear least squares converged; rate-noise-sensitivity=None; Jacobian condition=10.42691569606329; log RMSE=0.0007400672113192191; exclusions=194 x not a local maximum of absolute sampled displacement (indices [1, 2, 3, 4, 5, 6, 7, 8] ...). estimate passes the declared observation and residual screens only; local peak selection assumes sampled extrema approximate envelope peaks; sparse cadence can bias peak timing and height
- `oscillator_cadence_2` / `absolute_sample_peak_log_linear`: converged=True; message=weighted linear least squares converged; rate-noise-sensitivity=None; Jacobian condition=10.618323171526605; log RMSE=0.005866879720281825; exclusions=19 x not a local maximum of absolute sampled displacement (indices [1, 2, 3, 5, 6, 7, 9, 10] ...). estimate passes the declared observation and residual screens only; local peak selection assumes sampled extrema approximate envelope peaks; sparse cadence can bias peak timing and height; only 8 samples per estimated oscillation period; sampled-peak rate is cadence-sensitive
- `oscillator_cadence_3` / `absolute_sample_peak_log_linear`: converged=True; message=weighted linear least squares converged; rate-noise-sensitivity=None; Jacobian condition=9.339231872691341; log RMSE=0.06436038589559011; exclusions=7 x not a local maximum of absolute sampled displacement (indices [1, 3, 5, 7, 9, 10, 12]). estimate passes the declared observation and residual screens only; local peak selection assumes sampled extrema approximate envelope peaks; sparse cadence can bias peak timing and height; only 4 samples per estimated oscillation period; sampled-peak rate is cadence-sensitive
- `oscillator_cadence_4` / `absolute_sample_peak_log_linear`: converged=False; message=fewer than four usable absolute local peaks; rate-noise-sensitivity=None; Jacobian condition=None; log RMSE=None; exclusions=6 x not a local maximum of absolute sampled displacement (indices [1, 2, 3, 4, 5, 6]). absolute-peak envelope fit requires at least four above-floor peaks
- `oscillator_duration_1` / `absolute_sample_peak_log_linear`: converged=False; message=fewer than four usable absolute local peaks; rate-noise-sensitivity=None; Jacobian condition=None; log RMSE=None; exclusions=29 x not a local maximum of absolute sampled displacement (indices [1, 2, 3, 4, 5, 6, 7, 8] ...). absolute-peak envelope fit requires at least four above-floor peaks
- `oscillator_duration_2` / `absolute_sample_peak_log_linear`: converged=False; message=fewer than four usable absolute local peaks; rate-noise-sensitivity=None; Jacobian condition=None; log RMSE=None; exclusions=58 x not a local maximum of absolute sampled displacement (indices [1, 2, 3, 4, 5, 6, 7, 8] ...). absolute-peak envelope fit requires at least four above-floor peaks
- `oscillator_duration_3` / `absolute_sample_peak_log_linear`: converged=True; message=weighted linear least squares converged; rate-noise-sensitivity=None; Jacobian condition=10.42691569606329; log RMSE=0.0007400672113192191; exclusions=194 x not a local maximum of absolute sampled displacement (indices [1, 2, 3, 4, 5, 6, 7, 8] ...). estimate passes the declared observation and residual screens only; local peak selection assumes sampled extrema approximate envelope peaks; sparse cadence can bias peak timing and height
- `heat_single_mode_1` / `log_linear`: converged=True; message=weighted linear least squares converged; rate-noise-sensitivity=None; Jacobian condition=1354.187445064628; log RMSE=2.0821562424230766e-16; exclusions=none. estimate passes the declared observation and residual screens only
- `heat_single_mode_4` / `log_linear`: converged=True; message=weighted linear least squares converged; rate-noise-sensitivity=None; Jacobian condition=338.5544136790694; log RMSE=5.753606209632755e-16; exclusions=none. estimate passes the declared observation and residual screens only
- `mix_separated_equal` / `log_linear`: converged=True; message=weighted linear least squares converged; rate-noise-sensitivity=None; Jacobian condition=1374.4697511337702; log RMSE=0.06314337371799489; exclusions=none. estimate passes the declared observation and residual screens only
- `mix_separated_equal` / `log_linear`: converged=True; message=weighted linear least squares converged; rate-noise-sensitivity=None; Jacobian condition=130.0173087163557; log RMSE=0.02178726788060218; exclusions=none. estimate passes the declared observation and residual screens only
- `mix_separated_equal` / `log_linear`: converged=True; message=weighted linear least squares converged; rate-noise-sensitivity=None; Jacobian condition=4706.9284703542335; log RMSE=1.1474873078381576e-05; exclusions=none. estimate passes the declared observation and residual screens only
- `mix_weak_fast` / `log_linear`: converged=True; message=weighted linear least squares converged; rate-noise-sensitivity=None; Jacobian condition=1374.4697511337702; log RMSE=0.0015068436170454897; exclusions=none. estimate passes the declared observation and residual screens only
- `mix_weak_fast` / `log_linear`: converged=True; message=weighted linear least squares converged; rate-noise-sensitivity=None; Jacobian condition=130.0173087163557; log RMSE=0.0006707122362567495; exclusions=none. estimate passes the declared observation and residual screens only
- `mix_weak_fast` / `log_linear`: converged=True; message=weighted linear least squares converged; rate-noise-sensitivity=None; Jacobian condition=4706.9284703542335; log RMSE=2.2950510874205678e-07; exclusions=none. estimate passes the declared observation and residual screens only
- `mix_strong_fast` / `log_linear`: converged=True; message=weighted linear least squares converged; rate-noise-sensitivity=None; Jacobian condition=1374.4697511337702; log RMSE=0.8128952140033106; exclusions=none. log residual exceeds the declared 0.15 model-mismatch screen
- `mix_strong_fast` / `log_linear`: converged=True; message=weighted linear least squares converged; rate-noise-sensitivity=None; Jacobian condition=130.0173087163557; log RMSE=0.003772862362113367; exclusions=none. estimate passes the declared observation and residual screens only
- `mix_strong_fast` / `log_linear`: converged=True; message=weighted linear least squares converged; rate-noise-sensitivity=None; Jacobian condition=4706.9284703542335; log RMSE=0.0022795823351744304; exclusions=none. estimate passes the declared observation and residual screens only
- `mix_close_rates` / `log_linear`: converged=True; message=weighted linear least squares converged; rate-noise-sensitivity=None; Jacobian condition=91.30796527637663; log RMSE=0.022454505281739118; exclusions=none. estimate passes the declared observation and residual screens only
- `mix_close_rates` / `log_linear`: converged=True; message=weighted linear least squares converged; rate-noise-sensitivity=None; Jacobian condition=8.72421122184656; log RMSE=0.0003586775477892618; exclusions=none. estimate passes the declared observation and residual screens only
- `mix_close_rates` / `log_linear`: converged=True; message=weighted linear least squares converged; rate-noise-sensitivity=None; Jacobian condition=309.45015874472085; log RMSE=0.004879026988482682; exclusions=none. estimate passes the declared observation and residual screens only
- `mix_fast_component_below_noise_floor` / `weighted_log_linear`: converged=True; message=weighted linear least squares converged; rate-noise-sensitivity=9.668512697733819e-05; Jacobian condition=51.31381844013245; log RMSE=1.6716260186880345; exclusions=38 x value is at or below the declared three-sigma noise floor (indices [15, 17, 18, 19, 20, 21, 22, 24] ...); 5 x nonpositive value cannot be log-transformed (indices [27, 39, 46, 50, 54]). log residual exceeds the declared 0.15 model-mismatch screen; more than half of input samples were excluded
- `mix_fast_component_below_noise_floor` / `weighted_log_linear`: converged=True; message=weighted linear least squares converged; rate-noise-sensitivity=0.00012322745185450518; Jacobian condition=44.50394897743344; log RMSE=0.017942576020521364; exclusions=none. estimate passes the declared observation and residual screens only
- `mix_fast_component_below_noise_floor` / `log_linear`: converged=False; message=fewer than three usable positive observations; rate-noise-sensitivity=None; Jacobian condition=None; log RMSE=None; exclusions=26 x value is at or below the declared three-sigma noise floor (indices [0, 1, 2, 3, 4, 5, 6, 8] ...); 4 x nonpositive value cannot be log-transformed (indices [9, 16, 20, 24]). insufficient positive above-floor samples

Synthetic additive Gaussian noise uses NumPy default_rng with fixed scenario seeds recorded in each observation series. Noise-free cases are analytical samples.

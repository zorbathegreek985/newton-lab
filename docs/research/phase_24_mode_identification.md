# Phase 24 - Robustness of modal estimation and multi-mode identification

## Research question

Under which tested observation conditions can a sum of two positive decaying
modes be distinguished from a single mode, and how much do the conclusions
change across noise realizations, irregular timing, rate separation, amplitude
ratio, duration, and cadence? The study uses controlled synthetic signals;
estimators receive observations and initial guesses only. Noise-free curves
and true parameters are held in separate generator-side fields and used only
for post-fit evaluation.

The reproducible implementation is
`newton_lab.modal_identification.run_phase24_research()`. The computed tables
below come from `render_phase24_results()`. No market or external measurement
data are used.

## Models, units, and assumptions

The one-mode reference is

\[
y(t)=A e^{-\lambda t},
\]

and the positive two-mode model is

\[
y(t)=A_1e^{-\lambda_1t}+A_2e^{-\lambda_2t},
\qquad A_1,A_2,\lambda_1,\lambda_2>0,
\quad \lambda_1\ne\lambda_2.
\]

Time and decay rates use seconds and s^-1. Amplitudes and observations use
arbitrary `signal units`. Noise sigma uses the same signal units. These are
synthetic mathematical signals motivated by the exponential heat-mode form
from Phase 19. They are not automatically measurements of a physical heat
system or any other real process.

For a two-mode truth, the dimensionless separation reported is

\[
\Delta=\frac{|\lambda_2-\lambda_1|}{\max(\lambda_1,\lambda_2)}.
\]

The amplitude ratio is the faster-mode amplitude divided by the slower-mode
amplitude. Both quantities describe these experiment conditions; neither is a
universal identification threshold. Rates are matched for evaluation by
sorting fitted and true rates in ascending order. This removes the label-swap
symmetry without using truth in the fit.

## Estimation and model comparison

The single-mode fit reuses the Phase 23 nonlinear exponential least-squares
estimator. The two-mode fit uses bounded nonlinear least squares on the
original observations, with positive amplitudes and rates. It uses five
deterministic starts derived from the observation span and signal scale; a
separate sensitivity condition supplies four explicit alternative starts.
Each start is recorded. The lowest training residual among converged starts
is retained. The amplitude upper bound is `10^6` times the maximum observed
absolute value; the rate bounds are based on the observed time span
(`max(1e-12, 1e-8/span)` through `max(50/span, 100*lower_bound)`). The solver
uses Jacobian scaling and a maximum of 5000 function evaluations. Fitted rates
are independently constrained positive with no minimum-gap constraint, then
sorted after fitting; swapping component labels leaves the curve unchanged.
The two-mode fit can therefore collapse to nearly coincident rates; the
condition diagnostic and explicit rate estimates expose that behavior.

Model selection uses BIC only when the observation record declares a positive
iid additive Gaussian noise standard deviation. Treating that variance as
known, the compared scores omit shared likelihood constants:

\[
\mathrm{BIC}=\frac{\mathrm{RSS}}{\sigma^2}+k\ln(n),
\]

where RSS is the sum of squared observation-space residuals, `n` is the
sample count, and `k` is the number of fitted mean parameters: two for one
mode, four for two modes. `BIC(single)-BIC(two) >= 6` selects two modes;
`<= -6` selects one; values between are inconclusive. Six is a predeclared
decision rule for this study, not a universal cutoff or confidence level.
The BIC interpretation depends on independent Gaussian errors with the
declared known variance. It is not applied to noise-free data. A two-mode
training residual can be lower simply because it has two more free
parameters; residual reduction alone never selects the two-mode model.

Post-fit “supported” recovery requires the selected mode count to agree with
synthetic truth and each matched rate to be within 10% relative error. A
selection that disagrees with known truth, or a selected two-mode fit whose
rates miss that tolerance, is recorded as inconsistent. The 10% tolerance is
descriptive for this experiment only. Other outcomes distinguish
inconclusive selection, insufficient samples, numerical failure, and an
ill-conditioned supported two-mode fit. The scaled Jacobian condition screen
at `1e8` is diagnostic, not a proof of mathematical identifiability.

## Experiment matrix and reproducibility

### Repeated noise realizations

Five conditions use identical truth, timestamps, noise distribution, sigma,
and fitting procedure within each group. Only the fixed noise seed changes.
Each group uses 12 realizations, seeds 2401 through 2412, at 61 regular
timestamps over 15 s, with iid Gaussian noise sigma 0.015 signal units. Each
observation is the exact sum plus `default_rng(seed).normal(0, sigma, n)`:

- well-separated rates (0.2, 0.8) s^-1, amplitudes (1.0, 0.6);
- single-mode control at 0.2 s^-1, amplitude 1.0, for observing false-positive
  two-mode selections;
- moderately separated rates (0.2, 0.4) s^-1, amplitudes (1.0, 0.8);
- close rates (0.2, 0.24) s^-1, equal amplitudes;
- weak fast mode: rates (0.2, 0.8) s^-1, amplitudes (1.0, 0.03).

The count summaries distinguish both fits converging, failures, unreliable,
inconclusive, supported, and inconsistent outcomes. Rate distributions are
reported after ascending-rate matching. Median, minimum, maximum, and
population standard deviation are shown for fitted rates, absolute errors,
and relative errors, plus mean signed rate error to distinguish systematic
offset from across-seed dispersion. These are descriptive summaries of the
specified seeds, not universal probabilities, confidence intervals, or
estimates of real-world noise variability. Two-mode parameter distributions include
converged fits even when BIC selected one mode, so they show optimizer
variability rather than only selected-model estimates.

### Regular and irregular sampling

The main two-mode timing comparison uses the same true signal, sigma 0.015,
and seed 2440 over a 15 s span. It includes a 61-point regular baseline,
mild jitter (uniform perturbations within 0.2 of nominal spacing), stronger
jitter (within 0.4), two observation clusters separated by a long gap, and
observations with the middle 40%-60% of the time span removed. Jitter uses
NumPy `default_rng` with fixed timing seeds. The clustered design places about
60% of observations in the first 20% of the span and the remainder in the
last 20%; the missing-middle design deterministically removes the central
interval. A regular grid with the same count
and endpoints accompanies each irregular design. The renderer reports actual
sample counts and minimum, median, and maximum interval lengths.

The oscillator timing check reuses the Phase 23 analytical Phase 18
underdamped oscillator and its absolute-peak envelope estimator. Its dense
timing designs have 81 samples over 10 s; the missing-middle design is paired
with a regular 64-sample grid. Sparse checks use 7 or 13 samples over 10 s,
with regular, jittered, and clustered timing. They report the estimated
envelope status and rate as well as interval coverage. The oscillator is
`m=1 kg`, `c=0.1 kg/s`, `k=4 N/m`, with initial displacement 0.1 m and zero
initial velocity; its analytical envelope rate is 0.05 s^-1. Irregular
sampling is not assumed to remove aliasing or restore missed peaks.

### Rate, amplitude, duration, cadence, and noise matrix

The compact matrix includes a single-mode control; well-separated,
moderately separated, and close two-mode rates; weak and strong fast-mode
amplitudes; a fast mode whose initial amplitude is above three sigma but falls
below that floor during the observation; short and long observation spans;
sparse and dense grids; low and moderate noise; a noise-free overfit check;
and an initial-guess comparison.
The generated table records each condition's actual true and fitted rates
and amplitudes, rate separation, amplitude ratio, sample count, time span,
noise sigma, both residuals, parameter counts, BIC values, model selection,
post-fit evaluation, and amplitude error where mode matching is defined.

## Hypotheses and outcomes

| Hypothesis | Outcome in the tested designs | Evidence and boundary |
|---|---|---|
| H1 - Noise realizations change estimates or identification. | Supported. | With fixed 15 s/61-sample designs, the well-separated condition supported two-mode recovery in 8/12 runs and was inconsistent in 4/12; the moderate-separation case supported recovery in 4/12 and was inconsistent in 8/12. The single-mode control produced 11 supported single-mode results and one inconclusive result (no false two-mode selection in these seeds). Close-rate and weak-fast groups had no supported recoveries in these seeds. The rates, amplitudes, sigma, and times were held fixed within each group. |
| H2 - Irregular timing changes performance, with direction depending on design. | Supported. | At 61 samples, jittered timing still selected two modes; the clustered design with a 9.086 s maximum gap and the missing-middle design had different BIC and post-fit recovery outcomes than count-matched regular grids. In the sparse oscillator study (period about 3.143 s), 7 samples over 10 s gave fewer than 2 samples per period on average and all regular/jittered fits were insufficient. With 13 samples (about 3.77 per period), regular sampling estimated 0.07648 s^-1 and strong jitter estimated 0.04033 s^-1; both remained unreliable. |
| H3 - Mode separation, amplitude ratio, duration, sampling, and noise affect distinguishability. | Supported for the tested matrix. | Well-separated and moderate cases were recoverable in fixed examples, while close rates and the weak fast mode were often not selected as two modes across repeated seeds. A strong fast mode with amplitude ratio 10:1 was recoverable in the fixed matrix. When the initial fast component was 0.1 signal units, above 3 sigma (0.045), its analytic contribution fell below that level after about 1.00 s; the fitted two-mode model was selected but missed the 10% rate-recovery tolerance. A 3 s short window and 0.05 noise case also selected two modes but failed rate recovery. No universal threshold follows. |
| H4 - A more complex model can lower residual without identifying extra modes. | Supported. | In the noise-free two-mode check, the two-mode training RMSE was no greater than the one-mode RMSE, but BIC was unavailable because sigma was zero; the result remained inconclusive. In the 3 s noisy case, BIC selected two modes while their fitted rates missed the predeclared recovery tolerance. |
| H5 - Estimation limits are condition-dependent. | Supported within this grid. | The same fitting/model-selection rules behaved differently across rate spacing, amplitude ratios, sample design, duration, and noise. Sparse oscillator timing remained unreliable despite jitter, and there was no single condition boundary that worked across all designs. |

These outcomes apply only to the specified signals, finite parameter grid,
fixed seed lists, fit bounds, start sets, and diagnostics.

## Interpretation and limitations

One recoverable example is the 61-sample well-separated case: rates 0.2 and
0.8 s^-1, amplitudes 1.0 and 0.6, sigma 0.015, with a large positive BIC
difference favoring two modes and both matched rate errors within 10%. The
fixed-seed repetitions show that even this condition did not pass the recovery
criterion in every realization.

One non-distinguishable example is the close-rate (0.2, 0.24) s^-1 case with
equal amplitudes. BIC favored a single mode under the tested design, despite
two modes being present in the generator; this is an incorrect mode-count
selection relative to synthetic truth, not proof of mathematical
non-identifiability. The weak-fast mode group likewise often preferred one
mode or remained inconclusive. The noisy moderate-noise and short-window
cases show another limitation: BIC can prefer the more flexible model while
the fitted rates are inconsistent with known parameters.

The two-mode optimizer uses multiple initializations but can converge to
nearly coincident rates or rates at the bounds. The Jacobian condition,
initialization records, BIC, residuals, and post-fit truth errors are kept
separate for this reason. Sorting rates solves label matching only; it does
not make the physical modes identifiable. The BIC assumes independent
Gaussian errors with known sigma, and is not a confidence score. The selected
10% relative-rate recovery tolerance, BIC difference of six, and condition
screen are study rules, not universal sufficiency criteria. Short records
also weaken asymptotic model-selection assumptions.

Noise is iid additive Gaussian generated by NumPy `default_rng`; it is not
claimed to represent sensor data. The fit residual is RMSE in observation
units on the same samples used for fitting. There are no missing-not-at-random
mechanisms, colored noise, sigma-estimation procedures, or numerical solver
errors because these signals are sampled from analytic formulas. Amplitudes
are positive and the study does not examine cancellation from signed modes.

The existing plotting APIs target ODE/BVP trajectories and experiment metrics;
the Phase 24 output uses tables so truth, fitted parameters, BIC assumptions,
and condition coverage remain explicit without suggesting unsupported
precision. No earlier-phase artifacts or equation records are changed.

This study does not validate a real physical measurement system, establish
physical equivalence between oscillator and diffusion modes, or provide
evidence for denoising, market prediction, trading profitability, or a
financial risk-control method. Financial time series would require separate
evidence for sampling processes, noise distributions, model adequacy, and
out-of-sample performance; this phase makes none of those claims.

## Reproduce and next research question

From the project root with the existing environment:

```powershell
.venv\Scripts\python.exe -c "from newton_lab.modal_identification import run_phase24_research, render_phase24_results; print(render_phase24_results(run_phase24_research()))"
```

The generated run tables and timing summaries follow. A useful next question
is whether profile-likelihood or held-out prediction diagnostics can separate
weak or close modes without turning a single synthetic design's BIC behavior
into a general identifiability rule.\n\n\n\n\n\n

# Phase 24 computed modal-identification results

## Repeated noise realizations

| Condition | Runs | Both fits successful | Failed | Unreliable | Inconclusive | Supported | Inconsistent | Mode estimate and abs/rel error summaries (median, min, max, population SD) |
|---|---:|---:|---:|---:|---:|---:|---:|---|
| single_mode_control | 12 | 12 | 0 | 0 | 1 | 11 | 0 | sorted_mode_1: rate 0.2006/0.197/0.204/0.00188; signed bias 5.05e-05; abs err 0.00121/0.000435/0.00401/0.00108; rel err 0.00607/0.00218/0.0201/0.00538 |
  Seeds: `(2401, 2402, 2403, 2404, 2405, 2406, 2407, 2408, 2409, 2410, 2411, 2412)`. Descriptive summaries of these fixed seeds and this design only; counts are not universal probabilities or confidence intervals.
| well_separated_comparable | 12 | 12 | 0 | 0 | 0 | 8 | 4 | sorted_mode_1: rate 0.1992/0.1936/0.2095/0.00559; signed bias 0.00107; abs err 0.0043/5.23e-06/0.00946/0.00296; rel err 0.0215/2.61e-05/0.0473/0.0148; sorted_mode_2: rate 0.8423/0.7323/0.9498/0.0662; signed bias 0.0307; abs err 0.0551/0.00882/0.15/0.0374; rel err 0.0688/0.011/0.187/0.0468 |
  Seeds: `(2401, 2402, 2403, 2404, 2405, 2406, 2407, 2408, 2409, 2410, 2411, 2412)`. Descriptive summaries of these fixed seeds and this design only; counts are not universal probabilities or confidence intervals.
| moderately_separated | 12 | 12 | 0 | 0 | 0 | 4 | 8 | sorted_mode_1: rate 0.2083/0.1529/0.2336/0.0247; signed bias 0.00136; abs err 0.0194/0.00023/0.0471/0.0133; rel err 0.0971/0.00115/0.236/0.0665; sorted_mode_2: rate 0.4383/0.3207/0.6856/0.0995; signed bias 0.0513; abs err 0.0648/0.0211/0.286/0.0725; rel err 0.162/0.0528/0.714/0.181 |
  Seeds: `(2401, 2402, 2403, 2404, 2405, 2406, 2407, 2408, 2409, 2410, 2411, 2412)`. Descriptive summaries of these fixed seeds and this design only; counts are not universal probabilities or confidence intervals.
| close_rates | 12 | 12 | 0 | 0 | 3 | 0 | 9 | sorted_mode_1: rate 0.2151/6.667e-10/0.2198/0.0951; signed bias -0.0489; abs err 0.0187/0.0126/0.2/0.0796; rel err 0.0936/0.0629/1/0.398; sorted_mode_2: rate 0.2222/0.218/3.333/1.13; signed bias 0.638; abs err 0.0217/0.0162/3.09/1.11; rel err 0.0902/0.0677/12.9/4.64 |
  Seeds: `(2401, 2402, 2403, 2404, 2405, 2406, 2407, 2408, 2409, 2410, 2411, 2412)`. Descriptive summaries of these fixed seeds and this design only; counts are not universal probabilities or confidence intervals.
| weak_fast_mode | 12 | 12 | 0 | 0 | 8 | 0 | 4 | sorted_mode_1: rate 0.1971/6.667e-10/0.2053/0.0737; signed bias -0.0364; abs err 0.00416/0.000646/0.2/0.0727; rel err 0.0208/0.00323/1/0.363; sorted_mode_2: rate 1.016/0.2059/3.333/1.31; signed bias 0.805; abs err 0.593/0.0587/2.53/1.01; rel err 0.741/0.0734/3.17/1.26 |
  Seeds: `(2401, 2402, 2403, 2404, 2405, 2406, 2407, 2408, 2409, 2410, 2411, 2412)`. Descriptive summaries of these fixed seeds and this design only; counts are not universal probabilities or confidence intervals.

## Mode and sampling matrix

| Design | True modes | N | Span (s) | Sigma | A true | Rates true (s^-1) | A/rate single-fit | A two-fit | Rates two-fit (s^-1) | Afast/Aslow | Delta rate | RMSE single | RMSE two | Scaled Jacobian cond. | k single/two | BIC single | BIC two | Delta BIC | Selection | Evaluation | Rate abs/rel errors | Amplitude abs errors |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| single_mode_control | 1 | 61 | 15 | 0.015 | (1) | (0.2) | 0.99218/0.19864 | (0.49125, 0.50093) | (0.19864, 0.19864) | n/a | n/a | 0.014486 | 0.014486 | 1.5199e+16 | 2/4 | 65.116 | 73.337 | -8.2217 | single_mode | supported_under_tested_conditions | 0.001362/0.006809 | 0.007824 |
| well_separated | 2 | 61 | 15 | 0.015 | (1, 0.6) | (0.2, 0.8) | 1.452/0.26787 | (0.99579, 0.59828) | (0.20025, 0.79194) | 0.6 | 0.75 | 0.041616 | 0.016148 | 22.58 | 2/4 | 477.76 | 87.134 | 390.63 | two_modes | supported_under_tested_conditions | 0.0002465/0.001232; 0.008057/0.01007 | 0.004208, 0.001721 |
| moderate_separation | 2 | 61 | 15 | 0.015 | (1, 0.8) | (0.2, 0.4) | 1.7629/0.25878 | (0.8064, 1.0002) | (0.18602, 0.37099) | 0.8 | 0.5 | 0.022768 | 0.015648 | 219.67 | 2/4 | 148.76 | 82.829 | 65.935 | two_modes | supported_under_tested_conditions | 0.01398/0.0699; 0.02901/0.07253 | 0.1936, 0.2002 |
| close_separation | 2 | 61 | 15 | 0.015 | (1, 1) | (0.2, 0.24) | 1.9892/0.21734 | (0.98462, 1.0046) | (0.21734, 0.21734) | 1 | 0.16667 | 0.017189 | 0.017189 | 2.1249e+16 | 2/4 | 88.326 | 96.548 | -8.2217 | single_mode | fitted_model_inconsistent_with_synthetic_truth | 0.01734/0.08669; 0.02266/0.09442 | 0.01538, 0.004592 |
| weak_fast_amplitude | 2 | 61 | 15 | 0.015 | (1, 0.03) | (0.2, 0.8) | 1.0129/0.20054 | (0.51114, 0.50175) | (0.20054, 0.20054) | 0.03 | 0.75 | 0.015216 | 0.015216 | 1.9165e+16 | 2/4 | 70.991 | 79.213 | -8.2217 | single_mode | fitted_model_inconsistent_with_synthetic_truth | 0.0005426/0.002713; 0.5995/0.7493 | 0.4889, 0.4718 |
| fast_mode_crosses_noise_floor | 2 | 61 | 15 | 0.015 | (1, 0.1) | (0.2, 0.8) | 1.0712/0.21056 | (0.99015, 0.12494) | (0.19711, 0.92921) | 0.1 | 0.75 | 0.017502 | 0.014134 | 16.35 | 2/4 | 91.265 | 70.606 | 20.659 | two_modes | fitted_model_inconsistent_with_synthetic_truth | 0.00289/0.01445; 0.1292/0.1615 | 0.009853, 0.02494 |
| strong_fast_amplitude | 2 | 61 | 15 | 0.015 | (0.1, 1) | (0.2, 0.8) | 1.0696/0.66835 | (0.10499, 0.98874) | (0.21706, 0.78982) | 10 | 0.75 | 0.020577 | 0.01664 | 25.83 | 2/4 | 123.01 | 91.515 | 31.499 | two_modes | supported_under_tested_conditions | 0.01706/0.08532; 0.01018/0.01273 | 0.00499, 0.01126 |
| short_window | 2 | 13 | 3 | 0.015 | (1, 0.6) | (0.2, 0.8) | 1.5535/0.33684 | (0.35246, 1.2278) | (0.04344, 0.48697) | 0.6 | 0.75 | 0.018786 | 0.0099025 | 693.68 | 2/4 | 25.521 | 15.925 | 9.5955 | two_modes | fitted_model_inconsistent_with_synthetic_truth | 0.1566/0.7828; 0.313/0.3913 | 0.6475, 0.6278 |
| too_few_samples | 2 | 7 | 3 | 0.015 | (1, 0.6) | (0.2, 0.8) | 1.556/0.33338 | (failed) | (failed) | 0.6 | 0.75 | 0.017019 | n/a | n/a | 2/4 | n/a | n/a | n/a | unavailable | insufficient_observations | n/a | n/a |
| sparse_sampling | 2 | 16 | 15 | 0.015 | (1, 0.6) | (0.2, 0.8) | 1.5435/0.2895 | (1.0021, 0.61132) | (0.21018, 0.71587) | 0.6 | 0.75 | 0.039209 | 0.011132 | 31.137 | 2/4 | 114.87 | 19.903 | 94.963 | two_modes | fitted_model_inconsistent_with_synthetic_truth | 0.01018/0.05088; 0.08413/0.1052 | 0.002066, 0.01132 |
| low_noise | 2 | 61 | 15 | 0.003 | (1, 0.6) | (0.2, 0.8) | 1.4531/0.26724 | (0.98562, 0.61226) | (0.19803, 0.78563) | 0.6 | 0.75 | 0.039356 | 0.0024613 | 22.579 | 2/4 | 10506 | 57.505 | 10449 | two_modes | supported_under_tested_conditions | 0.001965/0.009825; 0.01437/0.01796 | 0.01438, 0.01226 |
| moderate_noise | 2 | 61 | 15 | 0.05 | (1, 0.6) | (0.2, 0.8) | 1.4746/0.26356 | (1.2177, 0.53695) | (0.21934, 2.0112) | 0.6 | 0.75 | 0.069647 | 0.046097 | 6.8999 | 2/4 | 126.58 | 68.292 | 58.287 | two_modes | fitted_model_inconsistent_with_synthetic_truth | 0.01934/0.09671; 1.211/1.514 | 0.2177, 0.06305 |
| noise_free_overfit_check | 2 | 61 | 15 | 0 | (1, 0.6) | (0.2, 0.8) | 1.4551/0.26749 | (1, 0.6) | (0.2, 0.8) | 0.6 | 0.75 | 0.038925 | 1.6285e-17 | 22.076 | 2/4 | n/a | n/a | n/a | unavailable | mode_identification_inconclusive | 0/0; 1.11e-16/1.388e-16 | 0, 1.11e-16 |
| explicit alternative initializations | 2 | 101 | 25 | 0.01 | (1, 0.6) | (0.18, 0.65) | 1.4659/0.23814 | (1.0063, 0.59672) | (0.18034, 0.6732) | 0.6 | 0.72308 | 0.032357 | 0.010387 | 22.301 | 2/4 | 1066.7 | 127.43 | 939.26 | two_modes | supported_under_tested_conditions | 0.0003446/0.001915; 0.0232/0.03569 | 0.006252, 0.003284 |
| regular_baseline | 2 | 61 | 15 | 0.015 | (1, 0.6) | (0.2, 0.8) | 1.455/0.26582 | (0.97986, 0.6153) | (0.19662, 0.75988) | 0.6 | 0.75 | 0.040994 | 0.01389 | 23.953 | 2/4 | 463.82 | 68.746 | 395.07 | two_modes | supported_under_tested_conditions | 0.003378/0.01689; 0.04012/0.05015 | 0.02014, 0.0153 |
| regular_matched_mild_jitter | 2 | 61 | 15 | 0.015 | (1, 0.6) | (0.2, 0.8) | 1.455/0.26582 | (0.97986, 0.6153) | (0.19662, 0.75988) | 0.6 | 0.75 | 0.040994 | 0.01389 | 23.953 | 2/4 | 463.82 | 68.746 | 395.07 | two_modes | supported_under_tested_conditions | 0.003378/0.01689; 0.04012/0.05015 | 0.02014, 0.0153 |
| mild_jitter | 2 | 61 | 15 | 0.015 | (1, 0.6) | (0.2, 0.8) | 1.455/0.26576 | (0.97812, 0.61651) | (0.19643, 0.75701) | 0.6 | 0.75 | 0.040847 | 0.013882 | 24.21 | 2/4 | 460.57 | 68.686 | 391.88 | two_modes | supported_under_tested_conditions | 0.003572/0.01786; 0.04299/0.05373 | 0.02188, 0.01651 |
| regular_matched_strong_jitter | 2 | 61 | 15 | 0.015 | (1, 0.6) | (0.2, 0.8) | 1.455/0.26582 | (0.97986, 0.6153) | (0.19662, 0.75988) | 0.6 | 0.75 | 0.040994 | 0.01389 | 23.953 | 2/4 | 463.82 | 68.746 | 395.07 | two_modes | supported_under_tested_conditions | 0.003378/0.01689; 0.04012/0.05015 | 0.02014, 0.0153 |
| strong_jitter | 2 | 61 | 15 | 0.015 | (1, 0.6) | (0.2, 0.8) | 1.4578/0.26652 | (0.98008, 0.61523) | (0.19666, 0.75987) | 0.6 | 0.75 | 0.041214 | 0.013889 | 24.072 | 2/4 | 468.73 | 68.739 | 399.99 | two_modes | supported_under_tested_conditions | 0.003335/0.01668; 0.04013/0.05016 | 0.01992, 0.01523 |
| regular_matched_clustered_long_gap | 2 | 61 | 15 | 0.015 | (1, 0.6) | (0.2, 0.8) | 1.455/0.26582 | (0.97986, 0.6153) | (0.19662, 0.75988) | 0.6 | 0.75 | 0.040994 | 0.01389 | 23.953 | 2/4 | 463.82 | 68.746 | 395.07 | two_modes | supported_under_tested_conditions | 0.003378/0.01689; 0.04012/0.05015 | 0.02014, 0.0153 |
| clustered_long_gap | 2 | 61 | 15 | 0.015 | (1, 0.6) | (0.2, 0.8) | 1.5338/0.31994 | (1.0235, 0.5793) | (0.20101, 0.82862) | 0.6 | 0.75 | 0.03871 | 0.01389 | 33.779 | 2/4 | 414.47 | 68.753 | 345.72 | two_modes | supported_under_tested_conditions | 0.001008/0.005038; 0.02862/0.03577 | 0.02345, 0.0207 |
| regular_matched_missing_middle | 2 | 48 | 15 | 0.015 | (1, 0.6) | (0.2, 0.8) | 1.4599/0.26657 | (0.96098, 0.63177) | (0.19404, 0.73734) | 0.6 | 0.75 | 0.041665 | 0.013677 | 24.96 | 2/4 | 378.08 | 55.389 | 322.69 | two_modes | supported_under_tested_conditions | 0.005961/0.0298; 0.06266/0.07832 | 0.03902, 0.03177 |
| missing_middle | 2 | 48 | 15 | 0.015 | (1, 0.6) | (0.2, 0.8) | 1.4752/0.27744 | (0.94137, 0.65131) | (0.19192, 0.71837) | 0.6 | 0.75 | 0.042615 | 0.01361 | 26.065 | 2/4 | 395.17 | 55 | 340.17 | two_modes | fitted_model_inconsistent_with_synthetic_truth | 0.008076/0.04038; 0.08163/0.102 | 0.05863, 0.05131 |

## Sampling coverage

| Design | Count | Span (s) | Min dt (s) | Median dt (s) | Max dt (s) |
|---|---:|---:|---:|---:|---:|
| regular_baseline | 61 | 15 | 0.25 | 0.25 | 0.25 |
| regular_matched_mild_jitter | 61 | 15 | 0.25 | 0.25 | 0.25 |
| mild_jitter | 61 | 15 | 0.16249 | 0.24707 | 0.34124 |
| regular_matched_strong_jitter | 61 | 15 | 0.25 | 0.25 | 0.25 |
| strong_jitter | 61 | 15 | 0.11024 | 0.24974 | 0.4057 |
| regular_matched_clustered_long_gap | 61 | 15 | 0.25 | 0.25 | 0.25 |
| clustered_long_gap | 61 | 15 | 0.085714 | 0.085714 | 9.0857 |
| regular_matched_missing_middle | 48 | 15 | 0.31915 | 0.31915 | 0.31915 |
| missing_middle | 48 | 15 | 0.25 | 0.25 | 3.5 |
| oscillator_regular | 81 | 10 | 0.125 | 0.125 | 0.125 |
| oscillator_mild_jitter | 81 | 10 | 0.080928 | 0.1229 | 0.16732 |
| oscillator_clustered | 81 | 10 | 0.042553 | 0.042553 | 6.0426 |
| oscillator_regular_matched_missing_middle | 64 | 10 | 0.15873 | 0.15873 | 0.15873 |
| oscillator_missing_middle | 64 | 10 | 0.125 | 0.125 | 2.25 |
| oscillator_sparse_7_regular | 7 | 10 | 1.6667 | 1.6667 | 1.6667 |
| oscillator_sparse_7_mild_jitter | 7 | 10 | 1.1749 | 1.7009 | 1.9121 |
| oscillator_sparse_7_strong_jitter | 7 | 10 | 1.0139 | 1.4944 | 2.9756 |
| oscillator_sparse_7_clustered | 7 | 10 | 0.66667 | 0.66667 | 6.6667 |
| oscillator_sparse_13_regular | 13 | 10 | 0.83333 | 0.83333 | 0.83333 |
| oscillator_sparse_13_mild_jitter | 13 | 10 | 0.5919 | 0.83011 | 1.0865 |
| oscillator_sparse_13_strong_jitter | 13 | 10 | 0.58119 | 0.75432 | 1.4078 |
| oscillator_sparse_13_clustered | 13 | 10 | 0.28571 | 0.34286 | 6.2857 |

## Oscillator timing comparison

- `regular`: n=81, span=10 s, available_interpretable_under_tested_assumptions; estimated envelope rate 0.05047036 s^-1.
- `mild_jitter`: n=81, span=10 s, available_interpretable_under_tested_assumptions; estimated envelope rate 0.05077429 s^-1.
- `clustered`: n=81, span=10 s, insufficient_observations; estimated envelope rate unavailable.
- `regular_matched_missing_middle`: n=64, span=10 s, available_interpretable_under_tested_assumptions; estimated envelope rate 0.05161443 s^-1.
- `missing_middle`: n=64, span=10 s, available_interpretable_under_tested_assumptions; estimated envelope rate 0.05047141 s^-1.
- `sparse_7_regular`: n=7, span=10 s, insufficient_observations; estimated envelope rate unavailable.
- `sparse_7_mild_jitter`: n=7, span=10 s, insufficient_observations; estimated envelope rate unavailable.
- `sparse_7_strong_jitter`: n=7, span=10 s, insufficient_observations; estimated envelope rate unavailable.
- `sparse_7_clustered`: n=7, span=10 s, insufficient_observations; estimated envelope rate unavailable.
- `sparse_13_regular`: n=13, span=10 s, available_but_unreliable_under_declared_diagnostics; estimated envelope rate 0.07647902 s^-1.
- `sparse_13_mild_jitter`: n=13, span=10 s, available_but_unreliable_under_declared_diagnostics; estimated envelope rate 0.07840383 s^-1.
- `sparse_13_strong_jitter`: n=13, span=10 s, available_but_unreliable_under_declared_diagnostics; estimated envelope rate 0.0403265 s^-1.
- `sparse_13_clustered`: n=13, span=10 s, insufficient_observations; estimated envelope rate unavailable.

## Initial-guess sensitivity

- start A=(0.9, 0.7), rates=(0.12, 0.55) s^-1: converged=True; `ftol` termination condition is satisfied.; fitted rates 0.1803446, 0.6732013 s^-1.
- start A=(0.7, 0.9), rates=(0.55, 0.12) s^-1: converged=True; `ftol` termination condition is satisfied.; fitted rates 0.1803446, 0.6732013 s^-1.
- start A=(1, 0.5), rates=(0.3, 0.9) s^-1: converged=True; `ftol` termination condition is satisfied.; fitted rates 0.1803446, 0.6732013 s^-1.
- start A=(0.4, 1.1), rates=(0.02, 0.4) s^-1: converged=True; `ftol` termination condition is satisfied.; fitted rates 0.1803446, 0.6732013 s^-1.

## Selection assumptions

- BIC uses the iid additive Gaussian likelihood with known positive noise variance declared in ObservationSeries; only shared constants are omitted.
- Two-mode selection requires BIC(single)-BIC(two) >= 6; a difference <= -6 selects one mode; intermediate values are inconclusive.
- The per-rate recovery tolerance is 10% and is used only after fitting to summarize agreement with synthetic truth.
- Noise-free cases have no BIC selection; lower residual in a four-parameter fit is reported but cannot establish an additional mode.
- Rate-separation Delta=abs(r2-r1)/max(r1,r2) is descriptive, not a universal threshold.
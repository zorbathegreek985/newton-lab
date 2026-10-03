# Phase 30: Causal diffusion-inspired filtering

## Research question and preregistered hypotheses

Can a causal one-sided filter suppress high-frequency variation while retaining slower components, and what delay and transition costs appear? This is a synthetic signal-processing study, not a trading or prediction test.

H1: measured sinusoidal response matches the transfer-function magnitude and phase within 0.003. H2: outputs do not anticipate later samples. H3: increasing beta broadens the step rise and increases the DC group delay over the tested sweep; nonzero-frequency group delay is frequency-dependent and need not be monotonic in beta. H4: at least one beta reduces error on the designed mode-12 contamination. H5: diffusion and the causal filter have no automatic global parameter equivalence.

## Filter and analytical response

For 0 <= beta < 1, y[0]=x[0] and y[n]=(1-beta)x[n]+beta*y[n-1]. This recurrence reads no index later than n. Its z-transform transfer function is H(z)=(1-beta)/(1-beta z^-1). For omega in radians/sample, |H|=(1-beta)/sqrt(1+beta^2-2 beta cos(omega)); phase is -atan2(beta sin(omega),1-beta cos(omega)). Group delay is beta(cos(omega)-beta)/(1+beta^2-2 beta cos(omega)); zero-frequency group delay is beta/(1-beta). Phase delay is -phase/omega and is undefined at DC. Step threshold delay is a separate finite-record measurement.
With zero-state transform notation, Y(z)=(1-beta)X(z)+beta*z^-1*Y(z). Solving for Y/X gives the stated transfer function. On the unit circle, the denominator magnitude is sqrt((1-beta*cos(omega))^2 + (beta*sin(omega))^2), yielding the magnitude expression. The reported phase uses H(e^(i*omega)) and is negative for a lag.

Initialization y[0]=x[0] avoids a zero-state startup jump, but the first observation is passed through unchanged; subsequent outputs contain a decaying initialization transient. The filter is causal on its finite record by construction.

For step metrics, the transition is located halfway between the last pre-change sample and the first post-change sample. This makes the beta=0 threshold delay zero; pre-transition samples are checked separately.

## Numerical frequency response

Sinusoids use beta=0.9, 8192 samples, and a least-squares sine/cosine fit after discarding 1500 startup samples. These angular frequencies are discrete radians per sample.

| omega | measured gain | analytic gain | gain error | measured phase | analytic phase | phase error |
|---:|---:|---:|---:|---:|---:|---:|
| 0.15708 | 0.557616 | 0.557616 | 1.67e-15 | -0.902815 | -0.902815 | 8.88e-16 |
| 0.785398 | 0.136436 | 0.136436 | 3.33e-16 | -1.05171 | -1.05171 | 1.78e-15 |
| 2.35619 | 0.0569545 | 0.0569545 | 5.55e-17 | -0.370902 | -0.370902 | 3.11e-15 |

## Relationship to Phase 29 diffusion

Phase 29 applies a noncausal, whole-interval sine transform with exact modal factor exp(-(n*pi)^2*s) on the unit interval. It uses both sides of a finite record through its modal projection; it is not online filtering. Here the sampled spatial coordinate is xi=j/M with M=sample_count-1, so mode n maps to discrete angular frequency omega=n*pi/M. For the aligned modes below, the causal magnitude is matched separately to the diffusion factor. The differing beta values demonstrate that a single beta does not exactly reproduce these multiple modal attenuations. This mapping compares amplitudes at selected frequencies only; it does not match phase, boundary treatment, or causality.

| mode n | omega (rad/sample) | exposure s | diffusion gain | beta matching this mode |
|---:|---:|---:|---:|---:|
| 1 | 0.006283185 | 0.001 | 0.9901789 | 0.9564784 |
| 4 | 0.02513274 | 0.001 | 0.8539235 | 0.9596023 |
| 12 | 0.07539822 | 0.001 | 0.2414186 | 0.9814221 |

## Controlled contamination and beta sweep

The deterministic Phase 29 signal is truth sin(pi*xi) plus a mode-12 perturbation of amplitude 0.35, sampled at 501 points. No random noise is used. The same input realization is compared with the unchanged Phase 29 diffusion output at s=0.001, the causal recurrence, and a causal trailing moving average of width 21 samples (partial renormalized windows at startup). Raw RMSE to truth is 0.24724; Phase 29 diffusion RMSE is 0.0600902; trailing-mean RMSE is 0.222299. The sweep table reports exponential-filter errors. These are synthetic reconstruction errors, not a denoising guarantee. The moving average and exponential filter are not asserted to have equal bandwidth.

| index | beta | status | filtered RMSE to truth | RMSE to input | slow retention | fast retention | step 10–90 width (samples) | half delay (samples) | plateau | pre-transition max | failure |
|---:|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---|
| 0 | 0 | succeeded | 0.24724 | 0 | 1 | 1 | 0.8 | 0 | 1 | 0 | — |
| 1 | 0.5 | succeeded | 0.24542 | 0.018999 | 1.00017 | 0.989711 | 3.2 | 0.5 | 1 | 0 | — |
| 2 | 0.9 | succeeded | 0.201637 | 0.139974 | 0.998516 | 0.6723 | 20.8608 | 6.09162 | 0.999536 | 0 | — |
| 3 | 0.98 | succeeded | 0.180965 | 0.29777 | 0.934084 | 0.0426073 | 108.757 | 33.8118 | 0.845222 | 0 | — |
| 4 | 0.9 | succeeded | 0.201637 | 0.139974 | 0.998516 | 0.6723 | 20.8608 | 6.09162 | 0.999536 | 0 | — |
| 5 | -0.1 | failed | — | — | — | — | — | — | — | — | beta must be a finite number in [0, 1) |
| 6 | 1 | failed | — | — | — | — | — | — | — | — | beta must be a finite number in [0, 1) |

## Hypothesis outcomes

- H1: supported for these three frequencies; maximum magnitude error 1.67e-15, maximum phase error 3.11e-15.
- H2: supported by direct prefix non-anticipation tests and recurrence; scope is the declared discrete-time input model.
- H3: supported for the DC group-delay formula and measured step widths in the tested sweep; nonzero-frequency group delay is frequency-dependent and is not generally monotonic in beta. There is no universal one-number delay claim.
- H4: supported for this deterministic mode-12 perturbation; the finding is restricted to this signal and metric.
- H5: supported for the displayed aligned modes: separately matched beta values differ. Amplitude matching at isolated frequencies does not establish equivalence of operators.

## Interpretation, limitations, and reproducibility

A causal feature can in principle be computed from observations available through the current index. This establishes only information timing. It does not establish predictive value, a trading signal, profitability, risk-adjusted returns, regime robustness, or performance after costs. No market data, financial outcomes, or backtest are used.

The experiment uses deterministic synthetic values, a finite grid, least-squares sinusoid fits, and one structured contamination. Diffusion uses a unit-interval spatial-mode mapping and whole-series modal projection; filter initialization and finite-sample transients affect early response. The numerical comparisons do not cover arbitrary signals or boundary conditions. The trailing moving average is a separately defined causal baseline, but its bandwidth is not matched to either smoother.

Reproduce the study in the project environment:

```python
from newton_lab.causal_diffusion_filter import run_causal_filter_study

study = run_causal_filter_study()
print(study.report_markdown)
```

The study is self-contained and is not registered as a general Newton Lab experiment record; integration would require a new typed result model.

### Next research question

Under a predeclared stationary signal family, can a trailing moving average and the exponential filter be compared at matched measured attenuation and latency, with repeated noise realizations and explicit edge handling? That would isolate operator choice without implying financial utility.

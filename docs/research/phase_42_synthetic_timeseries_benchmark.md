# Phase 42 — Synthetic Time-Series Research Benchmark

> **All observations in this experiment are synthetic. The results do not
> establish real-market forecasting performance, trading utility, or
> profitability.** No NSE data was read or used.

## Research question and scope

Can simple descriptive diagnostics identify known behaviours in deterministic,
controlled signals, and where do they confuse distinct generating processes?
This benchmark evaluates a lag-one autocorrelation flag, a two-half mean
difference, a two-half variance ratio, and the existing causal exponential
filter. It tests recognition of designed properties, not prediction of future
observations.

## Synthetic data-generating processes

All signals contain 512 samples. Each configuration has three fixed random
seeds. There are 16 configurations and 48 runs. Seeds follow
`42000 + configuration_index * 10 + replicate_index`; NumPy's
`default_rng` produces the Gaussian draws. The generator is
`newton_lab.synthetic_timeseries_benchmark.generate_signal`.

| Process | Definition and varied parameter values |
|---|---|
| Stationary noise | iid Gaussian observations, mean 0 and standard deviation 1; three seeds. |
| Trend plus noise | latent `0.002 * sample_index` plus iid Gaussian noise; noise SD 0.25, 0.75, or 1.5. |
| Autocorrelated | stationary AR(1), `z[t] = 0.8 z[t-1] + epsilon[t]`, unit innovation SD, 200-step burn-in; independent observation-noise SD 0.1, 0.5, or 1.0. |
| Damped oscillation | `exp(-0.002 t) sin(2 pi * 0.10 t)` plus iid Gaussian noise; noise SD 0.05, 0.2, or 0.5. |
| Variance shift | zero-mean Gaussian noise with SD 1 before the midpoint and SD multiplier 1.1, 1.4, or 2.0 after it. |
| Mean shift | latent mean 0 before and 0.5, 1.0, or 2.0 after the midpoint, with iid unit-SD Gaussian noise. |

For the filter error comparison only, the known deterministic latent component
is retained as synthetic truth. For the stationary noise and variance-shift
processes that truth is zero. This is not an available reference in ordinary
observational data.

## Methods and evaluation procedure

The lag-one statistic centers the full observed record and computes the lag-one
product sum divided by the full centered sum of squares. The diagnostic flag is
`abs(lag-1 autocorrelation) >= 0.25`.

The mean-shift diagnostic compares the first and second halves, which are split
at the known design midpoint. Its statistic is the difference of half means
divided by a standard error formed from the two within-half sample variances;
the flag is `abs(statistic) >= 3`. It is a fixed-location diagnostic, not a
general change-point search. Its standard error does not adjust for serial
dependence.

The variance diagnostic divides the larger within-half sample variance by the
smaller and flags ratios at least 2. Neither change diagnostic searches for an
unknown transition location. Thresholds are transparent demonstration rules,
not calibrated significance tests.

The filter comparison reuses
`newton_lab.causal_diffusion_filter.exponential_filter` with beta 0.8 and the
documented recurrence `y[0]=x[0]`,
`y[t]=(1-beta)*x[t] + beta*y[t-1]`. For synthetic runs only, input and filtered
RMSE are measured against the known latent component. An improved RMSE is a
post-hoc reconstruction metric, not forecasting accuracy or proof that the
filter denoises arbitrary signals.

The fixed grid is run once and not tuned against a separate evaluation set.
Multiple configurations and replicates expose seed/parameter sensitivity, but
three seeds are only a small illustration, not a statistical power analysis.
Run with the project environment from the workspace root:

```powershell
.venv\Scripts\python.exe -m newton_lab.synthetic_timeseries_benchmark
```

The command writes `reports/phase_42_synthetic_benchmark/results.csv`,
`metadata.json`, and the figures listed below. The CSV retains one row per run
and records all diagnostic statistics and flags.

## Results

The intended AR(1) case crossed the autocorrelation threshold in 9/9 runs.
However, the same flag also occurred in 9/9 damped-oscillation runs, 3/9 trend
runs, and 5/9 mean-shift runs. A large lag-one statistic therefore did not
uniquely identify an AR-generating process in this benchmark.

The two-half mean diagnostic flagged all 9/9 designed mean-shift runs, but also
8/9 trend runs and 4/9 AR(1) runs. A global drift, dependence, and a discrete
mean change can all produce a large two-half mean statistic. The statistic
detects a difference between halves, not its cause.

The variance ratio flagged 5/9 variance-shift runs overall: 0/3 at post-change
SD multiplier 1.1, 2/3 at 1.4, and 3/3 at 2.0. It also flagged 5/9 damped
oscillations: a decaying deterministic amplitude can appear as changing
half-record variance even without a random variance regime change. This
illustrates both severity dependence and a false-positive mechanism.

At beta 0.8, filtering reduced truth RMSE in 3/3 stationary-noise runs, 9/9
trend runs, 0/9 AR(1) runs, 1/9 damped-oscillation runs, 9/9 variance-shift
runs, and 9/9 mean-shift runs. In mean-shift and variance-shift cases, this
single aggregate RMSE can conceal transition smearing; it does not establish
good event timing. The AR signal's genuine structure was also attenuated, so
removing input variation is not always improvement against the designed
signal.

These results are observations from the declared finite grid. They are not
universal properties of the methods. The configuration-level values and every
seed result are in `reports/phase_42_synthetic_benchmark/results.csv`.

## Figures

- [Detection rates by process](../../reports/phase_42_synthetic_benchmark/detection_rates.png)
  shows how intended detections and false positives vary across signal kinds.
- [Diagnostic overlap](../../reports/phase_42_synthetic_benchmark/diagnostic_overlap.png)
  plots lag-one autocorrelation against the two-half mean statistic. Similar
  coordinates can arise from different generating processes.
- [Filter reconstruction error](../../reports/phase_42_synthetic_benchmark/filter_error.png)
  compares observed and filtered RMSE to known synthetic truth. It is not a
  market-data or forecasting metric.

## Failures, sensitivity, and interpretation

The small variance multiplier was missed in every replicate, while stronger
variance changes were increasingly detected. Mean changes had complete
detection at the chosen settings, but the same diagnostic frequently flagged
trends and autocorrelated signals. The lag-one threshold separated the white
noise baseline from AR(1) in these seeds, yet damped oscillations always crossed
it. A statistic that flags an observed pattern does not identify the unique
underlying process.

Filtering improved the measured reference error in the white-noise baseline,
but the chosen beta and RMSE aggregated across all samples. Alternative beta,
sample lengths, transition positions, amplitudes, and noise laws could change
the result. Filtering is causal by implementation, but causal availability
alone does not imply forecasting value.

## Limitations

- The models are deliberately simple Gaussian constructions on equally spaced
  samples. They omit missing data, outliers, heavy tails, seasonality, and
  irregular observation timing.
- Three seeds per setting cannot characterize rare outcomes or provide precise
  false-positive probabilities.
- Fixed midpoint comparisons have the transition location supplied by the
  generator. They are not unknown-time change-point estimators.
- The mean statistic's standard error assumes independent halves; its flags are
  especially difficult to interpret for autocorrelated observations.
- The thresholds were chosen as simple fixed demonstration criteria, not
  calibrated on an independent null distribution.
- Synthetic latent truth enables reconstruction-error comparisons; real
  market observations do not supply such a known clean reference.
- Nothing here resolves the separate data permission questions. Phase 32 stays
  `BLOCKED_AUTHORIZATION` as recorded in
  [the readiness gate](phase_32_study_readiness.md).

## Conclusions supported by this experiment

The selected diagnostics recover several intentionally designed patterns at
these settings, but they also confuse trends, oscillations, shifts, and
dependence. The detector outputs are descriptive evidence about these generated
records, not identification of a unique system. Filter performance depends on
the signal and metric; smoothing and preserving structure can conflict.
No result supports real-market forecasting, trading, investment, or
profitability claims.

## Concrete next research question

With a larger, preregistered synthetic grid and independent calibration seeds,
how do simple autocorrelation and fixed-window mean/variance diagnostics trade
off missed weak changes against false alarms under serial dependence and
decaying oscillations?

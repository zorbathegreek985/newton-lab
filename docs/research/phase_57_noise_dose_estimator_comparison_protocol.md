# Phase 57 â€” Noise-Dose and Estimator-Form Comparison Protocol

**Protocol version:** 1.0
**Status:** Frozen before implementation and held-out evaluation.
**Executable configuration:** [`protocol.json`](../../reports/phase_57_noise_dose_estimator_comparison/protocol.json)

## Question and estimands

For the Phase 53 fold, transcritical, and positive supercritical-pitchfork branches, how do additive observation-noise dose and estimator form affect recovery-rate error, and do rate-estimation differences translate into recovery-rate-exponent identifiability?

The primary estimands are model- and amplitude-specific (a) the mean signed and absolute relative error of the recovery-rate estimates, retaining invalid-fit proportions; (b) paired differences in absolute relative error between two observation-only estimator forms at each noise dose; and (c) their difference from the zero-noise estimator contrast (the noise-by-estimator interaction). Exponent identifiability is a distinct secondary outcome, estimated only for the two observation-only methods. The correct-normal-form fit is a separately labelled model-informed rate diagnostic and is excluded from the primary estimator interaction and exponent analysis.

No estimator receives the true equilibrium, margin, bifurcation boundary, recovery rate, or noise standard deviation. The named model is supplied only to the diagnostic normal-form fit; that fit is not treated as an equal-information competitor.

## Frozen generator and settings

The equations, branches, margins and theoretical rates are those of Phase 53, for \(\mu>0\):

| Model | Selected stable branch | Local recovery rate | (q,p) in x* ~ Î¼^q, Ï ~ Î¼^p |
|---|---|---|---|
| Fold: xdot = Î¼ âˆ’ xÂ² | x* = +sqrt(Î¼) | Ï = 2sqrt(Î¼) | (1/2, 1/2) |
| Transcritical: xdot = Î¼x âˆ’ xÂ² | x* = Î¼ | Ï = Î¼ | (1, 1) |
| Supercritical pitchfork: xdot = Î¼x âˆ’ xÂ³ | x* = +sqrt(Î¼) | Ï = 2Î¼ | (1/2, 1) |

Use margins Î¼ = (0.04, 0.09, 0.16, 0.25, 0.36), settings u = 1 + Î¼, and initial state x(0) = x*(1+a), with fixed amplitudes a = (0.0025, 0.04). These are the smallest and largest Phase 55 sub-5% amplitudes. Dynamics are deterministic and use the existing exact trajectories. Add independent Gaussian measurement noise to each sampled state, with Ïƒ = Î· x*, at dose fractions Î· = (0, 0.0025, 0.01, 0.05). Zero is a noiseless diagnostic; 0.0025 makes Ïƒ equal the initial displacement at the smallest amplitude; 0.01 is Phase 55's dose; 0.05 is Phase 54's high dose. No dose changes the dynamics.

### Fixed observation times (leakage correction)

Use the same absolute grid t = 0, 0.1, ..., 75.0 seconds (751 points) for every model, margin, dose, amplitude, and seed. Phase 55 and the Phase 56 sketch sampled at intervals proportional to the true recovery time. Since timestamps are estimator inputs, that schedule can reveal the true rate from the sample spacing. This protocol deliberately replaces the proposed 3Ï„ grid with a fixed grid to avoid that hidden information channel. The 75-second horizon spans 3Ï„ for the slowest included rate (transcritical Î¼=0.04, Ï=0.04); the 0.1-second interval gives over eight samples per Ï„ at the fastest included rate (fold Î¼=0.36, Ï=1.2). Relative window length consequently differs among settings. This is a known design limitation, reported by model and margin rather than hidden by pooling.

For each seed/model/margin, generate one standardized normal innovation vector and reuse it across amplitudes and all nonzero dose levels. For zero dose, the observed array is exactly the clean trajectory; seed-labelled zero-dose records are identical. Estimators receive only the shared timestamps and that condition's observed values, plus an observation-derived positive numerical initialization scale. The scale is max(0.001 Ã— observed range, 1eâˆ’12); it is not a noise estimate and cannot encode generator truth. The model-informed diagnostic also receives the declared model name.

## Estimators

1. **Free-offset exponential (reference):** Phase 53/55 fit `b + A exp(âˆ’Ït)`, with free offset, positive amplitude and positive rate. It estimates the equilibrium as `b`. The initialization scale is derived from the observed range by the common rule above. A valid result requires optimizer convergence and a finite positive rate; equilibrium and rate estimates are additionally recorded for parameter-recovery checks.
2. **Tail-offset log-linear (observation-only alternative):** estimate the offset as the arithmetic mean of the final 10% of samples (76 points); form residuals `y âˆ’ offset`; fit ordinary least squares to `log(residual)` against time using positive residual entries. Validity requires at least five positive residuals, a finite strictly negative slope, and a finite positive rate `âˆ’slope`. Record the offset and number of included points. The positive-residual selection is part of this estimator, not data cleaning.
3. **Correct-normal-form nonlinear fit (separate diagnostic):** use the existing Phase 53 named-model trajectory fit to estimate positive margin from the observed trajectory and then derive equilibrium and local rate from that fitted margin. It receives the model family, but no true state/parameter. This model-informed method is reported only as a diagnostic because its structural prior differs from the two model-agnostic fits. It does not enter exponent-identifiability criteria.

No estimator is selected or tuned on evaluation results. Calibration seeds are used only for execution, input-validation, and output-schema checks; no thresholds or fit settings are calibrated.

## Replicates and pairing

Calibration seeds are 58000â€“58007; held-out evaluation seeds are 58050â€“58099. Sets are disjoint. Pair by seed, model, margin, and amplitude across noise doses; for each observed trace, all applicable estimators receive the identical arrays. The zero-dose output is deterministic and repeated in seed-indexed rows; do not interpret those duplicate observations as independent zero-noise sampling uncertainty. The seed bootstrap at zero dose is omitted. Paired contrasts against zero dose bootstrap the nonzero seed-level paired differences.

Per partition, the design has 3 models Ã— 5 margins Ã— 2 amplitudes Ã— 4 doses Ã— seed count = 120 Ã— seed-count traces. Thus it generates 960 calibration traces and 6,000 evaluation traces. Each trace is sent to three rate estimators, giving 2,880 calibration and 18,000 evaluation fit attempts. The two observation-only estimators produce 2 Ã— 3 Ã— 2 Ã— 4 Ã— 50 = 2,400 evaluation exponent groups (each group uses five margins). Each group is one seed-level replicate; margins and observations within it are not treated as independent replicates.

## Validity, outcomes, uncertainty

Retain one result row for every attempted estimator fit. Record convergence, validity, and a reason code for invalidity; never replace invalid errors with zero or remove them from validity denominators. For each model/amplitude/dose/method, report attempted and valid counts, validity, mean signed relative rate error, mean absolute relative error, median and 90th percentile absolute relative error, rate-error distribution summaries, equilibrium relative error where estimated, and 95% seed-bootstrap intervals for mean errors. Error summaries are conditional on valid fits and always accompany validity counts. The exact theoretical Phase 53 rate is used only after fitting to score estimates.

Use 2,000 deterministic percentile-bootstrap resamples of evaluation seeds for condition means, paired rate-error differences and estimator-interaction differences. Resample seeds as clusters. Paired estimator contrasts use seeds valid for both methods and report the matched denominator; noise contrasts use within-seed differences and require both compared fits valid. At zero dose, report exact deterministic point estimates without a sampling interval. Intervals are descriptive and unadjusted across the prespecified model/condition family; no p-value or familywise claim is made. Do not interpret a missing-zero interval as evidence of certainty.

For each model/amplitude/dose and each observation-only estimator, form one exponent group per evaluation seed from its five estimated equilibria and five rates using the existing unknown-boundary scaling fit. Group validity requires five valid positive rates, five valid positive equilibria, and a converged finite scaling fit. Report attempted/valid groups, q and p estimates, estimated boundary, mean bias/absolute error, 95% seed-bootstrap intervals, interval width, and theoretical-exponent coverage. Preserve invalid groups in the denominator. Apply the Phase 55 p-specific criterion unchanged: at least 80% valid groups, mean absolute p error â‰¤ 0.20, and the 95% interval includes theoretical p. Report q separately; it is not part of this p-specific criterion. No criterion is applied to the model-informed diagnostic because its named equation supplies the scaling structure.

The descriptive SNR is RMS of the clean sampled excursion from true equilibrium divided by Ïƒ; it is infinite at zero noise. It is generator-only, not an estimator input, receives no interval, and is not interpreted causally.

## Comparisons, stopping and success rules

Report model-specific dose contrasts (each nonzero dose minus zero) for each observation-only estimator and amplitude. The primary estimator contrast is tail-offset minus free-offset mean absolute relative error at each dose. The interaction contrast is the estimator contrast at dose Î· minus that contrast at zero dose. Use paired seed bootstrap intervals and matched denominators. A contrast is **supported in the observed direction** only when its 95% interval excludes zero; otherwise call it inconclusive/not detected under this design. This is descriptive interval logic, with no multiplicity-adjusted significance claim. Do not label failure to detect an effect as equivalence. Do not declare an estimator universally superior.

The p-identifiability criterion above can pass or fail separately by condition. Phase 57 is a completed-method benchmark if all frozen cells are run, all attempted records and failure reasons are preserved, the outputs validate and reproduce, and the report distinguishes supported from unresolved contrasts. It does not require a preferred estimator or a successful p criterion. Stop without scientific conclusions if serialization, protocol validation, or deterministic reproduction fails; repair implementation defects without changing the frozen design, then rerun all outputs.

## Reproducibility and known limits

The frozen JSON is hashed into metadata. Run from repository root with `.venv\Scripts\python.exe -c "from pathlib import Path; from newton_lab.noise_dose_estimator_comparison import run_study; print(run_study(Path('reports/phase_57_noise_dose_estimator_comparison')))"`. Run twice from a clean output directory and compare hashes for every prescribed output; metadata records software versions but no wall-clock timestamp. The manifest lists every output except itself. Tests use synthetic inputs only.

Claims are conditional on exact deterministic normal forms, their selected positive branches, five margins, two amplitudes, Gaussian iid observation noise, and the fixed 75-second schedule. Different effective numbers of recovery times across margins can affect fit precision; that is intentionally visible. No process noise, drift, sensor response, missingness, model mismatch, real-system observations, application validation, market data, or trading outcome is studied. Exponents are encoded in the generator and are not discovered here.

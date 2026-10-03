# Phase 32: Study readiness and data authorization gate

## Purpose and scope

This document is the fillable preflight specification for one initial financial
time-series study using the Phase 30 causal exponential filter. Its purpose is
to expose unresolved decisions and establish a documented go/no-go state before
any real data is ingested.

**No real financial data was accessed, downloaded, scraped, or processed in this
phase.** It did not contact a provider, use credentials, access IEX endpoints,
run a backtest, or establish that any permission has been granted. Do not infer
rights from public accessibility or technical access.

The current readiness state is **`BLOCKED_AUTHORIZATION`**. No provider, terms,
permission evidence, or documented review of intended use is present. Study
design choices are also incomplete; these are recorded separately below.

`READY_FOR_DATA_INGESTION` would mean only that the specified authorization
evidence and pre-ingestion decisions have been recorded and reviewed for this
study. It would not certify legal sufficiency or establish data quality,
scientific validity, forecasting value, or financial performance.

## Relationship to prior phases

Phase 31 established the general protocol: chronological training, validation,
and untouched test periods; an information-availability rule; train-only fitted
preprocessing; a small baseline set; separate signal-processing and forecasting
tracks; and a compact experiment record. It left provider, dataset, instrument,
frequency, target, dates, filter settings, split boundaries, and metrics
unresolved.

Phase 30's actual candidate API is
`exponential_filter(values, beta) -> NDArray[np.float64]`. It initializes with
`y[0] = x[0]` and rejects non-finite inputs. It does not interpret timestamps,
missing market observations, data revisions, or provider permissions. Phase 29's
whole-record diffusion transform is noncausal and is not an eligible real-time
feature. Neither phase supplies financial-data evidence.

## Current readiness state

**`BLOCKED_AUTHORIZATION`**

Reason: authorization is unresolved because the provider, dataset, applicable
terms, intended-use scope, and evidence reviewed have not been supplied. The
study is also incomplete: instrument, market, frequency, timestamps, date range,
primary objective, target (if forecasting), evaluation boundaries, metrics,
baseline configurations, data-quality policies, and experiment identifier are
not specified.

No provider terms or permission evidence are asserted here. In particular, this
document does not assume that IEX access or any particular IEX use is authorized.

## Fillable specification for one initial study

Replace each `[UNRESOLVED]` only with a documented decision. If a field is not
applicable, enter `NOT APPLICABLE` with a short reason. Do not leave ambiguity
hidden behind a Boolean such as `permission_granted: true`.

### A. Study identity and research question

| Field | Entry |
|---|---|
| Experiment ID | `[UNRESOLVED]` |
| Specification version/date | `[UNRESOLVED]` |
| Responsible researcher/reviewer | `[UNRESOLVED]` |
| Study status | `BLOCKED_AUTHORIZATION` (current) |
| Research question | Does the Phase 30 causal filter provide measurable value for the declared objective, compared with the declared causal baselines, using only information available at each timestamp? `[REFINE BEFORE EVALUATION]` |
| One primary objective | `[UNRESOLVED: select at most one—signal-processing properties OR forecasting; economic performance is a separate later study]` |
| Hypothesis and evidence that would contradict it | `[UNRESOLVED; specify before outcomes are examined]` |
| Exploratory or confirmatory | `[UNRESOLVED]` |

The three claim levels remain separate:

1. **Signal processing:** measurable transformation behavior. Real market data
   generally has no observed noise-free reference signal; do not call a smoother
   its own ground truth.
2. **Forecasting:** a defined target and horizon evaluated out of sample. A
   smoother-looking series does not establish forecasting improvement.
3. **Economic/trading performance:** requires a separate decision/execution
   design, realistic costs, and risk evaluation. It is outside this initial
   readiness gate; no return or profitability claim is permitted here.

If the primary objective is forecasting, complete the target section below. If
it is signal processing, state exactly which observable properties can be
measured without claiming an unavailable true signal. The primary objective is
currently unresolved and must not be selected on the final test results.

### B. Data source and authorization evidence

| Required field | Entry/evidence reference |
|---|---|
| Provider name | `[UNRESOLVED]` |
| Dataset name and identifier/version | `[UNRESOLVED]` |
| Official terms/license URL and version/date, if available | `[UNRESOLVED; a URL alone is not proof of permission]` |
| Exact intended use to be reviewed | `[UNRESOLVED: local research, transformation, retention, sharing/publication, etc.]` |
| Evidence that the intended use is permitted | `[UNRESOLVED: cite the specific terms section, written provider response, license, or other reviewed record]` |
| Permission status | `UNRESOLVED` |
| Local analysis permitted? | `[UNRESOLVED]` |
| Derived results may be retained? | `[UNRESOLVED]` |
| Derived results may be published/shared? | `[UNRESOLVED]` |
| Raw data may be retained, and for how long? | `[UNRESOLVED]` |
| Raw data may be redistributed? | `[UNRESOLVED]` |
| Required attribution | `[UNRESOLVED; state none only when supported by reviewed terms]` |
| Automated-access restrictions | `[UNRESOLVED; no automated access is authorized by this phase]` |
| Expiration, renewal, account, or territory limits | `[UNRESOLVED]` |
| Terms reviewed by/date | `[UNRESOLVED]` |
| Review method and evidence location/reference | `[UNRESOLVED]` |
| Separate provider permission required? | `[UNRESOLVED; if not required, cite the reviewed terms supporting that conclusion]` |

Mark permission `UNRESOLVED` if evidence is missing, vague, outdated, does not
cover the intended use, or conflicts with another applicable restriction. Public
availability is not by itself evidence that local analysis, retention,
publication, or redistribution is authorized. This record is an audit aid, not
legal advice or an independent certification that terms are legally sufficient.

### C. Dataset and observation definition

| Required field | Entry |
|---|---|
| Instrument or explicitly defined universe | `[UNRESOLVED; user/researcher decision required]` |
| Market and trading venue | `[UNRESOLVED]` |
| Universe membership source and point-in-time rule, if applicable | `[UNRESOLVED / NOT APPLICABLE with reason]` |
| Observation frequency | `[UNRESOLVED]` |
| Timestamp meaning (event, bar open/close, publication, etc.) | `[UNRESOLVED]` |
| Time zone and daylight-saving convention | `[UNRESOLVED]` |
| Trading calendar, sessions, holidays, overnight boundaries | `[UNRESOLVED]` |
| Historical date range and rationale | `[UNRESOLVED; dates must be fixed before final-test inspection]` |
| Observed fields, units, adjustment/vintage status | `[UNRESOLVED]` |
| Availability/publication time of each input field | `[UNRESOLVED]` |
| Missing-observation policy | `[UNRESOLVED; must not use future values]` |
| Duplicate timestamp policy | `[UNRESOLVED]` |
| Corporate-action handling, if relevant | `[UNRESOLVED / NOT APPLICABLE with reason]` |
| Futures/contract-roll handling, if relevant | `[UNRESOLVED / NOT APPLICABLE with reason]` |
| Revisions/corrections and point-in-time vintages | `[UNRESOLVED]` |
| Coverage limitations, delistings, survivorship | `[UNRESOLVED]` |
| Planned data-quality exclusions and rule fixed before outcomes | `[UNRESOLVED]` |

Do not select the provider or instrument on the user's behalf. Choices that
depend on a dataset's actual timestamps, revisions, calendar, or missingness can
be finalized only after the dataset is identified and its metadata/terms are
reviewed, but before evaluation outcomes are inspected.

### D. Target and horizon (conditional)

If the selected objective is **forecasting**, all of these fields are required
before readiness:

| Required field | Entry |
|---|---|
| Exact target variable and units | `[UNRESOLVED]` |
| Target construction formula | `[UNRESOLVED]` |
| Forecast origin and information cutoff | `[UNRESOLVED]` |
| Target observation/realization timestamp | `[UNRESOLVED]` |
| Horizon and convention | `[UNRESOLVED]` |
| Availability time for every input | `[UNRESOLVED]` |
| Whether current-bar input is available at prediction time | `[UNRESOLVED]` |
| Overlapping target horizons and purge/embargo handling | `[UNRESOLVED]` |
| Prediction-to-target timestamp alignment rule | `[UNRESOLVED]` |

If the selected objective is **signal processing**, enter `NOT APPLICABLE` for
forecast target and horizon, explain the choice, and define the measurable
properties (for example, attenuation, delay, response to abrupt changes, or
sampling-frequency sensitivity). Do not invent a clean reference for observed
financial data. Reconstruction error requires an independently known reference,
such as a synthetic signal, and is not available merely by smoothing a market
series.

### E. Evaluation plan

| Required field | Entry |
|---|---|
| Training dates/boundaries | `[UNRESOLVED]` |
| Validation dates/boundaries | `[UNRESOLVED]` |
| Final untouched test dates/boundaries | `[UNRESOLVED]` |
| Boundary convention and observation ordering | `[UNRESOLVED]` |
| Embargo/gap and target-overlap rationale | `[UNRESOLVED; enter none only with a stated justification]` |
| Recalibration/retraining schedule | `[UNRESOLVED; fixed beta or explicit schedule]` |
| Primary metric and exact calculation | `[UNRESOLVED]` |
| Secondary metrics and exact calculations | `[UNRESOLVED]` |
| Aggregation across timestamps/windows/instruments | `[UNRESOLVED]` |
| Uncertainty/serial-dependence reporting | `[UNRESOLVED]` |
| No-smoothing baseline | `[PRESPECIFY]` |
| Causal trailing-mean window and edge handling | `[UNRESOLVED]` |
| Persistence baseline, if forecasting | `[UNRESOLVED / NOT APPLICABLE with reason]` |
| Phase 30 beta value(s) and selection procedure | `[UNRESOLVED; fixed in advance or selected on training/validation only]` |
| Warm-up and filter-state initialization/reset rule | `[UNRESOLVED; Phase 30 initializes y[0]=x[0]]` |
| Missing-data exclusions and state behavior | `[UNRESOLVED]` |
| Multiple-testing and reporting-window safeguards | `[UNRESOLVED; record all configurations and failures]` |
| Test-set protection and one-time use rule | `[UNRESOLVED; must be explicit before readiness]` |

Use chronological splits, not random splitting. Training may fit preprocessing;
validation may select among the declared small candidate set; the final test
must not select parameters, metrics, models, exclusions, or reporting windows.
Chronological ordering alone does not prevent leakage through revised data,
global preprocessing, overlapping targets, or repeated selection. Do not claim
filters are matched by bandwidth or delay unless a declared procedure measures
that match without using final-test results.

### F. Data custody and reproducibility

| Required field | Entry |
|---|---|
| Manual file handoff/source identifier | `[UNRESOLVED; no retrieval mechanism is selected]` |
| Retrieval/provision date and method | `[UNRESOLVED]` |
| Integrity digest/checksum and file inventory | `[UNRESOLVED]` |
| Permitted local storage location/retention period | `[UNRESOLVED]` |
| Raw-data deletion/renewal requirement | `[UNRESOLVED]` |
| Experiment ID and configuration record location | `[UNRESOLVED]` |
| Code/version identifier (release or source-tree digest) | `[UNRESOLVED; no Git operation required]` |
| Protocol deviations, exclusions, and failed runs log | `[UNRESOLVED]` |

Do not place restricted raw data in a report or experiment log. Preserve a
reproducible record of configuration, dates, metrics, exclusions, failures, and
protocol deviations within the provider's retention restrictions.

## Readiness states and conservative decision rule

Use exactly these state labels:

| State | Meaning |
|---|---|
| `NOT_STARTED` | No specific study proposal has been submitted for readiness review. |
| `BLOCKED_AUTHORIZATION` | Permission for the documented intended use is missing, ambiguous, expired, insufficiently evidenced, or not reviewed. |
| `BLOCKED_SPECIFICATION` | Authorization evidence has been documented and reviewed for the intended use, but required study-design or pre-ingestion fields remain unresolved. |
| `READY_FOR_DATA_INGESTION` | Authorization evidence and every applicable pre-ingestion decision are documented and reviewed. This permits beginning the separately chosen, authorized ingestion step only. |
| `ON_HOLD` | The study has been deliberately paused and the hold is recorded. |

Apply this precedence: an explicitly recorded hold yields `ON_HOLD`; when a
study proposal is under review, unresolved or ambiguous authorization yields
`BLOCKED_AUTHORIZATION`; before any specific proposal is submitted, use
`NOT_STARTED`; only after authorization is documented and reviewed may
incomplete design yield `BLOCKED_SPECIFICATION`; readiness is allowed only
after every applicable checklist item below passes. If authorization is
unresolved, never report ready, even if every study-design field is filled.

No Boolean field such as `permission_granted: true` is sufficient. Readiness
requires a traceable record of provider and dataset, terms or permission
evidence, exact intended use, evidence references/review date/method, applicable
retention and publication restrictions, and all applicable study decisions.
A URL or a user's assertion without the evidence record is not enough. This
document-only gate does not automatically fetch, interpret, or certify terms.
The state is a conservative project preflight status, not legal advice.

## Pre-ingestion checklist

Do not load real observations until every applicable item is checked and the
reviewed specification is updated. A blank, ambiguous, or unsupported critical
item blocks readiness.

- [ ] Provider, dataset, version, and source identifier are recorded.
- [ ] Intended research use is covered by reviewed, documented terms or
  permission evidence.
- [ ] Any required provider response or separate permission is recorded; if not
  required, the reviewed terms supporting that conclusion are referenced.
- [ ] Raw-data access, storage, retention, redistribution, derived-result
  retention/publication, attribution, restrictions, and expiry are understood.
- [ ] Instrument/universe, market/venue, frequency, timestamp, time zone, and
  session/calendar conventions are fixed.
- [ ] Target, forecast horizon, origin/realization timestamps, and formula are
  defined if forecasting is selected; otherwise the signal-processing objective
  and its available measurements are explicit.
- [ ] Input availability timing, revisions, and point-in-time treatment are
  understood.
- [ ] Training, validation, untouched test boundaries, and any gap/embargo are
  specified chronologically with rationale.
- [ ] Primary/secondary metrics, aggregation, and uncertainty reporting are
  prespecified.
- [ ] Baselines, beta/window candidates, and parameter-selection rules are
  prespecified; the final test is excluded from selection.
- [ ] Warm-up, missing/duplicate observations, market-calendar effects,
  corporate actions/rolls, and exclusions are documented as applicable.
- [ ] A reproducible experiment ID, configuration record, provenance/integrity
  metadata, and failure/deviation log are specified.

Conditional items: target/horizon/persistence requirements apply only to a
forecasting objective; corporate actions and contract rolls apply only where
relevant; separate provider responses apply only if required by the applicable
terms. For a conditional item that does not apply, record `NOT APPLICABLE` and
the reason rather than silently skipping it.

## User decisions, evidence, and dataset-dependent choices

**User/researcher decisions required:** provider and dataset; intended use;
permission evidence; instrument/universe; market; primary objective; target and
horizon if forecasting; observation period; final-test policy; and approval to
proceed once prerequisites are met. None is selected on the user's behalf here.

**Evidence to obtain/review:** the applicable official license/terms and version,
written provider permission if needed, exact research and derived-output uses,
retention/publication/redistribution rights, attribution, access restrictions,
expiration, and a dated record of how the terms were reviewed. Do not imply that
the Phase 31 protocol or this checklist supplies such evidence.

**Choices dependent on identifying a dataset:** timestamp semantics and
availability, timezone/DST, sessions/calendar, historical coverage, revisions,
missingness/duplicates, point-in-time universe coverage, adjustments/rolls,
and source-specific data-quality exclusions. These must be resolved before
evaluation outcomes are viewed and before readiness is granted.

**Established by Phase 31:** chronological evaluation, a protected final test,
causal information timing, the baseline families, train-only preprocessing, and
the separation of signal processing, forecasting, and economic evaluation.
Phase 32 makes these fields auditable for a single study but does not fill
unknown empirical choices.

## Conditions before Phase 33

Phase 33 must not begin real-data ingestion or analysis until:

1. Provider and dataset are identified, and the exact intended use is written.
2. Applicable terms and permission evidence are reviewed and recorded; any
   required provider response is present. Ambiguity means
   `BLOCKED_AUTHORIZATION` and stop.
3. Raw/derived data rights, storage, retention, attribution, restrictions, and
   expiry are documented for the proposed handling.
4. The instrument, market, frequency, timestamp/availability convention,
   observation period, and relevant calendar/corporate-action/revision policies
   are fixed.
5. Exactly one primary objective is selected. Any forecast target, horizon,
   timestamp alignment, and overlap rule are complete if applicable.
6. Chronological train/validation/test boundaries, purge/gap rationale,
   metrics, aggregation, uncertainty, baselines, filter configuration/selection,
   warm-up, missing-data, and multiple-testing policies are prespecified.
7. The experiment ID, configuration, source/integrity metadata, and log for
   exclusions, failures, and deviations are ready.
8. Every applicable pre-ingestion checklist item is checked and the readiness
   state is reviewed as `READY_FOR_DATA_INGESTION`.

This final state is not a claim that data is clean, the design is scientifically
valid, the filter predicts, or any strategy is profitable. The subsequent
ingestion method—including whether a manually supplied file is used—must be
decided only after this gate is satisfied and must remain within the documented
authorization.

## Implementation choice and verification

Phase 32 is document-only. The project has Pydantic configuration and validation
conventions, but a machine-readable authorization validator would not establish
that evidence is authentic, current, or applicable to the intended use. A
Boolean or URL-based validator could create false assurance, while a more
elaborate evidence workflow would exceed this focused phase. The auditable
specification therefore records evidence references and requires conservative
human review; no data-ingestion capability is added.

No real data was accessed, downloaded, scraped, or processed. No legal
authorization, forecast value, or financial performance is claimed.

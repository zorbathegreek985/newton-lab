# Phase 40 — Historical Coverage and Schema Plan

**Prepared:** 2026-10-01

## 1. Executive summary

The project currently has two user-provided NSE-style CM UDiFF Bhavcopy final
files: one report date, 2026-09-28, in the earlier CSV, and one report date,
2026-10-01, in the Phase 39 ZIP and its extracted CSV. This establishes two
observed report dates across separate local artifacts, but it does not yet
establish a consistent, permitted, instrument-level historical series. No larger
archive or continuous historical coverage has been verified.

The Phase 39 ZIP is **consistent with** NSE's named CM-UDiFF Common Bhavcopy
Final report family: the official reports page lists that product; the attached
filename follows the documented UDiFF BhavCopy naming pattern; and its CSV
contains `Sgmt=CM` and `Src=NSE`. Its immediate provenance is nevertheless the
user attachment: no original NSE download receipt or acquisition record is in the
local data tree. This is a format/report-family identification, not independent
proof of custody or authorization.

NSE's official Forms & Formats page links a CM Bhavcopy format workbook and
versioned UDiFF catalogues. The available page reader could not inspect the
linked XLSX field catalogue. The Phase 35 dictionary therefore remains correct:
field-level meanings, units, and field-specific validation are not verified.
In particular, this report does **not** claim that `TradDt` definitively means
trading date, `BizDt` definitively means business date, or that column labels
alone prove identifier, price, or quantity semantics.

For planning, collect at least **252 distinct valid daily observations for one
preselected instrument/series** (approximately one trading year) to run an
initial chronological workflow with rough 50/25/25 train/validation/test
partitions. This is a planning threshold for pipeline exploration, not a
guarantee of statistical adequacy or forecasting power. First verify the actual
field catalogue and applicable use terms. Phase 32 remains
**`BLOCKED_AUTHORIZATION`**.

## 2. Current verified dataset inventory

The local files were re-read on 2026-10-01 using the existing CSV inspector.
The two CSV schemas have 34 headers with the same names and order. This observed
match is not a substitute for verifying the applicable catalogue version.

| Local artifact | Immediate provenance | Observed report date(s) | Shape / integrity | Audit result |
|---|---|---|---|---|
| `data/raw/user_supplied_market_data.csv` | User-provided file from Phase 34; no source receipt/manifest found | `TradDt` and `BizDt`: 2026-09-28 only | 634,897 bytes; SHA-256 `3822761536596adb38171e51d3b1422d99004830c605172e4bb6b648f64858b5`; 3,682 rows × 34 columns | 0 exact duplicate rows; 5 missing `SctySrs`; 12 entirely blank fields; inspector found 0 invalid `TradDt` and 3,681 repeated date values (one date per row). |
| `data/raw/BhavCopy_NSE_CM_0_0_0_20261001_F_0000.csv.zip` | User attachment copied in Phase 39; original NSE retrieval method not documented | ZIP member contains `TradDt` and `BizDt`: 2026-10-01 only | 208,712 bytes; SHA-256 `ccc5fb27872716bbcc99d2d87e522ab304f6620e11c3a5045c30e1ff25cbfb73`; one member, `BhavCopy_NSE_CM_0_0_0_20261001_F_0000.csv`, 640,497 bytes | ZIP CRC valid; member CSV has 3,712 rows × 34 columns; 0 malformed-width rows; 0 exact duplicate rows; 6 missing `SctySrs`; 12 entirely blank fields; inspector found 0 invalid `TradDt` and 3,711 repeated date values. |
| `data/derived/phase_39_audit/BhavCopy_NSE_CM_0_0_0_20261001_F_0000.csv` | Extracted locally from the above ZIP | Same member content: 2026-10-01 only | 640,497 bytes; SHA-256 `5d5ac5591c32ca5d53539df02cfdee3795a2657183570a42eb1d714ab1822354` | Kept separate from raw ZIP. |

“Repeated date values” counts are not duplicate-row counts: all rows in each
cross-section share the same date. The exact duplicate-row count is zero in both
CSV files. Dates above are direct observations in CSV contents, not dates inferred
only from filenames.

The Phase 39 CSV has 3,712 populated `FinInstrmId` values and 3,712 distinct
values; 3,712 `ISIN` values and 3,709 distinct values; and 3,712 `TckrSymb`
values and 3,709 distinct values. These are within-file counts, not proof of
instrument semantics or a stable cross-date key. The Phase 35/36/39 reports
document the exact observed schema and additional missingness details.

No standalone Phase 34 or Phase 38 research document was found in
`docs/research/`; this inventory does not infer their contents from their phase
names.

## 3. Authoritative schema evidence

Official NSE pages were consulted on **2026-10-01**:

| Source | Verified support | Boundary of the evidence |
|---|---|---|
| [NSE All Reports](https://www.nseindia.com/all-reports), accessed 2026-10-01 | Lists **CM-UDiFF Common Bhavcopy Final (zip)** and says the older CM Bhavcopy CSV and CM Common Bhavcopy CSV reports were discontinued effective 2024-07-08. | Names the report family; does not establish the custody, historical completeness, access rights, or fields of the attached artifact. |
| [NSE Forms & Formats](https://www.nseindia.com/static/resources/forms-formats-members), page updated 2026-06-30, accessed 2026-10-01 | Describes UDiFF and lists a **Bhavcopy File** format for segment CM, with linked format and test files. It lists UDiFF catalog versions 1.0–4.0. | The field-format workbook is XLSX and was not readable through the available page reader. No specific catalogue version was matched to the local file. |
| [SEBI UDiFF Guidance Note V2023.10](https://nsearchives.nseindia.com/web/sites/default/files/inline-files/UDiFF%20guidance%20document_Ver1.0.pdf), accessed 2026-10-01 | Gives generic UDiFF rules: header row followed by data; comma-separated CSV; ISO 8601 date formatting; blank values can apply to non-applicable fields; generic filename components include format, entity, segment/market/member placeholders, date, provisional/final marker, and time marker. Its example identifies `F_0000` as a final file pattern. | Generic format guidance does not define the specific columns `TradDt`, `BizDt`, `FinInstrmId`, price-like names, or quantity-like names. The note directs readers to the UDiFF catalogue for field details. |

Together, the official report listing, filename pattern, and local values support
the cautious description **“user-supplied file consistent with NSE CM UDiFF
Common Bhavcopy Final.”** The official filename example is a format convention;
it is not provenance for this supplied file.

## 4. Unresolved field definitions and units

The Phase 35 dictionary remains the controlling field-level record. It states
that the linked official catalogue workbook could not be inspected and does not
assert individual field meanings or units. The same header names occur in the
two local CSVs, but repeated headers do not verify semantics.

| Observed header(s) | Cautious status | What must be verified in the applicable official CM Bhavcopy catalogue |
|---|---|---|
| `TradDt`, `BizDt` | Each is a date-like CSV field; their meanings are unverified. Both happen to equal the observed report date in these two files. | Whether each means exchange trading date, business/processing date, settlement date, or another date; rules for legitimate differences and applicable calendar. |
| `FinInstrmId`, `ISIN`, `TckrSymb` | Candidate instrument/security identifiers based on their labels and populated values; field semantics, uniqueness scope, and lifecycle are unverified. `ISIN` and `TckrSymb` repeat within the 2026-10-01 file. | Authoritative definition, scope, stability, permitted nulls, mapping over identifier changes, and whether date plus `FinInstrmId` is a valid key. |
| `OpnPric`, `HghPric`, `LwPric`, `ClsPric`, `LastPric`, `PrvsClsgPric`, `SttlmPric` | Decimal-like populated columns with price-like labels; not treated as verified prices. | Exact calculation/event, currency/unit, adjustments, timestamp/availability, validation range, and applicability by instrument type. |
| `TtlTradgVol`, `TtlTrfVal`, `TtlNbOfTxsExctd`, `NewBrdLotQty` | Numeric-like columns with quantity/value/count-like labels; not assigned units or economic meanings. Exact local field `TtlTrfVal` is present; `TtlTradgVal` is not. | Unit and scale, aggregation interval/scope, zero/blank rules, valid range, and exact definition of each field. |
| `Sgmt`, `Src`, `FinInstrmTp`, `SctySrs`, `SsnId` | Observed codes/text; only the file's values are direct evidence. | Catalogue standard values and intended use, including applicability and cross-date consistency. |
| Other blank fields | Empty in the observed counts listed in the inventory; meaning/applicability is not inferred. | Whether each field is optional for this report/instrument class, and the expected blank/null representation. |

Before constructing a series, use the official catalogue linked from the NSE
Forms & Formats page to match the relevant report, segment, exact headers, field
definitions, standard values, units, and catalogue effective version. If that
catalogue is still inaccessible, record the unanswered items and do not label
the candidate fields as verified targets.

## 5. Data-use and authorization status

The [NSE Data Sharing & Usage Policy](https://www.nseindia.com/static/market-data/nse-data-policy)
(page says updated 2025-12-11; accessed 2026-10-01) states that subscriber
access entails a request and a relevant agreement describing intended use and
handling. It says redistribution is not allowed except as agreed in that
agreement. Its research provisions route researcher/institution requests
through NSE's Economic Policy Research department; offered data has use limits
in underlying documentation, and research recipients sign confidentiality
declarations. It also allows that fee reductions/waivers may be considered for
non-commercial users; this is not a promise of free access.

The [NSE EOD/Historical Data page](https://www.nseindia.com/static/market-data/eod-historical-data-subscription)
(accessed 2026-10-01) describes subscribed EOD data and historical-data request
routes. The [NSE Research Initiatives page](https://www.nseindia.com/static/research/research-initiatives)
(accessed 2026-10-01) states research data under that initiative is limited to
academic and non-commercial use. These pages do not establish that the attached
public-report file is licensed to this project or that this project qualifies
under the research initiative.

The attachment's immediate source is known only as user provision. No terms,
license, written permission, retrieval receipt, retention rule, or publication
permission is recorded locally. Therefore **Phase 32 stays
`BLOCKED_AUTHORIZATION`**. No file acquisition from an official source and no
permission check with NSE occurred during this phase.

## 6. Historical coverage assessment

Two distinct dates are present across two separate files: 2026-09-28 and
2026-10-01. This is **not verified continuous coverage**. The files do not by
themselves show availability for intervening market sessions, whether the
reports use the same specification revision, or whether an individual target
instrument has a valid row on every date. No further historical dates are
claimed.

A folder of daily cross-sections is not yet a time series for a particular
instrument. The rows must first be matched using an officially verified,
stable identifier, one defined observation per intended date, and consistent
field semantics. Symbol text alone can change or repeat. Duplicate instrument
rows require an explicit documented rule; they must not be silently collapsed.
Cross-date missingness may indicate non-trading, instrument lifecycle, report
coverage, or data problems, and these explanations must not be conflated.

Trading holidays are not missing trading-day observations: compare dates to an
authoritative exchange calendar and preserve the distinction. Missing expected
report dates, missing instrument rows, date-field disagreement, corrections,
schema changes, and changes in identifier or instrument availability should be
logged and investigated. Do not forward-fill market values or delete records
without a prespecified, justified policy. If provisional and final reports both
exist, select the appropriate version under a documented rule while preserving
each source file unchanged.

One day's cross-section has no within-instrument temporal ordering beyond a
single observation. It cannot establish a historical forecasting baseline,
temporal validation, or out-of-sample performance.

## 7. Proposed minimum dataset specification

For an initial **daily exploratory time-series pipeline exercise**:

- One preselected instrument or explicitly defined series; do not choose it
  after looking at comparative outcomes.
- At least **252 distinct valid observation dates** for that one series, subject
  to available and authorized coverage. This roughly corresponds to one year of
  daily sessions as a planning approximation, not an exact calendar assertion.
- If the intended observation is daily, one documented observation per valid
  trading session and a documented rule for non-trading holidays and exceptional
  closures. Verify whether `TradDt` is the appropriate date field before use.
- One officially verified target field, unit, adjustment convention, and
  information-availability time, with the same meaning across the observation
  period. Establish whether data are corrected/revised and which vintage is
  used.
- A stable verified instrument key; a documented resolution for symbol/ID
  changes, duplicate keys, corporate actions, suspensions, delistings, and
  absent rows, as relevant to the selected instrument.
- A frozen chronological split plan, illustratively 50% train, 25% validation,
  25% final test (about 126/63/63 dates at the 252-date floor). These counts are
  a planning example, not proof that any segment is statistically adequate.
- Source records, checksums, schema fingerprint, applicable terms, intended use,
  retention and output restrictions recorded before analysis.

The target and horizon affect how much data is needed, how boundaries should be
placed, and whether gap/embargo periods are required. A 252-date planning floor
does not establish forecast skill or adequate statistical power. A later
prespecified design may require substantially more history.

## 8. Manual ingestion and validation procedure

No automatic downloader or new utility is introduced. For each file obtained
manually through an authorized route:

1. **Before retrieval**, confirm the selected product, exact date, format/version,
   and applicable written terms for the intended use. If any required right is
   unclear, stop and leave Phase 32 blocked.
2. **Preserve raw source:** copy the original download unchanged into
   `data/raw/`. Use the original filename when unused; on collision, choose a
   distinct suffix such as `_copy02` and never overwrite. Record original
   filename, bytes, SHA-256, retrieval date, claimed report date, exact source
   page or manual retrieval method, and a user/source receipt if available.
3. **Validate a ZIP before extracting.** With the project's existing Python:

   ```powershell
   .venv\Scripts\python.exe -c "from zipfile import ZipFile; z=ZipFile(r'data\raw\FILE.zip'); print(z.namelist()); print('bad member:', z.testzip())"
   ```

   Continue only if the archive opens, the expected member list is understood,
   and `testzip()` returns `None`. Record the output. Never extract over an
   existing directory or file.
4. **Extract separately** under a non-raw directory such as
   `data/derived/nse_cm_udiff/<claimed-report-date>/<archive-basename>/` after
   checking the destination does not exist. Keep the ZIP as the raw record.
5. **Inspect without repair** using
   `newton_lab.data_inspection` and `newton_lab.market_data_exploration` on the
   extracted CSV. Record exact headers, row/column counts, per-column missing
   counts, candidate date coverage, duplicate full rows, and within-file
   candidate identifier counts. The repeated-date count is not a duplicate-row
   count. The utilities' lexical types and candidate identifiers do not prove
   field meaning.
6. **Validate structure and dates:** compare every record width with header
   width; reject/record malformed rows rather than skipping them. Compare observed
   `TradDt`/`BizDt` values with one another and with the claimed date, preserving
   any disagreement for investigation. Confirm date parsing and expected
   report-date coverage. Do not correct cells in the raw file.
7. **Fingerprint schema:** compute SHA-256 over a canonical serialization of the
   ordered header list (for example, UTF-8 JSON with compact separators). Record
   the fingerprint and exact schema version for every file. Compare fingerprints
   across dates. A change blocks combination until the official specification
   and field mapping establish compatibility; never rename or drop columns
   silently.
8. **Keep representations separate:** `data/raw/` contains unchanged received
   files. Extraction and any later normalized/matched series belong in separate
   derived directories, with a mapping/provenance record. Do not commit raw or
   derived observations to source code or documentation.

For date consistency checks, a mismatch is a review condition, not an automatic
license to rewrite the date. For exact-row duplicate checks, report both the
duplicate count and the repeated timestamp/date count separately. Retain source
records and log all exclusions or transformations in a derived-data manifest.

## 9. Acceptance criteria for beginning a time-series experiment

The **dataset readiness** review may proceed only when all are true:

- Authorized product and intended use have written, reviewed evidence covering
  local access, storage/retention, derived results, attribution, publication and
  redistribution as applicable.
- The official CM UDiFF catalogue/version and exact field definitions, units,
  code lists, date meanings, identifier rules, and validations are recorded.
- At least 252 distinct valid dates for the same chosen series are present, or a
  different sample floor is justified and prespecified for the defined objective.
- Date coverage, expected holidays, missing observations, duplicate keys,
  identifier history, corrections, corporate actions, and schema changes have
  been assessed without silently altering raw files.
- The exact target, units, horizon, forecast origin/information cutoff, and
  chronological train/validation/test boundaries are documented before outcome
  inspection; the final test is protected from selection.
- A reproducible manifest records provenance, sizes, hashes, source versions,
  configurations, parsing results, schema fingerprints, deviations, and failures.

Even if these dataset conditions pass, starting an experiment additionally
requires the separate Phase 32 gate to satisfy every applicable pre-ingestion
item and be reviewed as `READY_FOR_DATA_INGESTION`. A dataset passing validation
does not change authorization status.

## 10. Risks, limitations, and next action

Key risks are unproven attachment provenance and usage rights; a field catalogue
that is linked but not yet inspected; a potentially changing data schema; short
or gapped coverage; unstable identifiers; and ambiguity in how dates, values,
adjustments, and missing rows should be interpreted. Two files with matching
headers do not resolve those risks.

**Next action:** manually verify the exact CM Bhavcopy field-format workbook and
its effective version against the attached report, then decide one target series
and date range. Separately review and document written permission/terms for that
intended use. Only after those decisions should the user manually obtain
additional authorized daily files, beginning with enough dates to evaluate the
proposed 252-date planning floor. No request has been sent and no new data has
been acquired in Phase 40.

### Related project records

- [Phase 32 — Study readiness and authorization gate](phase_32_study_readiness.md)
- [Phase 33 — Official data source audit](phase_33_data_source_audit.md)
- [Phase 34 — local CSV inspection implementation and tests](../../src/newton_lab/data_inspection.py) and [tests](../../tests/test_data_inspection.py); no separate Phase 34 research report found.
- [Phase 35 — Market data dictionary](phase_35_market_data_dictionary.md)
- [Phase 36 — Cross-sectional analysis](phase_36_cross_sectional_analysis.md)
- [Phase 37 — Historical acquisition plan](phase_37_historical_data_acquisition_plan.md)
- Phase 38 — no report found in `docs/research/`.
- [Phase 39 — Historical dataset audit](phase_39_historical_dataset_audit.md)

### Verification performed

The existing CSV inspector was run on the Phase 39 extracted CSV and the earlier
2026-09-28 CSV. It reported the row counts, missingness, zero exact duplicate
rows, date coverage, and repeated-date diagnostics summarized above. The Phase
39 audit already records ZIP integrity, SHA-256 and CSV row-width checks. This
phase changed documentation only; no files were added to or extracted from the
data directories, no existing report or raw file was modified, and no tests were
run. Official NSE pages listed above were checked on 2026-10-01; no file or data
was downloaded from them.

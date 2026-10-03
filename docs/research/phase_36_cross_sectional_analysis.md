# Phase 36 — Single-Day Market Data Exploration

## Research question and scope

What does the supplied single-report CSV contain, and what numerical patterns and
data-quality checks can be described without assuming unverified field meanings?
This is a descriptive, cross-sectional exploration. It does not forecast, construct
signals, recommend investments, or test profitability.

## Dataset identity and limitations

- File: `data/raw/user_supplied_market_data.csv` (634,897 bytes)
- Shape: 3,682 rows × 34 columns
- SHA-256: `3822761536596adb38171e51d3b1422d99004830c605172e4bb6b648f64858b5`
- Both `TradDt` and `BizDt` have one observed value: `2026-09-28`.
- The Phase 35 field dictionary records that field definitions and units remain
  unverified. No economic meaning or unit is assigned here. The filename and
  abbreviated headers alone do not establish an official report format.
- Provenance, licence, terms, and permitted research use remain unresolved. The
  file's presence does not clear Phase 32 authorization; readiness remains
  `BLOCKED_AUTHORIZATION`.

The file contains one report date, not a temporal series. Its rows may describe
different records on that date, but this analysis does not establish that each row
is one company or one unique security.

## Methods and assumptions

The local Python utility reads CSV text without changing it, dropping rows, filling
values, or deduplicating by either date field. It classifies values lexically as
date-like, integer-like, decimal-like, categorical text, mixed, or all missing.
These are observed formats, not semantic types. Blank strings and the tokens `NA`,
`N/A`, `NaN`, and `null` are treated as missing. A malformed row whose field count
does not match its header is rejected.

The utility counts exact duplicate rows; reports cardinality and missingness for
every column; calculates descriptive statistics for numeric-like columns; and
checks the within-file uniqueness of `FinInstrmId`, `ISIN`, and `TckrSymb` as
candidate identifiers. Pearson correlations use row-aligned complete pairs among
numeric-like columns, excluding those three candidate identifiers. Correlations
are numerical associations and have no causal interpretation.

For decimal-like columns except the candidate identifiers, the utility reports
values outside Tukey's 1.5-IQR fences as screening flags. These are not errors or
confirmed outliers. The thresholds are generic and no field-specific validation
rules or units have been verified. No numeric columns were treated as prices,
returns, volumes, or monetary amounts for analysis.

## Direct observations

### Structure and missingness

- Exact duplicate rows: **0**.
- `SctySrs`: **5** missing values; 54 distinct observed values among 3,677
  populated rows.
- The following 12 columns are entirely missing (3,682 cells each):
  `XpryDt`, `FininstrmActlXpryDt`, `StrkPric`, `OptnTp`, `UndrlygPric`,
  `OpnIntrst`, `ChngInOpnIntrst`, `Rmks`, `Rsvd1`, `Rsvd2`, `Rsvd3`, `Rsvd4`.
- Date coverage is one date in each of `TradDt` and `BizDt`; all 3,682 values in
  each column parsed as dates. A frequency estimate is not meaningful for a single
  distinct date.
- Categorical-text cardinalities: `Sgmt` 1, `Src` 1, `FinInstrmTp` 1,
  `ISIN` 3,681, `TckrSymb` 3,680, `SctySrs` 54, `FinInstrmNm` 3,671, and `SsnId` 1.
  These counts do not verify what any code or label means.
- Candidate identifier observations: `FinInstrmId` is populated and unique in all
  3,682 rows; `ISIN` has 3,681 distinct values (one repeated occurrence);
  `TckrSymb` has 3,680 distinct values (two repeated occurrences). This tests
  strings within this file only and does not establish identifier semantics or
  record identity.

### Numeric-like columns

The table gives median, observed range, sign counts, and (for decimal-like fields)
the number outside the generic IQR fences. Units and semantic meanings are
unverified. Values are retained at their observed row level, including repeated
identifiers.

| Exact column | Lexical type | Median | Minimum–maximum | Zero / negative / positive | IQR flags |
|---|---:|---:|---:|---:|---:|
| `FinInstrmId` | integer-like | 17,397.5 | 1–766,238 | 0 / 0 / 3,682 | — |
| `OpnPric` | decimal-like | 172 | 0.15–124,005 | 0 / 0 / 3,682 | 391 |
| `HghPric` | decimal-like | 174.355 | 0.15–124,035 | 0 / 0 / 3,682 | 389 |
| `LwPric` | decimal-like | 168.15 | 0.14–123,195 | 0 / 0 / 3,682 | 389 |
| `ClsPric` | decimal-like | 169.925 | 0.14–123,425 | 0 / 0 / 3,682 | 386 |
| `LastPric` | decimal-like | 170 | 0.15–123,420 | 0 / 0 / 3,682 | 386 |
| `PrvsClsgPric` | decimal-like | 171.54 | 0.14–124,035 | 0 / 0 / 3,682 | 390 |
| `SttlmPric` | decimal-like | 169.925 | 0.14–123,426.68 | 0 / 0 / 3,682 | 386 |
| `TtlTradgVol` | integer-like | 60,060.5 | 1–500,546,584 | 0 / 0 / 3,682 | — |
| `TtlTrfVal` | decimal-like | 9,390,085.375 | 49.14–20,479,635,124.5 | 0 / 0 / 3,682 | 590 |
| `TtlNbOfTxsExctd` | integer-like | 933.5 | 1–1,559,174 | 0 / 0 / 3,682 | — |
| `NewBrdLotQty` | integer-like | 1 | 1–999,999,999 | 0 / 0 / 3,682 | — |

Other summary statistics (quartiles, mean, population standard deviation, and
distinct counts) are preserved in the machine-readable summary. For example,
`TtlTrfVal` has mean 273,045,095.638 and 3,682 distinct observed values. These
descriptions do not establish that the column is a currency value or specify its
scale.

Among the computed complete-pair Pearson associations, the largest absolute
coefficient is `ClsPric` with `SttlmPric` (r ≈ 0.99999984, 3,682 row pairs).
Several other columns have similarly high numerical associations. These are
single-date, unadjusted cross-sectional values; correlated fields may be related by
their source construction. Without verified definitions, units, and data rules,
these results should not be interpreted as a financial relationship or discovery.

## Figures

Figures and the machine-readable profile are in
[`reports/phase_36_single_day_final/`](../../reports/phase_36_single_day_final/):

- `figures/missingness_by_column.png` — missing-cell counts by exact header.
- `figures/decimal_column_distributions.png` — observed distributions of the eight
  decimal-like columns; axes use exact names and state that units are unknown.
- `figures/sctysrs_code_counts.png` — frequent raw `SctySrs` values, explicitly
  labelled as unverified codes.
- `figures/numeric_column_correlations.png` — Pearson r among numeric-like columns
  (excluding the three candidate identifier fields); this is descriptive only.
- `analysis_summary.json` — all column profiles, numeric summaries, uniqueness
  checks, date coverage, pairwise statistics, warnings, environment versions,
  command, and source hash.

The histograms have long right tails and some extreme observed values compress the
bulk of each distribution near the left edge. They display raw counts and values;
no log transform or fitted distribution was applied. IQR flags identify candidates
for follow-up against official field rules, not erroneous observations.

## Interpretation and open questions

**Direct observations:** the shape, exact text and numeric values, missingness,
duplicate-row count, within-file distinct counts, date coverage, and calculated
statistics above.

**Interpretations requiring caution:** the file appears cross-sectional because
both date columns have one distinct date. The name `FinInstrmId` suggests an
identifier, but neither its official definition nor its scope was verified in this
phase. No other field meaning, code list, unit, valid range, or calculation rule is
asserted.

The Phase 35 data dictionary and official-source audit
([`phase_35_market_data_dictionary.md`](phase_35_market_data_dictionary.md),
[`phase_33_data_source_audit.md`](phase_33_data_source_audit.md)) record unresolved
specification and terms. Before field-level or research interpretation, obtain and
verify an authoritative specification covering at least:

- The report identity and exact schema/version corresponding to this file.
- Definitions, data types, units/scales, code lists, and null conventions for all
  fields, particularly the price-like labels, `TtlTrfVal`, count-like columns,
  `NewBrdLotQty`, and date fields.
- The intended record key, treatment of multiple rows sharing `ISIN` or
  `TckrSymb`, validation/range rules, and the meaning of the all-missing columns.
- Provenance, retrieval details, licence, terms, and permitted use for this file.

One report date cannot provide a sequence of observations across time. It therefore
cannot support temporal forecasting, trend or volatility estimation, or historical
validation by itself. Those questions require a properly authorized, consistently
defined multi-date dataset and a separate research protocol. Phase 32 remains
`BLOCKED_AUTHORIZATION`.

## Reproduction

From the workspace root, with the existing project environment:

```powershell
.venv\Scripts\python.exe -m newton_lab.market_data_exploration data\raw\user_supplied_market_data.csv --output-dir reports\phase_36_single_day_final
```

The recorded environment is Python 3.12.10 and Matplotlib 3.11.2. The profile
module itself uses the Python standard library; Matplotlib is used for PNG figures.
The analysis command reported that it could not write the user's Matplotlib font
cache due to a permission error, but it completed and saved all four figures and
the JSON summary. No data was downloaded or uploaded, and no dependencies or raw
data were changed.

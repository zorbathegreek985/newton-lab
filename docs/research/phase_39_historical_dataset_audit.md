# Phase 39 — Historical Dataset Audit

**Audit date:** 2026-10-01

## Result

**`NEEDS_MORE_DATA`** — the attached archive is readable and its CSV is
structurally valid, but its contents cover only one report date. It is not a
multi-date time series and is insufficient for historical forecasting.

Phase 32 remains **`BLOCKED_AUTHORIZATION`**. The user's provision of a file
does not establish the provider, applicable terms, licence, or permission for
the intended research use.

## Attachment copy and integrity

The attached source was available at
`E:\BhavCopy_NSE_CM_0_0_0_20261001_F_0000.csv.zip`. It was copied, not moved, to:

`data/raw/BhavCopy_NSE_CM_0_0_0_20261001_F_0000.csv.zip`

The destination filename did not already exist. Source and destination sizes are
both **208,712 bytes**, and their SHA-256 digests match:

`ccc5fb27872716bbcc99d2d87e522ab304f6620e11c3a5045c30e1ff25cbfb73`

The ZIP contains one member, `BhavCopy_NSE_CM_0_0_0_20261001_F_0000.csv`, with
an uncompressed size of **640,497 bytes**. ZIP CRC validation completed without
reporting a bad member. The CSV was extracted separately to
`data/derived/phase_39_audit/`; no existing CSV or Phase 36 artifact was
overwritten.

Extracted CSV SHA-256:

`5d5ac5591c32ca5d53539df02cfdee3795a2657183570a42eb1d714ab1822354`

## Dataset contents

The report date is based on the extracted CSV values, not the download date or
filename: every row has `TradDt = 2026-10-01` and `BizDt = 2026-10-01`. The CSV
has **3,712 data rows and 34 columns**. Its exact headers, in file order, are:

```text
TradDt, BizDt, Sgmt, Src, FinInstrmTp, FinInstrmId, ISIN, TckrSymb, SctySrs,
XpryDt, FininstrmActlXpryDt, StrkPric, OptnTp, FinInstrmNm, OpnPric, HghPric,
LwPric, ClsPric, LastPric, PrvsClsgPric, UndrlygPric, SttlmPric, OpnIntrst,
ChngInOpnIntrst, TtlTradgVol, TtlTrfVal, TtlNbOfTxsExctd, SsnId, NewBrdLotQty,
Rmks, Rsvd1, Rsvd2, Rsvd3, Rsvd4
```

Observed date and candidate identifier diagnostics:

- `TradDt`: 3,712 populated rows; one distinct date, `2026-10-01`.
- `BizDt`: 3,712 populated rows; one distinct date, `2026-10-01`.
- `FinInstrmId`: 3,712 populated values and 3,712 distinct values.
- `ISIN`: 3,712 populated values and 3,709 distinct values.
- `TckrSymb`: 3,712 populated values and 3,709 distinct values.

These are observed strings and cardinalities only. The candidate identifier
columns have not been treated as proof that each row represents one unique
company or security. Field definitions and units remain subject to the
unresolved Phase 35 specification questions.

## Structural and data-quality observations

- The existing `newton_lab.data_inspection` utility read the CSV successfully:
  3,712 rows, 34 headers, no invalid `TradDt` values, zero exact duplicate rows,
  and one distinct date.
- Its diagnostic count of duplicate `TradDt` values is 3,711 because every row
  shares the same date. This is not the exact duplicate-row count and does not
  show that the rows themselves are duplicate records.
- An independent CSV-width check found **zero malformed rows**: every parsed data
  row has 34 fields.
- Missing cells: `SctySrs` has **6**. These 12 columns are entirely empty
  (3,712 each): `XpryDt`, `FininstrmActlXpryDt`, `StrkPric`, `OptnTp`,
  `UndrlygPric`, `OpnIntrst`, `ChngInOpnIntrst`, `Rmks`, `Rsvd1`, `Rsvd2`,
  `Rsvd3`, and `Rsvd4`.
- No rows were dropped, filled, or deduplicated. Numeric-looking field values
  were not assigned meanings or units in this audit.

## Cross-section versus time series

This is a **single-day cross-section**: 3,712 rows share one `TradDt` and one
`BizDt`. The repeated date in a row-oriented report does not create temporal
observations. It cannot support historical forecasting, temporal trend or
volatility estimation, or chronological train/validation/test partitions.

A future forecasting study requires multiple distinct dates for a clearly
identified instrument or explicitly defined series; an authoritative and
version-matched schema; documented target-field meaning and units; consistent
identifier and adjustment conventions; sufficient observations for prespecified
chronological splits; and recorded source provenance, integrity, terms, and
permitted use. Do not concatenate this CSV with later files until the report
identity, specification version, identifiers, field definitions, and consistency
across dates are checked.

## Provenance and unresolved questions

The immediate source provenance recorded here is **user attachment**, supplied
on 2026-10-01 and copied locally. The archive filename and content date are
consistent, but they do not independently establish the original exchange
download route or custody. No receipt, acquisition record, licence, terms review,
or permission evidence accompanied it in the inspected `data` tree. The Phase 37
acquisition plan's NSE manual-history and formal subscription routes remain
options to investigate, not evidence that this particular file or intended use is
authorized.

Still unresolved:

- Original source URL/account and retrieval method, beyond the user-provided
  attachment.
- Applicable NSE data terms, licence, retention, derived-result, publication,
  redistribution, attribution, and intended-use permissions.
- Official confirmation that this exact archive and CSV version are the named
  CM-UDiFF report, and the applicable field definitions, units, null rules, code
  lists, valid ranges, and identifier semantics.
- Historical coverage and consistency for any additional dates.

No network access, downloads, model training, or external transfer occurred. No
dependencies or project code were changed.

## Reproduction and checks

Using the existing project environment, the inspection command was:

```powershell
.venv\Scripts\python.exe -m newton_lab.data_inspection data\derived\phase_39_audit\BhavCopy_NSE_CM_0_0_0_20261001_F_0000.csv --timestamp-column TradDt
```

The CSV archive and content audit used Python's standard-library `zipfile`, `csv`,
and `hashlib` facilities. ZIP CRC, member list, file sizes, SHA-256 hashes,
header/row widths, date values, missingness, exact duplicate rows, and observed
identifier cardinalities were checked. No automated test suite was run: no source
code or tests were changed.

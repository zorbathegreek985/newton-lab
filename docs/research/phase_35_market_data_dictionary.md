# Phase 35: Market Data Field Verification

**Reviewed:** 2026-09-29
**Local file:** `data/raw/user_supplied_market_data.csv`
**Scope:** verify report family, exact schema, observed value types, and which
interpretations are supported by accessible official documentation. No market
data were downloaded, changed, or uploaded for this phase.

## Identification and evidence limits

The CSV has 34 headers and 3,682 data rows. Its filename
`BhavCopy_NSE_CM_0_0_0_20260928_F_0000.csv`, the observed `Sgmt=CM` and
`Src=NSE` values, and its schema are consistent with an NSE Capital Market
UDiFF Bhavcopy file. NSE's [All Reports page][nse-reports] lists the
â€œCM-UDiFF Common Bhavcopy Final (zip)â€ report family. NSE's [Forms & Formats
page][nse-formats] identifies a UDiFF Bhavcopy format for the CM segment and
links a format workbook.

This supports the **format-family identification**, not proof of the attached
file's chain of custody. The user supplied the file; no retrieval record or
provider authorization evidence is attached. The Phase 32 state therefore
remains **`BLOCKED_AUTHORIZATION`**.

The official format workbook linked from the NSE Forms & Formats page is an
Excel file. The available page reader returned an unsupported-content-type
error when that workbook was opened. It was not downloaded or imported. Thus,
the workbook's field master, standard value lists, field lengths, and
field-specific units/validation rules could not be inspected. The official
[UDiFF guidance note][udiff-guidance] explains that the catalogue's Field
Master organizes fields by format, but is not itself the field-by-field
catalogue. **No individual field meaning or unit below is asserted unless the
accessible official documentation explicitly supports it.**

The guidance note does verify format-wide rules: a UDiFF CSV has one table,
headers in row one, ISO tags, comma-separated fields, ISO 8601 date/time
formatting, and values may be absent where fields are not applicable. It also
describes final-file naming conventions. These generic rules do not supply the
meaning of each header, the allowed values for each field, or a per-column
validation rule. The observed date strings conform to `YYYY-MM-DD`; that is an
observation, not a claim that both date fields have a verified business
meaning.

## Observed file profile

- Size: 634,897 bytes.
- Rows: 3,682; all rows have 34 fields, matching the 34 headers.
- Duplicate whole rows: 0.
- `TradDt`: 3,682 populated values, all `2026-09-28`.
- `BizDt`: 3,682 populated values, all `2026-09-28`.
- Candidate date columns: `TradDt`, `BizDt`; neither contains a time of day.
- Inferred timestamp frequency: unavailable. There is only one distinct date.
- The inspector reports 3,681 repeated `TradDt` values because the report date
  repeats across the 3,682 rows. This is not a duplicate-row count.
- `Sgmt` is `CM` on every row; `Src` is `NSE` on every row; `FinInstrmTp` is
  `STK` on every row; `SsnId` is `F1` on every row. Except for the documented
  `CM` segment code, these values remain uninterpreted codes here.
- `SctySrs` has 5 blank cells; its other 3,677 values are text codes.
- Twelve fields are blank on all rows: `XpryDt`, `FininstrmActlXpryDt`,
  `StrkPric`, `OptnTp`, `UndrlygPric`, `OpnIntrst`, `ChngInOpnIntrst`, `Rmks`,
  `Rsvd1`, `Rsvd2`, `Rsvd3`, and `Rsvd4`.

Types below are **value-based lexical inferences**, not official schema types.
CSV cells are text on disk. Integer-like identifiers must remain strings until
their official type and semantics are confirmed. â€œDecimal-likeâ€ does not
establish currency or a unit.

## Field dictionary

Citation keys refer to the official sources listed below. For all fields,
field-specific validations and standard-value rules remain unresolved because
the linked catalogue workbook could not be read. Generic UDiFF structure/date
rules are not substitutes for that field master.

| Exact CSV column | Verified meaning | Observed/inferred data type | Units | Values in this file | Officially supported validation / source | Unresolved questions |
|---|---|---|---|---|---|---|
| `TradDt` | Not verified at field level. | Date-like text (`YYYY-MM-DD`) | Not verified | Populated: 3,682; one distinct value, `2026-09-28` | Observed date format is consistent with generic ISO date rule [C]. | Does this tag mean trade date for this Bhavcopy format? What timezone/calendar applies? |
| `BizDt` | Not verified at field level. | Date-like text (`YYYY-MM-DD`) | Not verified | Populated: 3,682; one distinct value, `2026-09-28` | Observed date format is consistent with generic ISO date rule [C]. | How does business date differ from `TradDt` for this report? |
| `Sgmt` | The row value `CM` is identified as Capital Market by NSE's format page. [B] | Text/code | Not applicable | Populated: 3,682; all `CM` | `CM` is listed as the Bhavcopy segment on the official formats page [B]. Per-row allowed values unavailable. | Is the segment value set fixed, and are there exclusions within CM? |
| `Src` | Not verified. Observed value is `NSE`. | Text/code | Not applicable | Populated: 3,682; all `NSE` | No field-specific validation verified [A][B][C]. | Does this identify exchange, source system, or another origin field? |
| `FinInstrmTp` | Not verified. Observed value is `STK`. | Text/code | Not applicable | Populated: 3,682; all `STK` | No field-specific standard-value list verified [A][B][C]. | What does `STK` denote in this exact file format? |
| `FinInstrmId` | Not verified. | Integer-like text | Not applicable | Populated: 3,682 | No field-specific validation verified [A][B][C]. | Is this identifier unique within date/segment? Is it stable across dates? |
| `ISIN` | Not verified at field level. | Text | Not applicable | Populated: 3,682 | No field-specific length/check-digit rule verified [A][B][C]. | Are all rows expected to carry ISINs? How are non-ISIN instruments represented? |
| `TckrSymb` | Not verified at field level. | Text | Not applicable | Populated: 3,682 | No field-specific validation verified [A][B][C]. | Is this a trading symbol, and what are its symbol-change rules? |
| `SctySrs` | Not verified. | Text/code | Not applicable | Mixed: 3,677 populated; 5 blank | No allowed-code or blankability rule verified [A][B][C]. | What do series values mean? Are blanks valid for the affected rows? |
| `XpryDt` | Not verified. | Undetermined (all blank) | Not verified | Blank: 3,682 | Generic UDiFF guidance permits non-applicable fields to be blank [C]; applicability here is not established. | Is expiry date inapplicable for all rows in this CM file, or is a data issue present? |
| `FininstrmActlXpryDt` | Not verified. | Undetermined (all blank) | Not verified | Blank: 3,682 | Generic UDiFF blank-if-not-applicable rule only [C]. | Definition, expected applicability, date format, and reason for blanks. |
| `StrkPric` | Not verified. | Undetermined (all blank) | Not verified | Blank: 3,682 | Generic UDiFF blank-if-not-applicable rule only [C]. | Definition, unit, precision, and whether blank is expected in these rows. |
| `OptnTp` | Not verified. | Undetermined (all blank) | Not verified | Blank: 3,682 | Generic UDiFF blank-if-not-applicable rule only [C]. | Definition, permitted values, and applicability. |
| `FinInstrmNm` | Not verified at field level. | Text | Not applicable | Populated: 3,682 | No field-specific validation verified [A][B][C]. | Is this a formal instrument name, and what naming/version rules apply? |
| `OpnPric` | Not verified. | Decimal-like text | Not verified | Populated: 3,682 | No field-specific range, precision, or unit verified [A][B][C]. | Exact field meaning, currency/unit, adjustment basis, and price construction. |
| `HghPric` | Not verified. | Decimal-like text | Not verified | Populated: 3,682 | No field-specific range, precision, or unit verified [A][B][C]. | Exact field meaning, currency/unit, adjustment basis, and price construction. |
| `LwPric` | Not verified. | Decimal-like text | Not verified | Populated: 3,682 | No field-specific range, precision, or unit verified [A][B][C]. | Exact field meaning, currency/unit, adjustment basis, and price construction. |
| `ClsPric` | Not verified. | Decimal-like text | Not verified | Populated: 3,682 | No field-specific range, precision, or unit verified [A][B][C]. | Exact field meaning, currency/unit, adjustment basis, and price construction. |
| `LastPric` | Not verified. | Decimal-like text | Not verified | Populated: 3,682 | No field-specific range, precision, or unit verified [A][B][C]. | Exact field meaning, currency/unit, and what event or interval it represents. |
| `PrvsClsgPric` | Not verified. | Decimal-like text | Not verified | Populated: 3,682 | No field-specific range, precision, or unit verified [A][B][C]. | Which prior session/date and adjustment convention does this refer to? |
| `UndrlygPric` | Not verified. | Undetermined (all blank) | Not verified | Blank: 3,682 | Generic UDiFF blank-if-not-applicable rule only [C]. | Definition, applicable instrument classes, and expected blankness. |
| `SttlmPric` | Not verified. | Decimal-like text | Not verified | Populated: 3,682 | No field-specific range, precision, or unit verified [A][B][C]. | Exact settlement-price definition, currency/unit, and applicability. |
| `OpnIntrst` | Not verified. | Undetermined (all blank) | Not verified | Blank: 3,682 | Generic UDiFF blank-if-not-applicable rule only [C]. | Definition, unit, applicability to these records, and reason for blanks. |
| `ChngInOpnIntrst` | Not verified. | Undetermined (all blank) | Not verified | Blank: 3,682 | Generic UDiFF blank-if-not-applicable rule only [C]. | Definition, unit, interval/baseline, and reason for blanks. |
| `TtlTradgVol` | Not verified. | Integer-like text | Not verified | Populated: 3,682 | No field-specific unit, range, or integer rule verified [A][B][C]. | Does this count shares, units, contracts, or another quantity? What is its aggregation scope? |
| `TtlTrfVal` | Not verified. This exact field is present; `TtlTradgVal` is not a CSV header. | Decimal-like text | Not verified | Populated: 3,682 | Generic guidance says documented amount fields use INR and up to 6 decimals [C], but this field's classification as an amount is unverified. | Is this traded/transfer value? What currency, unit scale, and calculation apply? |
| `TtlNbOfTxsExctd` | Not verified. | Integer-like text | Not verified | Populated: 3,682 | No field-specific count definition/range verified [A][B][C]. | What event counts as a transaction, and what is the aggregation scope? |
| `SsnId` | Not verified. Observed value is `F1`. | Text/code | Not applicable | Populated: 3,682; all `F1` | No field-specific standard-value list verified [A][B][C]. | What session does this code identify? |
| `NewBrdLotQty` | Not verified. | Integer-like text | Not verified | Populated: 3,682 | No field-specific unit, range, or validation verified [A][B][C]. | What does â€œboard lotâ€ mean in this field and which security rules apply? |
| `Rmks` | Not verified at field level. | Undetermined (all blank) | Not applicable | Blank: 3,682 | Guidance describes a format-wide remarks field, but its exact mapping/length in this CM catalogue is unverified [C]. | Confirm mapping, max length, and whether blank is valid. |
| `Rsvd1` | Not verified at field level. | Undetermined (all blank) | Not applicable | Blank: 3,682 | Guidance describes four format-wide dummy fields; exact mapping/length is unverified [C]. | Confirm reserved-field definition, max length, and future-use policy. |
| `Rsvd2` | Not verified at field level. | Undetermined (all blank) | Not applicable | Blank: 3,682 | Guidance describes four format-wide dummy fields; exact mapping/length is unverified [C]. | Confirm reserved-field definition, max length, and future-use policy. |
| `Rsvd3` | Not verified at field level. | Undetermined (all blank) | Not applicable | Blank: 3,682 | Guidance describes four format-wide dummy fields; exact mapping/length is unverified [C]. | Confirm reserved-field definition, max length, and future-use policy. |
| `Rsvd4` | Not verified at field level. | Undetermined (all blank) | Not applicable | Blank: 3,682 | Guidance describes four format-wide dummy fields; exact mapping/length is unverified [C]. | Confirm reserved-field definition, max length, and future-use policy. |

## Research-use assessment

### Directly supported by verified contents

- Reproducible structural checks: header set/order, uniform row width, row
  count, blankness, date consistency, and exact code-value counts.
- Descriptive comparisons of raw field values or code groups can be computed
  mechanically, but their economic interpretation is not verified by this
  dictionary. Do not label a numeric column â€œpriceâ€, â€œvolumeâ€, or â€œvalueâ€ in
  an analysis until the official field master confirms it.
- The file can serve as a single-day ingestion/schema-validation fixture.

### Requires additional metadata

- Interpreting any field as an instrument identifier, price, transaction
  measure, quantity, value, settlement, or session requires the CM Bhavcopy
  field master and standard-value list, with version applicability confirmed.
- Comparing instruments requires authoritative identifier, instrument type,
  security-series, currency/unit, market status, and corporate-action or
  adjustment metadata as relevant.

### Requires multiple dates

- Price/quantity changes, returns, volatility, trend, filtering, event
  comparisons, or any temporal model need multiple permitted dates with
  consistent field definitions, coverage, revisions, calendars, and a
  prespecified research design. One common date repeated across rows does not
  make this file a time series.

### Not supported by this file

- Intraday sequencing, order-book reconstruction, or trade-level event analysis:
  no event timestamps or order/trade rows are established here.
- Company fundamentals, point-in-time universe membership, or adjustment
  history.
- Forecasting, investment/trading recommendations, signals, or profitability
  claims.

## Provenance, authorization, and next decision

The file's name and report layout are consistent with NSE's official CM UDiFF
Bhavcopy report family, but this does not prove source custody or legal
permission. Record source location, retrieval date/method, terms reviewed,
permitted intended use, attribution, and retention conditions before using it
for research beyond local schema inspection. **Phase 32 remains
`BLOCKED_AUTHORIZATION`.**

The next technical decision is to obtain/inspect the official CM Bhavcopy field
catalogue and match its version to this report. Until then, treat meanings,
units, field-specific validation, and analytical interpretation as unresolved.
If the catalogue confirms price/quantity semantics and applicable rights, a
later explicitly scoped cross-sectional descriptive study could be considered;
temporal research would still require multiple dates.

## Official sources

- [A. NSE All Reports](https://www.nseindia.com/all-reports) â€” identifies the
  current CM-UDiFF Common Bhavcopy report family.
- [B. NSE Forms & Formats](https://www.nseindia.com/static/resources/forms-formats-members) â€”
  describes UDiFF formats, identifies CM Bhavcopy, and links its file-format
  workbook (updated on the page 2026-06-30; the workbook could not be read by
  the available page reader because it is XLSX).
- [C. UDiFF Guidance Note, Version 1.0 / V2023.10](https://nsearchives.nseindia.com/web/sites/default/files/inline-files/UDiFF%20guidance%20document_Ver1.0.pdf) â€”
  generic structure, date, blank-field, amount, and naming guidance. It
  explicitly directs users to the Excel UDiFF catalogue for field-level
  definitions; this note is not a replacement.
- [Linked UDiFF Trade and Bhavcopy File Formats Workbook](https://nsearchives.nseindia.com/web/mediaattachment/2026-06/UDiFF_trade_and_Bhavcopy_file_formats_20260630115803.xlsx) â€”
  official field-format resource identified on [B], not inspected because its
  XLSX content type was unsupported by the reader.

[nse-reports]: https://www.nseindia.com/all-reports
[nse-formats]: https://www.nseindia.com/static/resources/forms-formats-members
[udiff-guidance]: https://nsearchives.nseindia.com/web/sites/default/files/inline-files/UDiFF%20guidance%20document_Ver1.0.pdf

# Phase 33: Official Research Data Source Audit

**Reviewed:** 2026-09-29
**Scope:** source discovery and local file inspection only. No data were
downloaded, scraped, submitted, or uploaded for this audit. A public page is not
evidence that a particular dataset is available, licensed, or authorized for a
proposed study.

## Authorization and study readiness

Phase 32 remains **`BLOCKED_AUTHORIZATION`**. This catalogue does not identify a
selected provider/product, establish the applicable terms for a specific file,
or evidence that a stated intended use is permitted. A later study still needs
the Phase 32 source, product/version, intended-use scope, terms reviewed,
permission evidence, and research specification. Public reports do not satisfy
that gate by themselves.

## Keep these data types distinct

| Data type | What it represents | Audit caution |
|---|---|---|
| Daily security-level observations | Per-security end-of-day values or daily report rows | Verify exact fields, adjustment policy, identifiers, and calendar in the chosen file. |
| Intraday or order/trade-level market data | Timestamped transactions, orders, or intraday observations | EOD products can include trade details without being a complete intraday event feed. Verify sampling and event definitions. |
| Company fundamentals | Company/corporate financial or reference information | A market-data policy mentioning company data does not establish that a given fundamentals dataset is supplied. |
| Electricity delivery-block prices | Price observations for defined market delivery intervals | Do not infer block frequency or history from a chart or monthly report. Verify product specifications and files. |
| Monthly electricity-market summaries | Aggregated monthly market-monitoring statistics/reports | Not equivalent to delivery-block prices or raw trades. |
| Economic indicators | Macroeconomic, financial, public-finance, or social time series | Each series can have its own unit, frequency, revision, and coverage. |
| Household/distribution electricity consumption | Meter or client load measurements | These are not wholesale electricity prices or Indian grid-market observations. |

## Source catalogue

Facts below are limited to what the linked official pages exposed on the review
date. Dynamic pages and report-specific files may change. â€œNot verifiedâ€ means
this audit did not establish the fact; it does not mean the source lacks it.

### Indian securities-market data

#### NSE All Reports

- **Official page:** [NSE All Reports](https://www.nseindia.com/all-reports)
- **Domain/use:** Indian securities-market reports; useful for locating daily
  exchange report files, subject to confirming the selected file's fields.
- **Granularity/coverage verified:** The page separates current and historical
  reports, includes equity/index categories, and exposes report/date selection.
  Its current list included a daily â€œCM-UDiFF Common Bhavcopy Final (zip)â€
  report. The page notes that the older CM Bhavcopy/common Bhavcopy CSV was
  discontinued from 2024-07-08. This does not establish the full history or
  schema of the replacement file.
- **Access/price:** Report selection and file download are exposed on the page;
  no price was verified for each report by this page.
- **Type:** Downloadable exchange report files, including daily reports; not a
  general-purpose grant of rights for reuse.
- **Terms/unknowns:** The page itself did not establish a license for the
  selected file. Check the report specification, applicable NSE data policy,
  access conditions, field meanings, adjustments, and reuse rights before use.
- **Reviewed:** 2026-09-29.

#### NSE EOD / Historical Data subscription

- **Official page:** [NSE EOD / Historical Data](https://www.nseindia.com/static/market-data/eod-historical-data-subscription)
- **Domain/use:** End-of-day and historical securities-market products.
- **Granularity/coverage verified:** NSE describes end-of-day files generated
  after each trading day, with EOD bhavcopy/security/trade details for Capital
  Market, Wholesale Debt, and F&O. The page separately describes historical
  order/trade data through SFTP and historical O&T data through an online
  platform, with listed segments including CM, F&O, CD, and COM. Product
  specifications must be checked before equating these with intraday complete
  order books or every trade event.
- **Access/price:** Subscription products; the page links a tariff document.
  Exact current price was not verified in this audit. The page was marked
  updated 2026-09-02.
- **Type:** Paid/subscription access described for downloadable files and
  historical O&T products.
- **Terms/unknowns:** Confirm eligible periods, exact fields, delivery cadence,
  subscription fee, user/entity scope, retention, derived outputs, and
  redistribution rights in the current product documents and agreement.
- **Reviewed:** 2026-09-29.

#### NSE Research Initiative

- **Official page:** [NSE Research Initiative](https://www.nseindia.com/static/research/research-initiatives)
- **Domain/use:** Formal research access route for qualifying academic or
  noncommercial research, as described by NSE.
- **Granularity/coverage verified:** The page describes a data-seeking
  application and an approval-based, time-bound data-room route. It does not
  establish which dataset, fields, or date coverage a particular applicant
  will receive.
- **Access/price:** Formal request/application is described; the page directs
  applicants to send the data-seeking form to `nseri@nse.co.in`. No request was
  made and no price was verified.
- **Type:** Formal access request and, where approved, controlled data-room
  access; not a public raw-data download catalogue.
- **Terms/unknowns:** Eligibility, product, dates, permitted analysis, storage,
  publication, confidentiality, fees, and access term are case-specific and
  require written confirmation. Page language about academic/noncommercial
  use is not approval for any particular study.
- **Reviewed:** 2026-09-29.

#### NSE Data Sharing & Usage Policy

- **Official page:** [NSE Data Sharing & Usage Policy](https://www.nseindia.com/static/market-data/nse-data-policy)
- **Domain/use:** Data-rights and policy framework, not a dataset.
- **Coverage verified:** The policy discusses price, identifiers, volume/trade,
  and company/corporate data across equity, debt, and derivatives, with
  real-time, streaming, snapshot, delayed, EOD, historical, and corporate-data
  categories.
- **Access/price:** It describes subscriber requests and agreements that
  specify intended use/handling. Fees are described as arm's-length; the policy
  mentions possible noncommercial fee arrangements/waivers, not a guaranteed
  waiver. The page reports policy update 2025-11-12.
- **Terms:** Redistribution is restricted except as an agreement permits. The
  policy also describes research access via the research initiative and a
  confidentiality undertaking. The actual signed agreement/product terms
  govern a particular acquisition.
- **Type:** Rights/policy document; it supplies no raw observations.
- **Unknowns:** Whether a particular file/product and intended use are
  available/permitted, fee, permitted retention, publication, derived data, and
  redistribution must be confirmed in writing for the specific case.
- **Reviewed:** 2026-09-29.

### Indian electricity-market data

#### IEX India

- **Official pages:** [IEX India](https://www.iexindia.com/),
  [IEX site map](https://www.iexindia.com/site-map),
  [IEX Terms of Use](https://www.iexindia.com/terms-of-uses)
- **Domain/use:** Indian electricity exchange. The homepage presents DAM and
  RTM price trends in â‚¹/kWh and daily trade volume/price updates; the sitemap
  lists market snapshots and aggregate DAM/RTM demand/supply, plus TAM trade
  details, monthly trade details, trade summary, and snapshots.
- **Granularity/coverage verified:** The public pages establish that those
  product/report categories exist, not a complete block-level data series,
  exact interval length, fields, or historical coverage. The `/marketdata`
  page returned a â€œSomething went wrongâ€ message during this review.
- **Access/price:** Homepage offers daily updates by WhatsApp registration;
  sitemap lists public-facing market-data sections. No complete dataset
  download, paid tariff, or historical product price was verified.
- **Type:** Public overview/snapshot and aggregate report categories; product
  access for transaction-grade data remains unverified.
- **Terms:** Website terms state personal/noncommercial use of website content,
  restrict copying/extraction/distribution absent express IEX permission, and
  require attribution. These website-content terms do not establish the terms
  of a separately licensed market-data product. WhatsApp registration is not
  evidence of download or research permission.
- **Unknowns:** Delivery interval, raw versus aggregate fields, history, formal
  data product, fees, licensing, retention, derived outputs, and redistribution
  require confirmation from the relevant IEX product terms. No access request
  was sent.
- **Reviewed:** 2026-09-29.

#### CEA Market Monitoring Reports

- **Official page:** [CEA Market Monitoring Reports](https://cea.nic.in/market-monitoring-report/?lang=en)
- **Domain/use:** Central Electricity Authority reports about Indian
  electricity-market monitoring.
- **Granularity/coverage verified:** The page offers a month selector and
  report entries with PDF and Excel view/download options. It displayed a June
  2026 entry when reviewed (page updated 2026-09-24). Treat these as monthly
  monitoring reports/aggregates; this page does not verify delivery-block
  price observations.
- **Access/price:** Report links are publicly exposed on the official page; no
  price was shown.
- **Type:** Report/document and spreadsheet outputs, not verified raw
  transaction-level observations.
- **Terms/unknowns:** The landing page did not establish a dataset license or
  redistribution terms. Inspect the selected report's fields, time span,
  methodology, attribution, and terms before reuse.
- **Reviewed:** 2026-09-29.

#### CERC Market Monitoring Reports

- **Official page:** [CERC Market Monitoring Reports](https://cercind.gov.in/report_MM-2025.html)
- **Domain/use:** Central Electricity Regulatory Commission monthly reports
  concerning short-term electricity transactions in India.
- **Granularity/coverage verified:** The page lists monthly reports for
  Januaryâ€“December 2025 and provides report/PDF and Excel data links. This is
  monthly reporting; it does not itself establish block-level price series.
- **Access/price:** Official page links reports and Excel data; no price was
  shown.
- **Type:** Monthly reports and spreadsheet data, likely aggregated; inspect
  actual workbook fields before treating it as observations at any finer
  frequency.
- **Terms:** The page says material may be used freely with due
  acknowledgement to CERC and that the PDF is final in case of discrepancy.
  Preserve attribution and verify any report-specific conditions.
- **Unknowns:** Exact workbook fields, coverage beyond the listed period,
  revision policy, and granularity of individual tables require inspection.
- **Reviewed:** 2026-09-29.

### Indian economic data

#### RBI Database on Indian Economy (DBIE)

- **User-specified official URL:** [RBI DBIE legacy URL](https://dbieold.rbi.org.in/DBIE/)
- **Official notice/current portal reference:** [RBI Annual Report 2024â€“25](https://rbi.org.in/scripts/AnnualReportPublications.aspx?Id=1440)
  reports the DBIE URL change to [data.rbi.org.in](https://data.rbi.org.in/).
- **Domain/use:** Macroeconomic and financial time series. Official DBIE
  descriptions span real, corporate, financial/financial-market, external,
  public-finance, and socio-economic statistics, including time-series and
  other data.
- **Granularity/coverage verified:** Varies by series. No particular series,
  frequency, units, historical coverage, or revision policy was verified for
  this audit. The legacy URL timed out during review; the current portal did not
  load successfully in the audit browser, so current interface details remain
  unverified.
- **Access/price:** Official RBI material describes a data-query portal. No
  series was queried or downloaded and no price was verified.
- **Type:** Statistical indicators/time series, not security-level market
  trades or electricity delivery-block prices.
- **Terms/unknowns:** Check series metadata, definitions, units, release and
  revision dates, export method, attribution, reuse, and redistribution terms
  for the selected indicator. Do not assume all DBIE series share one
  frequency or license.
- **Reviewed:** 2026-09-29.

### Open research datasets

#### UCI Individual Household Electric Power Consumption

- **Official page:** [UCI dataset 235](https://archive.ics.uci.edu/dataset/235/individual%2Bhousehold%2Belectric%2Bpower%2Bconsumption)
- **Domain/use:** Household electricity consumption for signal processing and
  consumption-pattern research; not Indian wholesale-market data.
- **Granularity/coverage verified:** One household in Sceaux, France; one
  measurement per minute from December 2006 through November 2010; 2,075,259
  measurements and 9 features. UCI reports missing values in about 1.25% of
  rows. Fields include date/time, active/reactive power, voltage, current, and
  sub-metering.
- **Access/price:** UCI page exposes a download link labeled 19.7 MB and an
  optional `ucimlrepo` access method; its dataset-files table lists a 126.8 MB
  artifact. The page does not explain this size difference. No price was
  displayed. No client was used and no data were fetched.
- **Type:** Household meter observations in a downloadable dataset.
- **License/attribution:** CC BY 4.0 as stated on the dataset page; attribution
  is required. Preserve dataset citation and review the license before
  redistribution or publication.
- **Unknowns:** Confirm file/version, missing-data handling, timestamp
  interpretation, local use, and any downstream publication requirements.
- **Reviewed:** 2026-09-29.

#### UCI Electricity Load Diagrams 2011â€“2014

- **Official page:** [UCI dataset 321](https://archive.ics.uci.edu/dataset/321/electricity)
- **Domain/use:** Multi-client electricity-load time series; useful for
  distribution/load signal analysis, not Indian wholesale prices.
- **Granularity/coverage verified:** 370 client series, 15-minute values in
  kW, named for 2011â€“2014. The page reports no missing values and documents
  daylight-saving conventions: 96 samples per day, a zero-valued spring
  1â€“2-hour interval, and a fall 1â€“2-hour interval aggregating two hours; later
  client entries are zeroed as described on the page.
- **Access/price:** Semicolon-separated file, 249.2 MB, and optional
  `ucimlrepo` access are described on the page. No data were fetched.
- **Type:** Multi-series load observations, not transaction-level market data.
- **License/attribution:** CC BY 4.0 per the dataset page; attribute the
  dataset and review license conditions before redistribution/publication.
- **Unknowns:** Verify exact file version, client anonymization/zero conventions,
  timestamp/calendar interpretation, and suitability for any specific
  comparison.
- **Reviewed:** 2026-09-29.

#### UCI power-consumption dataset catalogue

- **Official page:** [UCI catalogue search: power consumption](https://archive.ics.uci.edu/datasets?Keywords=power+consumption)
- **Domain/use:** Discovery index for datasets matching a keyword, not a
  unified dataset or license.
- **Granularity/coverage verified:** Results include heterogeneous datasets;
  the page surfaced the household dataset above and other regional/power
  datasets. No common sampling frequency or historical coverage applies.
- **Access/price:** Catalogue links to individual dataset pages; access method
  and price depend on the selected dataset.
- **Type:** Search/index page only; it contains no common raw data file.
- **Terms/unknowns:** Review the chosen dataset's page for its actual fields,
  frequency, coverage, license, attribution, and file version. Do not infer
  those properties from a catalogue keyword match.
- **Reviewed:** 2026-09-29.

## Manual local-data workflow

1. First obtain a CSV or report through a provider's permitted method and
   resolve intended-use permission. This phase does not provide or request
   access.
2. Place original files under `data/raw/` and keep transformed tables under
   `data/derived/`. These paths are ignored by the repository rules; do not add
   raw or derived datasets to source control. Keep any provider-mandated
   retention/deletion restrictions in force.
3. Maintain a separate provenance note recording provider/source URL,
   product and version, retrieval date, filename, SHA-256 hash, terms/license
   reviewed and date, intended/permitted use, and attribution requirements.
   A file hash identifies bytes; it does not establish provenance or permission.
4. Inspect a local CSV with the read-only helper (installed project):

   ```powershell
   .venv\Scripts\python.exe -m newton_lab.data_inspection `
     data\raw\permitted_file.csv --timestamp-column timestamp
   ```

   Or call `inspect_local_csv(path, timestamp_column="timestamp")` from Python.
   Output includes columns, row count, null counts, duplicate whole rows,
   invalid and repeated timestamps, date bounds, and a diagnostic inferred
   frequency. Empty cells and the text markers `NA`, `N/A`, `NaN`, and `null`
   (case-insensitive) count as missing. Timestamp parsing expects
   ISO-formatted values. The helper only
   reads a local CSV; it does not fetch URLs, infer provider semantics, modify
   the file, or establish that the data are authorized. Frequency inference
   can be `null` for short or irregular data; a constant interval is reported
   as a simple duration such as `1h` or `15min`.
5. Inspect report formats manually with an appropriate local reader before
   converting them. Preserve the untouched source file and record each
   transformation separately. Never treat a spreadsheet export as raw event
   data without checking its fields and aggregation.

## Research limits

Source discovery is not authorization or a study. A mathematical
correspondence, synthetic result, or data-inspection report does not establish
financial validity, forecasting value, trading profitability, or a validated
application in another domain. For NSE/IEX and other licensed data, follow the
actual product agreement and document intended use before analysis or sharing.

# Phase 37 — Historical Market Data Acquisition Plan

**Research reviewed:** 2026-09-29 (Asia/Kolkata)

## Objective and present status

Identify a practical route to a multi-date Indian dataset for a future time-series
study. The existing `data/raw/user_supplied_market_data.csv` has 3,682 rows for
one date, 2026-09-28. It is a single-date cross-section and cannot support a
temporal study by itself.

The existing file has no adjacent provenance manifest, receipt, terms review, or
retrieval record in `data/raw/` or the data-tree records inspected. Its filesystem
timestamp is not evidence of source or permission. Phase 35 says its field
definitions, units, provenance, and permitted use remain unverified. Phase 32
continues to state **`BLOCKED_AUTHORIZATION`**; this plan does not change it.

No additional historical dataset was obtained, and no download or access request
was made in Phase 37.

## Candidate routes, ranked by practical fit

This ranks next steps for investigation, not exchanges or financial choices.
“Access shown” describes what the cited page documents; it does not establish a
right to use, retain, publish, or redistribute data.

| Rank | Source and candidate | What official pages establish | Access, cost, and limitations | Fit for this study |
|---:|---|---|---|---|
| 1 | NSE All Reports — **CM-UDiFF Common Bhavcopy Final (zip)** | The page has Current Reports and Historical Reports sections, date/report controls, and lists this report. It says the older CM Bhavcopy CSV and CM Common Bhavcopy CSV were discontinued from 2024-07-08 and directs users to the UDiFF Common Bhavcopy Final ZIP. NSE's page also presents EOD reports as daily reports. | A manual historical-report selector is visible, but the public page does not document its earliest date, completeness, exact date-by-date availability, free status, or use rights. No current terms specific to downloading this report for this user's intended use were established. Do not treat visible access as permission. | Best first availability check because the existing file purports to be a single-day NSE-style CM report, though its provenance and schema match are not verified. A selected instrument and consistent date series still need to be confirmed from the official specification. |
| 2 | NSE paid EOD/historical data; separately, NSE research data/Data Room | The paid-data page says EOD files are generated at the end of each trading day and include bhavcopy plus security/trade details for listed segments. It describes binary files, technical specifications, sample files, and an enquiry route by email or phone. The NSE policy requires a subscriber request and a relevant agreement that sets intended use and handling. For research, requests are routed through NSE's Economic Policy & Research department; conditions are documented case by case. NSE's current Research Initiatives page says research data is limited to academic and non-commercial use, with a request form, and describes proposal review for Data Room access. | NSE publishes an EOD tariff PDF; its displayed amount is **₹100,000 domestic / US$5,500 international for Capital Market**, with the tariff notes framed per display medium. This is not a verified quote for an individual local research license or this project's use. Research/data-room eligibility, exact fields and dates, fees, turnaround, retention, publication, and any permission needed for non-academic/private work remain to be confirmed. | Most direct official route if the self-service archive does not supply a sufficient and permitted series. The Data Room route may not fit a private/non-academic project and is not guaranteed approval. |
| 3 | BSE Historical Data / daily EOD products | BSE's official Information Products tariff sheet describes daily EOD Bhav Copy and historical stock-trading/corporate data, including availability to researchers and academicians through its information products. | The sheet describes access by subscription. Exact daily equity product, historical date span, manual file handoff, fees for this case, and research/redistribution terms were not verified. A product listing is not a free download license. | A possible alternative only if the study chooses a BSE-listed instrument and BSE can document a consistent multi-date series and terms. It is a separate exchange and must not be spliced into NSE records without an explicit harmonization study. |
| 4 | RBI Database on Indian Economy (DBIE) | RBI DBIE provides time-series publications across economic and financial subject areas. Official DBIE help material describes historical time-series data, Excel downloads, and access without login. | Frequency and coverage depend on the selected table/indicator; there is no one universal period. No security-level daily share observations are established. Public browsing/download does not by itself settle reuse terms. | Practical for a macroeconomic time-series question, not a substitute for a daily individual-security dataset. The Phase 32 research question and target would need to be revised first. |
| 5 | CEA / CERC market-monitoring publications | CEA lists monthly Market Monitoring Reports with PDF and Excel options; its Regulatory Affairs description says it prepares brief/detailed monthly and annual reports on short-term electricity-market transactions. CERC's 2025 page lists monthly reports for all twelve months; a CERC report says this monthly report series has been published since 2008. | These are reports and summaries, not demonstrated daily instrument-level observations. Exact archive coverage, underlying workbook fields, permissions, and reuse terms must be checked. The Excel link may still contain aggregates rather than a raw time series. | Useful only for a separately defined electricity-market question at the frequency and granularity actually present. Not suitable for an NSE security-level study. |

### Evidence links checked

- NSE [All Reports](https://www.nseindia.com/all-reports) — Current/Historical Reports sections and named CM report family.
- NSE [Paid End of the day/Historical Data](https://www.nseindia.com/static/market-data/eod-historical-data-subscription) — EOD descriptions, subscription enquiry, specifications and tariff link. Page showed an update date of 2026-09-02 when reviewed.
- NSE [Data Sharing & Usage Policy](https://www.nseindia.com/static/market-data/nse-data-policy) — subscriber agreement and research access provisions.
- NSE [Research Initiatives](https://www.nseindia.com/static/research/research-initiatives) — data-seeking form, stated academic/non-commercial limit, and Data Room process.
- BSE [Information Products Domestic Tariff Sheet](https://www.bseindia.com/downloads1/Information_Products_Pricing_Sheet.pdf).
- RBI [DBIE](https://data.rbi.org.in/) and official [DBIE FAQ](https://dbieold.rbi.org.in/DBIE/doc/Frequently%20Asked%20Questions%20-%20DBIE.pdf).
- CEA [Market Monitoring Report page](https://cea.nic.in/market-monitoring-report/?lang=en) and [Regulatory Affairs Division](https://cea.nic.in/regulatory-affairs-division/?lang=en).
- CERC [Market Monitoring Reports 2025](https://cercind.gov.in/report_MM-2025.html) and its [June 2025 report](https://www.cercind.gov.in/2025/market_monitoring/MMC%20Report%20on%20Short%20term%20market%20for%20June%202025.pdf).

The NSE public report listing was reviewed on 2026-09-29. The EOD subscription
page says it was updated 2026-09-02. The BSE tariff sheet, RBI pages, CEA pages,
and CERC pages were reviewed on 2026-09-29; their contents do not establish a
single applicable license or complete date range for this project.

## Minimum viable dataset

Before any observations are acquired for analysis, specify one series with:

1. One chosen instrument or explicitly defined aggregate, with a documented,
   stable identifier and any symbol/identifier changes accounted for.
2. Multiple distinct dates and a documented observation frequency. For a daily
   series, the date must mean a specific trading date/session and use a defined
   trading calendar.
3. A documented target field, definition, unit, adjustment/vintage status, and
   availability time. Do not choose a column from its abbreviation alone.
4. Consistent schema and definitions across files; confirm report version and
   historical coverage before combining daily files.
5. Enough observations for prespecified chronological training, validation, and
   untouched test partitions. The sample length and boundaries depend on the
   eventual target and horizon and remain undecided.
6. Documented provenance, integrity hashes, applicable terms, intended-use
   permissions, retention, attribution, and any publication/sharing limits.

Do not concatenate the supplied one-day CSV with future report files until the
official report family, version, field definitions, identifier semantics, and
schema consistency are established.

## Recommended next route and manual actions

**First check the NSE Historical Reports interface for the exact report
“CM-UDiFF Common Bhavcopy Final (zip).”** It is the closest visible report family
to the already supplied single-day file, but a matching identity is not yet
verified.

1. Open [NSE All Reports](https://www.nseindia.com/all-reports), select
   **Historical Reports**, and look for **CM-UDiFF Common Bhavcopy Final (zip)**.
2. Before downloading, check whether the manual interface actually offers the
   date range needed for a feasible study, and note the first/last available
   dates, any gaps, download format, and report/version label. The official page
   reviewed here does not establish archive depth; if the UI does not expose it,
   mark coverage unknown.
3. Compare the report's official format specification with the supplied file and
   Phase 35 dictionary. Verify the target field, units, date convention, stable
   instrument key, corporate-action treatment, and whether each selected date has
   compatible definitions. Do not infer these from the filename or abbreviations.
4. Review the applicable NSE policy and report-specific conditions for the exact
   intended use, local analysis, raw-data storage/retention, derived results,
   attribution, publication, redistribution, and whether a separate agreement is
   required. If terms do not clearly authorize the intended use, stop and keep
   authorization unresolved.
5. If a sufficient self-service archive or its applicable permissions cannot be
   established, use the official NSE contact/request route to ask for a written
   quote and terms for a **small daily Capital Market series for one specified
   instrument and fixed date range**. Ask separately whether the academic,
   non-commercial research path is available to this project. Do not presume
   eligibility or submit a request as part of this phase.

This is a discovery route, not a finding that the series is freely downloadable
or licensed for this project's use. The public materials reviewed do **not** let
us verify a route that is both free and lawfully usable for the intended study.
The formal NSE route is actionable, but its eligibility, price, terms, coverage,
and permission must be confirmed before data access or analysis. Given the broader
project interest in quantitative finance/trading, explicitly disclose the
intended use when asking; do not assume that a non-commercial research license
would cover any later trading or commercial purpose.

## Data-use and provenance checklist

For a dataset obtained later, record this information before analysis, subject to
the provider's terms:

- Provider, exact product/report, segment, version, instrument identifier, and
  requested date range.
- Actual retrieval/provision date, manual method, filenames, file sizes, and
  SHA-256 for each unmodified raw file.
- Official field specification and effective version; date semantics, frequency,
  units, identifier mapping, revisions, corporate actions, and coverage gaps.
- Terms/license URL and version/date; exact intended use; review date and evidence;
  permission for local use, raw retention, derived-result retention, publication,
  attribution, and redistribution; restrictions, fees, expiry, and renewal.
- Keep raw files separate from derived files. Do not put data in source code,
  reports, or version control. This plan neither changes nor grants permission for
  the current CSV.

## Stopping condition

The stopping condition for **acquisition planning** is met when the user has a
specific product/report, a documented manual or formal access route, a verified
date range and schema suitable for the stated series, a price/access outcome, and
applicable use and retention conditions understood. It is **not met yet**:
archive depth, report/schema match, project eligibility, fees for this particular
use, and authorization remain unresolved. No dataset should be ingested under
Phase 32 until its authorization and applicable pre-ingestion requirements are
documented and reviewed.

## Project boundary and verification

This is a documentation-only phase. No source files, CSV, Phase 36 artifacts,
energy-demand work, or Phase 29/30 studies were changed. No packages were added,
no data were downloaded, and no access request or message was sent. Tests were
not run because only documentation was added.

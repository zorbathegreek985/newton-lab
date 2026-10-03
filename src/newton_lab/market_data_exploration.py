"""Schema-aware exploratory summaries for a local, single-day market CSV.

Field names and numeric shape do not establish economic meanings or units.
This module profiles raw values without cleaning, dropping, or imputing rows.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
import sys
from collections import Counter
from collections.abc import Sequence
from dataclasses import asdict, dataclass
from datetime import date
from pathlib import Path
from statistics import fmean, pstdev
from typing import TYPE_CHECKING

from newton_lab.exceptions import ScientificValidationError

if TYPE_CHECKING:
    from matplotlib.figure import Figure

_MISSING_MARKERS = {"", "na", "n/a", "nan", "null"}
_INTEGER_PATTERN = re.compile(r"^[+-]?\d+$")
_CANDIDATE_IDENTIFIER_COLUMNS = ("FinInstrmId", "ISIN", "TckrSymb")


@dataclass(frozen=True)
class ColumnProfile:
    """Observed values and inferred lexical type for one exact CSV column."""

    name: str
    inferred_type: str
    row_count: int
    missing_count: int
    non_missing_count: int
    distinct_count: int
    singleton_value_count: int


@dataclass(frozen=True)
class CategoryValueCount:
    """Observed raw text level and its row count."""

    value: str
    count: int


@dataclass(frozen=True)
class CategoricalSummary:
    """Frequency profile for a text-like column, without interpreting codes."""

    name: str
    missing_count: int
    non_missing_count: int
    distinct_count: int
    singleton_value_count: int
    value_counts: tuple[CategoryValueCount, ...]
    counts_truncated: bool


@dataclass(frozen=True)
class NumericSummary:
    """Descriptive statistics for values that are numeric-like in the file."""

    name: str
    inferred_type: str
    count: int
    missing_count: int
    distinct_count: int
    minimum: float
    percentile_25: float
    median: float
    percentile_75: float
    maximum: float
    mean: float
    population_standard_deviation: float
    zero_count: int
    negative_count: int
    positive_count: int
    tukey_fence_flag_count: int | None


@dataclass(frozen=True)
class DateCoverage:
    """Observed ISO date range for a date-like column."""

    name: str
    valid_count: int
    missing_count: int
    invalid_count: int
    distinct_count: int
    first_date: str | None
    last_date: str | None


@dataclass(frozen=True)
class IdentifierCandidate:
    """Uniqueness diagnostics for a possible identifier field."""

    name: str
    non_missing_count: int
    distinct_count: int
    repeated_value_count: int
    unique_within_file: bool


@dataclass(frozen=True)
class CorrelationPair:
    """Pearson correlation of two numeric-like value columns."""

    left_column: str
    right_column: str
    pair_count: int
    pearson_r: float | None


@dataclass(frozen=True)
class MarketDataExploration:
    """Reproducible descriptive profile of one CSV file."""

    row_count: int
    column_count: int
    duplicate_row_count: int
    columns: tuple[ColumnProfile, ...]
    categorical_summaries: tuple[CategoricalSummary, ...]
    numeric_summaries: tuple[NumericSummary, ...]
    date_coverage: tuple[DateCoverage, ...]
    identifier_candidates: tuple[IdentifierCandidate, ...]
    numeric_column_correlations: tuple[CorrelationPair, ...]
    warnings: tuple[str, ...]


def _is_missing(value: str | None) -> bool:
    return value is None or value.strip().casefold() in _MISSING_MARKERS


def _infer_type(values: Sequence[str]) -> str:
    present = [value.strip() for value in values if not _is_missing(value)]
    if not present:
        return "all_missing"
    if all(_is_iso_date(value) for value in present):
        return "date_like"
    if all(_INTEGER_PATTERN.fullmatch(value) for value in present):
        return "integer_like"
    numeric_count = sum(_is_finite_number(value) for value in present)
    if numeric_count == len(present):
        return "decimal_like"
    if numeric_count:
        return "mixed"
    return "categorical_text"


def _is_iso_date(value: str) -> bool:
    try:
        date.fromisoformat(value)
    except ValueError:
        return False
    return len(value) == 10


def _is_finite_number(value: str) -> bool:
    try:
        return math.isfinite(float(value))
    except ValueError:
        return False


def _quantile(sorted_values: Sequence[float], fraction: float) -> float:
    position = (len(sorted_values) - 1) * fraction
    lower_index = math.floor(position)
    upper_index = math.ceil(position)
    if lower_index == upper_index:
        return sorted_values[lower_index]
    weight = position - lower_index
    return (
        sorted_values[lower_index] * (1.0 - weight)
        + sorted_values[upper_index] * weight
    )


def _numeric_values(values: Sequence[str]) -> list[float]:
    return [float(value) for value in values if not _is_missing(value)]


def _numeric_summary(
    name: str, values: Sequence[str], inferred_type: str
) -> NumericSummary:
    numeric = _numeric_values(values)
    ordered = sorted(numeric)
    q1 = _quantile(ordered, 0.25)
    median = _quantile(ordered, 0.5)
    q3 = _quantile(ordered, 0.75)
    tukey_count: int | None = None
    if inferred_type == "decimal_like" and name not in _CANDIDATE_IDENTIFIER_COLUMNS:
        lower_fence = q1 - 1.5 * (q3 - q1)
        upper_fence = q3 + 1.5 * (q3 - q1)
        tukey_count = sum(
            value < lower_fence or value > upper_fence for value in numeric
        )
    return NumericSummary(
        name=name,
        inferred_type=inferred_type,
        count=len(numeric),
        missing_count=len(values) - len(numeric),
        distinct_count=len(set(numeric)),
        minimum=ordered[0],
        percentile_25=q1,
        median=median,
        percentile_75=q3,
        maximum=ordered[-1],
        mean=fmean(numeric),
        population_standard_deviation=pstdev(numeric),
        zero_count=sum(value == 0 for value in numeric),
        negative_count=sum(value < 0 for value in numeric),
        positive_count=sum(value > 0 for value in numeric),
        tukey_fence_flag_count=tukey_count,
    )


def _date_coverage(name: str, values: Sequence[str]) -> DateCoverage:
    valid_dates: list[date] = []
    missing_count = 0
    invalid_count = 0
    for value in values:
        if _is_missing(value):
            missing_count += 1
        else:
            try:
                valid_dates.append(date.fromisoformat(value.strip()))
            except ValueError:
                invalid_count += 1
    unique_dates = sorted(set(valid_dates))
    return DateCoverage(
        name=name,
        valid_count=len(valid_dates),
        missing_count=missing_count,
        invalid_count=invalid_count,
        distinct_count=len(unique_dates),
        first_date=unique_dates[0].isoformat() if unique_dates else None,
        last_date=unique_dates[-1].isoformat() if unique_dates else None,
    )


def _pearson(left: Sequence[float], right: Sequence[float]) -> float | None:
    if len(left) != len(right) or len(left) < 2:
        return None
    left_mean = fmean(left)
    right_mean = fmean(right)
    left_centered = [value - left_mean for value in left]
    right_centered = [value - right_mean for value in right]
    left_sum_sq = sum(value * value for value in left_centered)
    right_sum_sq = sum(value * value for value in right_centered)
    denominator = math.sqrt(left_sum_sq * right_sum_sq)
    if denominator == 0:
        return None
    return (
        sum(
            left_value * right_value
            for left_value, right_value in zip(
                left_centered, right_centered, strict=True
            )
        )
        / denominator
    )


def analyze_market_csv(path: str | Path) -> MarketDataExploration:
    """Read and summarize a local CSV without mutating or filtering its rows."""
    csv_path = Path(path)
    if csv_path.suffix.lower() != ".csv" or not csv_path.is_file():
        raise ScientificValidationError("path must refer to an existing local CSV file")
    try:
        with csv_path.open(encoding="utf-8-sig", newline="") as handle:
            reader = csv.DictReader(handle)
            if not reader.fieldnames:
                raise ScientificValidationError("CSV must have a header row")
            headers = tuple(reader.fieldnames)
            if len(set(headers)) != len(headers):
                raise ScientificValidationError("CSV headers must be unique")
            rows = []
            for row_number, row in enumerate(reader, start=2):
                if row.get(None) is not None or any(
                    row.get(name) is None for name in headers
                ):
                    raise ScientificValidationError(
                        f"CSV row {row_number} does not match the header width"
                    )
                rows.append(row)
    except ScientificValidationError:
        raise
    except (OSError, UnicodeError, csv.Error) as exc:
        raise ScientificValidationError(f"could not read local CSV: {exc}") from exc

    if not headers:
        raise ScientificValidationError("CSV must have at least one column")
    values_by_column = {
        name: [row.get(name, "") or "" for row in rows] for name in headers
    }
    columns: list[ColumnProfile] = []
    categorical: list[CategoricalSummary] = []
    numeric: list[NumericSummary] = []
    date_ranges: list[DateCoverage] = []
    values_by_numeric_column: dict[str, list[str]] = {}

    for name in headers:
        values = values_by_column[name]
        present = [value for value in values if not _is_missing(value)]
        kind = _infer_type(values)
        counts = Counter(present)
        distinct_count = len(counts)
        columns.append(
            ColumnProfile(
                name=name,
                inferred_type=kind,
                row_count=len(rows),
                missing_count=len(rows) - len(present),
                non_missing_count=len(present),
                distinct_count=distinct_count,
                singleton_value_count=sum(count == 1 for count in counts.values()),
            )
        )
        if kind in {"categorical_text", "mixed"}:
            ordered_counts = sorted(
                counts.items(), key=lambda pair: (-pair[1], pair[0])
            )
            limit = 60
            categorical.append(
                CategoricalSummary(
                    name=name,
                    missing_count=len(rows) - len(present),
                    non_missing_count=len(present),
                    distinct_count=distinct_count,
                    singleton_value_count=sum(count == 1 for count in counts.values()),
                    value_counts=tuple(
                        CategoryValueCount(value=value, count=count)
                        for value, count in ordered_counts[:limit]
                    ),
                    counts_truncated=distinct_count > limit,
                )
            )
        elif kind in {"integer_like", "decimal_like"}:
            numeric.append(_numeric_summary(name, values, kind))
            if name not in _CANDIDATE_IDENTIFIER_COLUMNS:
                values_by_numeric_column[name] = values
        elif kind == "date_like":
            date_ranges.append(_date_coverage(name, values))

    duplicate_count = len(rows) - len(
        {tuple(row.get(name, "") for name in headers) for row in rows}
    )
    candidates = tuple(
        IdentifierCandidate(
            name=name,
            non_missing_count=profile.non_missing_count,
            distinct_count=profile.distinct_count,
            repeated_value_count=profile.non_missing_count - profile.distinct_count,
            unique_within_file=(
                profile.non_missing_count > 0
                and profile.non_missing_count == profile.distinct_count
            ),
        )
        for name in _CANDIDATE_IDENTIFIER_COLUMNS
        for profile in columns
        if profile.name == name
    )
    correlations: list[CorrelationPair] = []
    decimal_names = tuple(values_by_numeric_column)
    for left_index, left_name in enumerate(decimal_names):
        for right_name in decimal_names[left_index + 1 :]:
            paired_values = [
                (float(left), float(right))
                for left, right in zip(
                    values_by_numeric_column[left_name],
                    values_by_numeric_column[right_name],
                    strict=True,
                )
                if not _is_missing(left) and not _is_missing(right)
            ]
            correlations.append(
                CorrelationPair(
                    left_column=left_name,
                    right_column=right_name,
                    pair_count=len(paired_values),
                    pearson_r=_pearson(
                        [left for left, _ in paired_values],
                        [right for _, right in paired_values],
                    ),
                )
            )

    warnings = (
        "Column meanings and units are unverified; numeric labels are not interpreted.",
        "Pearson correlations are numerical associations only and are not causal.",
        "Tukey 1.5-IQR flags are screening checks, not confirmed errors.",
        "Repeated date values are retained; no records are deduplicated by date.",
    )
    return MarketDataExploration(
        row_count=len(rows),
        column_count=len(headers),
        duplicate_row_count=duplicate_count,
        columns=tuple(columns),
        categorical_summaries=tuple(categorical),
        numeric_summaries=tuple(numeric),
        date_coverage=tuple(date_ranges),
        identifier_candidates=candidates,
        numeric_column_correlations=tuple(correlations),
        warnings=warnings,
    )


def _new_figure(*, figsize: tuple[float, float]) -> Figure:
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.figure import Figure as MatplotlibFigure

    figure = MatplotlibFigure(figsize=figsize, layout="constrained")
    FigureCanvasAgg(figure)
    return figure


def _save_figure(figure: Figure, path: Path) -> None:
    if path.exists():
        raise FileExistsError(f"refusing to overwrite existing figure: {path}")
    figure.savefig(path, dpi=140, format="png", metadata={"Software": "Newton Lab"})
    figure.clear()


def create_exploration_plots(
    profile: MarketDataExploration,
    output_directory: str | Path,
    *,
    source_path: str | Path,
) -> tuple[Path, ...]:
    """Save neutral cross-sectional plots using optional Matplotlib."""
    output_path = Path(output_directory)
    output_path.mkdir(parents=True, exist_ok=True)
    saved: list[Path] = []

    missing_profiles = [item for item in profile.columns if item.missing_count]
    if missing_profiles:
        figure = _new_figure(figsize=(max(9.0, len(missing_profiles) * 0.42), 5.2))
        axis = figure.add_subplot(1, 1, 1)
        axis.bar(
            [item.name for item in missing_profiles],
            [item.missing_count for item in missing_profiles],
            color="#4c78a8",
        )
        axis.set_title("Missing CSV cells by exact column name")
        axis.set_xlabel("CSV column")
        axis.set_ylabel("Missing cell count")
        axis.tick_params(axis="x", rotation=70)
        axis.grid(True, axis="y", alpha=0.25)
        target = output_path / "missingness_by_column.png"
        _save_figure(figure, target)
        saved.append(target)

    decimal_names = [
        item.name
        for item in profile.numeric_summaries
        if item.inferred_type == "decimal_like"
    ]
    if decimal_names:
        figure = _new_figure(figsize=(14.0, 7.0))
        axes = figure.subplots(2, 4).ravel()
        decimal_values = _read_decimal_columns(source_path, decimal_names)
        for axis, name in zip(axes, decimal_names, strict=False):
            axis.hist(
                decimal_values[name], bins="auto", color="#4c78a8", edgecolor="white"
            )
            axis.set_title(name)
            axis.set_xlabel(f"{name} (unit unknown)")
            axis.set_ylabel("Row count")
            axis.text(
                0.02,
                0.98,
                "Observed values; not a fitted distribution",
                transform=axis.transAxes,
                va="top",
                fontsize=7,
            )
            axis.grid(True, axis="y", alpha=0.2)
        for axis in axes[len(decimal_names) :]:
            axis.set_visible(False)
        figure.suptitle("Cross-sectional distributions of decimal-like CSV columns")
        target = output_path / "decimal_column_distributions.png"
        _save_figure(figure, target)
        saved.append(target)

    series_summary = next(
        (item for item in profile.categorical_summaries if item.name == "SctySrs"),
        None,
    )
    if series_summary is not None and series_summary.value_counts:
        levels = series_summary.value_counts[:20]
        figure = _new_figure(figsize=(10.0, max(4.0, len(levels) * 0.25)))
        axis = figure.add_subplot(1, 1, 1)
        axis.barh(
            [item.value for item in reversed(levels)],
            [item.count for item in reversed(levels)],
            color="#59a14f",
        )
        axis.set_title("Most frequent raw values in SctySrs (codes unverified)")
        axis.set_xlabel("Row count")
        axis.set_ylabel("SctySrs value")
        axis.grid(True, axis="x", alpha=0.25)
        target = output_path / "sctysrs_code_counts.png"
        _save_figure(figure, target)
        saved.append(target)

    correlation_names = tuple(
        dict.fromkeys(
            name
            for pair in profile.numeric_column_correlations
            for name in (pair.left_column, pair.right_column)
        )
    )
    if correlation_names:
        correlation_values = {
            (pair.left_column, pair.right_column): pair.pearson_r
            for pair in profile.numeric_column_correlations
        }
        figure = _new_figure(figsize=(max(8.0, len(correlation_names) * 0.9), 7.0))
        axis = figure.add_subplot(1, 1, 1)
        matrix: list[list[float]] = []
        for left in correlation_names:
            matrix_row: list[float] = []
            for right in correlation_names:
                if left == right:
                    matrix_row.append(1.0)
                else:
                    key = (
                        (left, right)
                        if (left, right) in correlation_values
                        else (right, left)
                    )
                    value = correlation_values.get(key)
                    matrix_row.append(float("nan") if value is None else value)
            matrix.append(matrix_row)
        image = axis.imshow(matrix, vmin=-1, vmax=1, cmap="coolwarm")
        axis.set_xticks(range(len(correlation_names)), correlation_names, rotation=70)
        axis.set_yticks(range(len(correlation_names)), correlation_names)
        axis.set_title("Pearson r among numeric-like columns (association only)")
        figure.colorbar(image, ax=axis, label="Pearson r")
        target = output_path / "numeric_column_correlations.png"
        _save_figure(figure, target)
        saved.append(target)

    return tuple(saved)


def _read_decimal_columns(
    source_path: str | Path, selected_columns: Sequence[str]
) -> dict[str, list[float]]:
    """Read selected decimal-like columns for plots without changing the CSV."""
    values: dict[str, list[float]] = {name: [] for name in selected_columns}
    with Path(source_path).open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None or any(
            name not in reader.fieldnames for name in selected_columns
        ):
            raise ScientificValidationError(
                "plot source no longer matches the analyzed CSV schema"
            )
        for row in reader:
            for name in selected_columns:
                raw_value = row.get(name)
                if _is_missing(raw_value):
                    continue
                if not isinstance(raw_value, str):
                    raise ScientificValidationError(
                        f"numeric column {name!r} contains an invalid CSV cell"
                    )
                try:
                    value = float(raw_value)
                except (TypeError, ValueError) as exc:
                    raise ScientificValidationError(
                        f"numeric column {name!r} changed during plot generation"
                    ) from exc
                if not math.isfinite(value):
                    raise ScientificValidationError(
                        f"numeric column {name!r} contains a non-finite value"
                    )
                values[name].append(value)
    return values


def write_analysis_outputs(
    profile: MarketDataExploration,
    *,
    output_directory: str | Path,
    source_path: str | Path,
    command: str,
) -> tuple[Path, tuple[Path, ...]]:
    """Write a deterministic JSON summary and plots without overwriting files."""
    import matplotlib

    output_path = Path(output_directory)
    output_path.mkdir(parents=True, exist_ok=True)
    summary_path = output_path / "analysis_summary.json"
    if summary_path.exists():
        raise FileExistsError(f"refusing to overwrite existing summary: {summary_path}")
    plots = create_exploration_plots(
        profile, output_path / "figures", source_path=source_path
    )
    payload = {
        "source_path": str(Path(source_path)),
        "source_sha256": hashlib.sha256(Path(source_path).read_bytes()).hexdigest(),
        "command": command,
        "python_version": sys.version.split()[0],
        "matplotlib_version": matplotlib.__version__,
        "profile": asdict(profile),
    }
    summary_path.write_text(
        json.dumps(payload, indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )
    return summary_path, plots


def main(argv: Sequence[str] | None = None) -> int:
    """Profile a local CSV and save JSON/PNG outputs to a dedicated directory."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("csv_path", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args(argv)
    try:
        profile = analyze_market_csv(args.csv_path)
        summary_path, plots = write_analysis_outputs(
            profile,
            output_directory=args.output_dir,
            source_path=args.csv_path,
            command=f'"{sys.executable}" -m newton_lab.market_data_exploration '
            f"{args.csv_path} --output-dir {args.output_dir}",
        )
    except ScientificValidationError as exc:
        parser.error(str(exc))
    print(f"Summary: {summary_path}")
    for plot_path in plots:
        print(f"Plot: {plot_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

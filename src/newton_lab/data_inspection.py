"""Read-only inspection of manually obtained local CSV files."""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from collections.abc import Sequence
from dataclasses import asdict, dataclass
from datetime import datetime, timedelta
from pathlib import Path

from newton_lab.exceptions import ScientificValidationError


@dataclass(frozen=True)
class LocalCsvInspection:
    """Summary of a local CSV without retaining or modifying its contents."""

    columns: tuple[str, ...]
    row_count: int
    missing_values_by_column: dict[str, int]
    duplicate_row_count: int
    timestamp_column: str
    invalid_timestamp_count: int
    duplicate_timestamp_count: int
    date_start: str | None
    date_end: str | None
    inferred_frequency: str | None


def inspect_local_csv(path: str | Path, *, timestamp_column: str) -> LocalCsvInspection:
    """Inspect a manually obtained local CSV; this function performs no network I/O.

    Frequency is inferred from sorted, unique, valid timestamps and is only a
    diagnostic. It does not establish the source's sampling convention.
    """
    csv_path = Path(path)
    if csv_path.suffix.lower() != ".csv" or not csv_path.is_file():
        raise ScientificValidationError("path must refer to an existing local CSV file")
    try:
        with csv_path.open(encoding="utf-8-sig", newline="") as handle:
            reader = csv.DictReader(handle)
            if reader.fieldnames is None:
                raise ScientificValidationError("CSV must have a header row")
            columns = tuple(reader.fieldnames)
            if timestamp_column not in columns:
                raise ScientificValidationError(
                    f"timestamp column {timestamp_column!r} is not present in the CSV"
                )
            missing = Counter({column: 0 for column in columns})
            seen_rows: set[tuple[str | None, ...]] = set()
            duplicate_rows = 0
            timestamp_values: list[datetime] = []
            invalid_timestamps = 0
            duplicate_timestamps = 0
            seen_timestamps: set[datetime] = set()
            row_count = 0
            for row in reader:
                row_count += 1
                values = tuple(row.get(column) for column in columns)
                for column, value in zip(columns, values, strict=True):
                    if _is_missing(value):
                        missing[column] += 1
                if values in seen_rows:
                    duplicate_rows += 1
                else:
                    seen_rows.add(values)
                raw_timestamp = row.get(timestamp_column)
                try:
                    timestamp = datetime.fromisoformat(
                        raw_timestamp.strip().replace("Z", "+00:00")
                        if raw_timestamp is not None
                        else ""
                    )
                except ValueError:
                    invalid_timestamps += 1
                    continue
                if timestamp in seen_timestamps:
                    duplicate_timestamps += 1
                else:
                    seen_timestamps.add(timestamp)
                    timestamp_values.append(timestamp)
    except ScientificValidationError:
        raise
    except (OSError, UnicodeError, csv.Error) as exc:
        raise ScientificValidationError(f"could not read local CSV: {exc}") from exc
    unique_sorted = sorted(timestamp_values)
    frequency: str | None = None
    if len(unique_sorted) >= 3:
        differences = {
            second - first
            for first, second in zip(unique_sorted, unique_sorted[1:], strict=False)
        }
        if len(differences) == 1:
            frequency = _format_frequency(differences.pop())
    start = unique_sorted[0].isoformat() if unique_sorted else None
    end = unique_sorted[-1].isoformat() if unique_sorted else None
    return LocalCsvInspection(
        columns=columns,
        row_count=row_count,
        missing_values_by_column=dict(missing),
        duplicate_row_count=duplicate_rows,
        timestamp_column=timestamp_column,
        invalid_timestamp_count=invalid_timestamps,
        duplicate_timestamp_count=duplicate_timestamps,
        date_start=start,
        date_end=end,
        inferred_frequency=frequency,
    )


def _is_missing(value: str | None) -> bool:
    """Recognize empty cells and common text null markers in CSV fields."""
    return value is None or value.strip().casefold() in {"", "na", "n/a", "nan", "null"}


def _format_frequency(delta: timedelta) -> str:
    """Format a constant positive interval using a compact SI-style unit."""
    seconds = delta.total_seconds()
    for scale, unit in ((86_400, "d"), (3_600, "h"), (60, "min"), (1, "s")):
        if seconds >= scale and seconds % scale == 0:
            return f"{int(seconds / scale)}{unit}"
    return f"{seconds:g}s"


def main(argv: Sequence[str] | None = None) -> int:
    """Print an inspection summary as JSON for a local CSV."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("csv_path", type=Path)
    parser.add_argument("--timestamp-column", required=True)
    args = parser.parse_args(argv)
    try:
        report = inspect_local_csv(
            args.csv_path, timestamp_column=args.timestamp_column
        )
    except ScientificValidationError as exc:
        parser.error(str(exc))
    print(json.dumps(asdict(report), indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

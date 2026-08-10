from __future__ import annotations

import warnings
from numbers import Real

import numpy as np


def create_sliding_phase_intervals(
    start_seconds: Real,
    end_seconds: Real,
    window_seconds: Real,
    *,
    overlap_seconds: Real = 0.0,
    label_prefix: object = "window",
    warn_on_trailing: bool = True,
) -> tuple[np.ndarray, list[str]]:
    """Return complete half-open sliding windows expressed in seconds."""

    start = _finite_number(start_seconds, "start_seconds")
    end = _finite_number(end_seconds, "end_seconds")
    window = _finite_number(window_seconds, "window_seconds")
    overlap = _finite_number(overlap_seconds, "overlap_seconds")
    if start < 0:
        raise ValueError("start_seconds must be nonnegative")
    if end <= start:
        raise ValueError("end_seconds must be greater than start_seconds")
    if window <= 0:
        raise ValueError("window_seconds must be positive")
    if overlap < 0:
        raise ValueError("overlap_seconds must be nonnegative")
    if overlap >= window:
        raise ValueError("overlap_seconds must be smaller than window_seconds")
    if not isinstance(warn_on_trailing, bool):
        raise TypeError("warn_on_trailing must be a boolean")

    prefix = str(label_prefix)
    if prefix == "":
        raise ValueError("label_prefix must convert to a nonempty string")

    step = window - overlap
    available = end - start
    tolerance = np.finfo(float).eps * max(1.0, abs(start), abs(end), abs(window)) * 8
    n_windows = int(np.floor((available - window + tolerance) / step)) + 1
    if n_windows < 1:
        raise ValueError("no complete window fits within the requested interval")

    indexes = np.arange(n_windows, dtype=float)
    starts = start + indexes * step
    stops = starts + window
    intervals = np.vstack((starts, stops)).astype(float, copy=False)
    trailing = end - float(stops[-1])
    if trailing > tolerance and warn_on_trailing:
        warnings.warn(
            f"dropped {trailing:g} trailing second(s) that do not fill a complete window",
            UserWarning,
            stacklevel=2,
        )
    labels = [f"{prefix}_{index + 1:03d}" for index in range(n_windows)]
    return intervals, labels


def _finite_number(value: Real, name: str) -> float:
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, Real):
        raise TypeError(f"{name} must be numeric")
    normalized = float(value)
    if not np.isfinite(normalized):
        raise ValueError(f"{name} must be finite")
    return normalized

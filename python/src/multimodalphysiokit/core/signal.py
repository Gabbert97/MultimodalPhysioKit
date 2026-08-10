from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass
from typing import Any

import numpy as np


@dataclass(frozen=True)
class Signal:
    """Sampled physiological signal aligned with a time vector.

    This class mirrors the terminology of MATLAB ``mpk.Signal`` while using
    Python-native conventions. Samples and time are stored as one-dimensional
    floating-point NumPy arrays.
    """

    samples: Any
    time: Any
    sampling_frequency: float
    label: str = ""
    units: str = ""
    metadata: dict[str, Any] | None = None

    def __post_init__(self) -> None:
        samples = np.asarray(self.samples, dtype=float).reshape(-1)
        time = np.asarray(self.time, dtype=float).reshape(-1)

        if samples.size == 0:
            raise ValueError("samples must be non-empty")
        if time.size == 0:
            raise ValueError("time must be non-empty")
        if samples.shape[0] != time.shape[0]:
            raise ValueError("samples and time must have the same length")
        if not np.all(np.isfinite(samples)):
            raise ValueError("samples must contain only finite values")
        if not np.all(np.isfinite(time)):
            raise ValueError("time must contain only finite values")
        if np.any(np.diff(time) <= 0):
            raise ValueError("time must be monotonically increasing")

        sampling_frequency = float(self.sampling_frequency)
        if not np.isfinite(sampling_frequency) or sampling_frequency <= 0:
            raise ValueError("sampling_frequency must be positive and finite")

        metadata = {} if self.metadata is None else dict(self.metadata)

        object.__setattr__(self, "samples", samples)
        object.__setattr__(self, "time", time)
        object.__setattr__(self, "sampling_frequency", sampling_frequency)
        object.__setattr__(self, "label", str(self.label))
        object.__setattr__(self, "units", str(self.units))
        object.__setattr__(self, "metadata", metadata)

    @property
    def n_samples(self) -> int:
        """Number of samples in the signal."""

        return int(self.samples.shape[0])

    @property
    def duration_seconds(self) -> float:
        """Elapsed time between the first and last sample."""

        return float(self.time[-1] - self.time[0])

    def copy_with_metadata(self, **metadata_updates: Any) -> "Signal":
        """Return a copy with updated non-sensitive metadata."""

        metadata = dict(self.metadata)
        metadata.update(metadata_updates)
        return Signal(
            samples=self.samples.copy(),
            time=self.time.copy(),
            sampling_frequency=self.sampling_frequency,
            label=self.label,
            units=self.units,
            metadata=metadata,
        )

    def remove_time_ranges(
        self,
        intervals_seconds: Iterable[tuple[float, float | None]],
        *,
        reset_time: bool = True,
    ) -> "Signal":
        """Return a copy with samples in half-open time intervals removed."""

        intervals = _validate_time_ranges(intervals_seconds)
        if not isinstance(reset_time, bool):
            raise TypeError("reset_time must be a boolean")

        keep = np.ones(self.n_samples, dtype=bool)
        for start_seconds, end_seconds in intervals:
            if end_seconds is None:
                remove = self.time >= start_seconds
            else:
                remove = (self.time >= start_seconds) & (self.time < end_seconds)
            keep &= ~remove
        if not np.any(keep):
            raise ValueError("removing time ranges would remove all samples")

        samples = self.samples[keep].copy()
        time = (
            np.arange(samples.size, dtype=float) / self.sampling_frequency
            if reset_time
            else self.time[keep].copy()
        )
        metadata = dict(self.metadata)
        metadata.update(
            {
                "time_ranges_removed_seconds": [
                    (start_seconds, end_seconds) for start_seconds, end_seconds in intervals
                ],
                "reset_time_after_removal": reset_time,
            }
        )
        return Signal(
            samples=samples,
            time=time,
            sampling_frequency=self.sampling_frequency,
            label=self.label,
            units=self.units,
            metadata=metadata,
        )


def _validate_time_ranges(intervals_seconds: Iterable[tuple[float, float | None]]) -> list[tuple[float, float | None]]:
    if isinstance(intervals_seconds, (str, bytes)):
        raise TypeError("intervals_seconds must be an iterable of (start, end) pairs")
    try:
        raw_intervals = list(intervals_seconds)
    except TypeError as exc:
        raise TypeError("intervals_seconds must be an iterable of (start, end) pairs") from exc
    if not raw_intervals:
        raise ValueError("intervals_seconds must be non-empty")

    intervals: list[tuple[float, float | None]] = []
    for interval in raw_intervals:
        if isinstance(interval, (str, bytes)):
            raise TypeError("each time interval must be a (start, end) pair")
        try:
            start_seconds, end_seconds = interval
        except (TypeError, ValueError) as exc:
            raise ValueError("each time interval must contain exactly two values") from exc

        start_value = float(start_seconds)
        if not np.isfinite(start_value) or start_value < 0:
            raise ValueError("interval start times must be finite and nonnegative")

        if end_seconds is None:
            end_value = None
        else:
            end_value = float(end_seconds)
            if not np.isfinite(end_value) or end_value <= 0:
                raise ValueError("interval end times must be None or finite positive values")
            if end_value <= start_value:
                raise ValueError("interval end times must be greater than start times")
        intervals.append((start_value, end_value))
    return intervals

from __future__ import annotations

from typing import Any

import numpy as np

from multimodalphysiokit.core import Recording, Signal


def recording_from_arrays(
    signals: dict[str, Any],
    *,
    sampling_frequency: float | None = None,
    time: Any | None = None,
    units: dict[str, str] | None = None,
    labels: dict[str, str] | None = None,
    name: str = "recording",
    metadata: dict[str, Any] | None = None,
) -> Recording:
    """Create a ``Recording`` from in-memory one-dimensional numeric arrays."""

    if not isinstance(signals, dict):
        raise TypeError("signals must be a dictionary")
    if not signals:
        raise ValueError("signals must be non-empty")

    has_sampling_frequency = sampling_frequency is not None
    has_time = time is not None
    if has_sampling_frequency == has_time:
        raise ValueError("exactly one of sampling_frequency or time must be provided")

    units = {} if units is None else dict(units)
    labels = {} if labels is None else dict(labels)

    sample_arrays = {
        signal_name: _validate_signal_samples(signal_name, samples)
        for signal_name, samples in signals.items()
    }

    if has_sampling_frequency:
        sampling_frequency = _validate_sampling_frequency(sampling_frequency)
        time_vectors = {
            signal_name: np.arange(samples.size, dtype=float) / sampling_frequency
            for signal_name, samples in sample_arrays.items()
        }
    else:
        explicit_time = _validate_time(time)
        time_vectors = {}
        for signal_name, samples in sample_arrays.items():
            if samples.size != explicit_time.size:
                raise ValueError(
                    f"signal '{signal_name}' length must match explicit time length"
                )
            time_vectors[signal_name] = explicit_time.copy()
        sampling_frequency = _sampling_frequency_from_time(explicit_time)

    recording = Recording(name=name, metadata=metadata)
    for signal_name, samples in sample_arrays.items():
        signal = Signal(
            samples=samples,
            time=time_vectors[signal_name],
            sampling_frequency=sampling_frequency,
            label=str(labels.get(signal_name, signal_name)),
            units=str(units.get(signal_name, "")),
            metadata={},
        )
        recording.add_signal(signal_name, signal)

    return recording


def _validate_signal_samples(signal_name: str, samples: Any) -> np.ndarray:
    if not isinstance(signal_name, str) or not signal_name:
        raise ValueError("signal names must be non-empty strings")

    array = np.asarray(samples)
    if array.ndim != 1:
        raise ValueError(f"signal '{signal_name}' samples must be one-dimensional")
    if array.size == 0:
        raise ValueError(f"signal '{signal_name}' samples must be non-empty")
    if not np.issubdtype(array.dtype, np.number):
        raise ValueError(f"signal '{signal_name}' samples must be numeric")

    array = array.astype(float)
    if not np.all(np.isfinite(array)):
        raise ValueError(f"signal '{signal_name}' samples must contain only finite values")
    return array


def _validate_sampling_frequency(sampling_frequency: float) -> float:
    value = float(sampling_frequency)
    if not np.isfinite(value) or value <= 0:
        raise ValueError("sampling_frequency must be positive and finite")
    return value


def _validate_time(time: Any) -> np.ndarray:
    array = np.asarray(time)
    if array.ndim != 1:
        raise ValueError("time must be one-dimensional")
    if array.size == 0:
        raise ValueError("time must be non-empty")
    if not np.issubdtype(array.dtype, np.number):
        raise ValueError("time must be numeric")

    array = array.astype(float)
    if not np.all(np.isfinite(array)):
        raise ValueError("time must contain only finite values")
    if np.any(np.diff(array) <= 0):
        raise ValueError("time must be monotonically increasing")
    return array


def _sampling_frequency_from_time(time: np.ndarray) -> float:
    if time.size == 1:
        return 1.0
    median_step = float(np.median(np.diff(time)))
    return 1.0 / median_step

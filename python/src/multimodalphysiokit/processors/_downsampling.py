from __future__ import annotations

import numpy as np
from scipy.signal import butter, sosfiltfilt

_ANTIALIAS_ORDER = 4
_ANTIALIAS_CUTOFF_FRACTION = 0.4


def antialias_downsample(
    samples: np.ndarray,
    factor: int,
    target_sampling_frequency: float,
) -> np.ndarray:
    """Low-pass filter, then keep every ``factor``-th sample.

    A zero-phase 4th-order Butterworth low-pass at 0.4 x the target sampling
    frequency removes content above the new Nyquist frequency before
    decimation, so it cannot alias into the analysis band. Sample selection
    (``samples[::factor]``) is unchanged, so timing is preserved. With
    ``factor == 1`` the samples are returned unfiltered.
    """

    values = np.asarray(samples, dtype=float)
    if factor == 1:
        return values.copy()
    original_sampling_frequency = target_sampling_frequency * factor
    sos = butter(
        _ANTIALIAS_ORDER,
        _ANTIALIAS_CUTOFF_FRACTION * target_sampling_frequency,
        btype="low",
        fs=original_sampling_frequency,
        output="sos",
    )
    try:
        filtered = sosfiltfilt(sos, values)
    except ValueError as exc:
        raise ValueError("signal is too short for anti-alias filtering") from exc
    return filtered[::factor]

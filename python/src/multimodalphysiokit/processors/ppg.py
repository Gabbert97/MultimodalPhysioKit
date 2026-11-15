from __future__ import annotations

from numbers import Real
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from numpy.lib.stride_tricks import sliding_window_view
from scipy.signal import butter, filtfilt, find_peaks

from multimodalphysiokit.core import ProcessingResult, Signal

from .ecg import (
    _ecg_time_domain_feature_values,
    _extract_ecg_frequency_domain_features,
    _plot_ecg_bpm_trend,
)

# Cleaning band: 0.5 Hz removes baseline drift, 8 Hz removes noise above the
# pulse harmonics. All pulse characteristics are measured on this signal.
_CLEAN_BAND_HZ = (0.5, 8.0)
# Beat-detection band: 0.5-3.0 Hz covers 30-180 bpm and removes breathing drift.
_BAND_HZ = (0.5, 3.0)
_FILTER_ORDER = 3
_MAX_BPM = 180.0
_MIN_BPM = 40.0
# Foot search window around each filtered trough, as fractions of the median
# trough distance: 10 % before and 30 % after.
_FOOT_SEARCH_BEFORE = 0.1
_FOOT_SEARCH_AFTER = 0.3
# A band-pass trough is a foot only if its prominence is at least 30 % of the
# 75th percentile of the prominences of up to 11 neighbouring troughs; the
# shallower dicrotic notch therefore does not split a pulse.
_FOOT_MIN_RELATIVE_PROMINENCE = 0.3
_FOOT_PROMINENCE_PERCENTILE = 75.0
# Systolic peak: the first local maximum of the pulse whose height above the
# foot-to-foot line is at least 70 % of the highest local maximum.
_PEAK_MIN_RELATIVE_HEIGHT = 0.7
# Pulse quality: duration within 0.8-1.2 x the running median duration of up
# to 11 neighbouring pulses, and systolic peak within the first
# max_peak_position fraction of the pulse (default: first half).
_DURATION_TOLERANCE = (0.8, 1.2)
_DEFAULT_MAX_PEAK_POSITION = 0.5
_RUNNING_MEDIAN_HALF_WIDTH = 5
_MIN_FREQUENCY_SECONDS = 60.0

_PULSE_FEATURE_NAMES = ("ppg_pulse_amplitude", "ppg_pulse_duration", "ppg_rise_time")


class PPGProcessor:
    """Detect PPG pulses and extract pulse-shape and heart-rate features.

    The supplied ``Signal`` is analysed as one piece, whatever its length
    (a recording, a phase, or a completed real-time window); all thresholds
    adapt to the pulses found in it.

    ``label_frequency_analysis`` selects HRV frequency analysis: ``1`` STFT,
    ``2`` PSD (both need at least 60 s), or ``None`` to skip it so short
    windows can be processed.

    ``max_peak_position`` is the latest accepted systolic-peak position as a
    fraction of the pulse (default 0.5, the first half). Rounded wrist PPG
    pulses may need a larger value, e.g. 0.65.
    """

    def __init__(
        self,
        *,
        label_frequency_analysis: int | None = 1,
        return_intermediates: bool = False,
        max_peak_position: float = _DEFAULT_MAX_PEAK_POSITION,
    ) -> None:
        if label_frequency_analysis not in (None, 1, 2):
            raise ValueError("label_frequency_analysis must be None, 1, or 2")
        if not isinstance(return_intermediates, bool):
            raise TypeError("return_intermediates must be a boolean")
        self.label_frequency_analysis = label_frequency_analysis
        self.return_intermediates = return_intermediates
        self.max_peak_position = _validate_max_peak_position(max_peak_position)

    def process(
        self,
        signal: Signal,
        *,
        show_plots: bool = False,
        save_plots: bool = False,
        output_dir: str | Path | None = None,
        output_label: str = "ppg",
    ) -> ProcessingResult:
        """Process exactly one supplied PPG signal."""

        _validate_inputs(signal, show_plots, save_plots, output_dir, output_label)
        if self.label_frequency_analysis is not None and (
            signal.n_samples / signal.sampling_frequency < _MIN_FREQUENCY_SECONDS
        ):
            raise ValueError(
                "PPG frequency analysis requires at least 60 seconds; "
                "use label_frequency_analysis=None for shorter signals"
            )

        pulses = detect_ppg_pulses(signal, max_peak_position=self.max_peak_position)
        ibi = _accepted_ibis(signal.time, pulses)
        if np.count_nonzero(pulses["valid"]) < 3 or ibi["values"].size < 3:
            raise ValueError("too few valid PPG pulses were detected")

        features = _pulse_features(pulses)
        time_domain = _ecg_time_domain_feature_values(
            60 / ibi["values"], ibi["values"], ibi["flags"]
        )
        features.update({_ppg_name(name): value for name, value in time_domain.items()})

        frequency_data: dict[str, object] = {}
        if self.label_frequency_analysis is not None:
            frequency_features, frequency_data = _extract_ecg_frequency_domain_features(
                ibi["values"],
                label_frequency_analysis=self.label_frequency_analysis,
                stft_context=_stft_context(signal, ibi)
                if self.label_frequency_analysis == 1
                else None,
                ibi_start_times=ibi["start_times"],
                show_plots=show_plots,
                save_plots=save_plots,
                output_dir=output_dir,
                output_label=output_label,
            )
            features.update(
                {_ppg_name(name): value for name, value in frequency_features.items()}
            )

        if show_plots or save_plots:
            _plot_ppg_pulses(
                signal, pulses, show_plots=show_plots, save_plots=save_plots,
                output_dir=output_dir, output_label=output_label,
            )
            _plot_ecg_bpm_trend(
                60 / ibi["values"], show_plots=show_plots, save_plots=save_plots,
                output_dir=output_dir, output_label=output_label,
            )

        intermediates: dict[str, object] = {}
        if self.return_intermediates:
            intermediates = {
                "cleaned_signal": pulses["cleaned_signal"],
                "filtered_signal": pulses["filtered_signal"],
                "pulse_onset_indices": pulses["onsets"],
                "pulse_peak_indices": pulses["peaks"],
                "pulse_peak_times": signal.time[pulses["peaks"]],
                "pulse_valid": pulses["valid"],
                "pulse_amplitudes": pulses["amplitudes"],
                "pulse_durations": pulses["durations"],
                "pulse_rise_times": pulses["rise_times"],
                "ibi": ibi["values"],
                "ibi_times": ibi["start_times"],
                "bpm": 60 / ibi["values"],
            }
            if self.label_frequency_analysis == 1:
                intermediates.update(
                    {
                        "stft_power": frequency_data["power"],
                        "stft_frequencies": frequency_data["frequencies"],
                        "stft_times": frequency_data["time"],
                    }
                )

        return ProcessingResult(
            modality="ppg",
            processor_name="process_ppg",
            features=features,
            signal_labels=[signal.label or "ppg"],
            metadata={
                "source": "process_ppg",
                "sampling_frequency": signal.sampling_frequency,
                "n_samples": signal.n_samples,
                "start_time_seconds": float(signal.time[0]),
                "end_time_seconds": float(signal.time[-1]),
                "frequency_analysis_method": self.label_frequency_analysis,
                "max_peak_position": self.max_peak_position,
                "n_pulses": int(pulses["valid"].size),
                "n_valid_pulses": int(np.count_nonzero(pulses["valid"])),
                "n_accepted_ibis": int(ibi["values"].size),
            },
            intermediates=intermediates,
        )


def process_ppg(
    signal: Signal,
    *,
    label_frequency_analysis: int | None = 1,
    return_intermediates: bool = False,
    max_peak_position: float = _DEFAULT_MAX_PEAK_POSITION,
    show_plots: bool = False,
    save_plots: bool = False,
    output_dir: str | Path | None = None,
    output_label: str = "ppg",
) -> ProcessingResult:
    """Functional wrapper for PPGProcessor.process."""

    return PPGProcessor(
        label_frequency_analysis=label_frequency_analysis,
        return_intermediates=return_intermediates,
        max_peak_position=max_peak_position,
    ).process(
        signal,
        show_plots=show_plots,
        save_plots=save_plots,
        output_dir=output_dir,
        output_label=output_label,
    )


def detect_ppg_pulses(
    signal: Signal,
    *,
    max_peak_position: float = _DEFAULT_MAX_PEAK_POSITION,
) -> dict[str, np.ndarray]:
    """Return pulse onsets (feet), systolic peaks, per-pulse features and validity.

    Beats are located on the 0.5-3.0 Hz band-passed signal, keeping only
    troughs deep enough not to be a dicrotic notch; all pulse characteristics
    (feet, peaks, amplitude, duration, rise time) are measured on the cleaned
    signal, i.e. the PPG band-passed at 0.5-8 Hz, which removes baseline drift
    and noise but keeps the pulse morphology. The systolic peak is the first
    local maximum reaching 70 % of the pulse's highest one. IBIs are measured
    between systolic peaks.

    Pulse ``k`` spans ``onsets[k]`` to ``onsets[k + 1]``; all index arrays are
    zero-based sample indexes into ``signal``.
    """

    fs = signal.sampling_frequency
    cleaned, filtered = _clean_and_bandpass(signal.samples, fs)

    troughs, properties = find_peaks(
        -filtered, distance=max(1, _matlab_round(fs * 60 / _MAX_BPM)), prominence=0
    )
    prominences = properties["prominences"]
    local_prominence = _centered_window_stat(
        prominences,
        _RUNNING_MEDIAN_HALF_WIDTH,
        lambda values, axis=None: np.percentile(values, _FOOT_PROMINENCE_PERCENTILE, axis=axis),
    )
    troughs = troughs[prominences >= _FOOT_MIN_RELATIVE_PROMINENCE * local_prominence]
    if troughs.size < 3:
        raise ValueError("too few PPG pulses were detected")
    onsets = _refine_feet(cleaned, troughs)
    if onsets.size < 3:
        raise ValueError("too few PPG pulses were detected")

    starts, stops = onsets[:-1], onsets[1:]
    peaks = _first_systolic_peaks(cleaned, starts, stops)
    max_peak_position = _validate_max_peak_position(max_peak_position)
    time = signal.time
    durations = time[stops] - time[starts]
    rise_times = time[peaks] - time[starts]
    baseline_at_peak = cleaned[starts] + (cleaned[stops] - cleaned[starts]) * rise_times / durations
    amplitudes = cleaned[peaks] - baseline_at_peak

    reference = _running_median(durations, _RUNNING_MEDIAN_HALF_WIDTH)
    valid = (
        (durations >= _DURATION_TOLERANCE[0] * reference)
        & (durations <= _DURATION_TOLERANCE[1] * reference)
        & (peaks > starts)
        & (peaks < stops)
        & (peaks - starts <= max_peak_position * (stops - starts))
    )
    return {
        "cleaned_signal": cleaned,
        "filtered_signal": filtered,
        "onsets": onsets,
        "peaks": peaks,
        "amplitudes": amplitudes,
        "durations": durations,
        "rise_times": rise_times,
        "valid": valid,
    }


def _clean_and_bandpass(samples: np.ndarray, fs: float) -> tuple[np.ndarray, np.ndarray]:
    """Return (cleaned, band-passed): 0.5 Hz high-pass, then 8.0 Hz / 3.0 Hz low-pass of it."""

    if _CLEAN_BAND_HZ[1] >= fs / 2:
        raise ValueError("sampling_frequency must exceed 16 Hz for 0.5-8.0 Hz PPG filtering")
    try:
        # padlen matches MATLAB filtfilt (3 * filter order) for cross-language parity.
        b, a = butter(_FILTER_ORDER, _CLEAN_BAND_HZ[0], btype="high", fs=fs)
        high_passed = filtfilt(b, a, samples, padlen=3 * _FILTER_ORDER)
        b, a = butter(_FILTER_ORDER, _CLEAN_BAND_HZ[1], btype="low", fs=fs)
        cleaned = filtfilt(b, a, high_passed, padlen=3 * _FILTER_ORDER)
        b, a = butter(_FILTER_ORDER, _BAND_HZ[1], btype="low", fs=fs)
        return cleaned, filtfilt(b, a, high_passed, padlen=3 * _FILTER_ORDER)
    except ValueError as exc:
        raise ValueError("PPG signal is too short for filtering") from exc


def _first_systolic_peaks(cleaned: np.ndarray, starts: np.ndarray, stops: np.ndarray) -> np.ndarray:
    """First local maximum of each pulse reaching 70 % of its highest one.

    Heights are measured above the foot-to-foot line. A local maximum is a
    sample higher than the previous one and not lower than the next one; a
    pulse without one falls back to its global maximum.
    """

    is_max = np.zeros(cleaned.size, dtype=bool)
    is_max[1:-1] = (cleaned[1:-1] > cleaned[:-2]) & (cleaned[1:-1] >= cleaned[2:])
    candidates = np.flatnonzero(is_max)
    first = np.searchsorted(candidates, starts, side="right")
    last = np.searchsorted(candidates, stops, side="left")
    peaks = np.empty(starts.size, dtype=int)
    for k, (start, stop) in enumerate(zip(starts, stops)):
        inside = candidates[first[k] : last[k]]
        if inside.size:
            baseline = cleaned[start] + (cleaned[stop] - cleaned[start]) * ((inside - start) / (stop - start))
            height = cleaned[inside] - baseline
            if height.max() > 0:
                peaks[k] = inside[np.argmax(height >= _PEAK_MIN_RELATIVE_HEIGHT * height.max())]
                continue
        peaks[k] = start + int(np.argmax(cleaned[start : stop + 1]))
    return peaks


def _refine_feet(cleaned: np.ndarray, troughs: np.ndarray) -> np.ndarray:
    """Move each filtered trough to the cleaned-signal minimum near it."""

    median_distance = float(np.median(np.diff(troughs)))
    before = _matlab_round(_FOOT_SEARCH_BEFORE * median_distance)
    after = _matlab_round(_FOOT_SEARCH_AFTER * median_distance)
    feet = []
    for trough in troughs:
        lo = max(0, trough - before)
        hi = min(cleaned.size - 1, trough + after)
        feet.append(lo + int(np.argmin(cleaned[lo : hi + 1])))
    return np.unique(feet)


def _matlab_round(value: float) -> int:
    """Round half away from zero, as MATLAB does (values here are positive)."""

    return int(np.floor(value + 0.5))


def _running_median(values: np.ndarray, half_width: int) -> np.ndarray:
    return _centered_window_stat(values, half_width, np.median)


def _centered_window_stat(values: np.ndarray, half_width: int, stat) -> np.ndarray:
    """``stat`` of each value's centred window, shrunk at the edges.

    Equivalent to ``stat(values[max(0, k - h) : k + h + 1])`` for every ``k``;
    full windows are computed in one vectorised call.
    """

    values = np.asarray(values, dtype=float)
    n, width = values.size, 2 * half_width + 1
    result = np.empty(n)
    if n >= width:
        result[half_width : n - half_width] = stat(sliding_window_view(values, width), axis=1)
        edges = [*range(half_width), *range(n - half_width, n)]
    else:
        edges = range(n)
    for k in edges:
        result[k] = stat(values[max(0, k - half_width) : k + half_width + 1])
    return result


def _accepted_ibis(time: np.ndarray, pulses: dict[str, np.ndarray]) -> dict[str, np.ndarray]:
    """Peak-to-peak IBIs between consecutive valid pulses with plausible rate.

    ``flags[k]`` refers to the interval between pulse ``k`` and ``k + 1``;
    rejected intervals are skipped, never imputed.
    """

    beat_times = time[pulses["peaks"]]
    intervals = np.diff(beat_times)
    valid = pulses["valid"]
    with np.errstate(divide="ignore"):
        bpm = 60 / intervals
    flags = (valid[:-1] & valid[1:] & (bpm >= _MIN_BPM) & (bpm <= _MAX_BPM)).astype(int)
    accepted = flags != 0
    return {
        "values": intervals[accepted],
        "start_times": beat_times[:-1][accepted],
        "end_times": beat_times[1:][accepted],
        "flags": flags,
    }


def _stft_context(signal: Signal, ibi: dict[str, np.ndarray]) -> dict[str, np.ndarray]:
    """Tachogram context for STFT mode: IBIs at their midpoints over the signal."""

    time_range = np.asarray([signal.time[0], signal.time[-1]], dtype=float)
    return {
        "ibi": ibi["values"],
        "ibi_time": (ibi["start_times"] + ibi["end_times"]) / 2,
        "time_range": time_range,
        "analysis_window": time_range,
    }


def _pulse_features(pulses: dict[str, np.ndarray]) -> dict[str, float]:
    valid = pulses["valid"]
    features: dict[str, float] = {}
    for name, key in zip(_PULSE_FEATURE_NAMES, ("amplitudes", "durations", "rise_times")):
        values = pulses[key][valid]
        features[f"{name}_mean"] = float(np.mean(values))
        features[f"{name}_std"] = float(np.std(values, ddof=1))
        features[f"{name}_median"] = float(np.median(values))
    return features


def _validate_max_peak_position(value) -> float:
    if isinstance(value, bool) or not isinstance(value, Real):
        raise TypeError("max_peak_position must be a number")
    value = float(value)
    if not np.isfinite(value) or not 0 < value <= 1:
        raise ValueError("max_peak_position must be in (0, 1]")
    return value


def _ppg_name(ecg_feature_name: str) -> str:
    return "ppg_" + ecg_feature_name.removeprefix("ecg_")


def _validate_inputs(signal, show_plots, save_plots, output_dir, output_label) -> None:
    if not isinstance(signal, Signal):
        raise TypeError("signal must be a Signal instance")
    if not isinstance(show_plots, bool):
        raise TypeError("show_plots must be a boolean")
    if not isinstance(save_plots, bool):
        raise TypeError("save_plots must be a boolean")
    if not isinstance(output_label, str):
        raise TypeError("output_label must be a string")
    if save_plots and output_dir is None:
        raise ValueError("output_dir must be provided when save_plots=True")


def _plot_ppg_pulses(
    signal: Signal,
    pulses: dict[str, np.ndarray],
    *,
    show_plots: bool,
    save_plots: bool,
    output_dir: str | Path | None,
    output_label: str,
) -> None:
    """Raw PPG, cleaned PPG with accepted (blue) / rejected (red) pulses, and band-pass."""

    time, cleaned = signal.time, pulses["cleaned_signal"]
    onsets, peaks, valid = pulses["onsets"], pulses["peaks"], pulses["valid"]
    figure, axes = plt.subplots(3, 1, sharex=True, figsize=(12, 8))
    axes[0].plot(time, signal.samples, color="0.4", linewidth=1)
    axes[0].set_ylabel("Raw PPG")
    axes[0].set_title(f"PPG pulses: {int(valid.sum())} accepted, {int((~valid).sum())} rejected")
    for k, (start, stop) in enumerate(zip(onsets[:-1], onsets[1:])):
        axes[1].plot(time[start : stop + 1], cleaned[start : stop + 1],
                     color="tab:blue" if valid[k] else "tab:red", linewidth=1.2)
    axes[1].plot(time[onsets], cleaned[onsets], "v", color="k", label="foot")
    axes[1].plot(time[peaks[valid]], cleaned[peaks[valid]], "o", color="tab:blue", label="accepted peak")
    axes[1].plot(time[peaks[~valid]], cleaned[peaks[~valid]], "x", color="tab:red", label="rejected pulse")
    axes[1].set_ylabel("Cleaned PPG (0.5-8 Hz)")
    axes[1].legend(loc="upper right")
    axes[2].plot(time, pulses["filtered_signal"], linewidth=1)
    axes[2].set_ylabel("Band-pass 0.5-3.0 Hz")
    axes[2].set_xlabel("Time (s)")
    for axis in axes:
        axis.grid(True)

    if save_plots:
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        figure.savefig(output_path / f"{output_label}_pulses.png", dpi=150)
    if show_plots:
        plt.show()
    else:
        plt.close(figure)

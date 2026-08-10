from __future__ import annotations

from collections.abc import Mapping
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import butter, filtfilt

from multimodalphysiokit.core import ProcessingResult, Signal


_FS_NEW = 100
_FEATURE_NAMES = (
    "respiration_breath_rate_mean",
    "respiration_breath_rate_std",
    "respiration_inspiratory_time_mean",
    "respiration_inspiratory_time_std",
    "respiration_expiratory_time_mean",
    "respiration_expiratory_time_std",
    "respiration_timing_ratio_mean",
    "respiration_timing_ratio_std",
    "respiration_amplitude_mean",
    "respiration_amplitude_std",
    "respiration_minute_ventilation_mean",
    "respiration_minute_ventilation_std",
    "respiration_coefficient_of_variation",
    "respiration_autocorrelation",
)


class RespirationProcessor:
    """Thin object-oriented wrapper for respiration processing."""

    def __init__(self) -> None:
        pass

    def extract_features(
        self,
        signal: Signal,
        *,
        show_plots: bool = False,
        save_plots: bool = False,
        output_dir: str | Path | None = None,
    ) -> dict[str, float]:
        """Return respiration features for the supplied signal."""

        return self.process(
            signal,
            show_plots=show_plots,
            save_plots=save_plots,
            output_dir=output_dir,
        ).features

    def process(
        self,
        signal: Signal,
        *,
        show_plots: bool = False,
        save_plots: bool = False,
        output_dir: str | Path | None = None,
    ) -> ProcessingResult:
        """Process exactly one supplied respiration signal."""

        return process_respiration(
            signal,
            show_plots=show_plots,
            save_plots=save_plots,
            output_dir=output_dir,
        )

def process_respiration(
    signal: Signal,
    *,
    show_plots: bool = False,
    save_plots: bool = False,
    output_dir: str | Path | None = None,
) -> ProcessingResult:
    """Extract MATLAB-compatible features from one respiration signal."""

    if not isinstance(signal, Signal):
        raise TypeError("signal must be a Signal instance")
    if not isinstance(show_plots, bool):
        raise TypeError("show_plots must be a boolean")
    if not isinstance(save_plots, bool):
        raise TypeError("save_plots must be a boolean")
    if save_plots and output_dir is None:
        raise ValueError("output_dir must be provided when save_plots=True")

    _validate_signal(signal)
    downsample_factor = _validate_sampling_frequency(signal.sampling_frequency)
    respiration_downsampled = signal.samples[::downsample_factor]
    time_downsampled = signal.time[::downsample_factor]
    filtered_signal = _filter_respiration_signal(respiration_downsampled)
    feature_values, plot_record = _extract_respiration_features(
        filtered_signal,
        time_downsampled,
    )
    features = dict(feature_values)

    if show_plots or save_plots:
        _plot_respiration_diagnostics(
            time_downsampled,
            respiration_downsampled,
            filtered_signal,
            *plot_record,
            show_plots=show_plots,
            save_plots=save_plots,
            output_dir=output_dir,
        )

    metadata = {
        "source": "process_respiration",
        "original_sampling_frequency": signal.sampling_frequency,
        "processed_sampling_frequency": float(_FS_NEW),
        "original_n_samples": signal.n_samples,
        "processed_n_samples": int(respiration_downsampled.size),
        "start_time_seconds": float(signal.time[0]),
        "end_time_seconds": float(signal.time[-1]),
    }
    return ProcessingResult(
        modality="respiration",
        processor_name="process_respiration",
        features=features,
        signal_labels=[signal.label or "respiration"],
        metadata=metadata,
    )


def _validate_signal(signal: Signal) -> None:
    if signal.duration_seconds == 0:
        raise ValueError("respiration signal duration must be nonzero")


def _validate_sampling_frequency(sampling_frequency: float) -> int:
    ratio = sampling_frequency / _FS_NEW
    if not np.isclose(ratio, round(ratio)) or ratio < 1:
        raise ValueError("sampling_frequency must be an integer multiple of 100 Hz")
    return int(round(ratio))


def _filter_respiration_signal(signal: np.ndarray) -> np.ndarray:
    b, a = butter(4, 1.1 / (_FS_NEW / 2), btype="low")
    try:
        return filtfilt(b, a, signal)
    except ValueError as exc:
        raise ValueError("signal is too short for filtering") from exc


def _find_peaks_respiration(signal: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    up = False
    down = False
    min_rel: list[int] = []
    max_rel: list[int] = []

    if signal.size < 2:
        return np.asarray(min_rel, dtype=int), np.asarray(max_rel, dtype=int)

    if signal[0] - signal[1] > 0:
        up = True
    else:
        down = True

    for index in range(signal.size - 1):
        if up and signal[index + 1] < signal[index]:
            down = True
            up = False
            max_rel.append(index)
        if down and signal[index + 1] > signal[index]:
            up = True
            down = False
            min_rel.append(index)

    return np.asarray(min_rel, dtype=int), np.asarray(max_rel, dtype=int)


def _extract_respiration_features(
    filtered_phase: np.ndarray,
    phase_time: np.ndarray,
) -> tuple[Mapping[str, float], tuple[np.ndarray, np.ndarray]]:
    i_min, i_max = _find_peaks_respiration(filtered_phase)
    if i_min.size < 3 or i_max.size < 2:
        raise ValueError("too few respiration peaks detected")

    if i_min[0] >= i_max[0]:
        i_max = i_max[1:]
    if i_min.size == 0 or i_max.size == 0:
        raise ValueError("too few respiration peaks detected")
    if i_min[-1] <= i_max[-1]:
        i_max = i_max[:-1]
    if i_min.size < 3 or i_max.size < 2:
        raise ValueError("too few respiration peaks detected")

    i_min_def, i_max_def = _definitive_breath_points(filtered_phase, phase_time, i_min, i_max)
    if i_min_def.size < 3 or i_max_def.size < 2:
        raise ValueError("too few valid breaths detected")

    feature_values = _features_from_breath_points(
        filtered_phase,
        phase_time,
        i_min_def,
        i_max_def,
    )
    return feature_values, (i_min_def, i_max_def)


def _definitive_breath_points(
    filtered_phase: np.ndarray,
    phase_time: np.ndarray,
    i_min: np.ndarray,
    i_max: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    minimum_time_distance = np.diff(phase_time[i_min])
    if minimum_time_distance.size == 0:
        raise ValueError("too few respiration minima detected")
    mean_time_distance = float(np.mean(minimum_time_distance))
    max_min_distance_down = np.abs(filtered_phase[i_min[:-1]] - filtered_phase[i_max])
    max_min_distance_up = np.abs(filtered_phase[i_min[1:]] - filtered_phase[i_max])
    mean_max_min_distance = float(
        np.mean(np.concatenate([max_min_distance_down, max_min_distance_up]))
    )
    if not np.isfinite(mean_time_distance) or not np.isfinite(mean_max_min_distance):
        raise ValueError("too few valid breaths detected")

    start_candidates = np.flatnonzero(max_min_distance_down >= 0.6 * mean_max_min_distance)
    if start_candidates.size == 0:
        raise ValueError("too few valid breaths detected")

    start_candidate = int(start_candidates[0])
    i_min_def = [int(i_min[start_candidate])]
    i_max_def: list[int] = []
    minimum_time_distance_def: list[float] = []
    max_min_distance_down_def: list[float] = []
    max_min_distance_up_def: list[float] = []
    count_distance_def = 1

    for index in range(start_candidate, minimum_time_distance.size):
        i_min_point_down = i_min_def[count_distance_def - 1]
        i_min_point_up = int(i_min[index + 1])
        segment = filtered_phase[i_min_point_down : i_min_point_up + 1]
        i_max_point = i_min_point_down + int(np.argmax(segment))

        minimum_time_distance_wave = phase_time[i_min_point_up] - phase_time[i_min_point_down]
        max_min_distance_down_wave = filtered_phase[i_max_point] - filtered_phase[i_min_point_down]
        max_min_distance_up_wave = filtered_phase[i_max_point] - filtered_phase[i_min_point_up]

        if index < minimum_time_distance.size - 1:
            max_min_distance_up_wave_wave = abs(
                filtered_phase[i_max[index + 1]] - filtered_phase[i_min_point_up]
            )
            accept = (
                minimum_time_distance_wave >= 0.1 * mean_time_distance
                and max_min_distance_down_wave > 0.2 * mean_max_min_distance
                and max_min_distance_up_wave > 0.2 * mean_max_min_distance
                and max_min_distance_up_wave_wave > 0.2 * mean_max_min_distance
            )
        else:
            accept = (
                minimum_time_distance_wave >= 0.5 * mean_time_distance
                and max_min_distance_down_wave > 0.4 * mean_max_min_distance
                and max_min_distance_up_wave > 0.4 * mean_max_min_distance
            )

        if accept:
            i_min_def.append(i_min_point_up)
            i_max_def.append(i_max_point)
            minimum_time_distance_def.append(float(minimum_time_distance_wave))
            max_min_distance_down_def.append(float(max_min_distance_down_wave))
            max_min_distance_up_def.append(float(max_min_distance_up_wave))

            if len(i_min_def) > 7:
                mean_interval = slice(count_distance_def - 7, count_distance_def)
                mean_time_distance = float(np.mean(minimum_time_distance_def[mean_interval]))
                mean_max_min_distance = float(
                    np.mean(
                        max_min_distance_down_def[mean_interval]
                        + max_min_distance_up_def[mean_interval]
                    )
                )
            elif index == minimum_time_distance.size - 1 and len(i_min_def) > 2:
                mean_interval = slice(count_distance_def - 2, count_distance_def)
                mean_time_distance = float(np.mean(minimum_time_distance_def[mean_interval]))
                mean_max_min_distance = float(
                    np.mean(
                        max_min_distance_down_def[mean_interval]
                        + max_min_distance_up_def[mean_interval]
                    )
                )
            count_distance_def += 1

    return np.asarray(i_min_def, dtype=int), np.asarray(i_max_def, dtype=int)


def _features_from_breath_points(
    filtered_phase: np.ndarray,
    phase_time: np.ndarray,
    i_min_def: np.ndarray,
    i_max_def: np.ndarray,
) -> Mapping[str, float]:
    diff_t_max = np.diff(phase_time[i_max_def])
    respiration_rate = 60 / diff_t_max
    indexes_rr_ok = np.flatnonzero((respiration_rate >= 5) & (respiration_rate <= 30))
    if indexes_rr_ok.size < 2:
        raise ValueError("too few valid breaths detected")

    i_max_rr_ok = np.concatenate(([i_max_def[0]], i_max_def[indexes_rr_ok + 1]))
    i_min_rr_down_ok = np.concatenate(([i_min_def[0]], i_min_def[indexes_rr_ok + 1]))
    i_min_rr_up_ok = np.concatenate(([i_min_def[1]], i_min_def[indexes_rr_ok + 2]))
    respiration_rate = respiration_rate[indexes_rr_ok]

    inspiratory_time = phase_time[i_max_rr_ok] - phase_time[i_min_rr_down_ok]
    expiratory_time = phase_time[i_min_rr_up_ok] - phase_time[i_max_rr_ok]
    timing_ratio = inspiratory_time / expiratory_time

    lower_point_amplitude = []
    for lower_start, lower_end, peak in zip(i_min_rr_down_ok, i_min_rr_up_ok, i_max_rr_ok):
        coefficients = np.polyfit(
            [phase_time[lower_start], phase_time[lower_end]],
            [filtered_phase[lower_start], filtered_phase[lower_end]],
            1,
        )
        lower_point_amplitude.append(np.polyval(coefficients, phase_time[peak]))
    lower_point_amplitude = np.asarray(lower_point_amplitude)
    tidal_volume = filtered_phase[i_max_rr_ok] - lower_point_amplitude
    minute_ventilation = respiration_rate * tidal_volume[1:]

    respiration_rate_1 = respiration_rate[:-1]
    respiration_rate_2 = respiration_rate[1:]
    if respiration_rate_1.size < 2:
        raise ValueError("too few valid breaths detected")
    autocorrelation = abs(np.corrcoef(respiration_rate_1, respiration_rate_2)[0, 1])

    values = {
        "respiration_breath_rate_mean": np.mean(respiration_rate),
        "respiration_breath_rate_std": np.std(respiration_rate, ddof=1),
        "respiration_inspiratory_time_mean": np.mean(inspiratory_time),
        "respiration_inspiratory_time_std": np.std(inspiratory_time, ddof=1),
        "respiration_expiratory_time_mean": np.mean(expiratory_time),
        "respiration_expiratory_time_std": np.std(expiratory_time, ddof=1),
        "respiration_timing_ratio_mean": np.mean(timing_ratio),
        "respiration_timing_ratio_std": np.std(timing_ratio, ddof=1),
        "respiration_amplitude_mean": np.mean(tidal_volume),
        "respiration_amplitude_std": np.std(tidal_volume, ddof=1),
        "respiration_minute_ventilation_mean": np.mean(minute_ventilation),
        "respiration_minute_ventilation_std": np.std(minute_ventilation, ddof=1),
        "respiration_coefficient_of_variation": (
            np.std(respiration_rate, ddof=1) / np.mean(respiration_rate)
        )
        * 100,
        "respiration_autocorrelation": autocorrelation,
    }

    if not all(np.isfinite(value) for value in values.values()):
        raise ValueError("too few valid breaths detected")
    return {feature_name: float(values[feature_name]) for feature_name in _FEATURE_NAMES}


def _plot_respiration_diagnostics(
    time_downsampled: np.ndarray,
    respiration_downsampled: np.ndarray,
    filtered_signal: np.ndarray,
    i_min: np.ndarray,
    i_max: np.ndarray,
    *,
    show_plots: bool,
    save_plots: bool,
    output_dir: str | Path | None,
) -> None:
    output_path = None
    if save_plots:
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)

    figures = [
        _plot_downsampled_respiration(time_downsampled, respiration_downsampled),
        _plot_filtered_respiration(time_downsampled, filtered_signal, i_min, i_max),
    ]

    if save_plots and output_path is not None:
        figures[0].savefig(output_path / "respiration_downsampled_signal.png", dpi=150)
        figures[1].savefig(output_path / "respiration_filtered_phase_peaks.png", dpi=150)

    if show_plots:
        plt.show()
    else:
        for figure in figures:
            plt.close(figure)


def _plot_downsampled_respiration(
    time_downsampled: np.ndarray,
    respiration_downsampled: np.ndarray,
):
    figure, axes = plt.subplots()
    axes.plot(time_downsampled, respiration_downsampled, linewidth=1.5)
    axes.set_title("Respiration downsampled signal")
    axes.set_xlabel("Time [s]")
    axes.set_ylabel("Respiration [s]")
    axes.grid(True)
    return figure


def _plot_filtered_respiration(
    signal_time: np.ndarray,
    filtered_signal: np.ndarray,
    i_min: np.ndarray,
    i_max: np.ndarray,
):
    figure, axes = plt.subplots()
    axes.plot(signal_time, filtered_signal, linewidth=1.5)
    axes.plot(signal_time[i_min], filtered_signal[i_min], "o", linewidth=1.5)
    axes.plot(signal_time[i_max], filtered_signal[i_max], "*", linewidth=1.5)
    axes.set_title("Filtered respiration signal with detected minima and maxima")
    axes.set_xlabel("Time [s]")
    axes.set_ylabel("Chest Displacement [%]")
    axes.grid(True)
    return figure

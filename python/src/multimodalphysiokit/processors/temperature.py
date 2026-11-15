from __future__ import annotations

from collections.abc import Mapping
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import butter, filtfilt

from multimodalphysiokit.core import ProcessingResult, Signal


# A window is an artifact candidate when an extremum deviates from the phase
# mean by more than this many degrees (degC or K; about 6 % of 33 degC).
_ARTIFACT_DEVIATION_DEGREES = 2.0

_FEATURE_NAMES = (
    "temperature_delta",
    "temperature_mean",
    "temperature_std",
    "temperature_delta_over_time",
    "temperature_linear_slope",
    "temperature_derivative_mean",
    "temperature_derivative_std",
    "temperature_derivative_linear_slope",
)


class TemperatureProcessor:
    """Thin object-oriented wrapper for temperature processing."""

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
        """Return temperature features for the supplied signal."""

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
        """Process exactly one supplied temperature signal."""

        return process_temperature(
            signal,
            show_plots=show_plots,
            save_plots=save_plots,
            output_dir=output_dir,
        )

def process_temperature(
    signal: Signal,
    *,
    show_plots: bool = False,
    save_plots: bool = False,
    output_dir: str | Path | None = None,
) -> ProcessingResult:
    """Extract MATLAB-compatible features from one temperature signal."""

    if not isinstance(signal, Signal):
        raise TypeError("signal must be a Signal instance")
    if not isinstance(show_plots, bool):
        raise TypeError("show_plots must be a boolean")
    if not isinstance(save_plots, bool):
        raise TypeError("save_plots must be a boolean")
    if save_plots and output_dir is None:
        raise ValueError("output_dir must be provided when save_plots=True")

    _validate_signal(signal)
    filtered_temperature = _filter_temperature(signal)
    features = dict(
        _extract_phase_features(
            filtered_temperature,
            signal.time,
            signal.sampling_frequency,
        )
    )

    if show_plots or save_plots:
        _plot_temperature_diagnostics(
            signal,
            filtered_temperature,
            show_plots=show_plots,
            save_plots=save_plots,
            output_dir=output_dir,
        )

    metadata = {
        "source": "process_temperature",
        "sampling_frequency": signal.sampling_frequency,
        "n_samples": signal.n_samples,
        "start_time_seconds": float(signal.time[0]),
        "end_time_seconds": float(signal.time[-1]),
    }

    return ProcessingResult(
        modality="temperature",
        processor_name="process_temperature",
        features=features,
        signal_labels=[signal.label or "temperature"],
        metadata=metadata,
    )


def _filter_temperature(signal: Signal) -> np.ndarray:
    nyquist = signal.sampling_frequency / 2
    cutoff = 1.0
    if cutoff >= nyquist:
        raise ValueError("filtering cannot be applied with sampling_frequency <= 2 Hz")

    b, a = butter(5, cutoff / nyquist, btype="low")
    try:
        return filtfilt(b, a, signal.samples)
    except ValueError as exc:
        raise ValueError("filtering cannot be applied because the signal is too short") from exc


def _validate_signal(signal: Signal) -> None:
    if signal.duration_seconds == 0:
        raise ValueError("temperature signal duration must be nonzero")
    if signal.n_samples < 3:
        raise ValueError("temperature signal must contain at least 3 samples")


def _extract_phase_features(
    filtered_phase: np.ndarray,
    time_phase: np.ndarray,
    sampling_frequency: float,
) -> Mapping[str, float]:
    filtered_phase, time_phase = _screen_temperature_phase_artifacts(
        filtered_phase,
        time_phase,
        sampling_frequency=sampling_frequency,
    )
    derivative = np.diff(filtered_phase) / np.diff(time_phase)
    phase_duration = time_phase[-1] - time_phase[0]
    if derivative.size == 0:
        raise ValueError("temperature derivative cannot be computed after artifact rejection")

    temperature_initial = filtered_phase[0]
    temperature_final = filtered_phase[-1]
    temperature_delta = temperature_final - temperature_initial

    values = {
        "temperature_delta": temperature_delta,
        "temperature_mean": np.mean(filtered_phase),
        "temperature_std": np.std(filtered_phase, ddof=1),
        "temperature_delta_over_time": temperature_delta / phase_duration,
        "temperature_linear_slope": np.polyfit(time_phase, filtered_phase, 1)[0],
        "temperature_derivative_mean": np.mean(derivative),
        "temperature_derivative_std": np.std(derivative, ddof=1),
        "temperature_derivative_linear_slope": np.polyfit(time_phase[1:], derivative, 1)[0],
    }
    return {feature_name: float(values[feature_name]) for feature_name in _FEATURE_NAMES}


def _screen_temperature_phase_artifacts(
    filtered_phase: np.ndarray,
    time_phase: np.ndarray,
    sampling_frequency: float,
) -> tuple[np.ndarray, np.ndarray]:
    length_phase = filtered_phase.size
    if length_phase < 2:
        raise ValueError("temperature phase must contain at least 2 samples after artifact rejection")

    mean_filtered_phase = np.mean(filtered_phase)
    # Robust baseline: artifacts cannot inflate the median of |diff|.
    median_abs_diff_filtered_phase = np.median(np.abs(np.diff(filtered_phase)))
    window_length = int(2 * sampling_frequency)
    if window_length <= 0:
        raise ValueError("temperature artifact screening window length must be positive")
    threshold_diff = 100

    if window_length < length_phase:
        indexes_windows_phase = list(range(1, length_phase + 1, window_length))
        if length_phase - indexes_windows_phase[-1] < window_length and length_phase - indexes_windows_phase[-1] != 0:
            indexes_windows_phase.append(length_phase)
    else:
        indexes_windows_phase = [1, length_phase]

    valid = np.zeros(length_phase, dtype=bool)
    for left_one_based, right_one_based in zip(indexes_windows_phase[:-1], indexes_windows_phase[1:], strict=True):
        start_index = left_one_based - 1
        end_index = right_one_based - 1
        window = filtered_phase[start_index : end_index + 1]
        diff_window = np.diff(window)
        mean_diff_window = np.mean(np.abs(diff_window)) if diff_window.size else 0.0
        maxima, minima = _find_max_min(window)
        max_min_indexes = np.concatenate((maxima, minima))
        if max_min_indexes.size == 0:
            max_min_indexes = np.asarray(
                [0, _matlab_round(window.size / 2) - 1, window.size - 1],
                dtype=int,
            )
            max_min_indexes = np.clip(max_min_indexes, 0, window.size - 1)

        has_amplitude_artifact = np.any(
            np.abs(window[max_min_indexes] - mean_filtered_phase)
            > _ARTIFACT_DEVIATION_DEGREES
        )
        has_diff_artifact = mean_diff_window > threshold_diff * median_abs_diff_filtered_phase
        valid[start_index : end_index + 1] = not (has_amplitude_artifact and has_diff_artifact)

    if not np.any(valid):
        raise ValueError("all temperature samples were rejected as artifacts")
    if np.count_nonzero(valid) < 2:
        raise ValueError("temperature phase must contain at least 2 valid samples after artifact rejection")
    return filtered_phase[valid], time_phase[valid]


def _find_max_min(signal: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    if signal.size < 3:
        return np.asarray([], dtype=int), np.asarray([], dtype=int)
    up = signal[1] - signal[0] > 0
    down = not up
    minima: list[int] = []
    maxima: list[int] = []
    for index in range(1, signal.size - 1):
        if up and signal[index + 1] < signal[index]:
            down = True
            up = False
            maxima.append(index)
        if down and signal[index + 1] > signal[index]:
            up = True
            down = False
            minima.append(index)
    return np.asarray(maxima, dtype=int), np.asarray(minima, dtype=int)


def _matlab_round(value: float) -> int:
    return int(np.floor(value + 0.5))


def _plot_temperature_diagnostics(
    signal: Signal,
    filtered_temperature: np.ndarray,
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
        _plot_temperature_vs_time(signal, filtered_temperature),
        _plot_temperature_first_derivative_vs_time(signal),
    ]

    if save_plots and output_path is not None:
        figures[0].savefig(output_path / "temperature_vs_time.png", dpi=150)
        figures[1].savefig(
            output_path / "temperature_first_derivative_vs_time.png",
            dpi=150,
        )

    if show_plots:
        plt.show()
    else:
        for figure in figures:
            plt.close(figure)


def _plot_temperature_vs_time(
    signal: Signal,
    filtered_temperature: np.ndarray,
):
    figure, axes = plt.subplots()
    axes.plot(signal.time, signal.samples)
    axes.plot(signal.time, filtered_temperature, linewidth=1.5)
    axes.set_title("Temperature vs time")
    axes.set_xlabel("Time (s)")
    axes.set_ylabel("Raw temperature (°C)")
    axes.grid(True)
    return figure


def _plot_temperature_first_derivative_vs_time(
    signal: Signal,
):
    figure, axes = plt.subplots()
    axes.plot(signal.time, signal.samples, linewidth=1.5)

    derivative_nonnegative = np.diff(signal.samples) >= 0
    derivative_time = signal.time[1:]
    derivative_temperature = signal.samples[1:]
    axes.plot(
        derivative_time[~derivative_nonnegative],
        derivative_temperature[~derivative_nonnegative],
        "*r",
        linewidth=0.5,
    )
    axes.plot(
        derivative_time[derivative_nonnegative],
        derivative_temperature[derivative_nonnegative],
        "*b",
        linewidth=0.5,
    )
    axes.set_title("Temperature First Derivative vs Time")
    axes.set_xlabel("Time (s)")
    axes.set_ylabel("Raw temperature (°C)")
    axes.grid(True)
    return figure

def _to_matlab_cell(values: list[np.ndarray]) -> np.ndarray:
    """Convert variable-length NumPy arrays to a MATLAB cell array."""

    matlab_cell = np.empty((1, len(values)), dtype=object)

    for index, value in enumerate(values):
        matlab_cell[0, index] = np.asarray(value)

    return matlab_cell

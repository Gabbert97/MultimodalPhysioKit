from __future__ import annotations

from numbers import Real
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import butter, buttord, sosfilt, welch
from scipy.stats import kurtosis, skew

from multimodalphysiokit.core import ProcessingResult, Signal


_FS_NEW = 10

_FNIRS_LOWPASS_PASSBAND_HZ = 0.2
_FNIRS_LOWPASS_STOPBAND_HZ = 0.4
_FNIRS_LOWPASS_PASSBAND_RIPPLE_DB = 0.5
_FNIRS_LOWPASS_STOPBAND_ATTENUATION_DB = 80.0

_FNIRS_HIGHPASS_FS = 10.0

_FNIRS_HIGHPASS_SOS = np.array(
    [
        [0.9995548354966819, -1.9991096709933638, 0.9995548354966819, 1.0, -1.999092594752352, 0.9991267472343754],
        [0.9986924652981852, -1.9973849305963705, 0.9986924652981852, 1.0, -1.9973678690879584, 0.9974019921047824],
        [0.997860512403017, -1.995721024806034, 0.997860512403017, 1.0, -1.9957039775105772, 0.9957380721014907],
        [0.9970773680414668, -1.9941547360829337, 0.9970773680414668, 1.0, -1.9941377021665947, 0.9941717699992727],
        [0.9963602606920263, -1.9927205213840526, 0.9963602606920263, 1.0, -1.9927034997186652, 0.9927375430494398],
        [0.9957248949196912, -1.9914497898393824, 0.9957248949196912, 1.0, -1.991432779028486, 0.9914668006502786],
        [0.9951851284005139, -1.9903702568010277, 0.9951851284005139, 1.0, -1.9903532552114198, 0.9903872583906359],
        [0.9947526918940314, -1.9895053837880627, 0.9947526918940314, 1.0, -1.9894883895861333, 0.9895223779899921],
        [0.994436955848053, -1.988873911696106, 0.994436955848053, 1.0, -1.9888569228881627, 0.9888909005040493],
        [0.9942447463787437, -1.9884894927574874, 0.9942447463787437, 1.0, -1.988472507233221, 0.9885064782817538],
        [0.9970858601800293, -0.9970858601800293, 0.0, 1.0, -0.9941717203600587, 0.0],
    ],
    dtype=float,
)

_FNIRS_FEATURE_NAMES = (
    "min",
    "max",
    "mean",
    "variance",
    "std",
    "skewness",
    "kurtosis",
    "power_band",
    "max_frequency",
    "median_frequency",
    "spectral_entropy",
    "mean_diff",
    "peak_to_peak",
    "slope",
    "zero_crossing",
    "polarity",
    "entropy",
    "auc",
    "total_energy",
)


class FNIRSProcessor:
    """Object-oriented fNIRS preprocessing and concentration conversion."""

    def __init__(
        self,
        *,
        age_years: Real,
    ) -> None:
        if not isinstance(age_years, Real):
            raise TypeError("age_years must be numeric")
        if age_years <= 0 or not np.isfinite(age_years):
            raise ValueError("age_years must be positive and finite")
        self.age_years = float(age_years)

    def downsample_currents(self, red_signal: Signal, infrared_signal: Signal) -> dict[str, object]:
        """Downsample red and infrared currents to 10 Hz with moving averages."""

        _validate_signal_pair(red_signal, infrared_signal)
        factor = int(round(red_signal.sampling_frequency / _FS_NEW))
        red, time = _downsample_signal(red_signal.samples, red_signal.time, factor)
        infrared, _ = _downsample_signal(
            infrared_signal.samples,
            infrared_signal.time,
            factor,
        )

        return {
            "red": red,
            "infrared": infrared,
            "time": time,
            "sampling_frequency": _FS_NEW,
        }

    def extract_downsampled_currents(self, red_signal: Signal, infrared_signal: Signal) -> dict[str, object]:
        """Return the downsampled fNIRS current stage."""

        return self.downsample_currents(red_signal, infrared_signal)

    def filter_currents(self, red_signal: Signal, infrared_signal: Signal) -> dict[str, object]:
        """Low-pass filter downsampled fNIRS currents."""

        downsampled = self.downsample_currents(red_signal, infrared_signal)
        red = _lowpass_filter_current(np.asarray(downsampled["red"], dtype=float))
        infrared = _lowpass_filter_current(np.asarray(downsampled["infrared"], dtype=float))
        time = downsampled["time"]

        return {
            "red": red,
            "infrared": infrared,
            "time": time,
            "sampling_frequency": downsampled["sampling_frequency"],
        }

    def extract_filtered_currents(self, red_signal: Signal, infrared_signal: Signal) -> dict[str, object]:
        """Return the low-pass filtered fNIRS current stage."""

        return self.filter_currents(red_signal, infrared_signal)

    def compute_baseline_values(
        self, red_signal: Signal, infrared_signal: Signal
    ) -> tuple[float, float]:
        """Compute filtered-current reference means after trimming 10 seconds per end."""

        filtered = self.filter_currents(red_signal, infrared_signal)
        return _baseline_values_from_filtered(
            np.asarray(filtered["red"], dtype=float),
            np.asarray(filtered["infrared"], dtype=float),
        )

    def compute_concentration_changes(
        self,
        red_signal: Signal,
        infrared_signal: Signal,
        *,
        baseline_red_value: object = None,
        baseline_infrared_value: object = None,
    ) -> dict[str, object]:
        """Convert filtered currents to oxy/deoxy concentration changes."""

        baseline_red_value, baseline_infrared_value, baseline_mode = _validate_baseline_values(
            baseline_red_value, baseline_infrared_value
        )
        filtered = self.filter_currents(red_signal, infrared_signal)
        red = np.asarray(filtered["red"], dtype=float)
        infrared = np.asarray(filtered["infrared"], dtype=float)
        if baseline_mode == "self":
            baseline_red_value, baseline_infrared_value = _baseline_values_from_filtered(red, infrared)
        oxy, deoxy = _concentration_changes(
            red, infrared, baseline_red_value, baseline_infrared_value, self.age_years
        )

        return {
            "oxy": oxy,
            "deoxy": deoxy,
            "time": filtered["time"],
            "sampling_frequency": filtered["sampling_frequency"],
            "baseline_mode": baseline_mode,
            "baseline_red_value": baseline_red_value,
            "baseline_infrared_value": baseline_infrared_value,
        }

    def extract_concentration_changes(
        self,
        red_signal: Signal,
        infrared_signal: Signal,
        *,
        baseline_red_value: object = None,
        baseline_infrared_value: object = None,
        show_plots: bool = False,
        save_plots: bool = False,
        output_dir: str | Path | None = None,
        output_label: str = "fnirs",
    ) -> dict[str, object]:
        """Return Lambert-Beer oxy/deoxy concentration changes."""

        _validate_plot_options(
            show_plots=show_plots,
            save_plots=save_plots,
            output_dir=output_dir,
            output_label=output_label,
        )
        result = self.compute_concentration_changes(
            red_signal,
            infrared_signal,
            baseline_red_value=baseline_red_value,
            baseline_infrared_value=baseline_infrared_value,
        )
        if show_plots or save_plots:
            _plot_fnirs_concentration_changes(
                np.asarray(result["oxy"], dtype=float),
                np.asarray(result["deoxy"], dtype=float),
                np.asarray(result["time"], dtype=float),
                show_plots=show_plots,
                save_plots=save_plots,
                output_dir=output_dir,
                output_label=output_label,
            )
        return result

    def extract_features(
        self,
        red_signal: Signal,
        infrared_signal: Signal,
        *,
        baseline_red_value: object = None,
        baseline_infrared_value: object = None,
        show_plots: bool = False,
        save_plots: bool = False,
        output_dir: str | Path | None = None,
        output_label: str = "fnirs",
    ) -> dict[str, float]:
        """Extract fNIRS features for oxy and deoxy signals."""

        _validate_plot_options(
            show_plots=show_plots,
            save_plots=save_plots,
            output_dir=output_dir,
            output_label=output_label,
        )
        concentration = self.extract_concentration_changes(
            red_signal,
            infrared_signal,
            baseline_red_value=baseline_red_value,
            baseline_infrared_value=baseline_infrared_value,
            show_plots=show_plots,
            save_plots=save_plots,
            output_dir=output_dir,
            output_label=output_label,
        )
        return _features_from_concentration(
            concentration,
            show_plots=show_plots,
            save_plots=save_plots,
            output_dir=output_dir,
            output_label=output_label,
        )

    def process(
        self,
        red_signal: Signal,
        infrared_signal: Signal,
        *,
        baseline_red_value: object = None,
        baseline_infrared_value: object = None,
        show_plots: bool = False,
        save_plots: bool = False,
        output_dir: str | Path | None = None,
        output_label: str = "fnirs",
    ) -> ProcessingResult:
        """Extract fNIRS features and return a ProcessingResult."""

        _validate_plot_options(
            show_plots=show_plots,
            save_plots=save_plots,
            output_dir=output_dir,
            output_label=output_label,
        )
        concentration = self.extract_concentration_changes(
            red_signal,
            infrared_signal,
            baseline_red_value=baseline_red_value,
            baseline_infrared_value=baseline_infrared_value,
            show_plots=show_plots,
            save_plots=save_plots,
            output_dir=output_dir,
            output_label=output_label,
        )
        features = _features_from_concentration(
            concentration,
            show_plots=show_plots,
            save_plots=save_plots,
            output_dir=output_dir,
            output_label=output_label,
        )
        return ProcessingResult(
            modality="fnirs",
            processor_name="process_fnirs",
            features=features,
            signal_labels=[red_signal.label or "fnirs_red", infrared_signal.label or "fnirs_infrared"],
            metadata={
                "source": "process_fnirs",
                "age_years": self.age_years,
                "original_sampling_frequency": red_signal.sampling_frequency,
                "processed_sampling_frequency": _FS_NEW,
                "original_n_samples": red_signal.n_samples,
                "processed_n_samples": _matlab_round(red_signal.n_samples / int(round(red_signal.sampling_frequency / _FS_NEW))),
                "start_time_seconds": float(red_signal.time[0]),
                "end_time_seconds": float(red_signal.time[-1]),
                "baseline_mode": concentration["baseline_mode"],
                "baseline_red_value": float(concentration["baseline_red_value"]),
                "baseline_infrared_value": float(concentration["baseline_infrared_value"]),
            },
        )


def process_fnirs(
    red_signal: Signal,
    infrared_signal: Signal,
    *,
    age_years: Real,
    baseline_red_value: object = None,
    baseline_infrared_value: object = None,
    show_plots: bool = False,
    save_plots: bool = False,
    output_dir: str | Path | None = None,
    output_label: str = "fnirs",
) -> ProcessingResult:
    """Functional wrapper for FNIRSProcessor.process."""

    return FNIRSProcessor(
        age_years=age_years,
    ).process(
        red_signal,
        infrared_signal,
        baseline_red_value=baseline_red_value,
        baseline_infrared_value=baseline_infrared_value,
        show_plots=show_plots,
        save_plots=save_plots,
        output_dir=output_dir,
        output_label=output_label,
    )


def _validate_signal_pair(red_signal: Signal, infrared_signal: Signal) -> None:
    if not isinstance(red_signal, Signal):
        raise TypeError("red_signal must be a Signal instance")
    if not isinstance(infrared_signal, Signal):
        raise TypeError("infrared_signal must be a Signal instance")
    if red_signal.n_samples != infrared_signal.n_samples:
        raise ValueError("red_signal and infrared_signal must have the same number of samples")
    if red_signal.sampling_frequency != infrared_signal.sampling_frequency:
        raise ValueError("red_signal and infrared_signal must have the same sampling_frequency")
    if red_signal.time.shape != infrared_signal.time.shape:
        raise ValueError("red_signal and infrared_signal must have the same time vector shape")
    if not np.array_equal(red_signal.time, infrared_signal.time):
        raise ValueError("red_signal and infrared_signal must have synchronized time vectors")
    if np.any(red_signal.samples <= 0) or np.any(infrared_signal.samples <= 0):
        raise ValueError("fNIRS current signals must contain only positive values")
    ratio = red_signal.sampling_frequency / _FS_NEW
    if not np.isclose(ratio, round(ratio)) or ratio < 1:
        raise ValueError("fNIRS sampling_frequency must be an integer multiple of 10 Hz")


def _validate_baseline_values(red_value: object, infrared_value: object) -> tuple[float | None, float | None, str]:
    if (red_value is None) != (infrared_value is None):
        raise ValueError("baseline_red_value and baseline_infrared_value must be supplied together")
    if red_value is None:
        return None, None, "self"
    return _validate_baseline_scalar(red_value, "baseline_red_value"), _validate_baseline_scalar(
        infrared_value, "baseline_infrared_value"
    ), "provided"


def _validate_baseline_scalar(value: object, name: str) -> float:
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, Real):
        raise TypeError(f"{name} must be a numeric scalar")
    normalized = float(value)
    if not np.isfinite(normalized) or normalized <= 0:
        raise ValueError(f"{name} must be positive and finite")
    return normalized


def _baseline_values_from_filtered(red: np.ndarray, infrared: np.ndarray) -> tuple[float, float]:
    margin = 10 * _FS_NEW
    if red.size <= 2 * margin or infrared.size <= 2 * margin:
        raise ValueError("fNIRS signal must exceed 20 seconds for baseline estimation")
    return float(np.mean(red[margin:-margin])), float(np.mean(infrared[margin:-margin]))


def _concentration_changes(
    red: np.ndarray,
    infrared: np.ndarray,
    baseline_red_value: float,
    baseline_infrared_value: float,
    age_years: float,
) -> tuple[np.ndarray, np.ndarray]:
    delta_od_red = np.log10(baseline_red_value / red)
    delta_od_infrared = np.log10(baseline_infrared_value / infrared)

    epsilon = np.asarray([3.4408, 0.3346, 0.7977, 1.2071], dtype=float)
    source_detector_distance_cm = 2
    lambda1 = 660
    lambda2 = 860
    dpf_lambda1 = _differential_pathlength_factor(lambda1, age_years)
    dpf_lambda2 = _differential_pathlength_factor(lambda2, age_years)
    denominator = (
        (epsilon[0] * epsilon[3]) - (epsilon[2] * epsilon[1])
    ) * source_detector_distance_cm

    oxy = (
        ((epsilon[0] * delta_od_infrared) / dpf_lambda2)
        - ((epsilon[2] * delta_od_red) / dpf_lambda1)
    ) / denominator
    deoxy = (
        ((epsilon[3] * delta_od_red) / dpf_lambda1)
        - ((epsilon[1] * delta_od_infrared) / dpf_lambda2)
    ) / denominator
    return oxy, deoxy


def _features_from_concentration(
    concentration: dict[str, object],
    *,
    show_plots: bool,
    save_plots: bool,
    output_dir: str | Path | None,
    output_label: str,
) -> dict[str, float]:
    common = {
        "sampling_frequency": float(concentration["sampling_frequency"]),
        "show_plots": show_plots,
        "save_plots": save_plots,
        "output_dir": output_dir,
        "output_label": output_label,
    }
    deoxy = _extract_fnirs_channel_features(
        np.asarray(concentration["deoxy"], dtype=float),
        np.asarray(concentration["time"], dtype=float),
        channel_prefix="fnirs_deoxy",
        **common,
    )
    oxy = _extract_fnirs_channel_features(
        np.asarray(concentration["oxy"], dtype=float),
        np.asarray(concentration["time"], dtype=float),
        channel_prefix="fnirs_oxy",
        **common,
    )
    return {**oxy, **deoxy}


def _extract_fnirs_channel_features(
    signal: np.ndarray,
    time: np.ndarray,
    sampling_frequency: float,
    *,
    channel_prefix: str,
    show_plots: bool = False,
    save_plots: bool = False,
    output_dir: str | Path | None = None,
    output_label: str = "fnirs",
) -> dict[str, float]:
    _validate_plot_options(
        show_plots=show_plots,
        save_plots=save_plots,
        output_dir=output_dir,
        output_label=output_label,
    )
    if not np.all(np.isfinite(signal)):
        raise ValueError("fNIRS concentration signal must contain finite values")
    if not np.all(np.isfinite(time)):
        raise ValueError("fNIRS time vector must contain finite values")
    if signal.size != time.size:
        raise ValueError("fNIRS signal and time must have the same length")

    filtered_signal = _highpass_filter_concentration(signal, sampling_frequency)
    if show_plots or save_plots:
        _plot_fnirs_highpass(
            signal,
            filtered_signal,
            time,
            channel_prefix=channel_prefix,
            show_plots=show_plots,
            save_plots=save_plots,
            output_dir=output_dir,
            output_label=output_label,
        )

    values = _fnirs_features(signal, filtered_signal, time, sampling_frequency)
    return {f"{channel_prefix}_{name}": float(values[name]) for name in _FNIRS_FEATURE_NAMES}


def _fnirs_features(
    signal: np.ndarray,
    filtered_signal: np.ndarray,
    time: np.ndarray,
    sampling_frequency: float,
) -> dict[str, float]:
    if signal.size < 4:
        raise ValueError("fNIRS signal is too short for feature extraction")
    duration = float(time[-1] - time[0])
    if duration == 0:
        raise ValueError("fNIRS signal duration must be nonzero")

    mean_feature = float(np.mean(signal))
    entropy_feature = _histogram_entropy(signal)
    auc = np.sum(0.5 * np.diff(time) * np.abs(signal[:-1] + signal[1:]))
    auc_feature = float(auc / duration)
    total_energy_feature = float(np.sum(signal**2 / duration))
    slope_feature = float(np.polyfit(time, signal, 1)[0])

    # Preserve the historical MATLAB feature convention explicitly: biased
    # finite-sample skewness and biased Pearson kurtosis.
    skewness_feature = float(skew(filtered_signal, bias=True))
    kurtosis_feature = float(kurtosis(filtered_signal, fisher=False, bias=True))
    if not np.isfinite(skewness_feature) or not np.isfinite(kurtosis_feature):
        raise ValueError("non-finite fNIRS skewness or kurtosis")

    spectral = _fnirs_spectral_features(filtered_signal, sampling_frequency)
    peak_to_peak = _fnirs_peak_to_peak(filtered_signal, sampling_frequency)
    min_value = float(np.min(filtered_signal))
    max_value = float(np.max(filtered_signal))
    if np.isclose(min_value, 0.0):
        raise ValueError("fNIRS polarity cannot be computed with zero minimum")

    values = {
        "min": min_value,
        "max": max_value,
        "mean": mean_feature,
        "variance": float(np.var(filtered_signal, ddof=1)),
        "std": float(np.std(filtered_signal, ddof=1)),
        "skewness": skewness_feature,
        "kurtosis": kurtosis_feature,
        "power_band": spectral["power_band"],
        "max_frequency": spectral["max_frequency"],
        "median_frequency": spectral["median_frequency"],
        "spectral_entropy": spectral["spectral_entropy"],
        "mean_diff": float(np.mean(np.diff(filtered_signal))),
        "peak_to_peak": peak_to_peak,
        "slope": slope_feature,
        "zero_crossing": _zero_crossing_rate(filtered_signal),
        "polarity": float(abs(max_value / min_value)),
        "entropy": entropy_feature,
        "auc": auc_feature,
        "total_energy": total_energy_feature,
    }
    if not all(np.isfinite(value) for value in values.values()):
        raise ValueError("non-finite fNIRS feature value")
    return values


def _highpass_filter_concentration(
    signal: np.ndarray,
    sampling_frequency: float,
    *,
    pad_seconds: float = 180.0,
) -> np.ndarray:
    """High-pass filter fNIRS concentration using MATLAB-derived SOS coefficients."""

    if not np.isclose(sampling_frequency, _FNIRS_HIGHPASS_FS):
        raise ValueError(
            "The stored fNIRS high-pass filter coefficients were designed for "
            f"{_FNIRS_HIGHPASS_FS} Hz, but received {sampling_frequency} Hz."
        )

    sos = np.ascontiguousarray(_FNIRS_HIGHPASS_SOS, dtype=float)

    signal = np.asarray(signal, dtype=float).squeeze()
    signal = np.ascontiguousarray(signal, dtype=float)

    if signal.ndim != 1:
        raise ValueError("fNIRS concentration signal must be one-dimensional.")

    if signal.size < 4:
        raise ValueError("fNIRS concentration signal is too short for high-pass filtering.")

    if not np.all(np.isfinite(signal)):
        raise ValueError("fNIRS concentration signal must contain only finite values.")

    pad_samples = min(
        int(round(pad_seconds * sampling_frequency)),
        signal.size - 2,
    )
    try:
        return _explicit_zero_phase_sos_filter(
            sos,
            signal,
            padding_length=pad_samples,
        )
    except ValueError as exc:
        raise ValueError(f"fNIRS high-pass filtering failed: {exc}") from exc


def _histogram_entropy(signal: np.ndarray) -> float:
    x = np.asarray(signal, dtype=float).squeeze()

    if x.size == 0:
        return 0.0

    if np.max(x) == np.min(x):
        return 0.0

    edges = np.linspace(np.min(x), np.max(x), 11)
    counts, _ = np.histogram(x, bins=edges)

    total = np.sum(counts)
    if total == 0:
        return 0.0

    probability = counts.astype(float) / total
    probability = probability[probability > 0]
    if probability.size == 0:
        return 0.0
    return float(-np.sum(probability * np.log2(probability)))


def _zero_crossing_rate(signal: np.ndarray) -> float:
    sample_sign = np.zeros(signal.size, dtype=int)
    sample_sign[signal > 0] = 1
    return float(np.count_nonzero(np.diff(sample_sign) != 0) / signal.size)


def _fnirs_spectral_features(
    filtered_signal: np.ndarray,
    sampling_frequency: float,
) -> dict[str, float]:
    signal_no_mean = filtered_signal - np.mean(filtered_signal)
    n = _next_power_of_two(signal_no_mean.size)
    frequencies_entropy = sampling_frequency * np.arange(n // 2 + 1) / n
    magnitude = np.abs(np.fft.fft(signal_no_mean, n))[: n // 2 + 1]
    power = magnitude**2

    if np.sum(power) == 0:
        spectral_entropy = 0.0
    else:
        probability = power / np.sum(power)
        probability = probability[probability != 0]
        spectral_entropy = (
            0.0
            if probability.size == 1
            else float(-np.sum(probability * np.log2(probability)) / np.log2(probability.size))
        )

    cumulative_magnitude = np.cumsum(magnitude)
    if cumulative_magnitude[-1] == 0:
        median_frequency = 0.0
    else:
        median_index = int(np.argmax(cumulative_magnitude > cumulative_magnitude[-1] * 0.50))
        median_frequency = float(frequencies_entropy[median_index])

    freq, psd_power = _fnirs_welch_psd(filtered_signal, sampling_frequency)
    if np.sum(psd_power) == 0:
        power_band = 0.0
        max_frequency = 0.0
    else:
        max_frequency = float(freq[int(np.argmax(psd_power))])
        lower_index, upper_index = _matlab_power_band_indices(psd_power)
        power_band = float(abs(freq[upper_index] - freq[lower_index]))

    return {
        "power_band": power_band,
        "max_frequency": max_frequency,
        "median_frequency": median_frequency,
        "spectral_entropy": spectral_entropy,
    }


def _fnirs_welch_psd(
    filtered_signal: np.ndarray,
    sampling_frequency: float,
) -> tuple[np.ndarray, np.ndarray]:
    window_length = 256
    overlap = window_length // 2
    if filtered_signal.size % window_length != 0:
        num_windows = int(filtered_signal.size // overlap)
    else:
        num_windows = int(filtered_signal.size // overlap) - 1
    if num_windows < 1:
        raise ValueError("fNIRS signal is too short for spectral feature extraction")

    pxx_values = []
    frequencies = None
    for window_index in range(num_windows):
        start_index = window_index * (window_length - overlap)
        if window_index == num_windows - 1 and filtered_signal.size % window_length != 0:
            end_index = filtered_signal.size
        else:
            end_index = start_index + window_length
        window_signal = filtered_signal[start_index:end_index]
        window_signal = window_signal - np.mean(window_signal)
        if window_signal.size < window_length:
            window_signal = np.pad(window_signal, (0, window_length - window_signal.size))
        frequencies, pxx = welch(
            window_signal,
            fs=sampling_frequency,
            window=np.hamming(window_signal.size),
            noverlap=0,
            nfft=4096,
        )
        pxx_values.append(pxx)

    if frequencies is None:
        raise ValueError("fNIRS signal is too short for spectral feature extraction")
    return frequencies, np.mean(np.column_stack(pxx_values), axis=1)


def _matlab_power_band_indices(psd_power: np.ndarray) -> tuple[int, int]:
    """Translate MATLAB's one-based 95% cumulative PSD indexes exactly."""

    cumulative_power = np.cumsum(psd_power)
    lower_index = int(np.argmax(cumulative_power >= cumulative_power[-1] * 0.95))

    cumulative_power_inverse = np.cumsum(psd_power[::-1])
    inverse_index = int(
        np.argmax(cumulative_power_inverse >= cumulative_power[-1] * 0.95)
    )
    inverse_index_one_based = inverse_index + 1
    upper_index_one_based = abs(
        inverse_index_one_based - psd_power.size + 1
    )
    if upper_index_one_based == 0:
        upper_index_one_based = 1
    upper_index = min(upper_index_one_based - 1, psd_power.size - 1)
    return lower_index, upper_index


def _fnirs_peak_to_peak(filtered_signal: np.ndarray, sampling_frequency: float) -> float:
    filter_order = 5
    sos = butter(
        filter_order,
        0.08,
        btype="low",
        fs=sampling_frequency,
        output="sos",
    )
    try:
        smoothed = _explicit_zero_phase_sos_filter(
            sos,
            filtered_signal,
            padding_length=3 * filter_order,
        )
    except ValueError as exc:
        raise ValueError("fNIRS signal is too short for peak-to-peak filtering") from exc

    derivative = np.diff(smoothed)
    maxima = np.flatnonzero((derivative[:-1] > 0) & (derivative[1:] <= 0)) + 1
    minima = np.flatnonzero((derivative[:-1] < 0) & (derivative[1:] >= 0)) + 1
    if maxima.size == 0 or minima.size == 0:
        return 0.0
    if minima.size == maxima.size == 1:
        if minima[0] < maxima[0]:
            return 0.0
        return float(abs(smoothed[maxima[0]] - smoothed[minima[0]]))

    if minima[0] < maxima[0]:
        minima = minima[1:]
    if maxima.size and minima.size and maxima[-1] > minima[-1]:
        maxima = maxima[:-1]
    if maxima.size == 0 or minima.size == 0:
        return 0.0
    n_pairs = min(maxima.size, minima.size)
    return float(np.max(np.abs(smoothed[maxima[:n_pairs]] - smoothed[minima[:n_pairs]])))


def _plot_fnirs_highpass(
    signal: np.ndarray,
    filtered_signal: np.ndarray,
    time: np.ndarray,
    *,
    channel_prefix: str,
    show_plots: bool,
    save_plots: bool,
    output_dir: str | Path | None,
    output_label: str,
) -> None:
    output_path = None
    if save_plots:
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)

    if channel_prefix == "fnirs_oxy":
        title = "fNIRS oxy concentration - Before and after HPF"
        signal_label = "Oxy concentration"
        filtered_label = "Oxy concentration after HPF"
    elif channel_prefix == "fnirs_deoxy":
        title = "fNIRS deoxy concentration - Before and after HPF"
        signal_label = "Deoxy concentration"
        filtered_label = "Deoxy concentration after HPF"
    else:
        title = "fNIRS concentration - Before and after HPF"
        signal_label = "Concentration"
        filtered_label = "Concentration after HPF"

    figure, axes = plt.subplots()
    axes.plot(time, signal, linewidth=1.1, label=signal_label)
    axes.plot(time, filtered_signal, linewidth=1.1, label=filtered_label)
    axes.set_ylabel("Concentration change [mM]")
    axes.set_xlabel("Time [s]")
    axes.set_title(title)
    axes.legend()
    axes.grid(True)

    if save_plots and output_path is not None:
        channel_suffix = channel_prefix.removeprefix("fnirs_")
        figure.savefig(output_path / f"{output_label}_{channel_suffix}_highpass.png", dpi=150)

    if show_plots:
        plt.show()
    else:
        plt.close(figure)


def _plot_fnirs_concentration_changes(
    oxy: np.ndarray,
    deoxy: np.ndarray,
    time: np.ndarray,
    *,
    show_plots: bool = False,
    save_plots: bool = False,
    output_dir: str | Path | None = None,
    output_label: str = "fnirs",
) -> None:
    _validate_plot_options(
        show_plots=show_plots,
        save_plots=save_plots,
        output_dir=output_dir,
        output_label=output_label,
    )

    output_path = None
    if save_plots:
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)

    figure, axes = plt.subplots()
    axes.plot(time, oxy, linewidth=1.2, label="Oxy")
    axes.plot(time, deoxy, linewidth=1.2, label="Deoxy")
    axes.set_title("fNIRS oxy/deoxy concentration changes")
    axes.set_xlabel("Time [s]")
    axes.set_ylabel("Concentration change [mM]")
    axes.legend()
    axes.grid(True)

    if save_plots and output_path is not None:
        figure.savefig(output_path / f"{output_label}_concentration_changes.png", dpi=150)

    if show_plots:
        plt.show()
    else:
        plt.close(figure)


def _validate_plot_options(
    *,
    show_plots: bool,
    save_plots: bool,
    output_dir: str | Path | None,
    output_label: str,
) -> None:
    if not isinstance(show_plots, bool):
        raise TypeError("show_plots must be a boolean")
    if not isinstance(save_plots, bool):
        raise TypeError("save_plots must be a boolean")
    if not isinstance(output_label, str):
        raise TypeError("output_label must be a string")
    if save_plots and output_dir is None:
        raise ValueError("output_dir must be provided when save_plots=True")


def _downsample_signal(
    samples: np.ndarray,
    time: np.ndarray,
    factor: int,
) -> tuple[np.ndarray, np.ndarray]:
    n_downsampled = _matlab_round(samples.size / factor)
    downsampled = np.zeros(n_downsampled, dtype=float)
    downsampled_time = np.zeros(n_downsampled, dtype=float)

    for index in range(n_downsampled):
        start_index = index * factor
        end_index = min(start_index + factor, samples.size)
        grouped = samples[start_index:end_index]
        downsampled[index] = np.mean(grouped)
        time_index = start_index + _matlab_round(grouped.size / 2) - 1
        downsampled_time[index] = time[time_index]

    return downsampled, downsampled_time


def _lowpass_filter_current(samples: np.ndarray) -> np.ndarray:
    sos, order, _ = _design_lowpass_current_sos()
    try:
        return _explicit_zero_phase_sos_filter(sos, samples, padding_length=3 * order)
    except ValueError as exc:
        raise ValueError("fNIRS signal is too short for low-pass filtering") from exc


def _design_lowpass_current_sos() -> tuple[np.ndarray, int, float]:
    """Design the canonical stopband-matched Butterworth current filter."""

    order, _ = buttord(
        _FNIRS_LOWPASS_PASSBAND_HZ,
        _FNIRS_LOWPASS_STOPBAND_HZ,
        _FNIRS_LOWPASS_PASSBAND_RIPPLE_DB,
        _FNIRS_LOWPASS_STOPBAND_ATTENUATION_DB,
        fs=_FS_NEW,
    )

    # MATLAB designfilt(..., DesignMethod="butter") matches the stopband
    # exactly. SciPy buttord returns a passband-matched cutoff, so calculate
    # the equivalent stopband-matched digital cutoff explicitly.
    warped_stopband = 2.0 * _FS_NEW * np.tan(
        np.pi * _FNIRS_LOWPASS_STOPBAND_HZ / _FS_NEW
    )
    warped_cutoff = warped_stopband / (
        10.0 ** (_FNIRS_LOWPASS_STOPBAND_ATTENUATION_DB / 10.0) - 1.0
    ) ** (1.0 / (2.0 * order))
    critical_frequency_hz = float(
        _FS_NEW / np.pi * np.arctan(warped_cutoff / (2.0 * _FS_NEW))
    )
    sos = butter(
        order,
        critical_frequency_hz,
        btype="low",
        fs=_FS_NEW,
        output="sos",
    )
    return np.ascontiguousarray(sos, dtype=float), int(order), critical_frequency_hz


def _explicit_zero_phase_sos_filter(
    sos: np.ndarray,
    samples: np.ndarray,
    *,
    padding_length: int,
) -> np.ndarray:
    """Apply SOS forward/backward filtering with explicit odd reflection."""

    values = np.ascontiguousarray(np.asarray(samples, dtype=float).squeeze())
    if values.ndim != 1 or values.size < 2:
        raise ValueError("zero-phase filtering requires at least two samples")
    if not np.all(np.isfinite(values)):
        raise ValueError("zero-phase filtering requires finite samples")
    if padding_length < 1:
        raise ValueError("padding_length must be positive")

    padding_length = min(int(padding_length), values.size - 1)
    pre_padding = 2.0 * values[0] - values[1 : padding_length + 1][::-1]
    post_padding = 2.0 * values[-1] - values[-padding_length - 1 : -1][::-1]
    padded = np.concatenate((pre_padding, values, post_padding))

    initial_conditions = _sos_steady_state_initial_conditions(sos)
    forward, _ = sosfilt(sos, padded, zi=initial_conditions * padded[0])
    backward, _ = sosfilt(
        sos,
        forward[::-1],
        zi=initial_conditions * forward[-1],
    )
    filtered = backward[::-1]
    return np.ascontiguousarray(
        filtered[padding_length : padding_length + values.size],
        dtype=float,
    )


def _sos_steady_state_initial_conditions(sos: np.ndarray) -> np.ndarray:
    """Return cascade-aware direct-form-II states for a unit constant input."""

    states = np.empty((sos.shape[0], 2), dtype=float)
    section_input = 1.0
    for index, section in enumerate(sos):
        b0, b1, b2, a0, a1, a2 = section
        section_output = (
            (b0 + b1 + b2) / (a0 + a1 + a2) * section_input
        )
        states[index] = (
            section_output - b0 * section_input,
            b2 * section_input - a2 * section_output,
        )
        section_input = section_output
    return states


def _differential_pathlength_factor(wavelength_nm: float, age_years: float) -> float:
    return (
        223.3
        + 0.05624 * age_years**0.8493
        - 5.723e-7 * wavelength_nm**3
        + 0.001245 * wavelength_nm**2
        - 0.9025 * wavelength_nm
    )


def _next_power_of_two(value: int) -> int:
    return 1 << (int(value) - 1).bit_length()


def _matlab_round(value: float) -> int:
    return int(np.floor(value + 0.5))


def _matlab_round_array(values: np.ndarray) -> np.ndarray:
    return np.floor(values + 0.5).astype(int)

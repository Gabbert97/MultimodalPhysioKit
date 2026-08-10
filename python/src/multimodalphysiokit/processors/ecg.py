from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.signal import butter, filtfilt, find_peaks

from multimodalphysiokit.core import ProcessingResult, Signal


_TIME_DOMAIN_FEATURE_NAMES = (
    "ecg_bpm_mean",
    "ecg_bpm_std",
    "ecg_bpm_median",
    "ecg_pnn50",
    "ecg_pnn50_max",
)

_FREQUENCY_DOMAIN_FEATURE_NAMES = (
    "ecg_plf_mean",
    "ecg_plf_std",
    "ecg_phf_mean",
    "ecg_phf_std",
    "ecg_lf_hf_mean",
    "ecg_lf_hf_std",
)


class ECGProcessor:
    """Object-oriented ECG processor interface."""

    def __init__(
        self,
        *,
        label_frequency_analysis: int = 1,
        return_intermediates: bool = False,
    ) -> None:
        if label_frequency_analysis not in (1, 2):
            raise ValueError("label_frequency_analysis must be 1 or 2")
        if not isinstance(return_intermediates, bool):
            raise TypeError("return_intermediates must be a boolean")
        self.label_frequency_analysis = label_frequency_analysis
        self.return_intermediates = return_intermediates

    def detect_r_peaks(
        self,
        signal: Signal,
        *,
        show_plots: bool = False,
        save_plots: bool = False,
        output_dir: str | Path | None = None,
    ) -> dict[str, np.ndarray]:
        """Detect R-peaks in the supplied signal."""

        return _detect_ecg_r_peaks_impl(
            signal,
            show_plots=show_plots,
            save_plots=save_plots,
            output_dir=output_dir,
        )

    def extract_time_domain_features(
        self,
        signal: Signal,
        *,
        show_plots: bool = False,
        save_plots: bool = False,
        output_dir: str | Path | None = None,
        output_label: str = "ecg",
    ) -> dict[str, object]:
        """Extract ECG time-domain features for the supplied signal."""

        _validate_ecg_plot_options(
            show_plots=show_plots,
            save_plots=save_plots,
            output_dir=output_dir,
            output_label=output_label,
        )
        result = _extract_ecg_time_domain_features(signal)
        if show_plots or save_plots:
            _plot_ecg_bpm_trend(
                result["bpm"],
                show_plots=show_plots,
                save_plots=save_plots,
                output_dir=output_dir,
                output_label=output_label,
            )
        return result

    def extract_frequency_domain_features(
        self,
        signal: Signal,
        *,
        show_plots: bool = False,
        save_plots: bool = False,
        output_dir: str | Path | None = None,
        output_label: str = "ecg",
    ) -> dict[str, float]:
        """Extract ECG frequency-domain HRV features for the supplied signal."""

        _validate_ecg_frequency_options(
            label_frequency_analysis=self.label_frequency_analysis,
            show_plots=show_plots,
            save_plots=save_plots,
            output_dir=output_dir,
            output_label=output_label,
        )
        time_domain = self.extract_time_domain_features(signal)
        stft_context = None
        if self.label_frequency_analysis == 1:
            _validate_ecg_stft_duration(signal)
            stft_context = _build_ecg_stft_context(signal)

        features, _ = _extract_ecg_frequency_domain_features(
            time_domain["ibi"],
            label_frequency_analysis=self.label_frequency_analysis,
            show_plots=show_plots,
            save_plots=save_plots,
            output_dir=output_dir,
            output_label=output_label,
            stft_context=stft_context,
        )
        return features

    def process(
        self,
        signal: Signal,
        *,
        show_plots: bool = False,
        save_plots: bool = False,
        output_dir: str | Path | None = None,
        output_label: str = "ecg",
    ) -> ProcessingResult:
        """Extract all MATLAB-compatible ECG and HRV features."""

        _validate_ecg_frequency_options(
            label_frequency_analysis=self.label_frequency_analysis,
            show_plots=show_plots,
            save_plots=save_plots,
            output_dir=output_dir,
            output_label=output_label,
        )
        time_domain = self.extract_time_domain_features(signal)
        if show_plots or save_plots:
            _plot_ecg_bpm_trend(
                time_domain["bpm"],
                show_plots=show_plots,
                save_plots=save_plots,
                output_dir=output_dir,
                output_label=output_label,
            )
        stft_context = None
        if self.label_frequency_analysis == 1:
            _validate_ecg_stft_duration(signal)
            stft_context = _build_ecg_stft_context(signal)

        frequency_features, frequency_data = _extract_ecg_frequency_domain_features(
            time_domain["ibi"],
            label_frequency_analysis=self.label_frequency_analysis,
            show_plots=show_plots,
            save_plots=save_plots,
            output_dir=output_dir,
            output_label=output_label,
            stft_context=stft_context,
        )

        features = dict(time_domain["features"])
        features.update(frequency_features)

        intermediates = {}
        if self.return_intermediates:
            intermediates = {
                "filtered_signal": time_domain["filtered_signal"],
                "r_peak_indices": time_domain["r_peak_indices"],
                "r_peak_times": signal.time[time_domain["r_peak_indices"]],
                "r_peak_amplitudes": time_domain["r_peak_amplitudes"],
                "accepted_r_peak_indices": time_domain["accepted_r_peak_indices"],
                "accepted_r_peak_times": signal.time[
                    time_domain["accepted_r_peak_indices"]
                ],
                "ibi": time_domain["ibi"],
                "ibi_times": time_domain["ibi_times"],
                "bpm": time_domain["bpm"],
                "bpm_times": time_domain["ibi_times"],
            }
            if self.label_frequency_analysis == 1:
                intermediates.update(
                    {
                        "stft_matrix": frequency_data["matrix"],
                        "stft_power": frequency_data["power"],
                        "stft_frequencies": frequency_data["frequencies"],
                        "stft_times": frequency_data["time"],
                    }
                )

        return ProcessingResult(
            modality="ecg",
            processor_name="process_ecg",
            features=features,
            signal_labels=[signal.label or "ecg"],
            metadata={
                "source": "process_ecg",
                "sampling_frequency": signal.sampling_frequency,
                "n_samples": signal.n_samples,
                "start_time_seconds": float(signal.time[0]),
                "end_time_seconds": float(signal.time[-1]),
                "frequency_analysis_method": self.label_frequency_analysis,
                "return_intermediates": self.return_intermediates,
            },
            intermediates=intermediates,
        )


def detect_ecg_r_peaks(
    signal: Signal,
    *,
    show_plots: bool = False,
    save_plots: bool = False,
    output_dir: str | Path | None = None,
) -> dict[str, np.ndarray]:
    """Detect R-peaks with an independently implemented Pan-Tompkins-style method.

    Methodological reference: Pan and Tompkins, IEEE Transactions on
    Biomedical Engineering, 1985, DOI 10.1109/TBME.1985.325532.
    """

    return ECGProcessor().detect_r_peaks(
        signal,
        show_plots=show_plots,
        save_plots=save_plots,
        output_dir=output_dir,
    )


def process_ecg(
    signal: Signal,
    *,
    label_frequency_analysis: int = 1,
    return_intermediates: bool = False,
    show_plots: bool = False,
    save_plots: bool = False,
    output_dir: str | Path | None = None,
    output_label: str = "ecg",
) -> ProcessingResult:
    """Extract MATLAB-compatible ECG and HRV features from one signal."""

    return ECGProcessor(
        label_frequency_analysis=label_frequency_analysis,
        return_intermediates=return_intermediates,
    ).process(
        signal,
        show_plots=show_plots,
        save_plots=save_plots,
        output_dir=output_dir,
        output_label=output_label,
    )


def _detect_ecg_r_peaks_impl(
    signal: Signal,
    *,
    show_plots: bool = False,
    save_plots: bool = False,
    output_dir: str | Path | None = None,
) -> dict[str, np.ndarray]:
    if not isinstance(signal, Signal):
        raise TypeError("signal must be a Signal instance")
    if not isinstance(show_plots, bool):
        raise TypeError("show_plots must be a boolean")
    if not isinstance(save_plots, bool):
        raise TypeError("save_plots must be a boolean")
    if save_plots and output_dir is None:
        raise ValueError("output_dir must be provided when save_plots=True")

    baseline_corrected_ecg = signal.samples - np.mean(signal.samples)
    filtered_ecg = _filter_ecg(baseline_corrected_ecg, signal.sampling_frequency)
    stages = _pan_tompkins_stages(filtered_ecg, signal.sampling_frequency)
    detections = _adaptive_threshold_detection(
        filtered_ecg,
        stages["integrated_ecg"],
        signal.sampling_frequency,
    )

    r_peak_indexes = detections["r_peak_indexes"]
    if r_peak_indexes.size == 0:
        raise ValueError("no R-peaks were detected")

    result = {
        "r_peak_amplitudes": filtered_ecg[r_peak_indexes],
        "r_peak_indexes": r_peak_indexes.astype(int),
        "filtered_ecg": filtered_ecg,
    }

    if show_plots or save_plots:
        _plot_ecg_diagnostics(
            signal,
            baseline_corrected_ecg,
            filtered_ecg,
            stages,
            detections,
            show_plots=show_plots,
            save_plots=save_plots,
            output_dir=output_dir,
        )

    return result


def _validate_ecg_frequency_options(
    *,
    label_frequency_analysis: int,
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
    if label_frequency_analysis not in (1, 2):
        raise ValueError("label_frequency_analysis must be 1 or 2")
    if save_plots and output_dir is None:
        raise ValueError("output_dir must be provided when save_plots=True")


def _validate_ecg_plot_options(
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


def _validate_ecg_stft_duration(signal: Signal) -> None:
    if signal.n_samples / signal.sampling_frequency < 60:
        raise ValueError(
            "STFT HRV analysis requires the input signal to be at least 60 seconds"
        )


def _build_ecg_stft_context(signal: Signal) -> dict[str, np.ndarray]:
    detections = detect_ecg_r_peaks(signal)
    r_peak_indexes = detections["r_peak_indexes"].astype(int)
    if r_peak_indexes.size < 3:
        raise ValueError("at least three R-peaks are required for global STFT analysis")

    ibi = np.diff(r_peak_indexes) / signal.sampling_frequency
    ibi_time = (signal.time[r_peak_indexes[:-1]] + signal.time[r_peak_indexes[1:]]) / 2

    valid = np.isfinite(ibi) & (ibi >= 60 / 160) & (ibi <= 60 / 40)
    if np.any(valid):
        median_ibi = np.median(ibi[valid])
        valid = valid & (ibi >= 0.5 * median_ibi) & (ibi <= 1.5 * median_ibi)

    ibi = ibi[valid]
    ibi_time = ibi_time[valid]
    if ibi.size < 3:
        raise ValueError("too few valid global IBIs remain for STFT analysis")

    return {
        "ibi": ibi,
        "ibi_time": ibi_time,
        "time_range": np.asarray([signal.time[0], signal.time[-1]], dtype=float),
        "analysis_window": np.asarray([signal.time[0], signal.time[-1]], dtype=float),
    }


def _filter_ecg(ecg: np.ndarray, sampling_frequency: float) -> np.ndarray:
    nyquist = sampling_frequency / 2
    if 16 >= nyquist:
        raise ValueError("sampling_frequency is too low for 16 Hz low-pass filtering")

    try:
        high_b, high_a = butter(2, 4 / nyquist, btype="high")
        high_passed = filtfilt(high_b, high_a, ecg)
        low_b, low_a = butter(2, 16 / nyquist, btype="low")
        return filtfilt(low_b, low_a, high_passed)
    except ValueError as exc:
        raise ValueError("signal is too short for filtering or initialization") from exc


def _pan_tompkins_stages(filtered_ecg: np.ndarray, sampling_frequency: float) -> dict[str, np.ndarray]:
    if filtered_ecg.size < int(2 * sampling_frequency):
        raise ValueError("signal is too short for filtering or initialization")

    differentiated_ecg = np.diff(filtered_ecg)
    max_derivative = np.max(differentiated_ecg)
    if max_derivative <= 0 or not np.isfinite(max_derivative):
        raise ValueError("no R-peaks were detected")

    normalized_derivative = differentiated_ecg / max_derivative
    squared_ecg = normalized_derivative**2
    integration_window_samples = int(round(0.15 * sampling_frequency))
    if integration_window_samples < 1:
        raise ValueError("signal is too short for filtering or initialization")
    integrated_ecg = np.convolve(
        squared_ecg,
        np.ones(integration_window_samples) / integration_window_samples,
        mode="full",
    )
    return {
        "differentiated_ecg": differentiated_ecg,
        "normalized_derivative": normalized_derivative,
        "squared_ecg": squared_ecg,
        "integrated_ecg": integrated_ecg,
        "integration_window_samples": np.asarray(integration_window_samples),
    }


def _adaptive_threshold_detection(
    filtered_ecg: np.ndarray,
    integrated_ecg: np.ndarray,
    sampling_frequency: float,
) -> dict[str, np.ndarray]:
    min_peak_distance = int(round(0.2 * sampling_frequency))
    candidate_amplitudes, candidate_indexes = find_peaks(
        integrated_ecg,
        distance=min_peak_distance,
    )
    candidate_values = integrated_ecg[candidate_amplitudes]
    if candidate_values.size == 0:
        raise ValueError("no R-peaks were detected")

    mean_candidate_value = np.mean(candidate_values)
    keep = candidate_values <= mean_candidate_value * 10
    locs = candidate_amplitudes[keep]
    pks = candidate_values[keep]
    if pks.size == 0:
        raise ValueError("no R-peaks were detected")

    initialization_samples = int(2 * sampling_frequency)
    if filtered_ecg.size < initialization_samples or integrated_ecg.size < initialization_samples:
        raise ValueError("signal is too short for filtering or initialization")

    integration_window_samples = int(round(0.15 * sampling_frequency))
    thr_sig = np.max(integrated_ecg[:initialization_samples]) * 0.25
    thr_noise = np.mean(integrated_ecg[:initialization_samples]) * 0.5
    sig_lev = thr_sig
    noise_lev = thr_noise
    thr_sig_filtered = np.max(filtered_ecg[:initialization_samples]) * 0.25
    thr_noise_filtered = np.mean(filtered_ecg[:initialization_samples]) * 0.5
    sig_lev_filtered = thr_sig_filtered
    noise_lev_filtered = thr_noise_filtered

    qrs_c: list[float] = []
    qrs_i: list[int] = []
    r_peak_indexes: list[int] = []
    r_peak_amplitudes: list[float] = []
    noise_values: list[float] = []
    noise_indexes: list[int] = []
    sig_buffer: list[float] = []
    noise_buffer: list[float] = []
    threshold_buffer: list[float] = []
    sig_filtered_buffer: list[float] = []
    noise_filtered_buffer: list[float] = []
    threshold_filtered_buffer: list[float] = []
    selected_rr_mean = 0.0
    mean_rr = 0.0

    for pks_i, loc_i in zip(pks, locs):
        y_i, x_i, ser_back = _locate_filtered_peak(
            filtered_ecg,
            int(loc_i),
            integration_window_samples,
        )
        skip = False

        if len(qrs_i) >= 9:
            diff_rr = np.diff(qrs_i[-9:])
            mean_rr = float(np.mean(diff_rr))
            latest_rr = qrs_i[-1] - qrs_i[-2]
            if latest_rr <= 0.92 * mean_rr or latest_rr >= 1.16 * mean_rr:
                thr_sig *= 0.5
                thr_sig_filtered *= 0.5
            else:
                selected_rr_mean = mean_rr

        test_m = selected_rr_mean or mean_rr
        if test_m and qrs_i and (loc_i - qrs_i[-1]) >= round(1.66 * test_m):
            search_start = qrs_i[-1] + int(round(0.2 * sampling_frequency))
            search_end = int(loc_i) - int(round(0.2 * sampling_frequency))
            if search_start < search_end:
                search_segment = integrated_ecg[search_start:search_end]
                if search_segment.size:
                    temp_offset = int(np.argmax(search_segment))
                    pks_temp = float(search_segment[temp_offset])
                    locs_temp = search_start + temp_offset
                    if pks_temp > thr_noise:
                        qrs_c.append(pks_temp)
                        qrs_i.append(locs_temp)
                        y_i_t, x_i_t, _ = _locate_filtered_peak(
                            filtered_ecg,
                            locs_temp,
                            integration_window_samples,
                        )
                        if y_i_t > thr_noise_filtered:
                            raw_index = max(
                                0,
                                locs_temp - integration_window_samples + x_i_t,
                            )
                            r_peak_indexes.append(raw_index)
                            r_peak_amplitudes.append(float(y_i_t))
                            sig_lev_filtered = 0.25 * y_i_t + 0.75 * sig_lev_filtered
                        sig_lev = 0.25 * pks_temp + 0.75 * sig_lev

        if pks_i >= thr_sig:
            if len(qrs_i) >= 3 and (loc_i - qrs_i[-1]) <= round(0.36 * sampling_frequency):
                slope_1 = _mean_recent_slope(integrated_ecg, int(loc_i), sampling_frequency)
                slope_2 = _mean_recent_slope(integrated_ecg, qrs_i[-1], sampling_frequency)
                if abs(slope_1) <= abs(0.5 * slope_2):
                    noise_values.append(float(pks_i))
                    noise_indexes.append(int(loc_i))
                    skip = True
                    noise_lev_filtered = 0.125 * y_i + 0.875 * noise_lev_filtered
                    noise_lev = 0.125 * pks_i + 0.875 * noise_lev

            if not skip:
                qrs_c.append(float(pks_i))
                qrs_i.append(int(loc_i))
                if y_i >= thr_sig_filtered:
                    if ser_back:
                        raw_index = x_i
                    else:
                        raw_index = int(loc_i) - integration_window_samples + x_i
                    raw_index = int(np.clip(raw_index, 0, filtered_ecg.size - 1))
                    r_peak_indexes.append(raw_index)
                    r_peak_amplitudes.append(float(y_i))
                    sig_lev_filtered = 0.125 * y_i + 0.875 * sig_lev_filtered
                sig_lev = 0.125 * pks_i + 0.875 * sig_lev
        elif thr_noise <= pks_i < thr_sig:
            noise_lev_filtered = 0.125 * y_i + 0.875 * noise_lev_filtered
            noise_lev = 0.125 * pks_i + 0.875 * noise_lev
        else:
            noise_values.append(float(pks_i))
            noise_indexes.append(int(loc_i))
            noise_lev_filtered = 0.125 * y_i + 0.875 * noise_lev_filtered
            noise_lev = 0.125 * pks_i + 0.875 * noise_lev

        if noise_lev != 0 or sig_lev != 0:
            thr_sig = noise_lev + 0.25 * abs(sig_lev - noise_lev)
            thr_noise = 0.5 * thr_sig
        if noise_lev_filtered != 0 or sig_lev_filtered != 0:
            thr_sig_filtered = noise_lev_filtered + 0.25 * abs(
                sig_lev_filtered - noise_lev_filtered
            )
            thr_noise_filtered = 0.5 * thr_sig_filtered

        sig_buffer.append(float(sig_lev))
        noise_buffer.append(float(noise_lev))
        threshold_buffer.append(float(thr_sig))
        sig_filtered_buffer.append(float(sig_lev_filtered))
        noise_filtered_buffer.append(float(noise_lev_filtered))
        threshold_filtered_buffer.append(float(thr_sig_filtered))

    return {
        "candidate_indexes": locs.astype(int),
        "candidate_values": pks.astype(float),
        "qrs_integrated_indexes": np.asarray(qrs_i, dtype=int),
        "qrs_integrated_amplitudes": np.asarray(qrs_c, dtype=float),
        "r_peak_indexes": _unique_sorted_indexes(r_peak_indexes, filtered_ecg.size),
        "r_peak_amplitudes": np.asarray(r_peak_amplitudes, dtype=float),
        "noise_indexes": np.asarray(noise_indexes, dtype=int),
        "noise_values": np.asarray(noise_values, dtype=float),
        "signal_level_buffer": np.asarray(sig_buffer, dtype=float),
        "noise_level_buffer": np.asarray(noise_buffer, dtype=float),
        "threshold_buffer": np.asarray(threshold_buffer, dtype=float),
        "signal_level_filtered_buffer": np.asarray(sig_filtered_buffer, dtype=float),
        "noise_level_filtered_buffer": np.asarray(noise_filtered_buffer, dtype=float),
        "threshold_filtered_buffer": np.asarray(threshold_filtered_buffer, dtype=float),
    }


def _locate_filtered_peak(
    filtered_ecg: np.ndarray,
    integrated_index: int,
    integration_window_samples: int,
) -> tuple[float, int, bool]:
    ser_back = False
    if integrated_index - integration_window_samples >= 0 and integrated_index < filtered_ecg.size:
        start = integrated_index - integration_window_samples
        stop = integrated_index + 1
    elif integrated_index < filtered_ecg.size:
        start = 0
        stop = integrated_index + 1
        ser_back = True
    else:
        start = max(0, integrated_index - integration_window_samples)
        stop = filtered_ecg.size

    segment = filtered_ecg[start:stop]
    if segment.size == 0:
        raise ValueError("no R-peaks were detected")
    offset = int(np.argmax(segment))
    return float(segment[offset]), offset if ser_back else offset, ser_back


def _mean_recent_slope(
    integrated_ecg: np.ndarray,
    index: int,
    sampling_frequency: float,
) -> float:
    window = int(round(0.075 * sampling_frequency))
    start = max(0, index - window)
    segment = integrated_ecg[start : index + 1]
    if segment.size < 2:
        return 0.0
    return float(np.mean(np.diff(segment)))


def _unique_sorted_indexes(indexes: list[int], signal_length: int) -> np.ndarray:
    if not indexes:
        return np.asarray([], dtype=int)
    return np.asarray(sorted(set(int(np.clip(index, 0, signal_length - 1)) for index in indexes)), dtype=int)


def _plot_ecg_diagnostics(
    signal: Signal,
    baseline_corrected_ecg: np.ndarray,
    filtered_ecg: np.ndarray,
    stages: dict[str, np.ndarray],
    detections: dict[str, np.ndarray],
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
        _plot_pan_tompkins_stages(signal, baseline_corrected_ecg, filtered_ecg, stages),
        _plot_qrs_detection_summary(signal, baseline_corrected_ecg, filtered_ecg, stages, detections),
    ]

    if save_plots and output_path is not None:
        figures[0].savefig(output_path / "ecg_pan_tompkins_stages.png", dpi=150)
        figures[1].savefig(output_path / "ecg_qrs_detection_summary.png", dpi=150)

    if show_plots:
        plt.show()
    else:
        for figure in figures:
            plt.close(figure)


def _plot_pan_tompkins_stages(
    signal: Signal,
    baseline_corrected_ecg: np.ndarray,
    filtered_ecg: np.ndarray,
    stages: dict[str, np.ndarray],
):
    figure, axes = plt.subplots(4, 1, sharex=False)
    figure.suptitle("Pan-Tompkins Algorithm")
    axes[0].plot(signal.time[1:], baseline_corrected_ecg[1:], linewidth=1.5)
    axes[0].plot(signal.time[1:], filtered_ecg[1:], linewidth=1.5)
    axes[0].set_ylabel("[mv]")
    axes[0].legend(["Raw ECG", "Filtered ECG"])
    axes[1].plot(signal.time[1:], stages["normalized_derivative"], linewidth=1.5)
    axes[1].set_ylabel("[-]")
    axes[1].legend(["Diff and normalized ECG"])
    axes[2].plot(signal.time[1:], stages["squared_ecg"], linewidth=1.5)
    axes[2].set_ylabel("[-]")
    axes[2].legend(["Squared ECG"])
    integrated = stages["integrated_ecg"][: stages["squared_ecg"].size]
    axes[3].plot(signal.time[1:], integrated, linewidth=1.5)
    axes[3].set_xlabel("time [s]")
    axes[3].set_ylabel("[-]")
    axes[3].legend(["Integrated ECG"])
    for axis in axes:
        axis.grid(True)
    return figure


def _plot_qrs_detection_summary(
    signal: Signal,
    baseline_corrected_ecg: np.ndarray,
    filtered_ecg: np.ndarray,
    stages: dict[str, np.ndarray],
    detections: dict[str, np.ndarray],
):
    figure, axes = plt.subplots(3, 1, sharex=False)
    r_peak_indexes = detections["r_peak_indexes"]
    axes[0].plot(filtered_ecg)
    axes[0].set_title("QRS on Filtered Signal")
    axes[0].scatter(r_peak_indexes, filtered_ecg[r_peak_indexes], c="m")
    candidate_indexes = detections["candidate_indexes"]
    axes[1].plot(stages["integrated_ecg"])
    axes[1].set_title("QRS on MVI signal and Noise level(black),Signal Level (red) and Adaptive Threshold(green)")
    axes[1].scatter(detections["qrs_integrated_indexes"], detections["qrs_integrated_amplitudes"], c="m")
    for values, style, color in [
        (detections["noise_level_buffer"], "--", "k"),
        (detections["signal_level_buffer"], "-.", "r"),
        (detections["threshold_buffer"], "-.", "g"),
    ]:
        if values.size == candidate_indexes.size:
            axes[1].plot(candidate_indexes, values, linestyle=style, color=color, linewidth=2)
    axes[2].plot(baseline_corrected_ecg)
    axes[2].set_title("Pulse train of the found QRS on ECG signal")
    if r_peak_indexes.size:
        axes[2].vlines(
            r_peak_indexes,
            np.min(baseline_corrected_ecg) / 2,
            np.max(baseline_corrected_ecg) / 2,
            colors="r",
            linestyles="-.",
            linewidth=2.5,
        )
    for axis in axes:
        axis.grid(True)
    return figure


def _plot_ecg_bpm_trend(
    bpm: np.ndarray,
    *,
    show_plots: bool = False,
    save_plots: bool = False,
    output_dir: str | Path | None = None,
    output_label: str = "ecg",
) -> None:
    _validate_ecg_plot_options(
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
    bpm_values = np.asarray(bpm, dtype=float).reshape(-1)
    axes.plot(np.arange(1, bpm_values.size + 1), bpm_values, linewidth=1.5)
    axes.set_xlabel("Accepted IBI sequence")
    axes.set_ylabel("BPM")
    axes.set_title("BPM trend")
    axes.grid(True)

    if save_plots and output_path is not None:
        figure.savefig(output_path / f"{output_label}_BPM_trend.png", dpi=150)

    if show_plots:
        plt.show()
    else:
        plt.close(figure)


def _extract_ecg_time_domain_features(
    signal: Signal,
) -> dict[str, object]:
    """Extract MATLAB-compatible ECG time-domain features."""

    if not isinstance(signal, Signal):
        raise TypeError("signal must be a Signal instance")
    detections = detect_ecg_r_peaks(signal)
    r_peak_indexes = detections["r_peak_indexes"]
    r_peak_amplitudes = detections["r_peak_amplitudes"]
    if r_peak_indexes.size < 3:
        raise ValueError("too few R-peaks were detected in the signal")

    accepted_peak_indexes = _screen_r_peak_amplitudes(
        r_peak_indexes,
        r_peak_amplitudes,
    )
    if accepted_peak_indexes.size < 3:
        raise ValueError("too few accepted R-peaks remain after amplitude screening")

    ibi_values, control_ibi_flag = _screen_ibi_values(
        accepted_peak_indexes,
        signal.sampling_frequency,
    )
    if ibi_values.size < 3:
        raise ValueError("too few valid IBIs remain after IBI screening")

    bpm = _replace_bpm_outliers(60 / ibi_values)
    features = _ecg_time_domain_feature_values(
        bpm,
        ibi_values,
        control_ibi_flag,
    )
    accepted_interval_starts = accepted_peak_indexes[:-1][control_ibi_flag != 0]
    ibi_times = signal.time[accepted_interval_starts[: ibi_values.size]]

    return {
        "features": features,
        "filtered_signal": detections["filtered_ecg"],
        "r_peak_indices": r_peak_indexes,
        "r_peak_amplitudes": r_peak_amplitudes,
        "accepted_r_peak_indices": accepted_peak_indexes,
        "ibi": ibi_values,
        "ibi_times": ibi_times,
        "bpm": bpm,
    }


def _screen_r_peak_amplitudes(
    r_peak_indexes: np.ndarray,
    r_peak_amplitudes: np.ndarray,
) -> np.ndarray:
    default_l_window_amp = 10
    l_window_amp = default_l_window_amp
    upper_threshold = np.percentile(r_peak_amplitudes, 99)
    lower_threshold = np.percentile(r_peak_amplitudes, 1)
    initial_mask = (r_peak_amplitudes < upper_threshold) & (
        r_peak_amplitudes > lower_threshold
    )
    initial_indexes = r_peak_indexes[initial_mask]
    initial_amplitudes = r_peak_amplitudes[initial_mask]

    if initial_indexes.size < l_window_amp:
        l_window_amp = round(l_window_amp / 2)
    if initial_indexes.size < l_window_amp or l_window_amp < 1:
        raise ValueError("too few accepted R-peaks remain after amplitude screening")

    accepted_indexes = list(initial_indexes[:l_window_amp])
    accepted_amplitudes = list(initial_amplitudes[:l_window_amp])
    count_matches = np.flatnonzero(r_peak_indexes == accepted_indexes[-1])
    if count_matches.size == 0:
        raise ValueError("too few accepted R-peaks remain after amplitude screening")
    count = int(count_matches[0])

    for index in range(count + 1, r_peak_indexes.size):
        recent_mean = np.mean(accepted_amplitudes[-l_window_amp:])
        if 0.6 * recent_mean <= r_peak_amplitudes[index] <= 1.5 * recent_mean:
            accepted_indexes.append(int(r_peak_indexes[index]))
            accepted_amplitudes.append(float(r_peak_amplitudes[index]))

    return np.asarray(accepted_indexes, dtype=int)


def _screen_ibi_values(
    accepted_peak_indexes: np.ndarray,
    sampling_frequency: float,
) -> tuple[np.ndarray, np.ndarray]:
    l_window_time = 5
    diff_t_r_peaks = np.diff(accepted_peak_indexes) / sampling_frequency
    if diff_t_r_peaks.size < l_window_time:
        raise ValueError("too few valid IBIs remain after IBI screening")

    control_ibi_flag = np.zeros(diff_t_r_peaks.size, dtype=int)
    mean_diff_t_peaks = np.mean(diff_t_r_peaks)
    preliminary_mask = (diff_t_r_peaks < 1.5 * mean_diff_t_peaks) & (
        diff_t_r_peaks > 0.5 * mean_diff_t_peaks
    )
    ibi_preliminary = diff_t_r_peaks[preliminary_mask]
    indexes_ibi_preliminary = accepted_peak_indexes[:-1][preliminary_mask]
    if ibi_preliminary.size < l_window_time:
        raise ValueError("too few valid IBIs remain after IBI screening")

    ibi_def_phase = list(ibi_preliminary[:l_window_time])
    indexes_ibi_def = indexes_ibi_preliminary[:l_window_time]
    for accepted_index in indexes_ibi_def:
        matches = np.flatnonzero(accepted_peak_indexes == accepted_index)
        if matches.size:
            flag_index = int(matches[0])
            if flag_index < control_ibi_flag.size:
                control_ibi_flag[flag_index] = 1

    count_matches = np.flatnonzero(accepted_peak_indexes == indexes_ibi_def[-1])
    if count_matches.size == 0:
        raise ValueError("too few valid IBIs remain after IBI screening")
    count_ibi = int(count_matches[0])

    for index in range(count_ibi + 1, diff_t_r_peaks.size):
        recent_mean = np.mean(ibi_def_phase[-l_window_time:])
        if 0.7 * recent_mean <= diff_t_r_peaks[index] <= 1.5 * recent_mean:
            ibi_def_phase.append(float(diff_t_r_peaks[index]))
            control_ibi_flag[index] = 1

    return np.asarray(ibi_def_phase, dtype=float), control_ibi_flag


def _replace_bpm_outliers(raw_bpm: np.ndarray) -> np.ndarray:
    bpm = raw_bpm.copy()
    outliers = np.flatnonzero((raw_bpm < 40) | (raw_bpm > 160))
    if outliers.size == 0:
        return bpm
    if raw_bpm.size < 3:
        raise ValueError("too few valid IBIs remain after IBI screening")

    for index in outliers:
        if index == 0:
            bpm[index] = np.mean(raw_bpm[1:3])
        elif index == raw_bpm.size - 1:
            bpm[index] = np.mean(raw_bpm[-3:-1])
        else:
            bpm[index] = (raw_bpm[index - 1] + raw_bpm[index + 1]) / 2
    return bpm


def _ecg_time_domain_feature_values(
    bpm: np.ndarray,
    ibi_values: np.ndarray,
    control_ibi_flag: np.ndarray,
) -> dict[str, float]:
    diff_ibi_phase = []
    ibi_counter = 0
    for index in range(max(control_ibi_flag.size - 1, 0)):
        if control_ibi_flag[index] != 0 and control_ibi_flag[index + 1] != 0:
            if ibi_counter + 1 < ibi_values.size:
                diff_ibi_phase.append(abs(ibi_values[ibi_counter + 1] - ibi_values[ibi_counter]))
                ibi_counter += 1
            else:
                diff_ibi_phase.append(0.0)
        else:
            diff_ibi_phase.append(0.0)

    length_adjacent_intervals = _adjacent_interval_lengths(np.asarray(diff_ibi_phase))
    mean_number_ibi_2_minutes = round(120 * ibi_values.size) / np.sum(ibi_values)
    if length_adjacent_intervals.size == 0:
        pnn50 = 0.0
        pnn50_max = 0.0
    else:
        pnn50 = float(np.sum(length_adjacent_intervals) / ibi_values.size)
        pnn50_max = float(np.max(length_adjacent_intervals) / mean_number_ibi_2_minutes)

    return {
        "ecg_bpm_mean": float(np.mean(bpm)),
        "ecg_bpm_std": float(np.std(bpm, ddof=1)),
        "ecg_bpm_median": float(np.median(bpm)),
        "ecg_pnn50": pnn50,
        "ecg_pnn50_max": pnn50_max,
    }


def _adjacent_interval_lengths(diff_ibi_phase: np.ndarray) -> np.ndarray:
    indexes = np.flatnonzero(diff_ibi_phase >= 0.05)
    if indexes.size < 2:
        return np.asarray([], dtype=float)

    lengths: list[int] = []
    current_index = 0
    for index in range(indexes.size - 1):
        if indexes[index + 1] - indexes[index] == 1 and len(lengths) < current_index + 1:
            lengths.append(2)
        elif indexes[index + 1] - indexes[index] == 1:
            lengths[current_index] += 1
        elif len(lengths) == current_index + 1:
            current_index += 1
    return np.asarray(lengths, dtype=float)


def _extract_ecg_frequency_domain_features(
    ibi: np.ndarray,
    *,
    label_frequency_analysis: int = 2,
    show_plots: bool = False,
    save_plots: bool = False,
    output_dir: str | Path | None = None,
    output_label: str = "ecg",
    stft_context: dict[str, object] | None = None,
) -> tuple[dict[str, float], dict[str, object]]:
    """Extract MATLAB-compatible ECG frequency-domain HRV features."""

    if not isinstance(show_plots, bool):
        raise TypeError("show_plots must be a boolean")
    if not isinstance(save_plots, bool):
        raise TypeError("save_plots must be a boolean")
    if not isinstance(output_label, str):
        raise TypeError("output_label must be a string")
    if save_plots and output_dir is None:
        raise ValueError("output_dir must be provided when save_plots=True")
    if label_frequency_analysis not in (1, 2):
        raise ValueError("label_frequency_analysis must be 1 or 2")

    fs_hrv = 4
    lf_band = (0.04, 0.15)
    hf_band = (0.15, 0.40)
    total_band = (0.04, 0.40)
    winlen_seconds = 60

    ibi_tot = _validate_ibi(ibi)

    beat_time = np.cumsum(ibi_tot)
    analysis_edges_seconds = np.asarray([0.0, beat_time[-1]], dtype=float)
    analysis_window_seconds = analysis_edges_seconds.reshape(1, 2)

    uniform_time, ibi_uniform = _uniform_ibi_tachogram_from_ibi(
        ibi_tot,
        fs_hrv,
    )

    if label_frequency_analysis == 1 and stft_context is not None:
        uniform_time, ibi_uniform, analysis_window_seconds = _stft_context_to_tachogram(
            stft_context,
            fs_hrv,
        )

    if label_frequency_analysis == 1:
        analysis_values, plot_data = _frequency_features_stft(
            uniform_time,
            ibi_uniform,
            analysis_window_seconds,
            fs_hrv,
            winlen_seconds,
            lf_band,
            hf_band,
            total_band,
        )
    else:
        analysis_values, plot_data = _frequency_features_psd(
            uniform_time,
            ibi_uniform,
            analysis_edges_seconds,
            fs_hrv,
            winlen_seconds,
            lf_band,
            hf_band,
            total_band,
        )

    features: dict[str, float] = {}
    for feature_name in _FREQUENCY_DOMAIN_FEATURE_NAMES:
        features[feature_name] = float(analysis_values[feature_name][0])

    if show_plots or save_plots:
        _plot_ecg_frequency_diagnostics(
            label_frequency_analysis,
            plot_data,
            show_plots=show_plots,
            save_plots=save_plots,
            output_dir=output_dir,
            output_label=output_label,
        )

    return features, plot_data


def _validate_ibi(ibi: np.ndarray) -> np.ndarray:
    array = np.asarray(ibi, dtype=float).reshape(-1)
    if array.size == 0:
        raise ValueError("IBI array must be non-empty")
    if not np.all(np.isfinite(array)):
        raise ValueError("IBI array must contain finite values")
    if np.any(array <= 0):
        raise ValueError("IBI array must contain positive values")
    return array


def _uniform_ibi_tachogram_from_ibi(
    ibi_tot: np.ndarray,
    fs_hrv: float,
) -> tuple[np.ndarray, np.ndarray]:
    if ibi_tot.size < 4:
        raise ValueError("too few IBI samples are available for frequency analysis")

    beat_time = np.cumsum(ibi_tot)
    ibi_time = np.concatenate(([0.0], beat_time))
    ibi_values = np.concatenate(([ibi_tot[0]], ibi_tot))
    uniform_time = np.arange(0, beat_time[-1] + 1 / fs_hrv / 2, 1 / fs_hrv)
    interpolator = PchipInterpolator(ibi_time, ibi_values)
    ibi_uniform = interpolator(uniform_time)
    ibi_uniform = ibi_uniform - np.mean(ibi_uniform)
    return uniform_time, ibi_uniform


def _stft_context_to_tachogram(
    stft_context: dict[str, object],
    fs_hrv: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    required_keys = {"ibi", "ibi_time", "time_range", "analysis_window"}
    if not isinstance(stft_context, dict) or not required_keys.issubset(stft_context):
        raise ValueError("stft_context is incomplete")

    global_ibi = np.asarray(stft_context["ibi"], dtype=float).reshape(-1)
    global_ibi_time = np.asarray(stft_context["ibi_time"], dtype=float).reshape(-1)
    time_range = np.asarray(stft_context["time_range"], dtype=float).reshape(-1)
    analysis_window = np.asarray(stft_context["analysis_window"], dtype=float).reshape(-1)

    if (
        global_ibi.size == 0
        or global_ibi_time.size != global_ibi.size
        or time_range.size != 2
        or analysis_window.size != 2
    ):
        raise ValueError("stft_context is incomplete")
    if np.any(~np.isfinite(global_ibi)) or np.any(global_ibi <= 0):
        raise ValueError("stft_context IBI values must be positive and finite")

    sort_index = np.argsort(global_ibi_time)
    global_ibi_time = global_ibi_time[sort_index]
    global_ibi = global_ibi[sort_index]

    interpolation_time = np.concatenate(([time_range[0]], global_ibi_time, [time_range[1]]))
    interpolation_values = np.concatenate(([global_ibi[0]], global_ibi, [global_ibi[-1]]))
    unique_time, unique_index = np.unique(interpolation_time, return_index=True)
    interpolation_values = interpolation_values[unique_index]

    step = 1 / fs_hrv
    sample_count = int(np.floor((time_range[1] - time_range[0]) / step)) + 1
    uniform_time = time_range[0] + np.arange(sample_count) * step
    interpolator = PchipInterpolator(unique_time, interpolation_values)
    ibi_uniform = interpolator(uniform_time)
    ibi_uniform = ibi_uniform - np.mean(ibi_uniform)
    return uniform_time, ibi_uniform, analysis_window.reshape(1, 2)


def _frequency_features_psd(
    uniform_time: np.ndarray,
    ibi_uniform: np.ndarray,
    analysis_edges_seconds: np.ndarray,
    fs_hrv: float,
    winlen_seconds: float,
    lf_band: tuple[float, float],
    hf_band: tuple[float, float],
    total_band: tuple[float, float],
) -> tuple[dict[str, np.ndarray], dict[str, object]]:
    winlen = max(16, min(round(winlen_seconds * fs_hrv), ibi_uniform.size))
    nfft = max(256, _next_power_of_two(winlen))
    analysis_values = _empty_frequency_values(1)
    analysis_mask = (uniform_time >= analysis_edges_seconds[0]) & (
        uniform_time <= analysis_edges_seconds[1]
    )
    x_analysis = ibi_uniform[analysis_mask]
    if x_analysis.size < 16:
        raise ValueError("too few IBI samples are available for frequency analysis")

    pxx, frequencies = _local_fft_psd(x_analysis, fs_hrv, nfft)
    lf_power = _local_band_power(pxx, frequencies, lf_band)
    hf_power = _local_band_power(pxx, frequencies, hf_band)
    total_power = _local_band_power(pxx, frequencies, total_band)

    analysis_values["ecg_plf_mean"][0] = lf_power / max(total_power, np.finfo(float).eps)
    analysis_values["ecg_phf_mean"][0] = hf_power / max(total_power, np.finfo(float).eps)
    analysis_values["ecg_lf_hf_mean"][0] = lf_power / max(hf_power, np.finfo(float).eps)
    analysis_values["ecg_plf_std"][0] = 0.0
    analysis_values["ecg_phf_std"][0] = 0.0
    analysis_values["ecg_lf_hf_std"][0] = 0.0

    _validate_frequency_values(analysis_values)
    return analysis_values, {"mode": 2, "lf_hf": analysis_values["ecg_lf_hf_mean"]}


def _frequency_features_stft(
    uniform_time: np.ndarray,
    ibi_uniform: np.ndarray,
    analysis_window_seconds: np.ndarray,
    fs_hrv: float,
    winlen_seconds: float,
    lf_band: tuple[float, float],
    hf_band: tuple[float, float],
    total_band: tuple[float, float],
) -> tuple[dict[str, np.ndarray], dict[str, object]]:
    winlen = round(winlen_seconds * fs_hrv)
    if ibi_uniform.size < winlen:
        raise ValueError("resampled tachogram is shorter than the 60-second STFT window")
    nfft = max(256, _next_power_of_two(winlen))
    stft_matrix, t_stft, frequencies = _centered_stft(
        ibi_uniform,
        nfft,
        fs_hrv,
        winlen,
    )
    t_stft = t_stft + uniform_time[0]
    power = np.abs(stft_matrix) ** 2

    lf_index = (frequencies >= lf_band[0]) & (frequencies <= lf_band[1])
    hf_index = (frequencies > hf_band[0]) & (frequencies <= hf_band[1])
    total_index = (frequencies >= total_band[0]) & (frequencies <= total_band[1])

    lf_power = np.sum(power[lf_index, :], axis=0)
    hf_power = np.sum(power[hf_index, :], axis=0)
    total_power = np.sum(power[total_index, :], axis=0)
    plf = lf_power / np.maximum(total_power, np.finfo(float).eps)
    phf = hf_power / np.maximum(total_power, np.finfo(float).eps)
    lf_hf = lf_power / np.maximum(hf_power, np.finfo(float).eps)

    analysis_values = _empty_frequency_values(1)
    analysis_mask = (t_stft >= analysis_window_seconds[0, 0]) & (
        t_stft <= analysis_window_seconds[0, 1]
    )
    if not np.any(analysis_mask):
        raise ValueError("input signal has no STFT samples")
    analysis_values["ecg_plf_mean"][0] = np.mean(plf[analysis_mask])
    analysis_values["ecg_plf_std"][0] = np.std(plf[analysis_mask], ddof=1)
    analysis_values["ecg_phf_mean"][0] = np.mean(phf[analysis_mask])
    analysis_values["ecg_phf_std"][0] = np.std(phf[analysis_mask], ddof=1)
    analysis_values["ecg_lf_hf_mean"][0] = np.mean(lf_hf[analysis_mask])
    analysis_values["ecg_lf_hf_std"][0] = np.std(lf_hf[analysis_mask], ddof=1)

    _validate_frequency_values(analysis_values)
    return analysis_values, {
        "mode": 1,
        "matrix": stft_matrix,
        "time": t_stft,
        "frequencies": frequencies,
        "power": power,
        "lf_hf": lf_hf,
    }


def _empty_frequency_values(count: int) -> dict[str, np.ndarray]:
    return {
        feature_name: np.full(count, np.nan, dtype=float)
        for feature_name in _FREQUENCY_DOMAIN_FEATURE_NAMES
    }


def _validate_frequency_values(analysis_values: dict[str, np.ndarray]) -> None:
    for values in analysis_values.values():
        if not np.all(np.isfinite(values)):
            raise ValueError("too few IBI samples are available for frequency analysis")


def _next_power_of_two(value: int) -> int:
    return 1 << (int(value) - 1).bit_length()


def _local_fft_psd(
    x: np.ndarray,
    fs: float,
    nfft: int,
) -> tuple[np.ndarray, np.ndarray]:
    x = np.asarray(x, dtype=float).reshape(-1)
    x = x - np.mean(x)
    window = np.hamming(x.size)
    xw = x * window
    spectrum = np.fft.fft(xw, nfft)
    n_freq = nfft // 2 + 1
    pxx = np.abs(spectrum[:n_freq]) ** 2 / (fs * np.sum(window**2))
    frequencies = np.arange(n_freq) * fs / nfft
    return pxx, frequencies


def _local_band_power(
    pxx: np.ndarray,
    frequencies: np.ndarray,
    band: tuple[float, float],
) -> float:
    index = (frequencies >= band[0]) & (frequencies <= band[1])
    if not np.any(index):
        return float("nan")
    df = float(np.mean(np.diff(frequencies)))
    return float(np.sum(pxx[index]) * df)


def _centered_stft(
    x: np.ndarray,
    nfft: int,
    fs: float,
    winlen: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    window = np.hamming(winlen)
    x = np.asarray(x, dtype=float).reshape(-1)
    if winlen % 2 == 0:
        padded = np.concatenate((np.zeros(winlen // 2 - 1), x, np.zeros(winlen // 2)))
    else:
        padded = np.concatenate((np.zeros((winlen - 1) // 2), x, np.zeros((winlen - 1) // 2)))

    segments = np.zeros((winlen, x.size), dtype=float)
    for index in range(x.size):
        segment = padded[index : index + winlen]
        segments[:, index] = window * (segment - np.mean(segment))

    spectrum = np.fft.fft(segments, nfft, axis=0)
    n_freq = nfft // 2 + 1
    frequencies = np.arange(n_freq) * fs / nfft
    times = np.arange(x.size) / fs
    return spectrum[:n_freq, :], times, frequencies


def _plot_ecg_frequency_diagnostics(
    label_frequency_analysis: int,
    plot_data: dict[str, object],
    *,
    show_plots: bool,
    save_plots: bool,
    output_dir: str | Path | None,
    output_label: str,
) -> None:
    output_path = None
    if save_plots:
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)

    if label_frequency_analysis == 1:
        figures = [
            _plot_stft_hrv_power(plot_data),
            _plot_stft_lf_hf(plot_data),
        ]
        filenames = [
            f"{output_label}_STFT_HRV.png",
            f"{output_label}_STFT_LFHF.png",
        ]
    else:
        figures = [_plot_lf_hf_summary(plot_data)]
        filenames = [f"{output_label}_PSD_LFHF.png"]

    if save_plots and output_path is not None:
        for figure, filename in zip(figures, filenames):
            figure.savefig(output_path / filename, dpi=150)

    if show_plots:
        plt.show()
    else:
        for figure in figures:
            plt.close(figure)


def _plot_stft_hrv_power(
    plot_data: dict[str, object],
):
    time = np.asarray(plot_data["time"])
    frequencies = np.asarray(plot_data["frequencies"])
    power = np.asarray(plot_data["power"])
    freq_index = (frequencies >= 0.03) & (frequencies <= 0.50)

    figure, axes = plt.subplots()
    axes.pcolormesh(time, frequencies[freq_index], power[freq_index, :], shading="auto")
    axes.axhline(0.15, color="w", linestyle="-")
    axes.set_xlabel("Time (s)")
    axes.set_ylabel("Frequency (Hz)")
    axes.set_title("STFT HRV power, window = 60.0 s")
    axes.set_ylim([0.04, 0.40])
    axes.grid(True)
    return figure


def _plot_stft_lf_hf(
    plot_data: dict[str, object],
):
    time = np.asarray(plot_data["time"])
    lf_hf = np.asarray(plot_data["lf_hf"])
    figure, axes = plt.subplots()
    axes.plot(time, lf_hf, linewidth=1.5)
    axes.set_xlabel("Time (s)")
    axes.set_ylabel("LF/HF (-)")
    axes.set_title("STFT LF/HF")
    axes.grid(True)
    return figure


def _plot_lf_hf_summary(plot_data: dict[str, object]):
    lf_hf = np.asarray(plot_data["lf_hf"])
    figure, axes = plt.subplots()
    axes.plot(np.arange(1, lf_hf.size + 1), lf_hf, "o-", linewidth=1.5)
    axes.set_xlabel("Analysis unit (-)")
    axes.set_ylabel("LF/HF (-)")
    axes.set_title("LF/HF summary")
    axes.grid(True)
    return figure

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pytest

from multimodalphysiokit import Signal
from multimodalphysiokit.processors import detect_ecg_r_peaks


def make_synthetic_ecg(duration_seconds=10, sampling_frequency=250):
    time = np.arange(0, duration_seconds, 1 / sampling_frequency)
    samples = 0.02 * np.sin(2 * np.pi * 0.5 * time)
    for beat_time in np.arange(0.8, duration_seconds, 1.0):
        samples += np.exp(-0.5 * ((time - beat_time) / 0.025) ** 2)
        samples += 0.15 * np.exp(-0.5 * ((time - (beat_time + 0.22)) / 0.06) ** 2)
    return Signal(samples, time, sampling_frequency, label="ecg", units="mV")


def test_valid_synthetic_ecg_returns_dictionary():
    result = detect_ecg_r_peaks(make_synthetic_ecg())

    assert isinstance(result, dict)


def test_returned_dictionary_has_exact_keys():
    result = detect_ecg_r_peaks(make_synthetic_ecg())

    assert set(result) == {"r_peak_amplitudes", "r_peak_indexes", "filtered_ecg"}


def test_r_peak_indexes_are_integer_indexes():
    result = detect_ecg_r_peaks(make_synthetic_ecg())

    assert np.issubdtype(result["r_peak_indexes"].dtype, np.integer)


def test_r_peak_indexes_are_within_signal_bounds():
    signal = make_synthetic_ecg()

    result = detect_ecg_r_peaks(signal)

    assert np.all(result["r_peak_indexes"] >= 0)
    assert np.all(result["r_peak_indexes"] < signal.n_samples)


def test_r_peak_amplitudes_length_matches_r_peak_indexes_length():
    result = detect_ecg_r_peaks(make_synthetic_ecg())

    assert result["r_peak_amplitudes"].shape == result["r_peak_indexes"].shape


def test_filtered_ecg_length_matches_signal_length():
    signal = make_synthetic_ecg()

    result = detect_ecg_r_peaks(signal)

    assert result["filtered_ecg"].shape == signal.samples.shape


def test_detects_approximately_expected_number_of_peaks():
    result = detect_ecg_r_peaks(make_synthetic_ecg(duration_seconds=10))

    assert 9 <= result["r_peak_indexes"].size <= 11


def test_invalid_signal_type_raises_type_error():
    with pytest.raises(TypeError, match="Signal"):
        detect_ecg_r_peaks([1, 2, 3])


@pytest.mark.parametrize(
    "option_name, option_value",
    [
        ("show_plots", "yes"),
        ("save_plots", "no"),
    ],
)
def test_invalid_show_plots_save_plots_raises_type_error(option_name, option_value):
    with pytest.raises(TypeError, match=option_name):
        detect_ecg_r_peaks(make_synthetic_ecg(), **{option_name: option_value})


def test_save_plots_true_without_output_dir_raises_value_error():
    with pytest.raises(ValueError, match="output_dir"):
        detect_ecg_r_peaks(make_synthetic_ecg(), save_plots=True)


def test_too_short_signal_raises_value_error():
    time = np.arange(50) / 250
    signal = Signal(np.sin(2 * np.pi * 5 * time), time, 250)

    with pytest.raises(ValueError, match="too short"):
        detect_ecg_r_peaks(signal)


def test_too_low_sampling_frequency_raises_value_error():
    time = np.arange(200) / 30
    signal = Signal(np.sin(2 * np.pi * 5 * time), time, 30)

    with pytest.raises(ValueError, match="16 Hz"):
        detect_ecg_r_peaks(signal)


def test_flat_no_peak_signal_raises_value_error():
    time = np.arange(1000) / 250
    signal = Signal(np.zeros_like(time), time, 250)

    with pytest.raises(ValueError, match="no R-peaks"):
        detect_ecg_r_peaks(signal)


def test_show_plots_false_does_not_leave_open_figures():
    plt.close("all")

    detect_ecg_r_peaks(make_synthetic_ecg(), show_plots=False, save_plots=False)

    assert plt.get_fignums() == []


def test_save_plots_true_writes_png_files_to_temporary_folder(tmp_path):
    plt.close("all")

    detect_ecg_r_peaks(make_synthetic_ecg(), save_plots=True, output_dir=tmp_path)

    assert (tmp_path / "ecg_pan_tompkins_stages.png").is_file()
    assert (tmp_path / "ecg_qrs_detection_summary.png").is_file()
    assert plt.get_fignums() == []

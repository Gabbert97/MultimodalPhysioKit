import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pytest

from multimodalphysiokit import ProcessingResult, Signal
from multimodalphysiokit.processors import process_temperature
from multimodalphysiokit.processors.temperature import _screen_temperature_phase_artifacts


FEATURE_NAMES = [
    "temperature_delta",
    "temperature_mean",
    "temperature_std",
    "temperature_delta_over_time",
    "temperature_linear_slope",
    "temperature_derivative_mean",
    "temperature_derivative_std",
    "temperature_derivative_linear_slope",
]


def make_temperature_signal(n_samples=240, sampling_frequency=10, label="temp"):
    time = np.arange(n_samples) / sampling_frequency
    samples = 32 + 0.02 * time + 0.1 * np.sin(2 * np.pi * 0.05 * time)
    return Signal(samples, time, sampling_frequency, label=label, units="degC")


def test_functional_wrapper_returns_processing_result():
    result = process_temperature(make_temperature_signal())

    assert isinstance(result, ProcessingResult)
    assert result.modality == "temperature"
    assert result.processor_name == "process_temperature"


def test_all_eight_exact_phase_independent_feature_names_are_returned():
    result = process_temperature(make_temperature_signal())

    assert list(result.features) == FEATURE_NAMES
    assert all("__" not in name for name in result.features)


def test_signal_label_is_used_in_processing_result():
    result = process_temperature(make_temperature_signal(label="skin_temperature"))

    assert result.signal_labels == ["skin_temperature"]


def test_empty_signal_label_uses_temperature_fallback():
    result = process_temperature(make_temperature_signal(label=""))

    assert result.signal_labels == ["temperature"]


def test_metadata_contains_signal_level_information_and_no_phase_information():
    signal = make_temperature_signal()

    result = process_temperature(signal)

    assert result.metadata == {
        "source": "process_temperature",
        "sampling_frequency": 10.0,
        "n_samples": 240,
        "start_time_seconds": 0.0,
        "end_time_seconds": 23.9,
    }
    assert all("phase" not in key for key in result.metadata)


def test_invalid_signal_type_raises_type_error():
    with pytest.raises(TypeError, match="Signal"):
        process_temperature([1, 2, 3])


def test_too_low_sampling_frequency_preserves_filtering_error():
    with pytest.raises(ValueError, match="sampling_frequency <= 2 Hz"):
        process_temperature(make_temperature_signal(sampling_frequency=2))


def test_too_short_signal_raises_clear_error():
    signal = Signal([32, 33], [0, 0.1], 10)

    with pytest.raises(ValueError, match="at least 3 samples"):
        process_temperature(signal)


def test_zero_duration_signal_is_rejected():
    signal = Signal([32], [4], 10)

    with pytest.raises(ValueError, match="duration must be nonzero"):
        process_temperature(signal)


def test_filtering_too_short_signal_preserves_filtering_error():
    with pytest.raises(ValueError, match="signal is too short"):
        process_temperature(make_temperature_signal(n_samples=10))


def test_output_features_are_finite_for_smooth_synthetic_signal():
    result = process_temperature(make_temperature_signal())

    assert all(np.isfinite(value) for value in result.features.values())


def test_feature_names_and_has_feature_use_phase_independent_keys():
    result = process_temperature(make_temperature_signal())

    assert result.feature_names() == FEATURE_NAMES
    assert result.has_feature("temperature_mean")
    assert not result.has_feature("temperature_mean__baseline")


def test_show_plots_false_does_not_leave_open_figures():
    plt.close("all")

    process_temperature(make_temperature_signal(), show_plots=False, save_plots=False)

    assert plt.get_fignums() == []


def test_save_plots_true_writes_expected_files(tmp_path):
    plt.close("all")

    process_temperature(
        make_temperature_signal(),
        save_plots=True,
        output_dir=tmp_path,
    )

    assert (tmp_path / "temperature_vs_time.png").is_file()
    assert (tmp_path / "temperature_first_derivative_vs_time.png").is_file()
    assert plt.get_fignums() == []


def test_save_plots_true_without_output_dir_raises_value_error():
    with pytest.raises(ValueError, match="output_dir"):
        process_temperature(make_temperature_signal(), save_plots=True)


def test_plotting_does_not_change_extracted_feature_values(tmp_path):
    signal = make_temperature_signal()

    without_plots = process_temperature(signal)
    with_saved_plots = process_temperature(signal, save_plots=True, output_dir=tmp_path)

    assert with_saved_plots.features == without_plots.features


def test_artifact_rejected_windows_are_excluded_from_feature_extraction():
    sampling_frequency = 10
    time = np.arange(80) / sampling_frequency
    filtered_phase = 32 + 0.01 * time
    filtered_phase[20:40] = 80

    screened_temperature, screened_time = _screen_temperature_phase_artifacts(
        filtered_phase,
        time,
        sampling_frequency,
    )

    assert screened_temperature.size < filtered_phase.size
    assert not np.any((screened_time >= 2.0) & (screened_time < 4.0))
    assert np.all(screened_temperature < 40)


def test_smooth_drift_is_kept_whole_for_cooling_and_warming():
    # A smooth, artifact-free 6 degC drift must not be screened out in
    # either direction: the delta must be the physical +/-6 degC change.
    sampling_frequency = 10.0
    time = np.arange(0, 300, 1 / sampling_frequency)
    warming = Signal(30 + 6 * time / time[-1], time, sampling_frequency, units="degC")
    cooling = Signal(36 - 6 * time / time[-1], time, sampling_frequency, units="degC")

    assert process_temperature(warming).features["temperature_delta"] == pytest.approx(6.0, abs=1e-3)
    assert process_temperature(cooling).features["temperature_delta"] == pytest.approx(-6.0, abs=1e-3)


def test_artifact_screening_is_the_same_in_celsius_and_kelvin():
    # A 1 s, 5 degC sensor-detachment drop on a stable 33 degC signal.
    # Temperature differences do not depend on the absolute scale, so every
    # feature except the mean must be identical in degC and K.
    sampling_frequency = 10.0
    time = np.arange(0, 300, 1 / sampling_frequency)
    rng = np.random.default_rng(1)
    celsius = 33 + 0.01 * rng.standard_normal(time.size)
    celsius[(time >= 150) & (time < 151)] -= 5

    features_c = process_temperature(Signal(celsius, time, sampling_frequency, units="degC")).features
    features_k = process_temperature(Signal(celsius + 273.15, time, sampling_frequency, units="K")).features

    for name, value in features_c.items():
        expected = value + 273.15 if name == "temperature_mean" else value
        assert features_k[name] == pytest.approx(expected, rel=1e-6, abs=1e-9), name

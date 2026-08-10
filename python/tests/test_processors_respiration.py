import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pytest

from multimodalphysiokit import ProcessingResult, Signal
from multimodalphysiokit.processors import process_respiration
from multimodalphysiokit.processors.respiration import (
    _definitive_breath_points,
    _filter_respiration_signal,
    _find_peaks_respiration,
)


FEATURE_NAMES = [
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
]

EXPECTED_FEATURES = {
    "respiration_breath_rate_mean": 14.394645317725068,
    "respiration_breath_rate_std": 0.3941623575698864,
    "respiration_inspiratory_time_mean": 2.0714814814814835,
    "respiration_inspiratory_time_std": 0.13309059100141352,
    "respiration_expiratory_time_mean": 2.091111111111109,
    "respiration_expiratory_time_std": 0.14229798996724516,
    "respiration_timing_ratio_mean": 0.9978298900267016,
    "respiration_timing_ratio_std": 0.12086226552906491,
    "respiration_amplitude_mean": 2.013419562316465,
    "respiration_amplitude_std": 0.009736309440961361,
    "respiration_minute_ventilation_mean": 28.986685480630854,
    "respiration_minute_ventilation_std": 0.6972856869535435,
    "respiration_coefficient_of_variation": 2.7382568230738444,
    "respiration_autocorrelation": 0.08009763324674718,
}


def make_respiration_signal(duration_seconds=120, sampling_frequency=200, label="respiration"):
    time = np.arange(0, duration_seconds, 1 / sampling_frequency)
    angle = 2 * np.pi * (
        0.24 * time + 0.015 * np.sin(2 * np.pi * 0.02 * time)
    )
    samples = np.sin(angle) + 0.05 * np.sin(2 * np.pi * 0.9 * time)
    return Signal(samples, time, sampling_frequency, label=label)


def controlled_intermediates():
    signal = make_respiration_signal()
    downsampled = signal.samples[::2]
    downsampled_time = signal.time[::2]
    filtered = _filter_respiration_signal(downsampled)
    minima, maxima = _find_peaks_respiration(filtered)
    if minima[0] >= maxima[0]:
        maxima = maxima[1:]
    if minima[-1] <= maxima[-1]:
        maxima = maxima[:-1]
    definitive_minima, definitive_maxima = _definitive_breath_points(
        filtered,
        downsampled_time,
        minima,
        maxima,
    )
    return downsampled, downsampled_time, filtered, minima, maxima, definitive_minima, definitive_maxima


def test_functional_wrapper_returns_processing_result():
    result = process_respiration(make_respiration_signal())

    assert isinstance(result, ProcessingResult)
    assert result.modality == "respiration"
    assert result.processor_name == "process_respiration"


def test_exactly_fourteen_validated_base_feature_names_are_returned():
    result = process_respiration(make_respiration_signal())

    assert list(result.features) == FEATURE_NAMES
    assert len(result.features) == 14
    assert all("__" not in name for name in result.features)


def test_signal_label_and_empty_label_fallback():
    assert process_respiration(make_respiration_signal(label="rip")).signal_labels == ["rip"]
    assert process_respiration(make_respiration_signal(label="")).signal_labels == ["respiration"]


def test_metadata_contains_signal_level_information_and_no_phase_information():
    result = process_respiration(make_respiration_signal())

    assert result.metadata == {
        "source": "process_respiration",
        "original_sampling_frequency": 200.0,
        "processed_sampling_frequency": 100.0,
        "original_n_samples": 24000,
        "processed_n_samples": 12000,
        "start_time_seconds": 0.0,
        "end_time_seconds": 119.995,
    }
    assert all("phase" not in key for key in result.metadata)


def test_invalid_signal_type_raises_type_error():
    with pytest.raises(TypeError, match="Signal"):
        process_respiration([1, 2, 3])


@pytest.mark.parametrize(
    "option_name, option_value",
    [("show_plots", "yes"), ("save_plots", "no")],
)
def test_invalid_plot_options_raise_type_error(option_name, option_value):
    with pytest.raises(TypeError, match=option_name):
        process_respiration(
            make_respiration_signal(),
            **{option_name: option_value},
        )


def test_unsupported_sampling_frequency_preserves_existing_error():
    with pytest.raises(ValueError, match="integer multiple of 100"):
        process_respiration(make_respiration_signal(sampling_frequency=150))


def test_too_short_input_raises_clear_error():
    signal = Signal(np.arange(20), np.arange(20) / 200, 200)

    with pytest.raises(ValueError, match="signal is too short"):
        process_respiration(signal)


def test_zero_duration_input_is_rejected():
    signal = Signal([1], [4], 200)

    with pytest.raises(ValueError, match="duration must be nonzero"):
        process_respiration(signal)


def test_downsampled_samples_and_time_match_existing_implementation():
    downsampled, downsampled_time, *_ = controlled_intermediates()

    np.testing.assert_allclose(
        downsampled[:8],
        [
            0.0,
            0.01802342161017329,
            0.0360342985542393,
            0.054020115858762784,
            0.07196841784313382,
            0.08986683753497904,
            0.10770312580919629,
            0.1254651801598614,
        ],
        rtol=1e-13,
        atol=1e-13,
    )
    np.testing.assert_allclose(downsampled_time[:8], np.arange(8) / 100, atol=0, rtol=0)
    assert downsampled.size == downsampled_time.size == 12000


def test_filtered_values_match_existing_implementation():
    *_, filtered, minima, maxima, definitive_minima, definitive_maxima = controlled_intermediates()

    # Accommodate cross-platform SciPy/NumPy floating-point variation while
    # preserving tight numerical regression checking.
    np.testing.assert_allclose(
        filtered[:8],
        [
            0.0036618748297932873,
            0.018098389893553245,
            0.032766524383624866,
            0.04765012048436981,
            0.06273259775396016,
            0.07799700947566192,
            0.09342609987855974,
            0.10900236194192842,
        ],
        rtol=1e-9,
        atol=2e-11,
    )
    assert minima.size == 29
    assert maxima.size == 28
    np.testing.assert_array_equal(definitive_minima[:5], [308, 733, 1141, 1547, 1974])
    np.testing.assert_array_equal(definitive_maxima[:5], [508, 926, 1353, 1776, 2176])
    assert definitive_minima.size == 28
    assert definitive_maxima.size == 27


def test_all_feature_values_match_former_complete_signal_path():
    result = process_respiration(make_respiration_signal())

    assert result.features.keys() == EXPECTED_FEATURES.keys()
    for name, expected in EXPECTED_FEATURES.items():
        assert result.features[name] == pytest.approx(expected, rel=1e-9, abs=2e-11)


def test_all_feature_values_are_finite():
    result = process_respiration(make_respiration_signal())

    assert all(np.isfinite(value) for value in result.features.values())


def test_feature_names_and_has_feature_use_base_names():
    result = process_respiration(make_respiration_signal())

    assert result.feature_names() == FEATURE_NAMES
    assert result.has_feature("respiration_breath_rate_mean")
    assert not result.has_feature("respiration_breath_rate_mean__baseline")


def test_show_plots_false_does_not_leave_open_figures():
    plt.close("all")

    process_respiration(make_respiration_signal(), show_plots=False, save_plots=False)

    assert plt.get_fignums() == []


def test_save_plots_writes_expected_files(tmp_path):
    plt.close("all")

    process_respiration(make_respiration_signal(), save_plots=True, output_dir=tmp_path)

    assert (tmp_path / "respiration_downsampled_signal.png").is_file()
    assert (tmp_path / "respiration_filtered_phase_peaks.png").is_file()
    assert plt.get_fignums() == []


def test_save_plots_without_output_dir_raises_value_error():
    with pytest.raises(ValueError, match="output_dir"):
        process_respiration(make_respiration_signal(), save_plots=True)


def test_plotting_does_not_change_feature_values(tmp_path):
    signal = make_respiration_signal()

    without_plots = process_respiration(signal)
    with_plots = process_respiration(signal, save_plots=True, output_dir=tmp_path)

    assert with_plots.features == without_plots.features

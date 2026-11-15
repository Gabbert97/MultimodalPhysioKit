import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pytest

from multimodalphysiokit import ProcessingResult, Signal
from multimodalphysiokit.processors import process_ecg


FEATURE_NAMES = [
    "ecg_bpm_mean", "ecg_bpm_std", "ecg_bpm_median", "ecg_pnn50", "ecg_pnn50_max",
    "ecg_plf_mean", "ecg_plf_std", "ecg_phf_mean", "ecg_phf_std",
    "ecg_lf_hf_mean", "ecg_lf_hf_std",
]


def make_synthetic_ecg(duration_seconds=40, sampling_frequency=250, label="ecg"):
    time = np.arange(0, duration_seconds, 1 / sampling_frequency)
    samples = 0.02 * np.sin(2 * np.pi * 0.5 * time)
    for beat_time in np.arange(0.8, duration_seconds, 1.0):
        samples += np.exp(-0.5 * ((time - beat_time) / 0.025) ** 2)
        samples += 0.15 * np.exp(-0.5 * ((time - (beat_time + 0.22)) / 0.06) ** 2)
    return Signal(samples, time, sampling_frequency, label=label, units="mV")


def test_functional_wrapper_returns_processing_result():
    result = process_ecg(make_synthetic_ecg(), label_frequency_analysis=2)
    assert isinstance(result, ProcessingResult)
    assert result.modality == "ecg"
    assert result.processor_name == "process_ecg"


def test_exact_validated_base_feature_names_are_returned_without_suffixes():
    result = process_ecg(make_synthetic_ecg(), label_frequency_analysis=2)
    assert result.feature_names() == FEATURE_NAMES
    assert len(result.features) == 11
    assert all("__" not in name for name in result.features)


def test_signal_label_and_empty_label_fallback():
    assert process_ecg(make_synthetic_ecg(label="lead_ii"), label_frequency_analysis=2).signal_labels == ["lead_ii"]
    assert process_ecg(make_synthetic_ecg(label=""), label_frequency_analysis=2).signal_labels == ["ecg"]


def test_metadata_is_signal_level_and_contains_no_phase_information():
    signal = make_synthetic_ecg()
    result = process_ecg(signal, label_frequency_analysis=2)
    assert result.metadata == {
        "source": "process_ecg",
        "sampling_frequency": 250.0,
        "n_samples": 10000,
        "start_time_seconds": 0.0,
        "end_time_seconds": 39.996,
        "frequency_analysis_method": 2,
        "return_intermediates": False,
    }
    assert all("phase" not in key for key in result.metadata)


def test_invalid_signal_and_options_raise_clear_errors():
    with pytest.raises(TypeError, match="Signal"):
        process_ecg([1, 2, 3], label_frequency_analysis=2)
    with pytest.raises(ValueError, match="label_frequency_analysis"):
        process_ecg(make_synthetic_ecg(), label_frequency_analysis=3)
    with pytest.raises(TypeError, match="output_label"):
        process_ecg(make_synthetic_ecg(), label_frequency_analysis=2, output_label=3)


def test_stft_mode_validates_supplied_signal_duration():
    with pytest.raises(ValueError, match="input signal.*60 seconds"):
        process_ecg(make_synthetic_ecg(40), label_frequency_analysis=1)


def test_psd_features_match_existing_complete_signal_values():
    result = process_ecg(make_synthetic_ecg(), label_frequency_analysis=2)
    expected = {
        "ecg_bpm_mean": 60.0, "ecg_bpm_std": 0.0, "ecg_bpm_median": 60.0,
        "ecg_pnn50": 0.0, "ecg_pnn50_max": 0.0, "ecg_plf_mean": 0.0,
        "ecg_plf_std": 0.0, "ecg_phf_mean": 0.0, "ecg_phf_std": 0.0,
        "ecg_lf_hf_mean": 0.0, "ecg_lf_hf_std": 0.0,
    }
    assert result.features == expected


def test_default_result_has_no_intermediates():
    result = process_ecg(make_synthetic_ecg(), label_frequency_analysis=2)
    assert result.intermediates == {}


def test_intermediate_arrays_match_existing_complete_signal_path():
    result = process_ecg(
        make_synthetic_ecg(), label_frequency_analysis=2, return_intermediates=True
    )
    data = result.intermediates
    assert set(data) == {
        "filtered_signal", "r_peak_indices", "r_peak_times", "r_peak_amplitudes",
        "accepted_r_peak_indices", "accepted_r_peak_times", "ibi", "ibi_times",
        "bpm", "bpm_times",
    }
    np.testing.assert_array_equal(data["r_peak_indices"][:8], np.arange(200, 1951, 250))
    np.testing.assert_array_equal(data["accepted_r_peak_indices"][:8], np.arange(450, 2201, 250))
    np.testing.assert_allclose(data["ibi"], 1.0, atol=0, rtol=0)
    np.testing.assert_allclose(data["bpm"], 60.0, atol=0, rtol=0)
    np.testing.assert_allclose(
        data["filtered_signal"][:4],
        [-1.10783021e-4, -6.89429583e-5, -2.84604654e-5, 9.46275060e-6],
        rtol=1e-8, atol=1e-12,
    )


def test_stft_intermediates_are_consistent_and_do_not_change_features():
    signal = make_synthetic_ecg(70)
    without = process_ecg(signal, label_frequency_analysis=1)
    with_data = process_ecg(signal, label_frequency_analysis=1, return_intermediates=True)
    data = with_data.intermediates
    assert with_data.features == without.features
    assert np.iscomplexobj(data["stft_matrix"])
    np.testing.assert_allclose(data["stft_power"], np.abs(data["stft_matrix"]) ** 2)
    assert data["stft_matrix"].shape == data["stft_power"].shape
    assert data["stft_matrix"].shape[0] == data["stft_frequencies"].size
    assert data["stft_matrix"].shape[1] == data["stft_times"].size


def test_save_plots_requires_output_dir():
    with pytest.raises(ValueError, match="output_dir"):
        process_ecg(make_synthetic_ecg(), label_frequency_analysis=2, save_plots=True)


def test_psd_plotting_writes_expected_files_and_preserves_features(tmp_path):
    plt.close("all")
    signal = make_synthetic_ecg()
    without = process_ecg(signal, label_frequency_analysis=2)
    with_plots = process_ecg(
        signal, label_frequency_analysis=2, return_intermediates=True,
        save_plots=True, output_dir=tmp_path, output_label="ecg_test",
    )
    assert with_plots.features == without.features
    assert (tmp_path / "ecg_test_PSD_LFHF.png").is_file()
    assert (tmp_path / "ecg_test_BPM_trend.png").is_file()
    assert plt.get_fignums() == []


def test_stft_plotting_writes_expected_files_and_preserves_intermediates(tmp_path):
    plt.close("all")
    signal = make_synthetic_ecg(70)
    without = process_ecg(signal, label_frequency_analysis=1, return_intermediates=True)
    with_plots = process_ecg(
        signal, label_frequency_analysis=1, return_intermediates=True,
        save_plots=True, output_dir=tmp_path, output_label="stft_test",
    )
    assert with_plots.features == without.features
    for key in without.intermediates:
        np.testing.assert_array_equal(with_plots.intermediates[key], without.intermediates[key])
    assert (tmp_path / "stft_test_STFT_HRV.png").is_file()
    assert (tmp_path / "stft_test_STFT_LFHF.png").is_file()
    assert plt.get_fignums() == []


def test_pnn50_compares_only_truly_adjacent_ibis_after_a_rejected_interval():
    from multimodalphysiokit.processors.ecg import _ecg_time_domain_feature_values

    # Interval 3 was rejected, so the accepted IBIs form two contiguous runs:
    # [0.8, 0.8, 0.8] and [1.0, 1.0, 1.1, 1.2]. Successive differences within
    # the runs are [0, 0] and [0, 0.1, 0.1]; no difference may span the gap.
    # Textbook pNN50: 2 of the 5 evaluated differences exceed 50 ms.
    flags = np.array([1, 1, 1, 0, 1, 1, 1, 1])
    ibi = np.array([0.8, 0.8, 0.8, 1.0, 1.0, 1.1, 1.2])

    features = _ecg_time_domain_feature_values(60 / ibi, ibi, flags)

    assert features["ecg_pnn50"] == pytest.approx(2 / 5)
    expected_ibis_in_2_minutes = 120 * ibi.size / ibi.sum()
    assert features["ecg_pnn50_max"] == pytest.approx(2 / expected_ibis_in_2_minutes)


def test_pnn50_ignores_step_change_across_a_rejected_interval():
    from multimodalphysiokit.processors.ecg import _ecg_time_domain_feature_values

    # Each contiguous run is perfectly regular; the only large change lies
    # across the rejected interval and therefore is not a successive difference.
    flags = np.array([1, 1, 1, 0, 1, 1, 1])
    ibi = np.array([0.8, 0.8, 0.8, 1.0, 1.0, 1.0])

    features = _ecg_time_domain_feature_values(60 / ibi, ibi, flags)

    assert features["ecg_pnn50"] == 0.0


def test_pnn50_counts_isolated_differences_above_50_ms():
    from multimodalphysiokit.processors.ecg import _ecg_time_domain_feature_values

    # Task Force (1996) pNN50: share of successive NN differences > 50 ms.
    # Alternating 0.8/0.9 s intervals: 5 of 9 differences are 100 ms.
    ibi = np.array([0.8, 0.9, 0.9, 0.8, 0.8, 0.9, 0.9, 0.8, 0.8, 0.9])
    flags = np.ones(ibi.size, dtype=int)

    features = _ecg_time_domain_feature_values(60 / ibi, ibi, flags)

    assert features["ecg_pnn50"] == pytest.approx(5 / 9)
    # Longest run of consecutive qualifying differences is 1.
    assert features["ecg_pnn50_max"] == pytest.approx(1 / (120 * ibi.size / ibi.sum()))


def test_pnn50_excludes_differences_of_50_ms_or_less():
    from multimodalphysiokit.processors.ecg import _ecg_time_domain_feature_values

    ibi = np.array([0.80, 0.84, 0.80, 0.84, 0.80, 0.84])  # all differences 40 ms
    flags = np.ones(ibi.size, dtype=int)

    features = _ecg_time_domain_feature_values(60 / ibi, ibi, flags)

    assert features["ecg_pnn50"] == 0.0
    assert features["ecg_pnn50_max"] == 0.0


def _rsa_beats(duration_seconds=300.0, rsa_frequency=0.25):
    beat_times = [0.0]
    while beat_times[-1] < duration_seconds:
        last = beat_times[-1]
        beat_times.append(last + 0.9 + 0.05 * np.sin(2 * np.pi * rsa_frequency * last))
    beat_times = np.asarray(beat_times)
    return np.diff(beat_times), beat_times[:-1]


def test_psd_tachogram_keeps_real_time_across_rejected_intervals():
    from multimodalphysiokit.processors.ecg import _uniform_ibi_tachogram_from_ibi

    ibi, starts = _rsa_beats()
    keep = np.ones(ibi.size, dtype=bool)
    keep[np.arange(15, ibi.size - 15, 25)] = False  # rejected intervals
    true_span = (starts[keep][-1] + ibi[keep][-1]) - starts[keep][0]

    uniform_time, tachogram = _uniform_ibi_tachogram_from_ibi(
        ibi[keep], 4, ibi_start_times=starts[keep]
    )

    # The resampled tachogram must cover the real elapsed time ...
    assert uniform_time[-1] == pytest.approx(true_span, abs=0.25)
    # ... so respiratory sinus arrhythmia stays at its true 0.25 Hz.
    n_fft = 1 << 16
    spectrum = np.abs(np.fft.rfft(tachogram * np.hanning(tachogram.size), n_fft))
    frequencies = np.fft.rfftfreq(n_fft, d=1 / 4)
    assert frequencies[np.argmax(spectrum)] == pytest.approx(0.25, abs=0.002)


def test_psd_tachogram_without_rejections_matches_cumulative_timeline():
    from multimodalphysiokit.processors.ecg import _uniform_ibi_tachogram_from_ibi

    ibi, starts = _rsa_beats(duration_seconds=120.0)
    time_a, values_a = _uniform_ibi_tachogram_from_ibi(ibi, 4)
    time_b, values_b = _uniform_ibi_tachogram_from_ibi(ibi, 4, ibi_start_times=starts + 7.0)

    np.testing.assert_allclose(time_b, time_a, atol=1e-9)
    np.testing.assert_allclose(values_b, values_a, atol=1e-9)

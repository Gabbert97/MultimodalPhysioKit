import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pytest

from multimodalphysiokit.processors.ecg import (
    _extract_ecg_frequency_domain_features,
    _stft_context_to_tachogram,
)


FEATURE_NAMES = [
    "ecg_plf_mean", "ecg_plf_std", "ecg_phf_mean", "ecg_phf_std",
    "ecg_lf_hf_mean", "ecg_lf_hf_std",
]


def make_ibi(n_samples=90, base=0.82):
    index = np.arange(n_samples)
    return base + 0.03 * np.sin(2 * np.pi * 0.08 * index) + 0.015 * np.sin(2 * np.pi * 0.23 * index)


def make_stft_context(total_seconds=140):
    ibi_time = np.arange(1.0, total_seconds - 1.0, 0.8)
    index = np.arange(ibi_time.size)
    ibi = 0.8 + 0.03 * np.sin(2 * np.pi * 0.08 * index)
    return {
        "ibi": ibi,
        "ibi_time": ibi_time,
        "time_range": np.array([0.0, float(total_seconds)]),
        "analysis_window": np.array([10.0, 120.0]),
    }


def test_psd_mode_returns_six_ordered_finite_features():
    features, data = _extract_ecg_frequency_domain_features(make_ibi(), label_frequency_analysis=2)
    assert list(features) == FEATURE_NAMES
    assert all(np.isfinite(value) for value in features.values())
    assert features["ecg_plf_std"] == features["ecg_phf_std"] == features["ecg_lf_hf_std"] == 0.0
    assert data["mode"] == 2


def test_stft_mode_returns_existing_complex_matrix_and_power():
    features, data = _extract_ecg_frequency_domain_features(
        make_ibi(180), label_frequency_analysis=1, stft_context=make_stft_context()
    )
    assert list(features) == FEATURE_NAMES
    assert all(np.isfinite(value) for value in features.values())
    assert np.iscomplexobj(data["matrix"])
    np.testing.assert_allclose(data["power"], np.abs(data["matrix"]) ** 2)
    assert data["matrix"].shape == data["power"].shape
    assert data["matrix"].shape[0] == data["frequencies"].size
    assert data["matrix"].shape[1] == data["time"].size


def test_stft_tachogram_grid_does_not_exceed_final_physical_timestamp():
    sampling_frequency = 100
    duration_seconds = 60
    final_timestamp = duration_seconds - 1 / sampling_frequency
    context = make_stft_context(total_seconds=60)
    context["time_range"] = np.array([0.0, final_timestamp])
    context["analysis_window"] = np.array([0.0, final_timestamp])

    uniform_time, _, _ = _stft_context_to_tachogram(context, fs_hrv=4)

    assert uniform_time.size == 240
    assert uniform_time[-1] == 59.75
    assert not np.any(uniform_time == 60.0)


@pytest.mark.parametrize("ibi, message", [(np.array([]), "non-empty"), (np.array([.8, 0, .9]), "positive"), (np.array([.8, np.nan, .9]), "finite")])
def test_invalid_ibi_is_rejected(ibi, message):
    with pytest.raises(ValueError, match=message):
        _extract_ecg_frequency_domain_features(ibi, label_frequency_analysis=2)


def test_invalid_options_are_rejected():
    with pytest.raises(ValueError, match="label_frequency_analysis"):
        _extract_ecg_frequency_domain_features(make_ibi(), label_frequency_analysis=3)
    with pytest.raises(TypeError, match="show_plots"):
        _extract_ecg_frequency_domain_features(make_ibi(), show_plots="yes")
    with pytest.raises(TypeError, match="output_label"):
        _extract_ecg_frequency_domain_features(make_ibi(), output_label=3)
    with pytest.raises(ValueError, match="output_dir"):
        _extract_ecg_frequency_domain_features(make_ibi(), save_plots=True)


def test_too_short_tachogram_raises_value_error():
    with pytest.raises(ValueError, match="too few|shorter"):
        _extract_ecg_frequency_domain_features(make_ibi(4), label_frequency_analysis=2)


def test_psd_plotting_writes_file_and_preserves_values(tmp_path):
    plt.close("all")
    without, _ = _extract_ecg_frequency_domain_features(make_ibi(), label_frequency_analysis=2)
    with_plots, _ = _extract_ecg_frequency_domain_features(
        make_ibi(), label_frequency_analysis=2, save_plots=True,
        output_dir=tmp_path, output_label="test_ecg",
    )
    assert with_plots == without
    assert (tmp_path / "test_ecg_PSD_LFHF.png").is_file()
    assert plt.get_fignums() == []

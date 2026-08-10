from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pytest

from multimodalphysiokit.core import ProcessingResult, Signal
from multimodalphysiokit.processors import FNIRSProcessor, process_fnirs


NAMES = (
    "min", "max", "mean", "variance", "std", "skewness", "kurtosis",
    "power_band", "max_frequency", "median_frequency", "spectral_entropy",
    "mean_diff", "peak_to_peak", "slope", "zero_crossing", "polarity",
    "entropy", "auc", "total_energy",
)


def signals(*, duration=60.0, fs=100.0, red_label="red", infrared_label="infrared", offset=0.0):
    time = offset + np.arange(0, duration, 1 / fs)
    red = 1.6 + 0.03 * np.sin(2 * np.pi * 0.04 * time) + 0.01 * np.sin(2 * np.pi * 0.3 * time)
    infrared = 1.9 + 0.04 * np.sin(2 * np.pi * 0.035 * time + 0.2) + 0.01 * np.sin(2 * np.pi * 0.25 * time)
    return Signal(red, time, fs, label=red_label, units="A"), Signal(infrared, time, fs, label=infrared_label, units="A")


def processor():
    return FNIRSProcessor(age_years=35)


def test_constructs_without_signals():
    assert processor().age_years == 35.0


def test_process_and_functional_wrapper_return_processing_result():
    red, infrared = signals()
    direct = processor().process(red, infrared)
    wrapped = process_fnirs(red, infrared, age_years=35)
    assert isinstance(direct, ProcessingResult)
    assert wrapped.features == pytest.approx(direct.features)


def test_exact_feature_names_and_count():
    result = processor().process(*signals())
    expected = tuple(f"fnirs_oxy_{name}" for name in NAMES) + tuple(
        f"fnirs_deoxy_{name}" for name in NAMES
    )
    assert tuple(result.features) == expected
    assert len(result.features) == 38
    assert all("__" not in name for name in result.features)
    assert all(np.isfinite(value) for value in result.features.values())


def test_signal_labels_and_fallbacks():
    assert processor().process(*signals(red_label="R", infrared_label="IR")).signal_labels == ["R", "IR"]
    assert processor().process(*signals(red_label="", infrared_label="")).signal_labels == [
        "fnirs_red", "fnirs_infrared"
    ]


def test_self_and_provided_baseline_metadata():
    red, infrared = signals()
    self_result = processor().process(red, infrared)
    supplied = processor().process(
        red, infrared, baseline_red_value=1.5, baseline_infrared_value=np.float64(1.8)
    )
    assert self_result.metadata["baseline_mode"] == "self"
    assert supplied.metadata["baseline_mode"] == "provided"
    assert supplied.metadata["baseline_red_value"] == 1.5
    assert supplied.metadata["baseline_infrared_value"] == 1.8
    assert isinstance(self_result.metadata["baseline_red_value"], float)
    assert isinstance(self_result.metadata["baseline_infrared_value"], float)
    assert not any("phase" in key for key in self_result.metadata)
    assert self_result.intermediates == {}


def test_extract_features_matches_process():
    red, infrared = signals()
    assert processor().extract_features(red, infrared) == pytest.approx(processor().process(red, infrared).features)


@pytest.mark.parametrize("kwargs", [
    {"baseline_red_value": 1.0}, {"baseline_infrared_value": 1.0},
])
def test_partial_baseline_is_rejected(kwargs):
    with pytest.raises(ValueError, match="supplied together"):
        processor().process(*signals(), **kwargs)


@pytest.mark.parametrize("value", ["1", [1.0], {"x": 1}, np.asarray([1.0, 2.0])])
def test_nonnumeric_or_nonscalar_baseline_is_rejected(value):
    with pytest.raises(TypeError, match="numeric scalar"):
        processor().process(*signals(), baseline_red_value=value, baseline_infrared_value=1.0)


@pytest.mark.parametrize("value", [np.nan, np.inf, -np.inf, 0.0, -1.0])
def test_nonpositive_or_nonfinite_baseline_is_rejected(value):
    with pytest.raises(ValueError, match="positive and finite"):
        processor().process(*signals(), baseline_red_value=value, baseline_infrared_value=1.0)


def test_plotting_preserves_features_and_baselines_and_closes_figures(tmp_path):
    plt.close("all")
    red, infrared = signals()
    plain = processor().process(red, infrared)
    plotted = processor().process(red, infrared, save_plots=True, output_dir=tmp_path)
    assert plotted.features == pytest.approx(plain.features)
    assert plotted.metadata["baseline_red_value"] == pytest.approx(plain.metadata["baseline_red_value"])
    assert plotted.metadata["baseline_infrared_value"] == pytest.approx(plain.metadata["baseline_infrared_value"])
    assert {path.name for path in tmp_path.glob("*.png")} == {
        "fnirs_concentration_changes.png", "fnirs_oxy_highpass.png", "fnirs_deoxy_highpass.png"
    }
    assert plt.get_fignums() == []


def test_save_plot_requires_output_directory():
    with pytest.raises(ValueError, match="output_dir"):
        processor().process(*signals(), save_plots=True)

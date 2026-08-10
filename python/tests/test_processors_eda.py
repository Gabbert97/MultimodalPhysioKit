from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pytest

from multimodalphysiokit.core import ProcessingResult, Signal
from multimodalphysiokit.processors import EDAProcessor, process_eda


FEATURE_NAMES = (
    "eda_scl_mean", "eda_scl_std", "eda_scl_slope", "eda_scr_amplitude_mean",
    "eda_scr_amplitude_std", "eda_scr_rise_time_mean", "eda_scr_rise_time_std",
    "eda_scr_per_min",
)


def make_signal(*, label="eda", duration=70.0, sampling_frequency=1000.0):
    time = np.arange(0, duration, 1 / sampling_frequency)
    samples = 2.0 + 0.1 * np.sin(2 * np.pi * 0.02 * time)
    return Signal(samples, time, sampling_frequency, label=label, units="uS")


def fake_cvxeda(eda, delta):
    assert delta == pytest.approx(0.01)
    time = np.arange(len(eda)) / 100
    scr = sum((0.4 * np.exp(-0.5 * ((time - center) / 0.35) ** 2)
               for center in np.arange(2.0, time[-1], 4.0)), np.zeros(len(eda)))
    scl = 1.5 + 0.001 * time + 0.02 * np.sin(2 * np.pi * 0.01 * time)
    zeros = np.zeros(len(eda))
    return scr, zeros, scl, zeros, zeros, zeros, 0.0


def processor():
    return EDAProcessor(decomposition_callable=fake_cvxeda)


def test_constructs_without_signal():
    instance = processor()
    assert isinstance(instance, EDAProcessor)
    assert instance.minimum_scr_amplitude == pytest.approx(0.01)


@pytest.mark.parametrize("value", [-1, np.nan, np.inf])
def test_minimum_scr_amplitude_is_validated(value):
    with pytest.raises(ValueError, match="minimum_scr_amplitude"):
        EDAProcessor(decomposition_callable=fake_cvxeda, minimum_scr_amplitude=value)


def test_configured_minimum_scr_amplitude_reaches_feature_screening():
    result = EDAProcessor(
        decomposition_callable=fake_cvxeda, minimum_scr_amplitude=1.0
    ).process(make_signal())
    assert result.features["eda_scr_amplitude_mean"] == 0
    assert result.features["eda_scr_amplitude_std"] == 0
    assert result.features["eda_scr_rise_time_mean"] == 0
    assert result.features["eda_scr_rise_time_std"] == 0
    assert result.features["eda_scr_per_min"] == 0


def test_process_and_wrapper_return_processing_result():
    signal = make_signal()
    direct = processor().process(signal)
    wrapped = process_eda(signal, decomposition_callable=fake_cvxeda)
    assert isinstance(direct, ProcessingResult)
    assert wrapped.features == pytest.approx(direct.features)


def test_exact_feature_names_and_values_are_finite():
    result = processor().process(make_signal())
    assert tuple(result.features) == FEATURE_NAMES
    assert len(result.features) == 8
    assert all("__" not in key for key in result.features)
    assert all(np.isfinite(value) for value in result.features.values())


def test_signal_label_and_fallback():
    assert processor().process(make_signal(label="skin")).signal_labels == ["skin"]
    assert processor().process(make_signal(label="")).signal_labels == ["eda"]


def test_metadata_is_signal_level_only():
    result = processor().process(make_signal())
    assert set(result.metadata) == {
        "source", "original_sampling_frequency", "processed_sampling_frequency",
        "original_n_samples", "processed_n_samples", "start_time_seconds",
        "end_time_seconds", "cvxeda_source",
        "minimum_scr_amplitude_us",
    }
    assert not any("phase" in key for key in result.metadata)
    assert result.intermediates == {}
    assert result.metadata["minimum_scr_amplitude_us"] == pytest.approx(0.01)


def test_extract_features_matches_process():
    signal = make_signal()
    assert processor().extract_features(signal) == pytest.approx(processor().process(signal).features)


def test_invalid_signal_type_and_frequency():
    with pytest.raises(TypeError, match="Signal"):
        processor().process("bad")
    with pytest.raises(ValueError, match="integer multiple of 100"):
        processor().process(make_signal(sampling_frequency=250.0))


def test_zero_duration_is_rejected():
    with pytest.raises(ValueError, match="monotonically increasing"):
        Signal(np.ones(3), np.zeros(3), 1000.0)


def test_plotting_is_numerically_inert_and_closes_figures(tmp_path):
    plt.close("all")
    signal = make_signal()
    baseline = processor().process(signal)
    plotted = processor().process(signal, save_plots=True, output_dir=tmp_path)
    assert plotted.features == pytest.approx(baseline.features)
    assert (tmp_path / "eda_eda_decomposition.png").is_file()
    assert plt.get_fignums() == []


def test_save_plot_requires_output_directory():
    with pytest.raises(ValueError, match="output_dir"):
        processor().process(make_signal(), save_plots=True)

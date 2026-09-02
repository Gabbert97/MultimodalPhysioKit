import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pytest

from multimodalphysiokit import ProcessingResult, Signal
from multimodalphysiokit.processors import RespirationProcessor, process_respiration


def make_respiration_signal(duration_seconds=120, sampling_frequency=200):
    time = np.arange(0, duration_seconds, 1 / sampling_frequency)
    angle = 2 * np.pi * (
        0.24 * time + 0.015 * np.sin(2 * np.pi * 0.02 * time)
    )
    samples = np.sin(angle) + 0.05 * np.sin(2 * np.pi * 0.9 * time)
    return Signal(samples, time, sampling_frequency, label="respiration")


def test_respiration_processor_can_be_constructed_without_signal():
    assert isinstance(RespirationProcessor(), RespirationProcessor)


def test_default_breathing_rate_range_preserves_existing_behavior():
    processor = RespirationProcessor()

    assert processor.min_breathing_rate == 5.0
    assert processor.max_breathing_rate == 30.0


@pytest.mark.parametrize(
    "minimum, maximum, message",
    [
        (0, 30, "min_breathing_rate"),
        (-1, 30, "min_breathing_rate"),
        (np.nan, 30, "min_breathing_rate"),
        (5, 5, "max_breathing_rate"),
        (30, 5, "max_breathing_rate"),
        (5, np.inf, "max_breathing_rate"),
    ],
)
def test_invalid_breathing_rate_range_is_rejected(minimum, maximum, message):
    with pytest.raises(ValueError, match=message):
        RespirationProcessor(minimum, maximum)


def test_custom_breathing_rate_range_is_forwarded_to_processing(monkeypatch):
    observed = {}

    def fake_process(signal, **kwargs):
        observed.update(kwargs)
        return ProcessingResult("respiration", "fake", {}, ["respiration"])

    monkeypatch.setattr(
        "multimodalphysiokit.processors.respiration.process_respiration",
        fake_process,
    )
    processor = RespirationProcessor(min_breathing_rate=3, max_breathing_rate=30)
    processor.process(make_respiration_signal())

    assert observed["min_breathing_rate"] == 3.0
    assert observed["max_breathing_rate"] == 30.0


def test_three_to_thirty_range_accepts_coherent_slow_respiration():
    sampling_frequency = 200
    duration_seconds = 180
    time = np.arange(0, duration_seconds, 1 / sampling_frequency)
    signal = Signal(
        np.sin(2 * np.pi * (4.4 / 60) * time),
        time,
        sampling_frequency,
        label="respiration",
    )

    with pytest.raises(ValueError, match="too few valid breaths detected"):
        RespirationProcessor().process(signal)

    result = RespirationProcessor(
        min_breathing_rate=3.0,
        max_breathing_rate=30.0,
    ).process(signal)

    assert result.features["respiration_breath_rate_mean"] == pytest.approx(
        4.4, rel=1e-4
    )
    functional_result = process_respiration(
        signal,
        min_breathing_rate=3.0,
        max_breathing_rate=30.0,
    )
    assert functional_result.features == result.features


def test_process_signal_returns_processing_result():
    result = RespirationProcessor().process(make_respiration_signal())

    assert isinstance(result, ProcessingResult)


def test_extract_features_matches_process_features():
    signal = make_respiration_signal()
    processor = RespirationProcessor()

    assert processor.extract_features(signal) == processor.process(signal).features


def test_process_matches_functional_wrapper():
    signal = make_respiration_signal()

    assert (
        RespirationProcessor().process(signal).features
        == process_respiration(signal).features
    )


def test_invalid_signal_raises_type_error_through_processor():
    with pytest.raises(TypeError, match="Signal"):
        RespirationProcessor().process([1, 2, 3])


def test_save_plots_without_output_dir_raises_value_error():
    with pytest.raises(ValueError, match="output_dir"):
        RespirationProcessor().process(make_respiration_signal(), save_plots=True)


def test_save_plots_writes_existing_files(tmp_path):
    plt.close("all")

    RespirationProcessor().process(
        make_respiration_signal(),
        save_plots=True,
        output_dir=tmp_path,
    )

    assert (tmp_path / "respiration_downsampled_signal.png").is_file()
    assert (tmp_path / "respiration_filtered_phase_peaks.png").is_file()
    assert plt.get_fignums() == []


def test_plotting_does_not_change_features(tmp_path):
    signal = make_respiration_signal()
    processor = RespirationProcessor()

    without_plots = processor.process(signal)
    with_plots = processor.process(signal, save_plots=True, output_dir=tmp_path)

    assert with_plots.features == without_plots.features

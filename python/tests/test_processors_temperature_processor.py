import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pytest

from multimodalphysiokit import ProcessingResult, Signal
from multimodalphysiokit.processors import TemperatureProcessor, process_temperature


def make_temperature_signal(n_samples=240, sampling_frequency=10):
    time = np.arange(n_samples) / sampling_frequency
    samples = 32 + 0.02 * time + 0.1 * np.sin(2 * np.pi * 0.05 * time)
    return Signal(samples, time, sampling_frequency, label="temp", units="degC")


def test_temperature_processor_can_be_constructed_without_signal():
    processor = TemperatureProcessor()

    assert isinstance(processor, TemperatureProcessor)


def test_process_signal_returns_processing_result():
    result = TemperatureProcessor().process(make_temperature_signal())

    assert isinstance(result, ProcessingResult)


def test_extract_features_returns_same_dictionary_as_process():
    signal = make_temperature_signal()
    processor = TemperatureProcessor()

    assert processor.extract_features(signal) == processor.process(signal).features


def test_process_returns_same_features_as_functional_wrapper():
    signal = make_temperature_signal()

    result = TemperatureProcessor().process(signal)
    expected = process_temperature(signal)

    assert result.features == expected.features


def test_invalid_signal_raises_type_error_through_processor():
    with pytest.raises(TypeError, match="Signal"):
        TemperatureProcessor().process([1, 2, 3])


def test_save_plots_without_output_dir_raises_value_error():
    with pytest.raises(ValueError, match="output_dir"):
        TemperatureProcessor().process(make_temperature_signal(), save_plots=True)


def test_save_plots_writes_existing_temperature_png_files(tmp_path):
    plt.close("all")

    TemperatureProcessor().process(
        make_temperature_signal(),
        save_plots=True,
        output_dir=tmp_path,
    )

    assert (tmp_path / "temperature_vs_time.png").is_file()
    assert (tmp_path / "temperature_first_derivative_vs_time.png").is_file()
    assert plt.get_fignums() == []


def test_plotting_does_not_change_feature_values(tmp_path):
    signal = make_temperature_signal()
    processor = TemperatureProcessor()

    without_plots = processor.process(signal)
    with_plots = processor.process(signal, save_plots=True, output_dir=tmp_path)

    assert with_plots.features == without_plots.features

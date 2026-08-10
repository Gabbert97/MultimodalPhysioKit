import numpy as np
import pytest

from multimodalphysiokit import ProcessingResult, Signal
from multimodalphysiokit.processors import ECGProcessor, detect_ecg_r_peaks, process_ecg


def make_synthetic_ecg(duration_seconds=40, sampling_frequency=250):
    time = np.arange(0, duration_seconds, 1 / sampling_frequency)
    samples = 0.02 * np.sin(2 * np.pi * 0.5 * time)
    for beat_time in np.arange(0.8, duration_seconds, 1.0):
        samples += np.exp(-0.5 * ((time - beat_time) / 0.025) ** 2)
        samples += 0.15 * np.exp(-0.5 * ((time - (beat_time + 0.22)) / 0.06) ** 2)
    return Signal(samples, time, sampling_frequency, label="ecg", units="mV")


def test_processor_constructs_without_signal_and_stores_configuration():
    processor = ECGProcessor(label_frequency_analysis=2, return_intermediates=True)
    assert processor.label_frequency_analysis == 2
    assert processor.return_intermediates is True


def test_invalid_constructor_configuration_is_rejected():
    with pytest.raises(ValueError, match="label_frequency_analysis"):
        ECGProcessor(label_frequency_analysis=3)
    with pytest.raises(TypeError, match="return_intermediates"):
        ECGProcessor(return_intermediates=1)


def test_detect_r_peaks_matches_functional_api():
    signal = make_synthetic_ecg()
    class_result = ECGProcessor().detect_r_peaks(signal)
    function_result = detect_ecg_r_peaks(signal)
    assert set(class_result) == set(function_result)
    for key in class_result:
        np.testing.assert_array_equal(class_result[key], function_result[key])


def test_extract_time_domain_features_accepts_signal():
    result = ECGProcessor(label_frequency_analysis=2).extract_time_domain_features(
        make_synthetic_ecg()
    )
    assert len(result["features"]) == 5


def test_extract_frequency_domain_features_accepts_signal():
    features = ECGProcessor(label_frequency_analysis=2).extract_frequency_domain_features(
        make_synthetic_ecg()
    )
    assert len(features) == 6


def test_process_returns_processing_result_and_matches_wrapper():
    signal = make_synthetic_ecg()
    result = ECGProcessor(label_frequency_analysis=2).process(signal)
    expected = process_ecg(signal, label_frequency_analysis=2)
    assert isinstance(result, ProcessingResult)
    assert result.features == expected.features


def test_invalid_signal_is_rejected_by_process():
    with pytest.raises(TypeError, match="Signal"):
        ECGProcessor(label_frequency_analysis=2).process([1, 2, 3])

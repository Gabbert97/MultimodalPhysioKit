import numpy as np
import pytest

from multimodalphysiokit import ProcessingResult


def test_valid_construction():
    result = ProcessingResult(
        modality="ECG",
        processor_name="ECGProcessor",
        features={"bpm": 72},
        signal_labels="ecg",
        options={"show_plots": False},
        metadata={"legacy_function": "ecg_signal_processing"},
        intermediates={"r_peak_indexes": [1, 4, 8]},
    )

    assert result.modality == "ECG"
    assert result.processor_name == "ECGProcessor"
    assert result.features == {"bpm": 72}
    assert result.signal_labels == ["ecg"]
    assert result.options == {"show_plots": False}
    assert result.metadata == {"legacy_function": "ecg_signal_processing"}
    assert result.intermediates == {"r_peak_indexes": [1, 4, 8]}


def test_default_containers_are_empty():
    result = ProcessingResult("EDA", "EDAProcessor")

    assert result.features == {}
    assert result.signal_labels == []
    assert result.options == {}
    assert result.metadata == {}
    assert result.intermediates == {}


def test_explicit_intermediates_store_arrays_and_nested_dictionaries():
    filtered = np.arange(5, dtype=float)
    diagnostics = {"method": "synthetic", "details": {"valid": True}}

    result = ProcessingResult(
        "ECG",
        "ECGProcessor",
        intermediates={"filtered_signal": filtered, "diagnostics": diagnostics},
    )

    assert result.intermediates["filtered_signal"] is filtered
    assert result.intermediates["diagnostics"] is diagnostics


@pytest.mark.parametrize("intermediates", [None, [1, 2], "invalid", 3])
def test_invalid_intermediates_raise_type_error(intermediates):
    with pytest.raises(TypeError, match="intermediates must be a mapping"):
        ProcessingResult("ECG", "ECGProcessor", intermediates=intermediates)


def test_intermediates_top_level_dictionary_is_shallow_copied():
    array = np.arange(3)
    supplied = {"array": array}

    result = ProcessingResult("ECG", "ECGProcessor", intermediates=supplied)
    supplied["added_later"] = "value"

    assert set(result.intermediates) == {"array"}
    assert result.intermediates["array"] is array


def test_large_intermediate_array_is_not_expanded_in_repr():
    large = np.arange(10_000)

    result = ProcessingResult(
        "ECG",
        "ECGProcessor",
        intermediates={"large_array": large},
    )

    representation = repr(result)
    assert "large_array" not in representation
    assert str(large[-1]) not in representation


def test_feature_names_works():
    result = ProcessingResult("Temperature", "TemperatureProcessor", features={"mean": 1, "std": 2})

    assert result.feature_names() == ["mean", "std"]


def test_has_feature_works():
    result = ProcessingResult("Respiration", "RespirationProcessor", features={"rate": 12})

    assert result.has_feature("rate")
    assert not result.has_feature("missing")

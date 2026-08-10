import numpy as np
import pytest

from multimodalphysiokit.core import Recording, Signal
from multimodalphysiokit.io import recording_from_arrays


def test_valid_construction_with_sampling_frequency():
    recording = recording_from_arrays({"ecg": [1, 2, 3]}, sampling_frequency=2)

    assert isinstance(recording, Recording)
    signal = recording.get_signal("ecg")
    assert isinstance(signal, Signal)
    np.testing.assert_allclose(signal.samples, [1, 2, 3])
    np.testing.assert_allclose(signal.time, [0, 0.5, 1.0])
    assert signal.sampling_frequency == 2


def test_valid_construction_with_explicit_time():
    recording = recording_from_arrays({"eda": [1, 2, 3]}, time=[10, 11, 12])

    signal = recording.get_signal("eda")
    np.testing.assert_allclose(signal.time, [10, 11, 12])
    assert signal.sampling_frequency == pytest.approx(1.0)


def test_multiple_signals_can_have_different_lengths_with_sampling_frequency():
    recording = recording_from_arrays(
        {"ecg": [1, 2, 3], "eda": [4, 5]},
        sampling_frequency=10,
    )

    assert recording.signal_names() == ["ecg", "eda"]
    np.testing.assert_allclose(recording.get_signal("ecg").time, [0, 0.1, 0.2])
    np.testing.assert_allclose(recording.get_signal("eda").time, [0, 0.1])


def test_units_dictionary():
    recording = recording_from_arrays(
        {"ecg": [1, 2], "eda": [3, 4]},
        sampling_frequency=1,
        units={"ecg": "mV", "eda": "uS"},
    )

    assert recording.get_signal("ecg").units == "mV"
    assert recording.get_signal("eda").units == "uS"


def test_labels_dictionary():
    recording = recording_from_arrays(
        {"temperature": [32, 33]},
        sampling_frequency=1,
        labels={"temperature": "temp"},
    )

    assert recording.get_signal("temperature").label == "temp"


def test_missing_unit_defaults_to_empty_string():
    recording = recording_from_arrays(
        {"ecg": [1, 2], "eda": [3, 4]},
        sampling_frequency=1,
        units={"ecg": "mV"},
    )

    assert recording.get_signal("eda").units == ""


def test_missing_label_defaults_to_signal_name():
    recording = recording_from_arrays(
        {"ecg": [1, 2], "eda": [3, 4]},
        sampling_frequency=1,
        labels={"ecg": "lead_1"},
    )

    assert recording.get_signal("ecg").label == "lead_1"
    assert recording.get_signal("eda").label == "eda"


def test_metadata_is_passed_to_recording():
    recording = recording_from_arrays(
        {"ecg": [1, 2]},
        sampling_frequency=1,
        name="session_001",
        metadata={"source": "array"},
    )

    assert recording.name == "session_001"
    assert recording.metadata == {"source": "array"}


def test_neither_sampling_frequency_nor_time_raises():
    with pytest.raises(ValueError, match="exactly one"):
        recording_from_arrays({"ecg": [1, 2]})


def test_both_sampling_frequency_and_time_raises():
    with pytest.raises(ValueError, match="exactly one"):
        recording_from_arrays({"ecg": [1, 2]}, sampling_frequency=1, time=[0, 1])


def test_empty_signals_dictionary_raises():
    with pytest.raises(ValueError, match="non-empty"):
        recording_from_arrays({}, sampling_frequency=1)


def test_non_dictionary_signals_raises_type_error():
    with pytest.raises(TypeError, match="dictionary"):
        recording_from_arrays([1, 2, 3], sampling_frequency=1)


def test_empty_signal_name_raises():
    with pytest.raises(ValueError, match="signal names"):
        recording_from_arrays({"": [1, 2]}, sampling_frequency=1)


@pytest.mark.parametrize(
    "samples",
    [
        [[1, 2], [3, 4]],
        [],
        ["a", "b"],
    ],
)
def test_invalid_sample_array_raises(samples):
    with pytest.raises(ValueError):
        recording_from_arrays({"ecg": samples}, sampling_frequency=1)


def test_non_finite_samples_raise():
    with pytest.raises(ValueError, match="finite"):
        recording_from_arrays({"ecg": [1, np.nan]}, sampling_frequency=1)


@pytest.mark.parametrize("sampling_frequency", [0, -1, np.inf, np.nan])
def test_invalid_sampling_frequency_raises(sampling_frequency):
    with pytest.raises(ValueError, match="positive and finite"):
        recording_from_arrays({"ecg": [1, 2]}, sampling_frequency=sampling_frequency)


def test_explicit_time_length_mismatch_raises():
    with pytest.raises(ValueError, match="length"):
        recording_from_arrays({"ecg": [1, 2, 3]}, time=[0, 1])


def test_explicit_non_monotonic_time_raises():
    with pytest.raises(ValueError, match="monotonically increasing"):
        recording_from_arrays({"ecg": [1, 2, 3]}, time=[0, 1, 1])

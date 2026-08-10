from __future__ import annotations

import h5py
import numpy as np
import pytest

from multimodalphysiokit.io import read_generic_hdf5


def write_hdf5(path, writer):
    with h5py.File(path, "w") as hdf5_file:
        writer(hdf5_file)


def test_reads_one_signal_with_explicit_sampling_frequency(tmp_path):
    file_path = tmp_path / "sample.h5"
    write_hdf5(file_path, lambda h5: h5.create_dataset("/signals/ecg", data=np.arange(5)))

    recording = read_generic_hdf5(file_path, {"ecg": "/signals/ecg"}, sampling_frequency=100)

    signal = recording.get_signal("ecg")
    np.testing.assert_allclose(signal.samples, np.arange(5))
    np.testing.assert_allclose(signal.time, np.arange(5) / 100)
    assert signal.sampling_frequency == 100


def test_reads_multiple_signals_from_nested_paths(tmp_path):
    file_path = tmp_path / "nested.hdf5"

    def writer(h5):
        h5.create_dataset("/subject/session/ecg", data=np.arange(4))
        h5.create_dataset("/subject/session/eda", data=np.arange(4) + 10)

    write_hdf5(file_path, writer)

    recording = read_generic_hdf5(
        file_path,
        {"ecg": "/subject/session/ecg", "eda": "/subject/session/eda"},
        sampling_frequency=50,
    )

    assert set(recording.signal_names()) == {"ecg", "eda"}
    np.testing.assert_allclose(recording.get_signal("eda").samples, np.arange(4) + 10)


def test_reads_sampling_frequency_from_path(tmp_path):
    file_path = tmp_path / "fs.h5"

    def writer(h5):
        h5.create_dataset("/signals/temp", data=np.arange(3))
        h5.create_dataset("/metadata/fs", data=np.asarray([25.0]))

    write_hdf5(file_path, writer)

    recording = read_generic_hdf5(
        file_path,
        {"temp": "/signals/temp"},
        sampling_frequency_path="/metadata/fs",
    )

    assert recording.get_signal("temp").sampling_frequency == 25.0


def test_reads_explicit_time_vector_from_path(tmp_path):
    file_path = tmp_path / "time.h5"

    def writer(h5):
        h5.create_dataset("/signals/ecg", data=np.arange(4))
        h5.create_dataset("/time", data=np.asarray([0.0, 0.1, 0.2, 0.3]))

    write_hdf5(file_path, writer)

    recording = read_generic_hdf5(file_path, {"ecg": "/signals/ecg"}, time_path="/time")

    np.testing.assert_allclose(recording.get_signal("ecg").time, [0.0, 0.1, 0.2, 0.3])


def test_units_and_labels_are_applied(tmp_path):
    file_path = tmp_path / "labels.h5"
    write_hdf5(file_path, lambda h5: h5.create_dataset("/ecg", data=np.arange(3)))

    recording = read_generic_hdf5(
        file_path,
        {"ecg": "/ecg"},
        sampling_frequency=100,
        units={"ecg": "mV"},
        labels={"ecg": "ECG Lead I"},
    )

    signal = recording.get_signal("ecg")
    assert signal.units == "mV"
    assert signal.label == "ECG Lead I"


def test_recording_name_defaults_to_file_stem(tmp_path):
    file_path = tmp_path / "subject_1.hdf5"
    write_hdf5(file_path, lambda h5: h5.create_dataset("/ecg", data=np.arange(3)))

    recording = read_generic_hdf5(file_path, {"ecg": "/ecg"}, sampling_frequency=100)

    assert recording.name == "subject_1"


def test_custom_recording_name_is_respected(tmp_path):
    file_path = tmp_path / "subject_1.hdf5"
    write_hdf5(file_path, lambda h5: h5.create_dataset("/ecg", data=np.arange(3)))

    recording = read_generic_hdf5(file_path, {"ecg": "/ecg"}, sampling_frequency=100, name="custom")

    assert recording.name == "custom"


def test_metadata_contains_source_file_name_and_user_metadata(tmp_path):
    file_path = tmp_path / "meta.h5"
    write_hdf5(file_path, lambda h5: h5.create_dataset("/ecg", data=np.arange(3)))

    recording = read_generic_hdf5(
        file_path,
        {"ecg": "/ecg"},
        sampling_frequency=100,
        metadata={"dataset": "test"},
    )

    assert recording.metadata["source"] == "read_generic_hdf5"
    assert recording.metadata["file_name"] == "meta.h5"
    assert recording.metadata["dataset"] == "test"
    assert "file_path" not in recording.metadata


def test_missing_file_raises_file_not_found_error(tmp_path):
    with pytest.raises(FileNotFoundError):
        read_generic_hdf5(tmp_path / "missing.h5", {"ecg": "/ecg"}, sampling_frequency=100)


def test_invalid_signal_paths_type_raises_type_error(tmp_path):
    file_path = tmp_path / "invalid.h5"
    write_hdf5(file_path, lambda h5: h5.create_dataset("/ecg", data=np.arange(3)))

    with pytest.raises(TypeError, match="signal_paths"):
        read_generic_hdf5(file_path, ["/ecg"], sampling_frequency=100)


def test_empty_signal_paths_raises_value_error(tmp_path):
    file_path = tmp_path / "empty.h5"
    write_hdf5(file_path, lambda h5: h5.create_dataset("/ecg", data=np.arange(3)))

    with pytest.raises(ValueError, match="signal_paths"):
        read_generic_hdf5(file_path, {}, sampling_frequency=100)


def test_empty_signal_name_raises_value_error(tmp_path):
    file_path = tmp_path / "empty_name.h5"
    write_hdf5(file_path, lambda h5: h5.create_dataset("/ecg", data=np.arange(3)))

    with pytest.raises(ValueError, match="signal names"):
        read_generic_hdf5(file_path, {"": "/ecg"}, sampling_frequency=100)


def test_empty_dataset_path_raises_value_error(tmp_path):
    file_path = tmp_path / "empty_path.h5"
    write_hdf5(file_path, lambda h5: h5.create_dataset("/ecg", data=np.arange(3)))

    with pytest.raises(ValueError, match="dataset paths"):
        read_generic_hdf5(file_path, {"ecg": ""}, sampling_frequency=100)


def test_missing_signal_dataset_raises_key_error(tmp_path):
    file_path = tmp_path / "missing_dataset.h5"
    write_hdf5(file_path, lambda h5: h5.create_dataset("/eda", data=np.arange(3)))

    with pytest.raises(KeyError, match="missing HDF5 dataset"):
        read_generic_hdf5(file_path, {"ecg": "/ecg"}, sampling_frequency=100)


def test_non_1d_signal_dataset_raises_value_error(tmp_path):
    file_path = tmp_path / "matrix.h5"
    write_hdf5(file_path, lambda h5: h5.create_dataset("/ecg", data=np.ones((2, 2))))

    with pytest.raises(ValueError, match="one-dimensional"):
        read_generic_hdf5(file_path, {"ecg": "/ecg"}, sampling_frequency=100)


def test_non_numeric_signal_dataset_raises_value_error(tmp_path):
    file_path = tmp_path / "strings.h5"
    write_hdf5(file_path, lambda h5: h5.create_dataset("/ecg", data=np.asarray([b"a", b"b"])))

    with pytest.raises(ValueError, match="numeric"):
        read_generic_hdf5(file_path, {"ecg": "/ecg"}, sampling_frequency=100)


def test_non_finite_signal_dataset_raises_value_error(tmp_path):
    file_path = tmp_path / "nan.h5"
    write_hdf5(file_path, lambda h5: h5.create_dataset("/ecg", data=np.asarray([1.0, np.nan])))

    with pytest.raises(ValueError, match="finite"):
        read_generic_hdf5(file_path, {"ecg": "/ecg"}, sampling_frequency=100)


@pytest.mark.parametrize(
    "kwargs",
    [
        {},
        {"sampling_frequency": 100, "time_path": "/time"},
        {"sampling_frequency": 100, "sampling_frequency_path": "/fs"},
        {"sampling_frequency_path": "/fs", "time_path": "/time"},
    ],
)
def test_invalid_timing_selection_raises_value_error(tmp_path, kwargs):
    file_path = tmp_path / "timing.h5"

    def writer(h5):
        h5.create_dataset("/ecg", data=np.arange(3))
        h5.create_dataset("/time", data=np.arange(3))
        h5.create_dataset("/fs", data=100.0)

    write_hdf5(file_path, writer)

    with pytest.raises(ValueError, match="exactly one timing source"):
        read_generic_hdf5(file_path, {"ecg": "/ecg"}, **kwargs)


def test_invalid_sampling_frequency_value_raises_value_error(tmp_path):
    file_path = tmp_path / "bad_fs_value.h5"
    write_hdf5(file_path, lambda h5: h5.create_dataset("/ecg", data=np.arange(3)))

    with pytest.raises(ValueError, match="sampling_frequency"):
        read_generic_hdf5(file_path, {"ecg": "/ecg"}, sampling_frequency=0)


def test_missing_sampling_frequency_path_raises_key_error(tmp_path):
    file_path = tmp_path / "missing_fs.h5"
    write_hdf5(file_path, lambda h5: h5.create_dataset("/ecg", data=np.arange(3)))

    with pytest.raises(KeyError, match="sampling frequency"):
        read_generic_hdf5(file_path, {"ecg": "/ecg"}, sampling_frequency_path="/fs")


def test_non_scalar_sampling_frequency_path_raises_value_error(tmp_path):
    file_path = tmp_path / "vector_fs.h5"

    def writer(h5):
        h5.create_dataset("/ecg", data=np.arange(3))
        h5.create_dataset("/fs", data=np.asarray([100.0, 200.0]))

    write_hdf5(file_path, writer)

    with pytest.raises(ValueError, match="scalar"):
        read_generic_hdf5(file_path, {"ecg": "/ecg"}, sampling_frequency_path="/fs")


def test_non_numeric_sampling_frequency_path_raises_value_error(tmp_path):
    file_path = tmp_path / "string_fs.h5"

    def writer(h5):
        h5.create_dataset("/ecg", data=np.arange(3))
        h5.create_dataset("/fs", data=np.asarray([b"100"]))

    write_hdf5(file_path, writer)

    with pytest.raises(ValueError, match="numeric"):
        read_generic_hdf5(file_path, {"ecg": "/ecg"}, sampling_frequency_path="/fs")


def test_time_length_mismatch_raises_value_error(tmp_path):
    file_path = tmp_path / "time_mismatch.h5"

    def writer(h5):
        h5.create_dataset("/ecg", data=np.arange(4))
        h5.create_dataset("/time", data=np.arange(3))

    write_hdf5(file_path, writer)

    with pytest.raises(ValueError, match="length must match"):
        read_generic_hdf5(file_path, {"ecg": "/ecg"}, time_path="/time")

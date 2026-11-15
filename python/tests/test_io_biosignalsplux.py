from __future__ import annotations

import h5py
import numpy as np
import pytest

from multimodalphysiokit.io import read_biosignalsplux_hdf5


def write_biosignalsplux_file(
    file_path,
    *,
    channels=(1,),
    sensors=("ECG",),
    data=None,
    resolution=16,
    sampling_frequency=1000,
    fs_attr_name="sampling rate",
    first_group_name="recording_1",
):
    if data is None:
        data = {1: np.asarray([0, 32768, 65535], dtype=float)}
    with h5py.File(file_path, "w") as h5:
        group = h5.create_group(first_group_name)
        group.attrs["channels"] = np.asarray(channels)
        group.attrs["resolution"] = resolution
        group.attrs[fs_attr_name] = sampling_frequency
        raw_group = group.create_group("raw")
        for channel_id, sensor in zip(channels, sensors, strict=True):
            dataset = raw_group.create_dataset(f"channel_{channel_id}", data=np.asarray(data[channel_id]))
            dataset.attrs["sensor"] = sensor


def expected_ecg(raw, resolution=16):
    return ((np.asarray(raw) / (2**resolution) - 0.5) * 3 / 1019) * 1e3


def expected_eda(raw, resolution=16):
    return np.asarray(raw) * 3 / (0.12 * 2**resolution)


def expected_temp(raw, resolution=16):
    raw = np.asarray(raw)
    vcc = 3.0
    ntc_v = raw * vcc / (2**resolution)
    ntc_ohm = (ntc_v * 1e4) / (vcc - ntc_v)
    a0 = 1.12764514e-3
    a1 = 2.34282709e-4
    a2 = 8.77303013e-8
    tmp_k = 1.0 / (a0 + a1 * np.log(ntc_ohm) + a2 * np.log(ntc_ohm) ** 3)
    return tmp_k - 273.15


def expected_rip(raw, resolution=16):
    return (np.asarray(raw) / (2**resolution) - 0.5) * 100


def expected_respiration(raw, resolution=16):
    return (np.asarray(raw) / (2**resolution - 1) - 0.5) * 3


def expected_fnirs(raw, resolution=16):
    return 0.15 * np.asarray(raw) / (2**resolution)


def test_reads_first_top_level_recording_group(tmp_path):
    file_path = tmp_path / "first.h5"
    with h5py.File(file_path, "w") as h5:
        h5.create_dataset("not_a_group", data=np.arange(3))
        group = h5.create_group("recording_a")
        group.attrs["channels"] = np.asarray([1])
        group.attrs["resolution"] = np.asarray([16])
        group.attrs["sampling rate"] = 1000
        raw = group.create_group("raw")
        dataset = raw.create_dataset("channel_1", data=np.asarray([0, 32768, 65535]))
        dataset.attrs["sensor"] = "ECG"

    recording = read_biosignalsplux_hdf5(file_path)

    assert recording.metadata["dataset"] == "/recording_a"
    assert recording.has_signal("ecg")


def test_reads_sampling_frequency_from_sampling_rate(tmp_path):
    file_path = tmp_path / "fs_sampling_rate.h5"
    write_biosignalsplux_file(file_path, sampling_frequency=500)

    recording = read_biosignalsplux_hdf5(file_path)

    assert recording.get_signal("ecg").sampling_frequency == 500


def test_reads_sampling_frequency_from_fs_fallback(tmp_path):
    file_path = tmp_path / "fs_fallback.h5"
    write_biosignalsplux_file(file_path, sampling_frequency=250, fs_attr_name="fs")

    recording = read_biosignalsplux_hdf5(file_path)

    assert recording.get_signal("ecg").sampling_frequency == 250


def test_broadcasts_scalar_resolution_to_multiple_channels(tmp_path):
    file_path = tmp_path / "scalar_resolution.h5"
    data = {1: np.asarray([0, 32768]), 2: np.asarray([0, 1000])}
    write_biosignalsplux_file(
        file_path,
        channels=(1, 2),
        sensors=("ECG", "EDA"),
        data=data,
        resolution=16,
    )

    recording = read_biosignalsplux_hdf5(file_path)

    assert recording.get_signal("ecg").metadata["adc_resolution_bits"] == 16
    assert recording.get_signal("eda").metadata["adc_resolution_bits"] == 16


def test_creates_ecg_signal_from_ecg_channel(tmp_path):
    file_path = tmp_path / "ecg.h5"
    raw = np.asarray([0, 32768, 65535])
    write_biosignalsplux_file(file_path, data={1: raw})

    signal = read_biosignalsplux_hdf5(file_path).get_signal("ecg")

    np.testing.assert_allclose(signal.samples, expected_ecg(raw))
    assert signal.units == "mV"


def test_creates_eda_signal_from_eda_channel(tmp_path):
    file_path = tmp_path / "eda.h5"
    raw = np.asarray([0, 1000, 2000])
    write_biosignalsplux_file(file_path, sensors=("EDA",), data={1: raw})

    signal = read_biosignalsplux_hdf5(file_path).get_signal("eda")

    np.testing.assert_allclose(signal.samples, expected_eda(raw))
    assert signal.units == "uS"


def test_creates_temp_signal_from_temp_channel(tmp_path):
    file_path = tmp_path / "temp.h5"
    raw = np.asarray([20000, 22000, 24000])
    write_biosignalsplux_file(file_path, sensors=("TEMP",), data={1: raw})

    signal = read_biosignalsplux_hdf5(file_path).get_signal("temp")

    np.testing.assert_allclose(signal.samples, expected_temp(raw))
    assert signal.units == "degC"


def test_creates_rip_signal_from_rip_channel(tmp_path):
    file_path = tmp_path / "rip.h5"
    raw = np.asarray([0, 32768, 65535])
    write_biosignalsplux_file(file_path, sensors=("RIP",), data={1: raw})

    signal = read_biosignalsplux_hdf5(file_path).get_signal("rip")

    np.testing.assert_allclose(signal.samples, expected_rip(raw))
    assert signal.units == "percent full scale"


def test_creates_ppg_signal_from_bvp_channel(tmp_path):
    # BVP (blood volume pulse) is PLUX's PPG sensor. PPG amplitude has no
    # standard physical unit, so samples are percent of the ADC full scale,
    # centred at mid-scale: 0 -> -50 %, 2^15 -> 0 %, 2^16 -> +50 %.
    file_path = tmp_path / "bvp.h5"
    raw = np.asarray([0, 32768, 65536])
    write_biosignalsplux_file(file_path, sensors=("BVP",), data={1: raw})

    signal = read_biosignalsplux_hdf5(file_path).get_signal("ppg")

    np.testing.assert_allclose(signal.samples, [-50.0, 0.0, 50.0])
    assert signal.units == "percent full scale"


def test_creates_respiration_signal_from_respiration_channel(tmp_path):
    file_path = tmp_path / "respiration.h5"
    raw = np.asarray([0, 32768, 65535])
    write_biosignalsplux_file(file_path, sensors=("RESPIRATION",), data={1: raw})

    signal = read_biosignalsplux_hdf5(file_path).get_signal("respiration")

    np.testing.assert_allclose(signal.samples, expected_respiration(raw))
    assert signal.units == "V"


def test_creates_fnirs_red_and_infrared_from_hspo2_channels(tmp_path):
    file_path = tmp_path / "fnirs.h5"
    raw_1 = np.asarray([100, 200, 300])
    raw_2 = np.asarray([400, 500, 600])
    write_biosignalsplux_file(
        file_path,
        channels=(1, 2),
        sensors=("hSpO2", "hSpO2"),
        data={1: raw_1, 2: raw_2},
        resolution=np.asarray([16, 16]),
    )

    recording = read_biosignalsplux_hdf5(file_path)

    np.testing.assert_allclose(recording.get_signal("fnirs_red").samples, expected_fnirs(raw_1))
    np.testing.assert_allclose(recording.get_signal("fnirs_infrared").samples, expected_fnirs(raw_2))
    assert recording.get_signal("fnirs_red").units == "uA"
    assert recording.get_signal("fnirs_infrared").units == "uA"


def test_fnirs_column_order_controls_red_infrared_assignment(tmp_path):
    file_path = tmp_path / "fnirs_order.h5"
    raw_1 = np.asarray([100, 200, 300])
    raw_2 = np.asarray([400, 500, 600])
    write_biosignalsplux_file(
        file_path,
        channels=(1, 2),
        sensors=("hSpO2", "hSpO2"),
        data={1: raw_1, 2: raw_2},
    )

    recording = read_biosignalsplux_hdf5(file_path, fnirs_column_order=(1, 0))

    np.testing.assert_allclose(recording.get_signal("fnirs_red").samples, expected_fnirs(raw_2))
    np.testing.assert_allclose(recording.get_signal("fnirs_infrared").samples, expected_fnirs(raw_1))


def test_unsupported_channels_are_skipped(tmp_path):
    file_path = tmp_path / "unsupported.h5"
    data = {1: np.arange(3), 2: np.asarray([0, 32768, 65535])}
    write_biosignalsplux_file(
        file_path,
        channels=(1, 2),
        sensors=("UNKNOWN", "ECG"),
        data=data,
    )

    recording = read_biosignalsplux_hdf5(file_path)

    assert recording.signal_names() == ["ecg"]


def test_duplicate_non_fnirs_signal_labels_are_skipped_after_first(tmp_path):
    file_path = tmp_path / "duplicate.h5"
    data = {1: np.asarray([0, 32768]), 2: np.asarray([65535, 32768])}
    write_biosignalsplux_file(
        file_path,
        channels=(1, 2),
        sensors=("ECG", "ECG"),
        data=data,
    )

    recording = read_biosignalsplux_hdf5(file_path)

    assert recording.signal_names() == ["ecg"]
    np.testing.assert_allclose(recording.get_signal("ecg").samples, expected_ecg(data[1]))


def test_recording_metadata_contains_expected_fields(tmp_path):
    file_path = tmp_path / "metadata.h5"
    write_biosignalsplux_file(file_path)

    metadata = read_biosignalsplux_hdf5(file_path).metadata

    assert set(metadata) == {"source", "source_format", "file_name", "dataset"}
    assert metadata["source"] == "read_biosignalsplux_hdf5"
    assert metadata["source_format"] == "BioSignalsPlux HDF5"
    assert metadata["file_name"] == "metadata.h5"
    assert metadata["dataset"] == "/recording_1"


def test_signal_metadata_contains_expected_fields(tmp_path):
    file_path = tmp_path / "signal_metadata.h5"
    write_biosignalsplux_file(file_path)

    metadata = read_biosignalsplux_hdf5(file_path).get_signal("ecg").metadata

    assert set(metadata) == {
        "source_format",
        "dataset",
        "channel_id",
        "channel_path",
        "sensor_label",
        "adc_resolution_bits",
    }
    assert metadata["source_format"] == "BioSignalsPlux HDF5"
    assert metadata["dataset"] == "/recording_1"
    assert metadata["channel_id"] == 1
    assert metadata["channel_path"] == "/recording_1/raw/channel_1"
    assert metadata["sensor_label"] == "ECG"
    assert metadata["adc_resolution_bits"] == 16


def test_time_vector_length_matches_samples(tmp_path):
    file_path = tmp_path / "time.h5"
    write_biosignalsplux_file(file_path, sampling_frequency=100)

    signal = read_biosignalsplux_hdf5(file_path).get_signal("ecg")

    assert signal.time.shape == signal.samples.shape
    np.testing.assert_allclose(signal.time, np.arange(signal.n_samples) / 100)


def test_sensor_label_is_trimmed_and_nulls_removed(tmp_path):
    file_path = tmp_path / "label.h5"
    with h5py.File(file_path, "w") as h5:
        group = h5.create_group("recording_1")
        group.attrs["channels"] = np.asarray([1])
        group.attrs["resolution"] = 16
        group.attrs["sampling rate"] = 1000
        raw = group.create_group("raw")
        dataset = raw.create_dataset("channel_1", data=np.asarray([0, 32768, 65535]))
        dataset.attrs["sensor"] = np.frombuffer(b" ECG\x00 ", dtype=np.uint8)

    recording = read_biosignalsplux_hdf5(file_path)

    assert recording.has_signal("ecg")


def test_missing_file_raises_file_not_found_error(tmp_path):
    with pytest.raises(FileNotFoundError):
        read_biosignalsplux_hdf5(tmp_path / "missing.h5")


def test_no_top_level_group_raises_value_error(tmp_path):
    file_path = tmp_path / "no_group.h5"
    with h5py.File(file_path, "w") as h5:
        h5.create_dataset("dataset", data=np.arange(3))

    with pytest.raises(ValueError, match="top-level recording group"):
        read_biosignalsplux_hdf5(file_path)


def test_missing_channels_attribute_raises_key_error(tmp_path):
    file_path = tmp_path / "missing_channels.h5"
    with h5py.File(file_path, "w") as h5:
        group = h5.create_group("recording_1")
        group.attrs["resolution"] = 16
        group.attrs["sampling rate"] = 1000

    with pytest.raises(KeyError, match="channels"):
        read_biosignalsplux_hdf5(file_path)


def test_missing_resolution_attribute_raises_key_error(tmp_path):
    file_path = tmp_path / "missing_resolution.h5"
    with h5py.File(file_path, "w") as h5:
        group = h5.create_group("recording_1")
        group.attrs["channels"] = np.asarray([1])
        group.attrs["sampling rate"] = 1000

    with pytest.raises(KeyError, match="resolution"):
        read_biosignalsplux_hdf5(file_path)


def test_invalid_sampling_frequency_raises_value_error(tmp_path):
    file_path = tmp_path / "bad_fs.h5"
    write_biosignalsplux_file(file_path, sampling_frequency=0)

    with pytest.raises(ValueError, match="sampling frequency"):
        read_biosignalsplux_hdf5(file_path)


def test_resolution_channel_length_mismatch_raises_value_error(tmp_path):
    file_path = tmp_path / "mismatch.h5"
    write_biosignalsplux_file(
        file_path,
        channels=(1, 2, 3),
        sensors=("ECG", "EDA", "TEMP"),
        data={1: np.arange(3), 2: np.arange(3), 3: np.arange(3) + 20000},
        resolution=np.asarray([16, 16]),
    )

    with pytest.raises(ValueError, match="resolution"):
        read_biosignalsplux_hdf5(file_path)


def test_no_supported_channels_raises_value_error(tmp_path):
    file_path = tmp_path / "none.h5"
    write_biosignalsplux_file(file_path, sensors=("UNKNOWN",))

    with pytest.raises(ValueError, match="No supported|no supported"):
        read_biosignalsplux_hdf5(file_path)


@pytest.mark.parametrize("order", [(0,), (0, 1, 2), (-1, 0)])
def test_invalid_fnirs_column_order_value_raises_value_error(tmp_path, order):
    file_path = tmp_path / "bad_order_value.h5"
    write_biosignalsplux_file(file_path)

    with pytest.raises(ValueError, match="fnirs_column_order"):
        read_biosignalsplux_hdf5(file_path, fnirs_column_order=order)


@pytest.mark.parametrize("order", ["01", (0, "1"), (True, 1)])
def test_invalid_fnirs_column_order_type_raises_type_error(tmp_path, order):
    file_path = tmp_path / "bad_order_type.h5"
    write_biosignalsplux_file(file_path)

    with pytest.raises(TypeError, match="fnirs_column_order"):
        read_biosignalsplux_hdf5(file_path, fnirs_column_order=order)


def test_missing_requested_fnirs_channels_raises_value_error(tmp_path):
    file_path = tmp_path / "missing_fnirs.h5"
    write_biosignalsplux_file(file_path, sensors=("hSpO2",))

    with pytest.raises(ValueError, match="hSpO2"):
        read_biosignalsplux_hdf5(file_path, fnirs_column_order=(0, 1))


def test_selected_fnirs_length_mismatch_raises_value_error(tmp_path):
    file_path = tmp_path / "fnirs_length_mismatch.h5"
    write_biosignalsplux_file(
        file_path,
        channels=(1, 2),
        sensors=("hSpO2", "hSpO2"),
        data={1: np.asarray([1, 2, 3]), 2: np.asarray([1, 2])},
    )

    with pytest.raises(ValueError, match="equal lengths"):
        read_biosignalsplux_hdf5(file_path)


def test_temp_nominal_ntc_resistance_reads_25_degrees_celsius(tmp_path):
    # Mid-scale ADC -> 1.5 V across the divider -> 10 kOhm, the nominal
    # resistance of the 10 kOhm NTC, whose Steinhart-Hart coefficients give
    # 298.15 K. Kelvin to Celsius must subtract 273.15 (the SI definition).
    file_path = tmp_path / "temp_nominal.h5"
    write_biosignalsplux_file(file_path, sensors=("TEMP",), data={1: np.asarray([32768, 32768, 32768])})

    signal = read_biosignalsplux_hdf5(file_path).get_signal("temp")

    np.testing.assert_allclose(signal.samples, 25.0, atol=1e-3)

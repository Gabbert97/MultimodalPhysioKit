from __future__ import annotations

from pathlib import Path
from typing import Any

import h5py
import numpy as np

from multimodalphysiokit.core import Recording, Signal


_SOURCE_FORMAT = "BioSignalsPlux HDF5"
_UNITS = {
    "ecg": "mV",
    "eda": "uS",
    "temp": "degC",
    "rip": "percent full scale",
    "respiration": "V",
    "fnirs_red": "uA",
    "fnirs_infrared": "uA",
}


def read_biosignalsplux_hdf5(
    file_path: str | Path,
    *,
    fnirs_column_order: tuple[int, int] | list[int] = (0, 1),
    name: str = "biosignalsplux",
) -> Recording:
    """Read the supported BioSignalsPlux/OpenSignals HDF5 layout."""

    path = Path(file_path)
    if not path.is_file():
        raise FileNotFoundError(f"HDF5 file not found: {path}")
    order = _validate_fnirs_column_order(fnirs_column_order)

    with h5py.File(path, "r") as hdf5_file:
        dataset_path, group = _first_top_level_group(hdf5_file)
        channel_ids = _read_required_integer_attribute(group, "channels")
        resolution_bits = _read_resolution_bits(group, channel_ids.size)
        sampling_frequency = _read_sampling_frequency(group)

        recording = Recording(
            name=name,
            metadata={
                "source": "read_biosignalsplux_hdf5",
                "source_format": _SOURCE_FORMAT,
                "file_name": path.name,
                "dataset": dataset_path,
            },
        )
        fnirs_raw: list[np.ndarray] = []
        fnirs_channel_ids: list[int] = []
        fnirs_resolution_bits: list[float] = []
        fnirs_channel_paths: list[str] = []
        fnirs_sensor_labels: list[str] = []

        for channel_id, adc_resolution_bits in zip(channel_ids, resolution_bits, strict=True):
            channel_path = f"{dataset_path}/raw/channel_{int(channel_id)}"
            if channel_path not in hdf5_file or not isinstance(hdf5_file[channel_path], h5py.Dataset):
                raise KeyError(f"missing BioSignalsPlux raw channel dataset: {channel_path}")
            dataset = hdf5_file[channel_path]
            sensor_label = _sensor_label(dataset.attrs["sensor"])
            canonical_sensor = _canonical_sensor(sensor_label)
            if canonical_sensor is None:
                continue
            raw = _read_raw_channel(dataset, channel_path)

            if canonical_sensor == "hSpO2":
                fnirs_raw.append(raw)
                fnirs_channel_ids.append(int(channel_id))
                fnirs_resolution_bits.append(float(adc_resolution_bits))
                fnirs_channel_paths.append(channel_path)
                fnirs_sensor_labels.append(sensor_label)
                continue

            signal_name = _signal_name(canonical_sensor)
            if signal_name is None or recording.has_signal(signal_name):
                continue
            samples = _convert_biosignalsplux_channel(
                raw,
                canonical_sensor,
                float(adc_resolution_bits),
            )
            recording.add_signal(
                signal_name,
                _make_signal(
                    samples,
                    sampling_frequency,
                    signal_name,
                    dataset_path,
                    int(channel_id),
                    channel_path,
                    sensor_label,
                    float(adc_resolution_bits),
                ),
            )

        if fnirs_raw:
            _add_fnirs_signals(
                recording,
                fnirs_raw,
                fnirs_channel_ids,
                fnirs_resolution_bits,
                fnirs_channel_paths,
                fnirs_sensor_labels,
                order,
                dataset_path,
                sampling_frequency,
            )

    if len(recording) == 0:
        raise ValueError("no supported BioSignalsPlux channels were found")
    return recording


def _first_top_level_group(hdf5_file: h5py.File) -> tuple[str, h5py.Group]:
    for key in hdf5_file.keys():
        item = hdf5_file[key]
        if isinstance(item, h5py.Group):
            return f"/{key}", item
    raise ValueError("the HDF5 file has no top-level recording group")


def _read_required_integer_attribute(group: h5py.Group, name: str) -> np.ndarray:
    if name not in group.attrs:
        raise KeyError(f"missing BioSignalsPlux group attribute: {name}")
    values = np.asarray(group.attrs[name]).reshape(-1)
    if values.size == 0 or not np.issubdtype(values.dtype, np.number):
        raise ValueError(f"BioSignalsPlux group attribute '{name}' must be numeric")
    if not np.all(np.isfinite(values)):
        raise ValueError(f"BioSignalsPlux group attribute '{name}' must contain finite values")
    return values.astype(int)


def _read_resolution_bits(group: h5py.Group, n_channels: int) -> np.ndarray:
    if "resolution" not in group.attrs:
        raise KeyError("missing BioSignalsPlux group attribute: resolution")
    values = np.asarray(group.attrs["resolution"], dtype=float).reshape(-1)
    if values.size == 1 and n_channels > 1:
        values = np.repeat(values, n_channels)
    if values.size != n_channels:
        raise ValueError("channel and ADC resolution metadata have different lengths")
    if not np.all(np.isfinite(values)):
        raise ValueError("ADC resolution metadata must contain finite values")
    return values


def _read_sampling_frequency(group: h5py.Group) -> float:
    if "sampling rate" in group.attrs:
        value = group.attrs["sampling rate"]
    elif "fs" in group.attrs:
        value = group.attrs["fs"]
    else:
        raise KeyError("missing BioSignalsPlux sampling frequency attribute")
    array = np.asarray(value, dtype=float).reshape(-1)
    if array.size != 1:
        raise ValueError("sampling frequency must be a positive scalar")
    sampling_frequency = float(array[0])
    if not np.isfinite(sampling_frequency) or sampling_frequency <= 0:
        raise ValueError("sampling frequency must be a positive scalar")
    return sampling_frequency


def _sensor_label(value: Any) -> str:
    array = np.asarray(value)
    if array.dtype.kind in {"S", "U", "O"}:
        if array.shape:
            value = array.reshape(-1)[0]
        if isinstance(value, bytes):
            label = value.decode(errors="ignore")
        else:
            label = str(value)
    elif np.issubdtype(array.dtype, np.number):
        label = bytes(array.astype(np.uint8).reshape(-1)).decode(errors="ignore")
    else:
        label = str(value)
    return label.replace("\x00", "").strip()


def _canonical_sensor(sensor_label: str) -> str | None:
    if sensor_label.lower() == "ecg":
        return "ECG"
    if sensor_label.lower() == "eda":
        return "EDA"
    if sensor_label.lower() == "temp":
        return "TEMP"
    if sensor_label.lower() == "rip":
        return "RIP"
    if sensor_label.lower() == "respiration":
        return "RESPIRATION"
    if sensor_label.lower() == "hspo2":
        return "hSpO2"
    return None


def _signal_name(canonical_sensor: str) -> str | None:
    return {
        "ECG": "ecg",
        "EDA": "eda",
        "TEMP": "temp",
        "RIP": "rip",
        "RESPIRATION": "respiration",
    }.get(canonical_sensor)


def _read_raw_channel(dataset: h5py.Dataset, channel_path: str) -> np.ndarray:
    array = np.asarray(dataset[()])
    if array.ndim != 1:
        array = array.reshape(-1)
    if not np.issubdtype(array.dtype, np.number):
        raise ValueError(f"BioSignalsPlux channel dataset must be numeric: {channel_path}")
    array = array.astype(float)
    if not np.all(np.isfinite(array)):
        raise ValueError(f"BioSignalsPlux channel dataset must contain finite values: {channel_path}")
    return array


def _convert_biosignalsplux_channel(
    raw: np.ndarray,
    canonical_sensor: str,
    adc_resolution_bits: float,
) -> np.ndarray:
    vcc = 3.0
    denominator = 2.0**adc_resolution_bits
    if canonical_sensor == "TEMP":
        a0 = 1.12764514e-3
        a1 = 2.34282709e-4
        a2 = 8.77303013e-8
        ntc_v = raw * vcc / denominator
        ntc_ohm = (ntc_v * 1e4) / (vcc - ntc_v)
        tmp_k = 1.0 / (a0 + a1 * np.log(ntc_ohm) + a2 * np.log(ntc_ohm) ** 3)
        return tmp_k - 273.5
    if canonical_sensor == "RIP":
        return (raw / denominator - 0.5) * 100
    if canonical_sensor == "EDA":
        return raw * vcc / (0.12 * denominator)
    if canonical_sensor == "ECG":
        ecg_v = (raw / denominator - 0.5) * vcc / 1019
        return ecg_v * 1e3
    if canonical_sensor == "RESPIRATION":
        return (raw / (denominator - 1) - 0.5) * vcc
    raise ValueError(f"unsupported BioSignalsPlux sensor for conversion: {canonical_sensor}")


def _convert_biosignalsplux_fnirs(raw: np.ndarray, adc_resolution_bits: np.ndarray) -> np.ndarray:
    converted = np.zeros_like(raw, dtype=float)
    for index in range(raw.shape[1]):
        converted[:, index] = 0.15 * raw[:, index] / (2.0 ** adc_resolution_bits[index])
    return converted


def _add_fnirs_signals(
    recording: Recording,
    fnirs_raw: list[np.ndarray],
    channel_ids: list[int],
    resolution_bits: list[float],
    channel_paths: list[str],
    sensor_labels: list[str],
    order: tuple[int, int],
    dataset_path: str,
    sampling_frequency: float,
) -> None:
    if max(order) >= len(fnirs_raw):
        raise ValueError("fnirs_column_order requires more hSpO2 channels than the file contains")
    selected_raw = [fnirs_raw[index] for index in order]
    sample_counts = {samples.size for samples in selected_raw}
    if len(sample_counts) != 1:
        raise ValueError("selected hSpO2 channels must have equal lengths")
    raw = np.column_stack(selected_raw)
    selected_resolution_bits = np.asarray([resolution_bits[index] for index in order], dtype=float)
    converted = _convert_biosignalsplux_fnirs(raw, selected_resolution_bits)
    for output_column, signal_name in enumerate(("fnirs_red", "fnirs_infrared")):
        selected_index = order[output_column]
        recording.add_signal(
            signal_name,
            _make_signal(
                converted[:, output_column],
                sampling_frequency,
                signal_name,
                dataset_path,
                channel_ids[selected_index],
                channel_paths[selected_index],
                sensor_labels[selected_index],
                resolution_bits[selected_index],
            ),
        )


def _make_signal(
    samples: np.ndarray,
    sampling_frequency: float,
    signal_name: str,
    dataset_path: str,
    channel_id: int,
    channel_path: str,
    sensor_label: str,
    adc_resolution_bits: float,
) -> Signal:
    return Signal(
        samples=samples,
        time=np.arange(samples.size, dtype=float) / sampling_frequency,
        sampling_frequency=sampling_frequency,
        label=signal_name,
        units=_UNITS[signal_name],
        metadata={
            "source_format": _SOURCE_FORMAT,
            "dataset": dataset_path,
            "channel_id": int(channel_id),
            "channel_path": channel_path,
            "sensor_label": sensor_label,
            "adc_resolution_bits": float(adc_resolution_bits),
        },
    )


def _validate_fnirs_column_order(value: tuple[int, int] | list[int]) -> tuple[int, int]:
    if not isinstance(value, (tuple, list)):
        raise TypeError("fnirs_column_order must be a tuple or list of two zero-based integer indices")
    if len(value) != 2:
        raise ValueError("fnirs_column_order must contain exactly two indices")
    if not all(isinstance(index, int) and not isinstance(index, bool) for index in value):
        raise TypeError("fnirs_column_order indices must be integers")
    if any(index < 0 for index in value):
        raise ValueError("fnirs_column_order indices must be zero-based nonnegative integers")
    return int(value[0]), int(value[1])

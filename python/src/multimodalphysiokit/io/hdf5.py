from __future__ import annotations

from pathlib import Path
from typing import Any

import h5py
import numpy as np

from multimodalphysiokit.core import Recording

from .arrays import recording_from_arrays


def read_generic_hdf5(
    file_path: str | Path,
    signal_paths: dict[str, str],
    *,
    sampling_frequency: float | None = None,
    sampling_frequency_path: str | None = None,
    time_path: str | None = None,
    units: dict[str, str] | None = None,
    labels: dict[str, str] | None = None,
    name: str | None = None,
    metadata: dict[str, Any] | None = None,
) -> Recording:
    """Read explicitly selected one-dimensional numeric datasets from HDF5."""

    path = Path(file_path)
    if not path.is_file():
        raise FileNotFoundError(f"HDF5 file not found: {path}")
    _validate_signal_paths(signal_paths)
    _validate_timing_selection(
        sampling_frequency=sampling_frequency,
        sampling_frequency_path=sampling_frequency_path,
        time_path=time_path,
    )

    recording_metadata = {} if metadata is None else dict(metadata)
    recording_metadata.update(
        {
            "source": "read_generic_hdf5",
            "file_name": path.name,
        }
    )

    with h5py.File(path, "r") as hdf5_file:
        signals = {
            signal_name: _read_numeric_1d_dataset(hdf5_file, dataset_path, f"signal '{signal_name}'")
            for signal_name, dataset_path in signal_paths.items()
        }
        if sampling_frequency_path is not None:
            sampling_frequency = _read_sampling_frequency(hdf5_file, sampling_frequency_path)
            time = None
        elif time_path is not None:
            time = _read_numeric_1d_dataset(hdf5_file, time_path, "time")
        else:
            time = None

    return recording_from_arrays(
        signals,
        sampling_frequency=sampling_frequency,
        time=time,
        units=units,
        labels=labels,
        name=path.stem if name is None else name,
        metadata=recording_metadata,
    )


def _validate_signal_paths(signal_paths: dict[str, str]) -> None:
    if not isinstance(signal_paths, dict):
        raise TypeError("signal_paths must be a dictionary")
    if not signal_paths:
        raise ValueError("signal_paths must be non-empty")
    for signal_name, dataset_path in signal_paths.items():
        if not isinstance(signal_name, str) or not signal_name:
            raise ValueError("signal names must be non-empty strings")
        if not isinstance(dataset_path, str) or not dataset_path:
            raise ValueError("HDF5 dataset paths must be non-empty strings")


def _validate_timing_selection(
    *,
    sampling_frequency: float | None,
    sampling_frequency_path: str | None,
    time_path: str | None,
) -> None:
    selected = [
        sampling_frequency is not None,
        sampling_frequency_path is not None,
        time_path is not None,
    ]
    if sum(selected) != 1:
        raise ValueError("exactly one timing source must be provided")
    if sampling_frequency_path is not None and (
        not isinstance(sampling_frequency_path, str) or not sampling_frequency_path
    ):
        raise ValueError("sampling_frequency_path must be a non-empty string")
    if time_path is not None and (not isinstance(time_path, str) or not time_path):
        raise ValueError("time_path must be a non-empty string")


def _read_numeric_1d_dataset(hdf5_file: h5py.File, dataset_path: str, description: str) -> np.ndarray:
    if dataset_path not in hdf5_file:
        raise KeyError(f"missing HDF5 dataset for {description}: {dataset_path}")
    dataset = hdf5_file[dataset_path]
    if not isinstance(dataset, h5py.Dataset):
        raise KeyError(f"HDF5 path is not a dataset for {description}: {dataset_path}")
    array = np.asarray(dataset[()])
    if array.ndim != 1:
        raise ValueError(f"HDF5 dataset for {description} must be one-dimensional")
    if not np.issubdtype(array.dtype, np.number):
        raise ValueError(f"HDF5 dataset for {description} must be numeric")
    array = array.astype(float)
    if not np.all(np.isfinite(array)):
        raise ValueError(f"HDF5 dataset for {description} must contain finite values")
    return array


def _read_sampling_frequency(hdf5_file: h5py.File, dataset_path: str) -> float:
    if dataset_path not in hdf5_file:
        raise KeyError(f"missing sampling frequency dataset: {dataset_path}")
    dataset = hdf5_file[dataset_path]
    if not isinstance(dataset, h5py.Dataset):
        raise KeyError(f"sampling frequency path is not a dataset: {dataset_path}")
    array = np.asarray(dataset[()])
    if array.size != 1:
        raise ValueError("sampling_frequency_path must point to a scalar dataset")
    if not np.issubdtype(array.dtype, np.number):
        raise ValueError("sampling_frequency_path dataset must be numeric")
    value = float(array.reshape(-1)[0])
    if not np.isfinite(value) or value <= 0:
        raise ValueError("sampling_frequency_path dataset must contain a positive finite value")
    return value

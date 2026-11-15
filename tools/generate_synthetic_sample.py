#!/usr/bin/env python3
"""Generate the deterministic public MultimodalPhysioKit sample fixture."""

from __future__ import annotations

import argparse
import sys
import tempfile
from pathlib import Path

import h5py
import numpy as np


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
PYTHON_SOURCE = REPOSITORY_ROOT / "python" / "src"
OUTPUT_PATH = REPOSITORY_ROOT / "examples" / "sample_data" / "synthetic_multimodal_180s.hdf5"

SAMPLING_FREQUENCY = 100
DURATION_SECONDS = 180
N_SAMPLES = SAMPLING_FREQUENCY * DURATION_SECONDS
ADC_BITS = 16
ADC_LEVELS = 2**ADC_BITS
RANDOM_SEED = 20260808
GENERATOR_ID = "MultimodalPhysioKit synthetic sample generator v2"
LICENSE_STATEMENT = (
    "CC0-1.0: mathematically generated data; no human measurements. "
    "To the extent possible under law, the author waives all copyright and related rights."
)
CHANNELS = (
    (1, "ECG", "ecg"),
    (2, "EDA", "eda"),
    (3, "RIP", "rip"),
    (4, "TEMP", "temp"),
    (5, "BVP", "ppg"),
    (9, "hSpO2", "fnirs_red"),
    (10, "hSpO2", "fnirs_infrared"),
)
EXPECTED_UNITS = {
    "ecg": "mV",
    "eda": "uS",
    "rip": "percent full scale",
    "temp": "degC",
    "ppg": "percent full scale",
    "fnirs_red": "uA",
    "fnirs_infrared": "uA",
}


def _time() -> np.ndarray:
    return np.arange(N_SAMPLES, dtype=float) / SAMPLING_FREQUENCY


def _beat_times() -> list[float]:
    """Heartbeat times shared by the ECG and the PPG (72 bpm with slow variability)."""
    beat_times = []
    beat_time = 0.45
    while beat_time < DURATION_SECONDS + 0.5:
        beat_times.append(beat_time)
        heart_rate = 72.0 + 2.0 * np.sin(2 * np.pi * 0.10 * beat_time)
        heart_rate += 1.0 * np.sin(2 * np.pi * 0.025 * beat_time + 0.7)
        beat_time += 60.0 / heart_rate
    return beat_times


def _synthetic_ecg(time: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    signal = 0.025 * np.sin(2 * np.pi * 0.20 * time)
    signal += 0.008 * np.sin(2 * np.pi * 0.05 * time + 0.4)

    for beat_time in _beat_times():
        for offset, amplitude, width in (
            (-0.18, 0.08, 0.035),
            (-0.035, -0.12, 0.012),
            (0.0, 1.00, 0.014),
            (0.030, -0.22, 0.014),
            (0.26, 0.24, 0.070),
        ):
            signal += amplitude * np.exp(-0.5 * ((time - (beat_time + offset)) / width) ** 2)

    signal += rng.normal(0.0, 0.006, time.size)
    return signal


def _synthetic_ppg(time: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    """Fingertip-like PPG in percent of full scale, one pulse per ECG beat.

    Each pulse starts 0.20 s after the R-peak (pulse arrival time) and has a
    systolic wave plus a smaller dicrotic wave; breathing modulates the baseline.
    """

    signal = -10.0 + 1.5 * np.sin(2 * np.pi * 0.25 * time) + 0.8 * np.sin(2 * np.pi * time / 90.0)
    for beat_time in _beat_times():
        onset = beat_time + 0.20
        for offset, amplitude, width in ((0.13, 12.0, 0.055), (0.40, 4.0, 0.08)):
            signal += amplitude * np.exp(-0.5 * ((time - (onset + offset)) / width) ** 2)
    signal += rng.normal(0.0, 0.05, time.size)
    return signal


def _scr_kernel(time_since_event: np.ndarray) -> np.ndarray:
    response = np.zeros_like(time_since_event)
    active = time_since_event >= 0
    elapsed = time_since_event[active]
    response[active] = np.exp(-elapsed / 3.0) - np.exp(-elapsed / 0.7)
    return response


def _synthetic_eda(time: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    signal = 7.0 + 0.0025 * time + 0.10 * np.sin(2 * np.pi * time / 120.0)
    event_times = (12, 28, 46, 68, 83, 105, 128, 146, 166)
    event_amplitudes = (0.35, 0.55, 0.42, 0.60, 0.38, 0.50, 0.45, 0.62, 0.40)
    for event_time, amplitude in zip(event_times, event_amplitudes, strict=True):
        signal += amplitude * _scr_kernel(time - event_time)
    signal += rng.normal(0.0, 0.004, time.size)
    return signal


def _synthetic_respiration(time: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    instantaneous_frequency = 0.25 + 0.015 * np.sin(2 * np.pi * time / 75.0)
    phase = 2 * np.pi * np.cumsum(instantaneous_frequency) / SAMPLING_FREQUENCY
    amplitude = 22.0 * (1.0 + 0.08 * np.sin(2 * np.pi * time / 55.0 + 0.3))
    signal = amplitude * (np.sin(phase) + 0.12 * np.sin(2 * phase - 0.4))
    signal += 1.5 * np.sin(2 * np.pi * time / 90.0)
    signal += rng.normal(0.0, 0.12, time.size)
    return signal


def _synthetic_temperature(time: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    signal = 31.4 + 0.006 * time + 0.16 * np.sin(2 * np.pi * time / 140.0 + 0.2)
    signal += rng.normal(0.0, 0.003, time.size)
    return signal


def _smooth_step(time: np.ndarray, start: float, stop: float, edge_seconds: float = 3.0) -> np.ndarray:
    rise = 0.5 * (1.0 + np.tanh((time - start) / edge_seconds))
    fall = 0.5 * (1.0 + np.tanh((stop - time) / edge_seconds))
    return rise * fall


def _synthetic_fnirs(
    time: np.ndarray, rng: np.random.Generator
) -> tuple[np.ndarray, np.ndarray]:
    physiological_red = 0.00012 * np.sin(2 * np.pi * 0.10 * time)
    physiological_red += 0.00005 * np.sin(2 * np.pi * 0.25 * time + 0.4)
    physiological_infrared = 0.00016 * np.sin(2 * np.pi * 0.10 * time + 0.6)
    physiological_infrared += 0.00004 * np.sin(2 * np.pi * 0.25 * time + 1.1)

    task_1 = _smooth_step(time, 63.0, 117.0)
    task_2 = _smooth_step(time, 123.0, 177.0)
    red = 0.0150 + physiological_red - 0.00032 * task_1 + 0.00020 * task_2
    infrared = 0.0280 + physiological_infrared - 0.00012 * task_1 - 0.00038 * task_2
    red += rng.normal(0.0, 0.000008, time.size)
    infrared += rng.normal(0.0, 0.000010, time.size)
    return red, infrared


def generate_physical_signals() -> dict[str, np.ndarray]:
    rng = np.random.default_rng(RANDOM_SEED)
    time = _time()
    fnirs_red, fnirs_infrared = _synthetic_fnirs(time, rng)
    return {
        "ecg": _synthetic_ecg(time, rng),
        "eda": _synthetic_eda(time, rng),
        "rip": _synthetic_respiration(time, rng),
        "temp": _synthetic_temperature(time, rng),
        "ppg": _synthetic_ppg(time, rng),
        "fnirs_red": fnirs_red,
        "fnirs_infrared": fnirs_infrared,
    }


def _temperature_from_raw(raw: np.ndarray) -> np.ndarray:
    voltage = raw * 3.0 / ADC_LEVELS
    resistance = voltage * 1.0e4 / (3.0 - voltage)
    log_resistance = np.log(resistance)
    kelvin = 1.0 / (
        1.12764514e-3
        + 2.34282709e-4 * log_resistance
        + 8.77303013e-8 * log_resistance**3
    )
    return kelvin - 273.15


def _temperature_to_raw(temperature: np.ndarray) -> np.ndarray:
    low = np.ones(temperature.size, dtype=float)
    high = np.full(temperature.size, ADC_LEVELS - 1.0)
    for _ in range(32):
        middle = (low + high) / 2.0
        decoded = _temperature_from_raw(middle)
        lower_raw_needed = decoded < temperature
        high[lower_raw_needed] = middle[lower_raw_needed]
        low[~lower_raw_needed] = middle[~lower_raw_needed]
    return np.rint((low + high) / 2.0)


def convert_to_raw(signals: dict[str, np.ndarray]) -> dict[str, np.ndarray]:
    candidates = {
        "ecg": ((signals["ecg"] / 1000.0) * 1019.0 / 3.0 + 0.5) * ADC_LEVELS,
        "eda": signals["eda"] * (0.12 * ADC_LEVELS) / 3.0,
        "rip": (signals["rip"] / 100.0 + 0.5) * ADC_LEVELS,
        "temp": _temperature_to_raw(signals["temp"]),
        "ppg": (signals["ppg"] / 100.0 + 0.5) * ADC_LEVELS,
        "fnirs_red": signals["fnirs_red"] * ADC_LEVELS / 0.15,
        "fnirs_infrared": signals["fnirs_infrared"] * ADC_LEVELS / 0.15,
    }
    raw_signals: dict[str, np.ndarray] = {}
    for name, candidate in candidates.items():
        if not np.all(np.isfinite(candidate)):
            raise ValueError(f"{name} inverse conversion produced non-finite ADC values")
        rounded = np.rint(candidate)
        if np.any(rounded < 0) or np.any(rounded > ADC_LEVELS - 1):
            minimum = float(np.min(candidate))
            maximum = float(np.max(candidate))
            raise ValueError(f"{name} cannot be represented in 16 bits: range {minimum}..{maximum}")
        raw_signals[name] = rounded.astype(np.uint16)
    return raw_signals


def write_fixture(path: Path, raw_signals: dict[str, np.ndarray]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(path, "w") as hdf5_file:
        recording = hdf5_file.create_group("synthetic_recording")
        recording.attrs["channels"] = np.asarray([item[0] for item in CHANNELS], dtype=np.uint16)
        recording.attrs["resolution"] = np.full(len(CHANNELS), ADC_BITS, dtype=np.uint16)
        recording.attrs["sampling rate"] = np.uint16(SAMPLING_FREQUENCY)
        recording.attrs["duration seconds"] = np.uint16(DURATION_SECONDS)
        recording.attrs["nsamples"] = np.uint32(N_SAMPLES)
        recording.attrs["synthetic"] = np.bool_(True)
        recording.attrs["generator"] = GENERATOR_ID
        recording.attrs["license"] = LICENSE_STATEMENT
        raw_group = recording.create_group("raw")
        for channel_id, sensor, signal_name in CHANNELS:
            dataset = raw_group.create_dataset(
                f"channel_{channel_id}",
                data=raw_signals[signal_name].reshape(1, -1),
                dtype=np.uint16,
                compression="gzip",
                compression_opts=9,
                shuffle=True,
            )
            dataset.attrs["sensor"] = sensor


def _import_reader():
    sys.path.insert(0, str(PYTHON_SOURCE))
    from multimodalphysiokit.io import read_biosignalsplux_hdf5

    return read_biosignalsplux_hdf5


def validate_fixture(path: Path, physical_signals: dict[str, np.ndarray]) -> None:
    recording = _import_reader()(path)
    expected_names = ["ecg", "eda", "rip", "temp", "ppg", "fnirs_red", "fnirs_infrared"]
    if recording.signal_names() != expected_names:
        raise AssertionError(f"unexpected canonical signals: {recording.signal_names()}")
    forbidden_metadata_keys = {"source_path", "file_path", "path"}
    if forbidden_metadata_keys.intersection(recording.metadata):
        raise AssertionError("recording metadata leaks a filesystem path")

    tolerances = {
        "ecg": 2.3e-5,
        "eda": 2.0e-4,
        "rip": 8.0e-4,
        "temp": 1.0e-3,
        "ppg": 8.0e-4,
        "fnirs_red": 1.2e-6,
        "fnirs_infrared": 1.2e-6,
    }
    for name in expected_names:
        signal = recording.get_signal(name)
        if signal.sampling_frequency != SAMPLING_FREQUENCY:
            raise AssertionError(f"{name} has sampling frequency {signal.sampling_frequency}")
        if signal.units != EXPECTED_UNITS[name]:
            raise AssertionError(f"{name} has units {signal.units!r}")
        if signal.n_samples != N_SAMPLES:
            raise AssertionError(f"{name} has {signal.n_samples} samples")
        if forbidden_metadata_keys.intersection(signal.metadata):
            raise AssertionError(f"{name} metadata leaks a filesystem path")
        np.testing.assert_allclose(
            signal.samples,
            physical_signals[name],
            atol=tolerances[name],
            rtol=0.0,
            err_msg=f"{name} round-trip conversion exceeded ADC quantization tolerance",
        )


def _normalized_attribute(value: object) -> tuple[str, tuple[int, ...], bytes]:
    array = np.asarray(value)
    if array.dtype.kind in {"O", "S", "U"}:
        return "text", array.shape, str(value).encode("utf-8")
    return array.dtype.str, array.shape, array.tobytes()


def compare_logical_contents(expected_path: Path, actual_path: Path) -> None:
    with h5py.File(expected_path, "r") as expected, h5py.File(actual_path, "r") as actual:
        expected_paths: list[str] = []
        actual_paths: list[str] = []
        expected.visit(expected_paths.append)
        actual.visit(actual_paths.append)
        if expected_paths != actual_paths:
            raise AssertionError("tracked and regenerated HDF5 object paths differ")

        for object_path in ["/", *expected_paths]:
            expected_object = expected[object_path]
            actual_object = actual[object_path]
            if type(expected_object) is not type(actual_object):
                raise AssertionError(f"HDF5 object types differ at {object_path}")
            if set(expected_object.attrs) != set(actual_object.attrs):
                raise AssertionError(f"HDF5 attributes differ at {object_path}")
            for name in expected_object.attrs:
                if _normalized_attribute(expected_object.attrs[name]) != _normalized_attribute(
                    actual_object.attrs[name]
                ):
                    raise AssertionError(f"HDF5 attribute {name!r} differs at {object_path}")
            if isinstance(expected_object, h5py.Dataset):
                if expected_object.dtype != actual_object.dtype:
                    raise AssertionError(f"dataset dtype differs at {object_path}")
                if expected_object.shape != actual_object.shape:
                    raise AssertionError(f"dataset shape differs at {object_path}")
                if not np.array_equal(expected_object[...], actual_object[...]):
                    raise AssertionError(f"dataset values differ at {object_path}")


def generate(path: Path) -> None:
    physical_signals = generate_physical_signals()
    raw_signals = convert_to_raw(physical_signals)
    write_fixture(path, raw_signals)
    validate_fixture(path, physical_signals)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--check",
        action="store_true",
        help="regenerate temporarily and compare logical contents with the tracked fixture",
    )
    arguments = parser.parse_args()

    if arguments.check:
        if not OUTPUT_PATH.is_file():
            raise FileNotFoundError(f"tracked fixture is missing: {OUTPUT_PATH}")
        with tempfile.TemporaryDirectory(prefix="mpk-synthetic-sample-") as temporary_directory:
            regenerated_path = Path(temporary_directory) / OUTPUT_PATH.name
            generate(regenerated_path)
            compare_logical_contents(OUTPUT_PATH, regenerated_path)
        print(f"Synthetic fixture check passed: {OUTPUT_PATH.relative_to(REPOSITORY_ROOT)}")
    else:
        generate(OUTPUT_PATH)
        print(f"Generated and validated: {OUTPUT_PATH.relative_to(REPOSITORY_ROOT)}")


if __name__ == "__main__":
    main()

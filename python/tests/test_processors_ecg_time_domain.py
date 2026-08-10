import numpy as np
import pytest

from multimodalphysiokit import Signal
from multimodalphysiokit.processors.ecg import _extract_ecg_time_domain_features


FEATURE_NAMES = ["ecg_bpm_mean", "ecg_bpm_std", "ecg_bpm_median", "ecg_pnn50", "ecg_pnn50_max"]


def make_synthetic_ecg(duration_seconds=30, sampling_frequency=250):
    time = np.arange(0, duration_seconds, 1 / sampling_frequency)
    samples = 0.02 * np.sin(2 * np.pi * 0.5 * time)
    for beat_time in np.arange(0.8, duration_seconds, 1.0):
        samples += np.exp(-0.5 * ((time - beat_time) / 0.025) ** 2)
        samples += 0.15 * np.exp(-0.5 * ((time - (beat_time + 0.22)) / 0.06) ** 2)
    return Signal(samples, time, sampling_frequency, label="ecg", units="mV")


def test_helper_returns_single_signal_outputs():
    result = _extract_ecg_time_domain_features(make_synthetic_ecg())
    assert set(result) == {
        "features", "filtered_signal", "r_peak_indices", "r_peak_amplitudes",
        "accepted_r_peak_indices", "ibi", "ibi_times", "bpm",
    }
    assert list(result["features"]) == FEATURE_NAMES


def test_arrays_are_finite_and_scientifically_consistent():
    result = _extract_ecg_time_domain_features(make_synthetic_ecg())
    assert np.all(np.isfinite(result["filtered_signal"]))
    assert np.all(np.isfinite(result["r_peak_amplitudes"]))
    assert np.all(np.diff(result["r_peak_indices"]) > 0)
    assert np.all(np.diff(result["accepted_r_peak_indices"]) > 0)
    assert np.all(np.isfinite(result["ibi"])) and np.all(result["ibi"] > 0)
    assert np.all(np.isfinite(result["bpm"]))
    assert result["ibi"].shape == result["ibi_times"].shape == result["bpm"].shape


def test_time_domain_features_match_existing_complete_signal_values():
    result = _extract_ecg_time_domain_features(make_synthetic_ecg())
    assert result["features"] == {
        "ecg_bpm_mean": 60.0, "ecg_bpm_std": 0.0, "ecg_bpm_median": 60.0,
        "ecg_pnn50": 0.0, "ecg_pnn50_max": 0.0,
    }


def test_invalid_signal_and_flat_signal_raise_clear_errors():
    with pytest.raises(TypeError, match="Signal"):
        _extract_ecg_time_domain_features([1, 2, 3])
    time = np.arange(0, 10, 1 / 250)
    with pytest.raises(ValueError, match="no R-peaks|too few"):
        _extract_ecg_time_domain_features(Signal(np.zeros_like(time), time, 250))

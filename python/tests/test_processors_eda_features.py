from __future__ import annotations

import numpy as np
import pytest

import multimodalphysiokit.processors.eda as eda_module
from multimodalphysiokit.processors.eda import (
    _extract_eda_features,
    _find_max_min,
    _pair_scr_events,
    _scr_processing,
)


def test_extrema_detection_regression():
    maxima, minima = _find_max_min(np.asarray([0, 1, 0, -1, 0, 1, 0], dtype=float))
    np.testing.assert_array_equal(maxima, [1, 5])
    np.testing.assert_array_equal(minima, [3])


def test_scr_processing_is_deterministic():
    time = np.arange(0, 30, 0.01)
    scr = sum((np.exp(-0.5 * ((time - center) / 0.2) ** 2)
               for center in (3, 8, 13, 18, 23, 28)), np.zeros(time.size))
    first = _scr_processing(scr, time, 100.0)
    second = _scr_processing(scr, time, 100.0)
    np.testing.assert_allclose(first[0], second[0])
    np.testing.assert_allclose(first[1], second[1])
    assert first[2] == second[2]


def test_feature_formulas_regression():
    time = np.arange(0, 30, 0.01)
    scl = 1.5 + 0.001 * np.arange(time.size)
    scr = sum((0.4 * np.exp(-0.5 * ((time - center) / 0.2) ** 2)
               for center in (3, 8, 13, 18, 23, 28)), np.zeros(time.size))
    values = _extract_eda_features(scl, scr, time, 100.0)
    assert tuple(values) == (
        "scl_mean", "scl_std", "scl_slope", "scr_amplitude_mean",
        "scr_amplitude_std", "scr_rise_time_mean", "scr_rise_time_std", "scr_per_min",
    )
    assert values["scl_mean"] == pytest.approx(np.mean(scl))
    assert values["scl_std"] == pytest.approx(np.std(scl, ddof=1))
    assert values["scl_slope"] == pytest.approx(0.001)


def test_no_scr_events_produce_zero_scr_features_and_short_input_raises():
    values = _extract_eda_features(
        np.ones(1000), np.zeros(1000), np.arange(1000) / 100, 100.0
    )
    assert values["scr_amplitude_mean"] == 0
    assert values["scr_amplitude_std"] == 0
    assert values["scr_rise_time_mean"] == 0
    assert values["scr_rise_time_std"] == 0
    assert values["scr_per_min"] == 0
    with pytest.raises(ValueError, match="too short"):
        _extract_eda_features(np.ones(2), np.ones(2), np.arange(2), 100.0)


def test_event_pairing_uses_immediately_preceding_minimum():
    maxima, minima = _pair_scr_events(
        np.asarray([2, 5, 8]), np.asarray([1, 4, 6, 9])
    )
    np.testing.assert_array_equal(maxima, [2, 5, 8])
    np.testing.assert_array_equal(minima, [1, 4, 6])
    maxima, minima = _pair_scr_events(np.asarray([1, 4]), np.asarray([2, 3]))
    np.testing.assert_array_equal(maxima, [4])
    np.testing.assert_array_equal(minima, [3])


@pytest.mark.parametrize(
    "amplitude,expected_count",
    [(1e-8, 0), (0.009999, 0), (0.01, 1), (0.02, 1)],
)
def test_minimum_scr_amplitude_is_inclusive(monkeypatch, amplitude, expected_count):
    monkeypatch.setattr(eda_module, "filtfilt", lambda b, a, values: np.asarray(values))
    scr = np.asarray([0.0, -amplitude / 2, amplitude / 2, 0.0])
    accepted, rise_times, count = _scr_processing(
        scr, np.arange(scr.size) / 100, 100.0, minimum_scr_amplitude=0.01
    )
    assert count == expected_count
    if expected_count:
        np.testing.assert_allclose(accepted, [amplitude], rtol=0, atol=1e-15)
        np.testing.assert_allclose(rise_times, [0.01], rtol=0, atol=0)
    else:
        assert accepted.size == rise_times.size == 0


def test_multiple_and_single_valid_trough_to_peak_events(monkeypatch):
    monkeypatch.setattr(eda_module, "filtfilt", lambda b, a, values: np.asarray(values))
    multiple = np.asarray([0.0, -0.02, 0.02, 0.0, -0.03, 0.03, 0.0])
    amplitudes, rise_times, count = _scr_processing(
        multiple, np.arange(multiple.size) / 100, 100.0
    )
    np.testing.assert_allclose(amplitudes, [0.04, 0.06])
    np.testing.assert_allclose(rise_times, [0.01, 0.01])
    assert count == 2

    one = np.asarray([0.0, -0.005, 0.005, 0.0])
    amplitudes, rise_times, count = _scr_processing(
        one, np.arange(one.size) / 100, 100.0
    )
    np.testing.assert_allclose(amplitudes, [0.01])
    assert count == 1
    assert rise_times.tolist() == [0.01]

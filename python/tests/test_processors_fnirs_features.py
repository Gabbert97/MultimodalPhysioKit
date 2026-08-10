from __future__ import annotations

import numpy as np
import pytest

from multimodalphysiokit.processors.fnirs import (
    _extract_fnirs_channel_features, _fnirs_features, _highpass_filter_concentration,
    _matlab_power_band_indices,
)


def concentration_signal():
    time = np.arange(0, 60, 0.1)
    signal = 0.002 * np.sin(2 * np.pi * 0.04 * time) + 0.0005 * np.sin(2 * np.pi * 0.2 * time)
    return signal, time


def test_highpass_filter_is_finite_and_deterministic():
    signal, _ = concentration_signal()
    first = _highpass_filter_concentration(signal, 10.0)
    second = _highpass_filter_concentration(signal, 10.0)
    np.testing.assert_array_equal(first, second)
    assert np.all(np.isfinite(first))


def test_highpass_filter_suppresses_dc_and_trend_without_phase_shift():
    time = np.arange(600) / 10.0

    constant = _highpass_filter_concentration(np.ones(time.size), 10.0)
    trend = _highpass_filter_concentration(np.linspace(0.0, 1.0, time.size), 10.0)
    passband = np.sin(2 * np.pi * 0.1 * time)
    passband_filtered = _highpass_filter_concentration(passband, 10.0)

    assert constant.shape == time.shape
    assert np.all(np.isfinite(constant))
    np.testing.assert_allclose(constant, 0.0, rtol=0, atol=1e-13)
    assert np.std(trend) < 0.02 * np.std(np.linspace(0.0, 1.0, time.size))
    assert np.max(np.abs(passband_filtered)) < 1.2

    interior = slice(100, -100)
    retained_amplitude = np.std(passband_filtered[interior]) / np.std(passband[interior])
    assert 0.98 < retained_amplitude < 1.02
    correlation = np.correlate(passband_filtered[interior], passband[interior], mode="full")
    lag = int(np.argmax(correlation) - (passband[interior].size - 1))
    assert lag == 0


def test_all_nineteen_channel_features_remain_finite_and_ordered():
    signal, time = concentration_signal()
    features = _extract_fnirs_channel_features(
        signal, time, 10.0, channel_prefix="fnirs_oxy"
    )
    assert len(features) == 19
    assert all(key.startswith("fnirs_oxy_") and "__" not in key for key in features)
    assert all(np.isfinite(value) for value in features.values())


def test_feature_formulas_use_existing_complete_signal_calculation():
    signal, time = concentration_signal()
    filtered = _highpass_filter_concentration(signal, 10.0)
    values = _fnirs_features(signal, filtered, time, 10.0)
    assert values["mean"] == pytest.approx(np.mean(signal))
    assert values["variance"] == pytest.approx(np.var(filtered, ddof=1))
    assert values["std"] == pytest.approx(np.std(filtered, ddof=1))
    assert values["slope"] == pytest.approx(np.polyfit(time, signal, 1)[0])
    expected_auc = np.sum(0.5 * np.diff(time) * np.abs(signal[:-1] + signal[1:])) / (time[-1] - time[0])
    assert values["auc"] == pytest.approx(expected_auc)


def test_moment_features_use_explicit_historical_biased_convention():
    signal, time = concentration_signal()
    filtered = _highpass_filter_concentration(signal, 10.0)
    values = _fnirs_features(signal, filtered, time, 10.0)

    from scipy.stats import kurtosis, skew

    assert values["skewness"] == pytest.approx(skew(filtered, bias=True))
    assert values["kurtosis"] == pytest.approx(
        kurtosis(filtered, fisher=False, bias=True)
    )


def test_power_band_indices_match_matlab_one_based_translation():
    # The reverse 95% threshold is reached at MATLAB index 3. MATLAB maps
    # that to upper index 2 (one-based), hence Python index 1—not index 3.
    psd = np.array([0.0, 1.0, 2.0, 3.0, 4.0, 90.0])
    assert _matlab_power_band_indices(psd) == (5, 1)


def test_invalid_highpass_frequency_is_rejected():
    with pytest.raises(ValueError, match="designed for"):
        _highpass_filter_concentration(concentration_signal()[0], 20.0)

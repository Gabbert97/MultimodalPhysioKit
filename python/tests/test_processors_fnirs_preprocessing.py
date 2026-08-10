from __future__ import annotations

import numpy as np
import pytest
from scipy.signal import sosfreqz

from multimodalphysiokit.core import Signal
from multimodalphysiokit.processors import FNIRSProcessor
from multimodalphysiokit.processors.fnirs import (
    _baseline_values_from_filtered, _concentration_changes, _downsample_signal,
    _design_lowpass_current_sos, _explicit_zero_phase_sos_filter,
    _lowpass_filter_current,
)


def signals(*, duration=60.0, fs=100.0, offset=0.0):
    time = offset + np.arange(0, duration, 1 / fs)
    red = 1.5 + 0.02 * np.sin(2 * np.pi * 0.05 * time)
    infrared = 1.8 + 0.025 * np.sin(2 * np.pi * 0.04 * time + 0.1)
    return Signal(red, time, fs, label="red"), Signal(infrared, time, fs, label="infrared")


def processor():
    return FNIRSProcessor(age_years=35)


@pytest.mark.parametrize("age,error", [("35", TypeError), (0, ValueError), (-1, ValueError), (np.inf, ValueError)])
def test_invalid_age(age, error):
    with pytest.raises(error, match="age_years"):
        FNIRSProcessor(age_years=age)


def test_invalid_target_types():
    red, infrared = signals()
    with pytest.raises(TypeError, match="red_signal"):
        processor().process("red", infrared)
    with pytest.raises(TypeError, match="infrared_signal"):
        processor().process(red, "infrared")


def test_mismatched_lengths_frequencies_and_times():
    red, infrared = signals()
    shorter = Signal(infrared.samples[:-1], infrared.time[:-1], infrared.sampling_frequency)
    with pytest.raises(ValueError, match="same number"):
        processor().process(red, shorter)
    _, different_fs = signals(fs=200.0)
    with pytest.raises(ValueError, match="same number|sampling_frequency"):
        processor().process(red, different_fs)
    shifted = Signal(infrared.samples, infrared.time + 0.01, infrared.sampling_frequency)
    with pytest.raises(ValueError, match="synchronized"):
        processor().process(red, shifted)


def test_nonpositive_current_and_unsupported_frequency():
    red, infrared = signals()
    invalid = Signal(np.where(np.arange(red.n_samples) == 10, 0.0, red.samples), red.time, red.sampling_frequency)
    with pytest.raises(ValueError, match="positive"):
        processor().process(invalid, infrared)
    red_25, infrared_25 = signals(fs=25.0)
    with pytest.raises(ValueError, match="integer multiple of 10"):
        processor().process(red_25, infrared_25)


def test_downsampling_is_numerically_identical_to_validated_helper():
    red, infrared = signals(offset=8.0)
    result = processor().downsample_currents(red, infrared)
    expected_red, expected_time = _downsample_signal(red.samples, red.time, 10)
    expected_infrared, _ = _downsample_signal(infrared.samples, infrared.time, 10)
    np.testing.assert_array_equal(result["red"], expected_red)
    np.testing.assert_array_equal(result["infrared"], expected_infrared)
    np.testing.assert_array_equal(result["time"], expected_time)


def test_lowpass_filtering_is_numerically_identical():
    red, infrared = signals()
    down = processor().downsample_currents(red, infrared)
    filtered = processor().filter_currents(red, infrared)
    np.testing.assert_allclose(filtered["red"], _lowpass_filter_current(down["red"]), rtol=0, atol=0)
    np.testing.assert_allclose(filtered["infrared"], _lowpass_filter_current(down["infrared"]), rtol=0, atol=0)


def test_lowpass_filter_design_matches_canonical_butterworth_specification():
    sos, order, critical_frequency_hz = _design_lowpass_current_sos()

    assert order == 15
    assert critical_frequency_hz == pytest.approx(0.2172764098151067, rel=1e-13)
    assert sos.shape == (8, 6)
    assert np.all(np.isfinite(sos))
    assert all(np.all(np.abs(np.roots(section[3:])) < 1.0) for section in sos)

    frequencies = np.array([0.0, 0.2, 0.4, 1.0])
    _, response = sosfreqz(sos, worN=frequencies, fs=10.0)
    attenuation_db = -20.0 * np.log10(np.abs(response))

    assert abs(response[0]) == pytest.approx(1.0, abs=1e-12)
    assert attenuation_db[1] <= 0.5
    assert attenuation_db[2] == pytest.approx(80.0, abs=1e-9)
    assert attenuation_db[3] > 80.0


def test_explicit_zero_phase_filter_preserves_length_and_constant_signal():
    sos, order, _ = _design_lowpass_current_sos()
    samples = np.full(600, 2.75)

    filtered = _explicit_zero_phase_sos_filter(
        sos, samples, padding_length=3 * order
    )

    assert filtered.shape == samples.shape
    assert np.all(np.isfinite(filtered))
    np.testing.assert_allclose(filtered, samples, rtol=0, atol=5e-13)


def test_explicit_zero_phase_filter_is_deterministic_without_endpoint_explosion():
    sos, order, _ = _design_lowpass_current_sos()
    time = np.arange(600) / 10.0
    samples = np.sin(2 * np.pi * 0.1 * time)

    first = _explicit_zero_phase_sos_filter(
        sos, samples, padding_length=3 * order
    )
    second = _explicit_zero_phase_sos_filter(
        sos, samples, padding_length=3 * order
    )

    np.testing.assert_array_equal(first, second)
    assert np.max(np.abs(first)) < 1.1
    interior = slice(100, -100)
    correlation = np.correlate(first[interior], samples[interior], mode="full")
    lag = int(np.argmax(correlation) - (samples[interior].size - 1))
    assert lag == 0


def test_compute_baseline_values_uses_filtered_middle_exactly():
    red, infrared = signals()
    filtered = processor().filter_currents(red, infrared)
    expected = (float(np.mean(filtered["red"][100:-100])), float(np.mean(filtered["infrared"][100:-100])))
    actual = processor().compute_baseline_values(red, infrared)
    assert actual == pytest.approx(expected)
    assert all(isinstance(value, float) for value in actual)
    assert actual == pytest.approx(_baseline_values_from_filtered(filtered["red"], filtered["infrared"]))


def test_self_baseline_requires_more_than_twenty_seconds():
    with pytest.raises(ValueError, match="exceed 20 seconds"):
        processor().process(*signals(duration=20.0))


def test_self_process_uses_compute_baseline_values():
    red, infrared = signals()
    expected = processor().compute_baseline_values(red, infrared)
    result = processor().process(red, infrared)
    assert result.metadata["baseline_red_value"] == pytest.approx(expected[0])
    assert result.metadata["baseline_infrared_value"] == pytest.approx(expected[1])


def test_supplied_values_are_used_directly_and_skip_reference_calculation(monkeypatch):
    red, infrared = signals(duration=30.0)
    monkeypatch.setattr(
        "multimodalphysiokit.processors.fnirs._baseline_values_from_filtered",
        lambda *args: (_ for _ in ()).throw(AssertionError("must not calculate reference")),
    )
    result = processor().process(red, infrared, baseline_red_value=1.5, baseline_infrared_value=1.8)
    assert result.metadata["baseline_red_value"] == 1.5
    assert result.metadata["baseline_infrared_value"] == 1.8


def test_optical_density_and_concentrations_match_existing_equations():
    red, infrared = signals()
    filtered = processor().filter_currents(red, infrared)
    baseline = processor().compute_baseline_values(red, infrared)
    expected_oxy, expected_deoxy = _concentration_changes(
        filtered["red"], filtered["infrared"], baseline[0], baseline[1], 35.0
    )
    result = processor().compute_concentration_changes(
        red, infrared, baseline_red_value=baseline[0], baseline_infrared_value=baseline[1]
    )
    np.testing.assert_allclose(result["oxy"], expected_oxy, rtol=0, atol=0)
    np.testing.assert_allclose(result["deoxy"], expected_deoxy, rtol=0, atol=0)

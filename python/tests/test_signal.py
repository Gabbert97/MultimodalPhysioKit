import numpy as np
import pytest

from multimodalphysiokit import Signal


def test_valid_construction():
    signal = Signal([1, 2, 3], [0, 0.5, 1.0], 2, label="ecg", units="mV")

    assert isinstance(signal.samples, np.ndarray)
    assert isinstance(signal.time, np.ndarray)
    assert signal.samples.dtype == float
    assert signal.time.dtype == float
    assert signal.label == "ecg"
    assert signal.units == "mV"


def test_mismatched_samples_time_length_raises():
    with pytest.raises(ValueError, match="same length"):
        Signal([1, 2], [0], 1)


def test_empty_samples_time_raises():
    with pytest.raises(ValueError, match="samples must be non-empty"):
        Signal([], [], 1)

    with pytest.raises(ValueError, match="time must be non-empty"):
        Signal([1], [], 1)


def test_non_finite_samples_raise():
    with pytest.raises(ValueError, match="samples.*finite"):
        Signal([1, np.nan], [0, 1], 1)


def test_non_finite_time_raises():
    with pytest.raises(ValueError, match="time.*finite"):
        Signal([1, 2], [0, np.inf], 1)


def test_non_monotonic_time_raises():
    with pytest.raises(ValueError, match="monotonically increasing"):
        Signal([1, 2, 3], [0, 1, 1], 1)


def test_zero_or_negative_sampling_frequency_raises():
    with pytest.raises(ValueError, match="positive and finite"):
        Signal([1], [0], 0)

    with pytest.raises(ValueError, match="positive and finite"):
        Signal([1], [0], -1)


def test_n_samples_works():
    signal = Signal([1, 2, 3], [0, 1, 2], 1)

    assert signal.n_samples == 3


def test_duration_seconds_works():
    signal = Signal([1, 2, 3], [2, 2.5, 3.25], 2)

    assert signal.duration_seconds == pytest.approx(1.25)


def test_metadata_defaults_to_dict():
    signal = Signal([1], [0], 1)

    assert signal.metadata == {}
    assert isinstance(signal.metadata, dict)


def test_copy_with_metadata_works():
    signal = Signal([1, 2], [0, 1], 1, metadata={"source": "array"})

    copied = signal.copy_with_metadata(subject="anonymous")

    assert copied is not signal
    assert copied.metadata == {"source": "array", "subject": "anonymous"}
    assert signal.metadata == {"source": "array"}
    np.testing.assert_array_equal(copied.samples, signal.samples)
    np.testing.assert_array_equal(copied.time, signal.time)


def test_remove_time_ranges_removes_one_leading_interval():
    signal = Signal(np.arange(5), np.arange(5), 1)

    trimmed = signal.remove_time_ranges([(0, 2)])

    np.testing.assert_array_equal(trimmed.samples, [2, 3, 4])
    np.testing.assert_allclose(trimmed.time, [0, 1, 2])


def test_remove_time_ranges_removes_one_middle_interval():
    signal = Signal(np.arange(6), np.arange(6), 1)

    trimmed = signal.remove_time_ranges([(2, 4)])

    np.testing.assert_array_equal(trimmed.samples, [0, 1, 4, 5])
    np.testing.assert_allclose(trimmed.time, [0, 1, 2, 3])


def test_remove_time_ranges_removes_trailing_open_ended_interval():
    signal = Signal(np.arange(6), np.arange(6), 1)

    trimmed = signal.remove_time_ranges([(3, None)])

    np.testing.assert_array_equal(trimmed.samples, [0, 1, 2])
    np.testing.assert_allclose(trimmed.time, [0, 1, 2])


def test_remove_time_ranges_removes_multiple_intervals():
    signal = Signal(np.arange(8), np.arange(8), 1)

    trimmed = signal.remove_time_ranges([(1, 3), (5, 7)])

    np.testing.assert_array_equal(trimmed.samples, [0, 3, 4, 7])
    np.testing.assert_allclose(trimmed.time, [0, 1, 2, 3])


def test_remove_time_ranges_reset_time_false_preserves_timestamps():
    signal = Signal(np.arange(6), np.arange(6), 1)

    trimmed = signal.remove_time_ranges([(2, 4)], reset_time=False)

    np.testing.assert_array_equal(trimmed.samples, [0, 1, 4, 5])
    np.testing.assert_allclose(trimmed.time, [0, 1, 4, 5])


def test_remove_time_ranges_preserves_signal_fields():
    signal = Signal([1, 2, 3, 4], [0, 1, 2, 3], 1, label="ecg", units="mV")

    trimmed = signal.remove_time_ranges([(1, 2)])

    assert trimmed.label == "ecg"
    assert trimmed.units == "mV"
    assert trimmed.sampling_frequency == signal.sampling_frequency


def test_remove_time_ranges_metadata_is_preserved_and_augmented():
    signal = Signal([1, 2, 3, 4], [0, 1, 2, 3], 1, metadata={"source": "array"})

    trimmed = signal.remove_time_ranges([(1, 2)], reset_time=False)

    assert trimmed.metadata["source"] == "array"
    assert trimmed.metadata["time_ranges_removed_seconds"] == [(1.0, 2.0)]
    assert trimmed.metadata["reset_time_after_removal"] is False
    assert signal.metadata == {"source": "array"}


def test_remove_time_ranges_removing_all_samples_raises_value_error():
    signal = Signal([1, 2, 3], [0, 1, 2], 1)

    with pytest.raises(ValueError, match="all samples"):
        signal.remove_time_ranges([(0, None)])


@pytest.mark.parametrize(
    "intervals",
    [
        [],
        [(np.inf, None)],
        [(-1, 2)],
        [(1, 1)],
        [(2, 1)],
        [(1, np.nan)],
        [(1, 0)],
        [(1,)],
    ],
)
def test_remove_time_ranges_invalid_intervals_raise_value_error(intervals):
    signal = Signal([1, 2, 3], [0, 1, 2], 1)

    with pytest.raises(ValueError):
        signal.remove_time_ranges(intervals)


@pytest.mark.parametrize("intervals", [None, "0,1", [(0, "end")]])
def test_remove_time_ranges_invalid_intervals_raise_type_error(intervals):
    signal = Signal([1, 2, 3], [0, 1, 2], 1)

    with pytest.raises((TypeError, ValueError)):
        signal.remove_time_ranges(intervals)


def test_remove_time_ranges_invalid_reset_time_raises_type_error():
    signal = Signal([1, 2, 3], [0, 1, 2], 1)

    with pytest.raises(TypeError, match="reset_time"):
        signal.remove_time_ranges([(0, 1)], reset_time="yes")

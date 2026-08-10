from __future__ import annotations

import numpy as np
import pytest
import warnings

from multimodalphysiokit import Recording
from multimodalphysiokit.utils import create_sliding_phase_intervals


def test_one_valid_temporal_window():
    intervals, labels = create_sliding_phase_intervals(0, 10, 10)
    np.testing.assert_array_equal(intervals, [[0.0], [10.0]])
    assert labels == ["window_001"]


def test_multiple_nonoverlapping_windows_are_half_open():
    intervals, labels = create_sliding_phase_intervals(2, 8, 2)
    np.testing.assert_array_equal(intervals, [[2, 4, 6], [4, 6, 8]])
    assert intervals.shape == (2, 3)
    assert intervals.dtype.kind == "f"
    assert labels == ["window_001", "window_002", "window_003"]
    np.testing.assert_array_equal(intervals[1, :-1], intervals[0, 1:])


def test_overlapping_windows_and_custom_prefix():
    intervals, labels = create_sliding_phase_intervals(
        0, 5, 2, overlap_seconds=1, label_prefix="segment"
    )
    np.testing.assert_array_equal(intervals, [[0, 1, 2, 3], [2, 3, 4, 5]])
    assert labels == ["segment_001", "segment_002", "segment_003", "segment_004"]


def test_only_complete_windows_are_generated_and_trailing_warns():
    with pytest.warns(UserWarning, match="trailing"):
        intervals, _ = create_sliding_phase_intervals(0, 5.5, 2)
    np.testing.assert_array_equal(intervals, [[0, 2], [2, 4]])


def test_trailing_warning_can_be_disabled():
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        intervals, _ = create_sliding_phase_intervals(0, 5.5, 2, warn_on_trailing=False)
    assert caught == []
    np.testing.assert_array_equal(intervals, [[0, 2], [2, 4]])


def test_no_complete_window_raises():
    with pytest.raises(ValueError, match="no complete window"):
        create_sliding_phase_intervals(0, 1, 2)


@pytest.mark.parametrize("overlap", [2, 3])
def test_overlap_must_be_smaller_than_window(overlap):
    with pytest.raises(ValueError, match="smaller"):
        create_sliding_phase_intervals(0, 10, 2, overlap_seconds=overlap)


def test_negative_overlap_and_nonpositive_window_are_rejected():
    with pytest.raises(ValueError, match="nonnegative"):
        create_sliding_phase_intervals(0, 10, 2, overlap_seconds=-1)
    with pytest.raises(ValueError, match="positive"):
        create_sliding_phase_intervals(0, 10, 0)


@pytest.mark.parametrize("start,end", [(-1, 5), (2, 2), (3, 2)])
def test_invalid_start_end_range(start, end):
    with pytest.raises(ValueError):
        create_sliding_phase_intervals(start, end, 1)


@pytest.mark.parametrize("argument", ["start", "end", "window", "overlap"])
def test_nonfinite_inputs(argument):
    values = {"start": 0.0, "end": 10.0, "window": 2.0, "overlap": 0.0}
    values[argument] = np.inf
    with pytest.raises(ValueError, match="finite"):
        create_sliding_phase_intervals(
            values["start"], values["end"], values["window"],
            overlap_seconds=values["overlap"],
        )


def test_empty_converted_prefix_is_rejected():
    with pytest.raises(ValueError, match="nonempty"):
        create_sliding_phase_intervals(0, 4, 2, label_prefix="")


def test_generated_values_can_be_assigned_directly_to_recording():
    intervals, labels = create_sliding_phase_intervals(0, 6, 2)
    recording = Recording()
    recording.set_phases(intervals, labels)
    np.testing.assert_array_equal(recording.phase_intervals_seconds, intervals)
    assert recording.phase_labels == labels


def test_window_starts_do_not_accumulate_floating_point_addition_error():
    intervals, _ = create_sliding_phase_intervals(
        0.1, 1.0, 0.2, overlap_seconds=0.1, warn_on_trailing=False
    )
    expected = 0.1 + np.arange(intervals.shape[1]) * 0.1
    np.testing.assert_array_equal(intervals[0], expected)

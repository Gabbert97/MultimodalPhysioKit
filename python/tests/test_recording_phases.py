import numpy as np
import pytest

from multimodalphysiokit import Recording, Signal


def make_signal(*, sampling_frequency=1.0, start=0.0, count=6, label="ecg"):
    time = start + np.arange(count, dtype=float) / sampling_frequency
    return Signal(
        np.arange(count, dtype=float),
        time,
        sampling_frequency,
        label=label,
        units="mV",
        metadata={"channel": 2},
    )


def test_construction_without_phases_uses_empty_phase_values():
    recording = Recording()

    assert recording.phase_intervals_seconds.shape == (2, 0)
    assert recording.phase_intervals_seconds.dtype.kind == "f"
    assert recording.phase_labels == []
    assert recording.n_phases == 0


def test_construction_with_one_valid_phase():
    recording = Recording(phase_intervals_seconds=[[0], [2]], phase_labels=["baseline"])

    np.testing.assert_allclose(recording.phase_intervals_seconds, [[0], [2]])
    assert recording.phase_labels == ["baseline"]


def test_construction_with_multiple_valid_phases_preserves_order():
    recording = Recording(
        phase_intervals_seconds=[[2, 0, 1], [3, 1, 2]],
        phase_labels=["last", "first", "middle"],
    )

    np.testing.assert_allclose(recording.phase_intervals_seconds, [[2, 0, 1], [3, 1, 2]])
    assert recording.phase_labels == ["last", "first", "middle"]


def test_generated_default_labels():
    recording = Recording(phase_intervals_seconds=[[0, 1, 2], [1, 2, 3]])

    assert recording.phase_labels == ["phase_1", "phase_2", "phase_3"]


def test_set_phases_replaces_existing_phases_and_converts_labels_to_strings():
    recording = Recording(phase_intervals_seconds=[[0], [1]], phase_labels=["old"])

    recording.set_phases([[1, 3], [2, 5]], labels=[10, "new"])

    np.testing.assert_allclose(recording.phase_intervals_seconds, [[1, 3], [2, 5]])
    assert recording.phase_labels == ["10", "new"]


def test_clear_phases():
    recording = Recording(phase_intervals_seconds=[[0], [1]])

    recording.clear_phases()

    assert recording.phase_intervals_seconds.shape == (2, 0)
    assert recording.phase_intervals_seconds.dtype.kind == "f"
    assert recording.phase_labels == []


def test_n_phases_and_phase_durations_seconds():
    recording = Recording(phase_intervals_seconds=[[0, 2], [1.5, 5]])

    assert recording.n_phases == 2
    np.testing.assert_allclose(recording.phase_durations_seconds, [1.5, 3.0])


@pytest.mark.parametrize("intervals", [[0, 1], [[0, 1, 2]], np.zeros((3, 2))])
def test_invalid_interval_shape(intervals):
    with pytest.raises(ValueError, match=r"shape \(2, n_phases\)"):
        Recording(phase_intervals_seconds=intervals)


def test_nonnumeric_intervals():
    with pytest.raises(ValueError, match="numeric"):
        Recording(phase_intervals_seconds=[["start"], ["end"]])


@pytest.mark.parametrize("value", [np.nan, np.inf, -np.inf])
def test_nonfinite_interval_values(value):
    with pytest.raises(ValueError, match="finite"):
        Recording(phase_intervals_seconds=[[0], [value]])


def test_negative_interval_values():
    with pytest.raises(ValueError, match="nonnegative"):
        Recording(phase_intervals_seconds=[[-1], [1]])


@pytest.mark.parametrize("intervals", [[[1], [1]], [[2], [1]]])
def test_phase_start_must_be_smaller_than_end(intervals):
    with pytest.raises(ValueError, match="smaller"):
        Recording(phase_intervals_seconds=intervals)


def test_labels_without_intervals():
    with pytest.raises(ValueError, match="without phase intervals"):
        Recording(phase_labels=["baseline"])


def test_label_count_mismatch():
    with pytest.raises(ValueError, match="one label per phase"):
        Recording(phase_intervals_seconds=[[0, 1], [1, 2]], phase_labels=["one"])


def test_empty_label():
    with pytest.raises(ValueError, match="nonempty"):
        Recording(phase_intervals_seconds=[[0], [1]], phase_labels=[""])


def test_duplicate_labels():
    with pytest.raises(ValueError, match="unique"):
        Recording(
            phase_intervals_seconds=[[0, 1], [1, 2]],
            phase_labels=["same", "same"],
        )


def test_external_interval_array_mutation_does_not_modify_recording():
    intervals = np.array([[0.0], [2.0]])
    recording = Recording(phase_intervals_seconds=intervals)

    intervals[1, 0] = 99

    assert recording.phase_intervals_seconds[1, 0] == 2


def test_external_label_list_mutation_does_not_modify_recording():
    labels = ["baseline"]
    recording = Recording(phase_intervals_seconds=[[0], [2]], phase_labels=labels)

    labels[0] = "changed"

    assert recording.phase_labels == ["baseline"]


def test_same_temporal_interval_maps_to_different_sampling_frequencies():
    recording = Recording(
        signals={"slow": make_signal(count=5), "fast": make_signal(sampling_frequency=2, count=9)},
        phase_intervals_seconds=[[1], [3]],
    )

    np.testing.assert_array_equal(recording.phase_indices("slow"), [[1], [3]])
    np.testing.assert_array_equal(recording.phase_indices("fast"), [[2], [6]])


def test_boundaries_exactly_matching_sample_timestamps():
    signal = make_signal()
    recording = Recording(phase_intervals_seconds=[[1], [4]])

    np.testing.assert_array_equal(recording.phase_indices(signal), [[1], [4]])


def test_boundaries_between_sample_timestamps():
    signal = make_signal()
    recording = Recording(phase_intervals_seconds=[[1.2], [3.8]])

    np.testing.assert_array_equal(recording.phase_indices(signal), [[2], [4]])


def test_contiguous_phases_do_not_share_samples():
    signal = make_signal()
    recording = Recording(
        phase_intervals_seconds=[[0, 3], [3, 6]],
        phase_labels=["first", "second"],
    )

    indices = recording.phase_indices(signal)
    assert indices[1, 0] == indices[0, 1]
    assert set(range(*indices[:, 0])).isdisjoint(range(*indices[:, 1]))


def test_overlapping_phases_are_allowed():
    signal = make_signal()
    recording = Recording(
        phase_intervals_seconds=[[0, 2], [4, 5]],
        phase_labels=["wide", "overlap"],
    )

    np.testing.assert_array_equal(recording.phase_indices(signal), [[0, 2], [4, 5]])


def test_signal_with_nonzero_time_origin():
    signal = make_signal(start=10)
    recording = Recording(phase_intervals_seconds=[[11], [13]])

    np.testing.assert_array_equal(recording.phase_indices(signal), [[1], [3]])


def test_selecting_all_phases_preserves_recording_order():
    signal = make_signal()
    recording = Recording(
        phase_intervals_seconds=[[0, 2], [2, 4]],
        phase_labels=["first", "second"],
    )

    np.testing.assert_array_equal(recording.phase_indices(signal), [[0, 2], [2, 4]])


def test_selecting_one_phase():
    signal = make_signal()
    recording = Recording(
        phase_intervals_seconds=[[0, 2], [2, 4]],
        phase_labels=["first", "second"],
    )

    np.testing.assert_array_equal(recording.phase_indices(signal, "second"), [[2], [4]])


def test_selecting_multiple_phases_preserves_requested_order():
    signal = make_signal()
    recording = Recording(
        phase_intervals_seconds=[[0, 2], [2, 4]],
        phase_labels=["first", "second"],
    )

    np.testing.assert_array_equal(
        recording.phase_indices(signal, ["second", "first"]),
        [[2, 0], [4, 2]],
    )


def test_requesting_unknown_phase():
    recording = Recording(phase_intervals_seconds=[[0], [1]])

    with pytest.raises(ValueError, match="missing"):
        recording.phase_indices(make_signal(), "missing")


def test_requesting_unknown_signal_name():
    recording = Recording(phase_intervals_seconds=[[0], [1]])

    with pytest.raises(KeyError, match="missing"):
        recording.phase_indices("missing")


def test_phase_outside_signal_range_is_rejected():
    recording = Recording(phase_intervals_seconds=[[10], [11]], phase_labels=["outside"])

    with pytest.raises(ValueError, match="not fully represented"):
        recording.phase_indices(make_signal())


def test_partially_represented_phase_is_rejected():
    recording = Recording(phase_intervals_seconds=[[0], [7]], phase_labels=["partial"])

    with pytest.raises(ValueError, match="not fully represented"):
        recording.phase_indices(make_signal())


def test_phase_mapping_to_no_samples_is_rejected():
    recording = Recording(phase_intervals_seconds=[[0.1], [0.2]], phase_labels=["empty"])

    with pytest.raises(ValueError, match="maps to no signal samples"):
        recording.phase_indices(make_signal())


def test_phase_indices_rejects_recording_without_phases():
    with pytest.raises(ValueError, match="no phases"):
        Recording().phase_indices(make_signal())


def test_extract_phase_with_reset_time_true():
    recording = Recording(
        "session",
        signals={"ecg": make_signal()},
        phase_intervals_seconds=[[2], [5]],
        phase_labels=["task"],
    )

    extracted = recording.extract_phase("ecg", "task")

    np.testing.assert_array_equal(extracted.samples, [2, 3, 4])
    np.testing.assert_allclose(extracted.time, [0, 1, 2])


def test_extract_phase_with_reset_time_false():
    recording = Recording(
        signals={"ecg": make_signal(start=10)},
        phase_intervals_seconds=[[12], [15]],
        phase_labels=["task"],
    )

    extracted = recording.extract_phase("ecg", "task", reset_time=False)

    np.testing.assert_allclose(extracted.time, [12, 13, 14])


def test_extracted_signal_preserves_properties_and_augments_copied_metadata():
    source = make_signal()
    recording = Recording(
        "session",
        signals={"lead": source},
        phase_intervals_seconds=[[1], [3]],
        phase_labels=["baseline"],
    )

    extracted = recording.extract_phase("lead", "baseline")

    assert extracted.sampling_frequency == source.sampling_frequency
    assert extracted.label == source.label
    assert extracted.units == source.units
    assert extracted.metadata == {
        "channel": 2,
        "source_recording": "session",
        "phase_label": "baseline",
        "phase_start_seconds": 1.0,
        "phase_end_seconds": 3.0,
        "time_reset": True,
    }
    assert extracted.metadata is not source.metadata
    assert source.metadata == {"channel": 2}


def test_remove_time_ranges_still_works_without_phases():
    recording = Recording(signals={"ecg": make_signal()})

    trimmed = recording.remove_time_ranges([(2, 4)])

    np.testing.assert_array_equal(trimmed.get_signal("ecg").samples, [0, 1, 4, 5])
    assert trimmed.n_phases == 0


def test_remove_time_ranges_rejects_recording_with_phases():
    recording = Recording(
        signals={"ecg": make_signal()},
        phase_intervals_seconds=[[0], [2]],
    )

    with pytest.raises(ValueError, match="cleared or redefined"):
        recording.remove_time_ranges([(2, 4)])

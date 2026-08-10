import pytest
import numpy as np

from multimodalphysiokit import Recording, Signal


def make_signal(label="ecg"):
    return Signal([1, 2, 3], [0, 1, 2], 1, label=label)


def test_valid_construction_from_signal_objects():
    ecg = make_signal("ecg")
    recording = Recording("session", signals={"ecg": ecg}, metadata={"site": "lab"})

    assert recording.name == "session"
    assert recording.get_signal("ecg") is ecg
    assert recording.metadata == {"site": "lab"}


def test_add_signal_get_signal_has_signal():
    recording = Recording()
    signal = make_signal("eda")

    recording.add_signal("eda", signal)

    assert recording.get_signal("eda") is signal
    assert recording.has_signal("eda")
    assert not recording.has_signal("ecg")


def test_missing_signal_raises_key_error():
    recording = Recording()

    with pytest.raises(KeyError, match="missing"):
        recording.get_signal("missing")


def test_signal_names_returns_expected_names():
    recording = Recording(signals={"ecg": make_signal("ecg"), "eda": make_signal("eda")})

    assert recording.signal_names() == ["ecg", "eda"]


def test_len_recording_works():
    recording = Recording(signals={"ecg": make_signal("ecg"), "eda": make_signal("eda")})

    assert len(recording) == 2


def test_contains_works():
    recording = Recording(signals={"ecg": make_signal("ecg")})

    assert "ecg" in recording
    assert "eda" not in recording


def test_invalid_non_signal_value_raises():
    with pytest.raises(TypeError, match="Signal"):
        Recording(signals={"ecg": [1, 2, 3]})

    recording = Recording()
    with pytest.raises(TypeError, match="Signal"):
        recording.add_signal("ecg", [1, 2, 3])


def test_remove_time_ranges_applies_to_all_signals():
    recording = Recording(
        "session",
        signals={
            "ecg": Signal(np.arange(6), np.arange(6), 1, label="ecg"),
            "eda": Signal(np.arange(10, 16), np.arange(6), 1, label="eda"),
        },
    )

    trimmed = recording.remove_time_ranges([(2, 4)])

    np.testing.assert_array_equal(trimmed.get_signal("ecg").samples, [0, 1, 4, 5])
    np.testing.assert_array_equal(trimmed.get_signal("eda").samples, [10, 11, 14, 15])
    np.testing.assert_allclose(trimmed.get_signal("ecg").time, [0, 1, 2, 3])
    np.testing.assert_allclose(trimmed.get_signal("eda").time, [0, 1, 2, 3])


def test_remove_time_ranges_does_not_mutate_original_recording():
    recording = Recording(signals={"ecg": Signal(np.arange(5), np.arange(5), 1)})

    trimmed = recording.remove_time_ranges([(1, 3)])

    assert trimmed is not recording
    np.testing.assert_array_equal(recording.get_signal("ecg").samples, [0, 1, 2, 3, 4])
    np.testing.assert_allclose(recording.get_signal("ecg").time, [0, 1, 2, 3, 4])


def test_remove_time_ranges_preserves_recording_name_and_augments_metadata():
    recording = Recording(
        "session",
        signals={"ecg": Signal(np.arange(5), np.arange(5), 1)},
        metadata={"site": "lab"},
    )

    trimmed = recording.remove_time_ranges([(1, 3)], reset_time=False)

    assert trimmed.name == "session"
    assert trimmed.metadata["site"] == "lab"
    assert trimmed.metadata["time_ranges_removed_seconds"] == [(1.0, 3.0)]
    assert trimmed.metadata["reset_time_after_removal"] is False
    assert recording.metadata == {"site": "lab"}


def test_remove_time_ranges_reset_time_false_preserves_signal_timestamps():
    recording = Recording(signals={"ecg": Signal(np.arange(5), np.arange(5), 1)})

    trimmed = recording.remove_time_ranges([(1, 3)], reset_time=False)

    np.testing.assert_allclose(trimmed.get_signal("ecg").time, [0, 3, 4])

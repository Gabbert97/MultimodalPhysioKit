from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pytest

from multimodalphysiokit import Recording, Signal
from multimodalphysiokit.utils import plot_signal_with_phases


def recording(*, sampling_frequency=1.0):
    time = np.arange(0, 8, 1 / sampling_frequency)
    signal = Signal(np.sin(time), time, sampling_frequency, label="ecg", units="mV")
    return Recording(
        signals={"ecg": signal},
        phase_intervals_seconds=[[0, 4], [4, 8]],
        phase_labels=["baseline", "task"],
    )


def test_accepts_recording_and_signal_name():
    fig, ax = plot_signal_with_phases(recording(), "ecg")
    assert fig is ax.figure
    assert len(ax.lines) == 5
    assert [text.get_text() for text in ax.texts] == ["baseline", "task"]


def test_missing_phases_and_unknown_signal_raise():
    empty = Recording(signals={"ecg": Signal([1, 2], [0, 1], 1)})
    with pytest.raises(ValueError, match="no phases"):
        plot_signal_with_phases(empty, "ecg")
    with pytest.raises(KeyError, match="missing"):
        plot_signal_with_phases(recording(), "missing")


def test_uses_recording_phase_mapper(monkeypatch):
    value = recording()
    calls = []
    original = value.phase_indices
    def capture(name, phases=None):
        calls.append((name, phases))
        return original(name, phases)
    monkeypatch.setattr(value, "phase_indices", capture)
    plot_signal_with_phases(value, "ecg")
    assert calls == [("ecg", None)]


def test_different_sampling_frequencies_map_to_same_temporal_boundaries():
    slow = recording(sampling_frequency=1)
    fast = recording(sampling_frequency=4)
    _, slow_ax = plot_signal_with_phases(slow, "ecg")
    _, fast_ax = plot_signal_with_phases(fast, "ecg")
    slow_boundaries = [line.get_xdata()[0] for line in slow_ax.lines[1:]]
    fast_boundaries = [line.get_xdata()[0] for line in fast_ax.lines[1:]]
    assert slow_boundaries == fast_boundaries == [0.0, 4.0, 4.0, 8.0]


def test_show_false_closes_figure():
    plt.close("all")
    plot_signal_with_phases(recording(), "ecg", show=False)
    assert plt.get_fignums() == []


def test_save_behavior(tmp_path):
    plot_signal_with_phases(recording(), "ecg", save_plots=True, output_dir=tmp_path)
    assert (tmp_path / "ecg_phases.png").is_file()
    with pytest.raises(ValueError, match="output_dir"):
        plot_signal_with_phases(recording(), "ecg", save_plots=True)


def test_plotting_does_not_modify_recording_or_signal():
    value = recording()
    signal = value.get_signal("ecg")
    samples = signal.samples.copy()
    time = signal.time.copy()
    intervals = value.phase_intervals_seconds.copy()
    labels = list(value.phase_labels)
    plot_signal_with_phases(value, "ecg")
    np.testing.assert_array_equal(signal.samples, samples)
    np.testing.assert_array_equal(signal.time, time)
    np.testing.assert_array_equal(value.phase_intervals_seconds, intervals)
    assert value.phase_labels == labels


def test_existing_visual_options_and_validation():
    fig, ax = plt.subplots()
    returned_fig, returned_ax = plot_signal_with_phases(
        recording(), "ecg", ax=ax, show_phase_labels=False,
        signal_display_name="Lead I", title="Title",
    )
    assert returned_fig is fig
    assert returned_ax is ax
    assert ax.lines[0].get_label() == "Lead I"
    assert ax.get_title() == "Title"
    assert len(ax.texts) == 0
    with pytest.raises(ValueError, match="phase_alpha"):
        plot_signal_with_phases(recording(), "ecg", phase_alpha=2)

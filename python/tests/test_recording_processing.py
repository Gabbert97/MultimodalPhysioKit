from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from multimodalphysiokit import Recording, Signal
from multimodalphysiokit.processors import TemperatureProcessor


class CaptureProcessor:
    def __init__(self, result=None):
        self.calls = []
        self.result = object() if result is None else result

    def process(self, *signals, **kwargs):
        self.calls.append((signals, kwargs))
        return self.result


def make_recording(*, labels=("baseline", "task")):
    time = np.arange(8, dtype=float)
    return Recording(
        "session",
        signals={
            "first": Signal(np.arange(8), time, 1, label="first"),
            "second": Signal(np.arange(10, 18), time, 1, label="second"),
        },
        phase_intervals_seconds=[[0, 4], [4, 8]],
        phase_labels=list(labels),
    )


def test_process_signal_one_signal_returns_exact_result():
    recording = make_recording()
    expected = object()
    processor = CaptureProcessor(expected)
    result = recording.process_signal(processor, "first")
    assert result is expected
    assert processor.calls[0][0] == (recording.get_signal("first"),)


def test_process_signal_multiple_signals_preserves_requested_order():
    recording = make_recording()
    processor = CaptureProcessor()
    recording.process_signal(processor, ("second", "first"))
    assert processor.calls[0][0] == (
        recording.get_signal("second"), recording.get_signal("first")
    )


def test_process_signal_forwards_copied_kwargs():
    recording = make_recording()
    processor = CaptureProcessor()
    kwargs = {"flag": [1]}
    recording.process_signal(processor, "first", process_kwargs=kwargs)
    assert processor.calls[0][1] == kwargs
    assert processor.calls[0][1] is not kwargs
    assert kwargs == {"flag": [1]}


@pytest.mark.parametrize("processor", [object(), type("Bad", (), {"process": 1})()])
def test_invalid_processor_is_rejected(processor):
    with pytest.raises(TypeError, match="process"):
        make_recording().process_signal(processor, "first")


@pytest.mark.parametrize(
    "names,error,match",
    [([], ValueError, "empty"), (["first", 2], TypeError, "strings"),
     (["first", "first"], ValueError, "duplicates"), (["missing"], KeyError, "missing")],
)
def test_invalid_signal_names_are_rejected(names, error, match):
    with pytest.raises(error, match=match):
        make_recording().process_signal(CaptureProcessor(), names)


def test_invalid_process_kwargs_is_rejected():
    with pytest.raises(TypeError, match="process_kwargs"):
        make_recording().process_signal(CaptureProcessor(), "first", process_kwargs=[])


def test_process_phases_all_in_recording_order_and_once_each():
    recording = make_recording()
    processor = CaptureProcessor()
    results = recording.process_phases(processor, "first")
    assert list(results) == ["baseline", "task"]
    assert len(processor.calls) == 2
    assert all(result is processor.result for result in results.values())


def test_process_phases_accepts_string_and_requested_sequence_order():
    recording = make_recording()
    assert list(recording.process_phases(CaptureProcessor(), "first", phases="task")) == ["task"]
    assert list(
        recording.process_phases(CaptureProcessor(), "first", phases=["task", "baseline"])
    ) == ["task", "baseline"]


@pytest.mark.parametrize(
    "recording,selection,error,match",
    [(make_recording(), [], ValueError, "empty"),
     (make_recording(), ["missing"], ValueError, "missing"),
     (Recording(signals={"first": Signal([1], [0], 1)}), None, ValueError, "no phases"),
     (make_recording(), [1], TypeError, "phase-label")],
)
def test_invalid_phase_selection_is_rejected(recording, selection, error, match):
    with pytest.raises(error, match=match):
        recording.process_phases(CaptureProcessor(), "first", phases=selection)


def test_phase_calls_receive_extract_phase_membership_and_reset_times():
    recording = make_recording()
    processor = CaptureProcessor()
    recording.process_phases(processor, "first")
    for label, call in zip(recording.phase_labels, processor.calls, strict=True):
        received = call[0][0]
        expected = recording.extract_phase("first", label)
        np.testing.assert_array_equal(received.samples, expected.samples)
        np.testing.assert_array_equal(received.time, expected.time)
        assert received.metadata["phase_label"] == label


def test_reset_time_false_preserves_original_timestamps():
    recording = make_recording()
    processor = CaptureProcessor()
    recording.process_phases(processor, "first", phases="task", reset_time=False)
    np.testing.assert_array_equal(processor.calls[0][0][0].time, [4, 5, 6, 7])


def test_multimodal_phase_signals_preserve_argument_order():
    recording = make_recording()
    processor = CaptureProcessor()
    recording.process_phases(processor, ("second", "first"), phases="baseline")
    received = processor.calls[0][0]
    assert [signal.label for signal in received] == ["second", "first"]
    np.testing.assert_array_equal(received[0].samples, [10, 11, 12, 13])
    np.testing.assert_array_equal(received[1].samples, [0, 1, 2, 3])


def test_phase_results_and_feature_names_are_not_modified():
    result = type("Result", (), {"features": {"feature": 1.0}})()
    results = make_recording().process_phases(CaptureProcessor(result), "first")
    assert all(value is result for value in results.values())
    assert result.features == {"feature": 1.0}


def test_processor_exception_propagates_unchanged():
    expected = RuntimeError("scientific failure")
    class Failing:
        def process(self, signal):
            raise expected
    with pytest.raises(RuntimeError) as caught:
        make_recording().process_phases(Failing(), "first")
    assert caught.value is expected


def test_kwargs_are_forwarded_to_every_phase_without_mutation():
    processor = CaptureProcessor()
    kwargs = {"baseline_red_value": 1.25, "baseline_infrared_value": 2.5}
    make_recording().process_phases(processor, ("first", "second"), process_kwargs=kwargs)
    assert [call[1] for call in processor.calls] == [kwargs, kwargs]
    assert all(call[1] is not kwargs for call in processor.calls)
    assert kwargs == {"baseline_red_value": 1.25, "baseline_infrared_value": 2.5}


def test_phase_output_directories_are_created_and_forwarded(tmp_path):
    recording = make_recording(labels=("base line", "task/one"))
    processor = CaptureProcessor()
    kwargs = {"output_dir": tmp_path}
    recording.process_phases(processor, "first", process_kwargs=kwargs)
    directories = [call[1]["output_dir"] for call in processor.calls]
    assert directories == [tmp_path / "base_line", tmp_path / "task_one"]
    assert all(directory.is_dir() for directory in directories)
    assert kwargs == {"output_dir": tmp_path}


def test_sanitized_directory_collisions_are_unique_and_deterministic(tmp_path):
    recording = make_recording(labels=("a b", "a@b"))
    first = CaptureProcessor()
    second = CaptureProcessor()
    recording.process_phases(first, "first", process_kwargs={"output_dir": tmp_path})
    recording.process_phases(second, "first", process_kwargs={"output_dir": tmp_path})
    assert [call[1]["output_dir"].name for call in first.calls] == ["a_b", "a_b_2"]
    assert [call[1]["output_dir"].name for call in second.calls] == ["a_b", "a_b_2"]


def test_process_signal_forwards_output_directory_without_subdirectory(tmp_path):
    processor = CaptureProcessor()
    make_recording().process_signal(
        processor, "first", process_kwargs={"output_dir": tmp_path}
    )
    assert processor.calls[0][1]["output_dir"] == tmp_path
    assert list(tmp_path.iterdir()) == []


def test_phase_specific_intermediates_remain_separate():
    class IntermediateProcessor:
        def process(self, signal):
            return {"intermediates": {"stft_matrix": np.asarray([[signal.samples[0]]])}}
    results = make_recording().process_phases(IntermediateProcessor(), "first")
    assert results["baseline"] is not results["task"]
    assert results["baseline"]["intermediates"]["stft_matrix"][0, 0] == 0
    assert results["task"]["intermediates"]["stft_matrix"][0, 0] == 4


def test_temperature_processor_integration_returns_one_phase_independent_result_each():
    sampling_frequency = 10.0
    time = np.arange(0, 80, 1 / sampling_frequency)
    samples = 32 + 0.02 * np.sin(2 * np.pi * 0.02 * time)
    recording = Recording(
        signals={"temperature": Signal(samples, time, sampling_frequency, label="temperature")},
        phase_intervals_seconds=[[0, 40], [40, 80]],
        phase_labels=["first", "second"],
    )
    results = recording.process_phases(TemperatureProcessor(), "temperature")
    assert list(results) == ["first", "second"]
    assert all(len(result.features) == 8 for result in results.values())
    assert all("__" not in key for result in results.values() for key in result.features)

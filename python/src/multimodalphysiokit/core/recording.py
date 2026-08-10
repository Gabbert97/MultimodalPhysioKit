from __future__ import annotations

import re
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

from .signal import Signal, _validate_time_ranges


@dataclass
class Recording:
    """Named collection of physiological ``Signal`` objects."""

    name: str = "recording"
    signals: dict[str, Signal] | None = None
    metadata: dict[str, Any] | None = None
    phase_intervals_seconds: Any = None
    phase_labels: Sequence[str] | None = None

    def __post_init__(self) -> None:
        self.name = str(self.name)
        self.metadata = {} if self.metadata is None else dict(self.metadata)

        signals = {} if self.signals is None else dict(self.signals)
        for name, signal in signals.items():
            if not isinstance(signal, Signal):
                raise TypeError(f"signal '{name}' must be a Signal instance")
        self.signals = {str(name): signal for name, signal in signals.items()}

        labels = self.phase_labels
        intervals = self.phase_intervals_seconds
        self.clear_phases()
        if intervals is None:
            if labels is not None:
                raise ValueError("phase labels cannot be supplied without phase intervals")
        else:
            self.set_phases(intervals, labels)

    def set_phases(self, intervals_seconds, labels=None) -> None:
        """Validate and replace this recording's half-open phase intervals."""

        try:
            intervals = np.asarray(intervals_seconds, dtype=float)
        except (TypeError, ValueError) as exc:
            raise ValueError("phase intervals must be numeric") from exc
        if intervals.ndim != 2 or intervals.shape[0] != 2:
            raise ValueError("phase intervals must have shape (2, n_phases)")
        if not np.all(np.isfinite(intervals)):
            raise ValueError("phase intervals must contain only finite values")
        if np.any(intervals < 0):
            raise ValueError("phase intervals must be nonnegative")
        if np.any(intervals[0, :] >= intervals[1, :]):
            raise ValueError("each phase start must be smaller than its end")

        n_phases = int(intervals.shape[1])
        if labels is None:
            phase_labels = [f"phase_{index + 1}" for index in range(n_phases)]
        else:
            if isinstance(labels, (str, bytes)):
                raise ValueError("phase labels must contain one label per phase")
            try:
                phase_labels = [str(label) for label in labels]
            except TypeError as exc:
                raise ValueError("phase labels must contain one label per phase") from exc
            if len(phase_labels) != n_phases:
                raise ValueError("phase labels must contain one label per phase")
            if any(label == "" for label in phase_labels):
                raise ValueError("phase labels must be nonempty")
            if len(set(phase_labels)) != len(phase_labels):
                raise ValueError("phase labels must be unique")

        self.phase_intervals_seconds = intervals.copy()
        self.phase_labels = list(phase_labels)

    def clear_phases(self) -> None:
        """Remove all phase definitions from this recording."""

        self.phase_intervals_seconds = np.empty((2, 0), dtype=float)
        self.phase_labels = []

    @property
    def n_phases(self) -> int:
        """Return the number of defined phases."""

        return int(self.phase_intervals_seconds.shape[1])

    @property
    def phase_durations_seconds(self) -> np.ndarray:
        """Return phase durations in recording order."""

        return self.phase_intervals_seconds[1, :] - self.phase_intervals_seconds[0, :]

    def phase_indices(self, signal_or_name, phases=None) -> np.ndarray:
        """Map selected temporal phases to half-open indexes for a signal."""

        if self.n_phases == 0:
            raise ValueError("recording has no phases")
        signal = self._resolve_signal(signal_or_name)
        selected_indexes = self._select_phase_indexes(phases)
        mapped = np.empty((2, len(selected_indexes)), dtype=int)
        signal_start = float(signal.time[0])
        signal_stop = float(signal.time[-1] + 1.0 / signal.sampling_frequency)

        for output_index, phase_index in enumerate(selected_indexes):
            start_seconds = float(self.phase_intervals_seconds[0, phase_index])
            end_seconds = float(self.phase_intervals_seconds[1, phase_index])
            label = self.phase_labels[phase_index]
            if start_seconds < signal_start or end_seconds > signal_stop:
                raise ValueError(
                    f"phase '{label}' is not fully represented by the signal time range"
                )

            start_index = int(np.searchsorted(signal.time, start_seconds, side="left"))
            stop_index = int(np.searchsorted(signal.time, end_seconds, side="left"))
            if start_index >= stop_index:
                raise ValueError(f"phase '{label}' maps to no signal samples")
            mapped[:, output_index] = (start_index, stop_index)

        return mapped

    def extract_phase(
        self,
        signal_name: str,
        phase: str,
        reset_time: bool = True,
    ) -> Signal:
        """Return one named signal segment for a recording phase."""

        if not isinstance(phase, str):
            raise TypeError("phase must be a phase-label string")
        if not isinstance(reset_time, bool):
            raise TypeError("reset_time must be a boolean")
        signal = self.get_signal(signal_name)
        start_index, stop_index = self.phase_indices(signal, phase)[:, 0]
        selected_time = signal.time[start_index:stop_index].copy()
        if reset_time:
            selected_time -= selected_time[0]

        phase_index = self.phase_labels.index(phase)
        metadata = dict(signal.metadata)
        metadata.update(
            {
                "source_recording": self.name,
                "phase_label": phase,
                "phase_start_seconds": float(self.phase_intervals_seconds[0, phase_index]),
                "phase_end_seconds": float(self.phase_intervals_seconds[1, phase_index]),
                "time_reset": reset_time,
            }
        )
        return Signal(
            samples=signal.samples[start_index:stop_index].copy(),
            time=selected_time,
            sampling_frequency=signal.sampling_frequency,
            label=signal.label,
            units=signal.units,
            metadata=metadata,
        )

    def process_signal(
        self,
        processor,
        signal_names,
        *,
        process_kwargs=None,
    ):
        """Process one or more complete stored signals with a configured processor."""

        process = _validate_processor(processor)
        names = self._normalize_signal_names(signal_names)
        kwargs = _normalize_process_kwargs(process_kwargs)
        signals = tuple(self.get_signal(name) for name in names)
        return process(*signals, **kwargs)

    def process_phases(
        self,
        processor,
        signal_names,
        *,
        phases=None,
        reset_time: bool = True,
        process_kwargs=None,
    ) -> dict[str, object]:
        """Extract and independently process selected recording phases."""

        process = _validate_processor(processor)
        names = self._normalize_signal_names(signal_names)
        if not isinstance(reset_time, bool):
            raise TypeError("reset_time must be a boolean")
        if self.n_phases == 0:
            raise ValueError("recording has no phases")
        selected_indexes = self._select_phase_indexes(phases)
        selected_labels = tuple(self.phase_labels[index] for index in selected_indexes)
        base_kwargs = _normalize_process_kwargs(process_kwargs)
        output_directories = _phase_output_directories(
            selected_labels,
            base_kwargs.get("output_dir"),
        )

        results: dict[str, object] = {}
        for index, label in enumerate(selected_labels):
            signals = tuple(
                self.extract_phase(name, label, reset_time=reset_time) for name in names
            )
            kwargs = dict(base_kwargs)
            if output_directories is not None:
                kwargs["output_dir"] = output_directories[index]
            results[label] = process(*signals, **kwargs)
        return results

    def add_signal(self, name: str, signal: Signal) -> None:
        """Add or replace a named signal."""

        if not isinstance(signal, Signal):
            raise TypeError("signal must be a Signal instance")
        self.signals[str(name)] = signal

    def get_signal(self, name: str) -> Signal:
        """Return a named signal or raise ``KeyError`` if it is missing."""

        key = str(name)
        try:
            return self.signals[key]
        except KeyError as exc:
            raise KeyError(f"recording does not contain signal '{key}'") from exc

    def has_signal(self, name: str) -> bool:
        """Return whether a named signal exists."""

        return str(name) in self.signals

    def signal_names(self) -> list[str]:
        """Return signal names in insertion order."""

        return list(self.signals.keys())

    def remove_time_ranges(
        self,
        intervals_seconds,
        *,
        reset_time: bool = True,
    ) -> "Recording":
        """Return a new recording with time ranges removed from every signal."""

        if self.n_phases:
            raise ValueError(
                "phases must be cleared or redefined before removing time ranges"
            )

        intervals = _validate_time_ranges(intervals_seconds)
        if not isinstance(reset_time, bool):
            raise TypeError("reset_time must be a boolean")
        signals = {
            name: signal.remove_time_ranges(intervals, reset_time=reset_time)
            for name, signal in self.signals.items()
        }
        metadata = dict(self.metadata)
        metadata.update(
            {
                "time_ranges_removed_seconds": [
                    (start_seconds, end_seconds) for start_seconds, end_seconds in intervals
                ],
                "reset_time_after_removal": reset_time,
            }
        )
        return Recording(name=self.name, signals=signals, metadata=metadata)

    def _resolve_signal(self, signal_or_name) -> Signal:
        if isinstance(signal_or_name, Signal):
            return signal_or_name
        if isinstance(signal_or_name, str):
            return self.get_signal(signal_or_name)
        raise TypeError("signal_or_name must be a Signal or stored signal name")

    def _select_phase_indexes(self, phases) -> list[int]:
        if phases is None:
            labels = list(self.phase_labels)
        elif isinstance(phases, str):
            labels = [phases]
        else:
            if isinstance(phases, (bytes, bytearray)):
                raise TypeError("phases must be phase-label strings")
            try:
                labels = list(phases)
            except TypeError as exc:
                raise TypeError("phases must be phase-label strings") from exc
            if any(not isinstance(label, str) for label in labels):
                raise TypeError("phases must be phase-label strings")
            if not labels:
                raise ValueError("phases must not be empty")

        indexes = []
        for label in labels:
            try:
                indexes.append(self.phase_labels.index(label))
            except ValueError as exc:
                raise ValueError(f"recording does not contain phase '{label}'") from exc
        return indexes

    def _normalize_signal_names(self, signal_names) -> tuple[str, ...]:
        if isinstance(signal_names, str):
            names = (signal_names,)
        else:
            if isinstance(signal_names, (bytes, bytearray)):
                raise TypeError("signal_names must contain strings")
            try:
                names = tuple(signal_names)
            except TypeError as exc:
                raise TypeError("signal_names must be a string or sequence of strings") from exc
            if not names:
                raise ValueError("signal_names must not be empty")
            if any(not isinstance(name, str) for name in names):
                raise TypeError("signal_names must contain strings")
        if len(set(names)) != len(names):
            raise ValueError("signal_names must not contain duplicates")
        for name in names:
            if not self.has_signal(name):
                raise KeyError(f"recording does not contain signal '{name}'")
        return names

    def __contains__(self, name: object) -> bool:
        return str(name) in self.signals

    def __len__(self) -> int:
        return len(self.signals)


def _validate_processor(processor):
    try:
        process = processor.process
    except AttributeError as exc:
        raise TypeError("processor must expose a callable process attribute") from exc
    if not callable(process):
        raise TypeError("processor.process must be callable")
    return process


def _normalize_process_kwargs(process_kwargs) -> dict[str, object]:
    if process_kwargs is None:
        return {}
    if not isinstance(process_kwargs, dict):
        raise TypeError("process_kwargs must be a dictionary or None")
    return dict(process_kwargs)


def _phase_output_directories(
    labels: tuple[str, ...],
    parent_directory,
) -> tuple[Path, ...] | None:
    if parent_directory is None:
        return None
    parent = Path(parent_directory)
    used: set[str] = set()
    directories = []
    for index, label in enumerate(labels, start=1):
        base = re.sub(r"[^A-Za-z0-9_-]+", "_", label).strip("_-") or f"phase_{index}"
        candidate = base
        suffix = 2
        while candidate in used:
            candidate = f"{base}_{suffix}"
            suffix += 1
        used.add(candidate)
        directory = parent / candidate
        directory.mkdir(parents=True, exist_ok=True)
        directories.append(directory)
    return tuple(directories)

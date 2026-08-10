from __future__ import annotations

import numbers
from pathlib import Path

import matplotlib.pyplot as plt

from multimodalphysiokit.core import Recording


def plot_signal_with_phases(
    recording: Recording,
    signal_name: str,
    *,
    title: str | None = None,
    show_phase_labels: bool = True,
    phase_alpha: float = 0.08,
    ax=None,
    signal_display_name: str | None = None,
    max_phase_labels: int = 25,
    show: bool = False,
    save_plots: bool = False,
    output_dir: str | Path | None = None,
):
    """Plot one stored signal with the recording's temporal phase overlays."""

    if not isinstance(recording, Recording):
        raise TypeError("recording must be a Recording instance")
    if not isinstance(signal_name, str):
        raise TypeError("signal_name must be a string")
    signal = recording.get_signal(signal_name)
    if recording.n_phases == 0:
        raise ValueError("recording has no phases")
    mapped = recording.phase_indices(signal_name)

    if not isinstance(show_phase_labels, bool):
        raise TypeError("show_phase_labels must be a boolean")
    if not isinstance(phase_alpha, numbers.Real) or not 0 <= phase_alpha <= 1:
        raise ValueError("phase_alpha must be between 0 and 1")
    if (
        isinstance(max_phase_labels, bool)
        or not isinstance(max_phase_labels, int)
        or max_phase_labels < 0
    ):
        raise ValueError("max_phase_labels must be a nonnegative integer")
    if not isinstance(show, bool):
        raise TypeError("show must be a boolean")
    if not isinstance(save_plots, bool):
        raise TypeError("save_plots must be a boolean")
    if save_plots and output_dir is None:
        raise ValueError("output_dir must be provided when save_plots=True")

    if ax is None:
        fig, ax = plt.subplots()
    else:
        fig = ax.figure

    display_name = signal_display_name or signal.label or "Signal"
    ax.plot(signal.time, signal.samples, label=display_name)

    _, y_max = ax.get_ylim()
    for index, label in enumerate(recording.phase_labels):
        start_index = int(mapped[0, index])
        stop_index = int(mapped[1, index])
        start_time = float(signal.time[start_index])
        end_time = (
            float(signal.time[stop_index])
            if stop_index < signal.n_samples
            else float(recording.phase_intervals_seconds[1, index])
        )
        ax.axvspan(start_time, end_time, alpha=phase_alpha)
        ax.axvline(start_time, linestyle="--", linewidth=0.8)
        ax.axvline(end_time, linestyle=":", linewidth=0.8)

        if show_phase_labels and recording.n_phases <= max_phase_labels:
            ax.text(
                (start_time + end_time) / 2,
                y_max,
                label,
                ha="center",
                va="top",
                fontsize=8,
            )

    ax.set_title(title or signal.label or "Signal with phases")
    ax.set_xlabel("Time (s)")
    ax.set_ylabel(f"Signal ({signal.units})" if signal.units else "Signal")
    ax.grid(True)

    if save_plots:
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        fig.savefig(output_path / f"{signal_name}_phases.png", dpi=150, bbox_inches="tight")
    if show:
        plt.show()
    else:
        plt.close(fig)
    return fig, ax

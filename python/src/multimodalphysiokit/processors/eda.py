from __future__ import annotations

import importlib
import importlib.util
import hashlib
from collections.abc import Callable
from pathlib import Path
from types import ModuleType

import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import butter, filtfilt

from multimodalphysiokit.core import ProcessingResult, Signal


_FS_NEW = 100
_DEFAULT_MINIMUM_SCR_AMPLITUDE = 0.01
_CVXEDA_PARAMETERS = {
    "tau0": 2.0,
    "tau1": 0.7,
    "delta_knot": 10.0,
    "alpha": 8e-4,
    "gamma": 1e-2,
}
_MISSING_CVXEDA_MESSAGE = (
    "EDA decomposition requires external cvxEDA v1.1.0. Install it with "
    "`pip install cvxeda==1.1.0`, "
    "provide cvxeda_path, or provide decomposition_callable."
)
_EDA_FEATURE_NAMES = (
    "scl_mean",
    "scl_std",
    "scl_slope",
    "scr_amplitude_mean",
    "scr_amplitude_std",
    "scr_rise_time_mean",
    "scr_rise_time_std",
    "scr_per_min",
)


class EDAProcessor:
    """Object-oriented EDA preprocessing and external cvxEDA decomposition."""

    def __init__(
        self,
        *,
        cvxeda_path: str | Path | None = None,
        decomposition_callable: Callable | None = None,
        cache_directory: str | Path | None = None,
        output_label: str = "eda",
        minimum_scr_amplitude: float = _DEFAULT_MINIMUM_SCR_AMPLITUDE,
    ) -> None:
        if not isinstance(output_label, str):
            raise TypeError("output_label must be a string")
        if cvxeda_path is not None and not isinstance(cvxeda_path, str | Path):
            raise TypeError("cvxeda_path must be None, a string, or a Path")
        if cache_directory is not None and not isinstance(cache_directory, str | Path):
            raise TypeError("cache_directory must be None, a string, or a Path")
        if decomposition_callable is not None and not callable(decomposition_callable):
            raise TypeError("decomposition_callable must be callable")
        if (
            not isinstance(minimum_scr_amplitude, (int, float, np.integer, np.floating))
            or not np.isfinite(minimum_scr_amplitude)
            or minimum_scr_amplitude < 0
        ):
            raise ValueError("minimum_scr_amplitude must be a finite value greater than or equal to 0")

        self.cvxeda_path = None if cvxeda_path is None else Path(cvxeda_path)
        self.decomposition_callable = decomposition_callable
        self.cache_directory = None if cache_directory is None else Path(cache_directory)
        self.output_label = output_label
        self.minimum_scr_amplitude = float(minimum_scr_amplitude)

    def downsample_signal(self, signal: Signal) -> dict[str, object]:
        """Downsample EDA to 100 Hz."""

        _validate_signal(signal)
        down = int(round(signal.sampling_frequency / _FS_NEW))
        eda_downsampled = signal.samples[::down].copy()
        time_downsampled = signal.time[::down].copy()

        return {
            "eda": eda_downsampled,
            "time": time_downsampled,
            "sampling_frequency": _FS_NEW,
        }

    def filter_signal(self, signal: Signal) -> dict[str, object]:
        """Low-pass filter downsampled EDA."""

        downsampled = self.downsample_signal(signal)
        eda_downsampled = np.asarray(downsampled["eda"], dtype=float)
        eda_filtered = _filter_eda(eda_downsampled)

        return {
            "eda_filtered": eda_filtered,
            "eda_downsampled": eda_downsampled,
            "time": downsampled["time"],
            "sampling_frequency": downsampled["sampling_frequency"],
        }

    def decompose(
        self,
        signal: Signal,
        *,
        show_plots: bool = False,
        save_plots: bool = False,
        output_dir: str | Path | None = None,
        output_label: str | None = None,
    ) -> dict[str, object]:
        """Run or load external cvxEDA decomposition."""

        _validate_plot_options(show_plots=show_plots, save_plots=save_plots, output_dir=output_dir)
        if output_label is not None and not isinstance(output_label, str):
            raise TypeError("output_label must be a string")
        plot_output_label = self.output_label if output_label is None else output_label
        filtered = self.filter_signal(signal)
        eda_filtered = np.asarray(filtered["eda_filtered"], dtype=float)
        cache_path = self._cache_path(signal)
        if cache_path is not None:
            cached = _load_decomposition_cache(cache_path)
            if cached is not None:
                result = {
                    **cached,
                    "eda_filtered": eda_filtered,
                    "eda_downsampled": np.asarray(filtered["eda_downsampled"], dtype=float),
                    "time": filtered["time"],
                    "sampling_frequency": filtered["sampling_frequency"],
                }
                _plot_decomposition_if_requested(
                    result,
                    minimum_scr_amplitude=self.minimum_scr_amplitude,
                    output_label=plot_output_label,
                    show_plots=show_plots,
                    save_plots=save_plots,
                    output_dir=output_dir,
                )
                return result

        cvxeda_callable = _resolve_cvxeda_callable(
            decomposition_callable=self.decomposition_callable,
            cvxeda_path=self.cvxeda_path,
        )
        sampling_frequency = float(filtered["sampling_frequency"])
        if self.decomposition_callable is not None:
            output = list(cvxeda_callable(eda_filtered, 1 / sampling_frequency))
        else:
            output = list(
                cvxeda_callable(
                    eda_filtered,
                    1 / sampling_frequency,
                    **_CVXEDA_PARAMETERS,
                )
            )
        scr, scl, residuals = _normalize_decomposition_output(output, eda_filtered.size)

        if cache_path is not None:
            _save_decomposition_cache(cache_path, scr=scr, scl=scl, residuals=residuals)

        result = {
            "scr": scr,
            "scl": scl,
            "residuals": residuals,
            "eda_filtered": eda_filtered,
            "eda_downsampled": np.asarray(filtered["eda_downsampled"], dtype=float),
            "time": filtered["time"],
            "sampling_frequency": filtered["sampling_frequency"],
        }
        _plot_decomposition_if_requested(
            result,
            minimum_scr_amplitude=self.minimum_scr_amplitude,
            output_label=plot_output_label,
            show_plots=show_plots,
            save_plots=save_plots,
            output_dir=output_dir,
        )
        return result

    def extract_features(
        self,
        signal: Signal,
        *,
        show_plots: bool = False,
        save_plots: bool = False,
        output_dir: str | Path | None = None,
        output_label: str | None = None,
    ) -> dict[str, float]:
        """Extract MATLAB-derived EDA features."""

        if output_label is not None and not isinstance(output_label, str):
            raise TypeError("output_label must be a string")
        decomposition = self.decompose(
            signal,
            show_plots=show_plots,
            save_plots=save_plots,
            output_dir=output_dir,
            output_label=output_label,
        )

        values = _extract_eda_features(
            np.asarray(decomposition["scl"], dtype=float),
            np.asarray(decomposition["scr"], dtype=float),
            np.asarray(decomposition["time"], dtype=float),
            float(decomposition["sampling_frequency"]),
            minimum_scr_amplitude=self.minimum_scr_amplitude,
        )
        features = {f"eda_{name}": float(values[name]) for name in _EDA_FEATURE_NAMES}

        return features

    def process(
        self,
        signal: Signal,
        *,
        show_plots: bool = False,
        save_plots: bool = False,
        output_dir: str | Path | None = None,
        output_label: str | None = None,
    ) -> ProcessingResult:
        """Extract EDA features and return a ProcessingResult."""

        features = self.extract_features(
            signal,
            show_plots=show_plots,
            save_plots=save_plots,
            output_dir=output_dir,
            output_label=output_label,
        )
        return ProcessingResult(
            modality="eda",
            processor_name="process_eda",
            features=features,
            signal_labels=[signal.label or "eda"],
            metadata={
                "source": "process_eda",
                "original_sampling_frequency": signal.sampling_frequency,
                "processed_sampling_frequency": _FS_NEW,
                "original_n_samples": signal.samples.size,
                "processed_n_samples": signal.samples[:: int(round(signal.sampling_frequency / _FS_NEW))].size,
                "start_time_seconds": float(signal.time[0]),
                "end_time_seconds": float(signal.time[-1]),
                "cvxeda_source": _cvxeda_source(self.decomposition_callable, self.cvxeda_path),
                "minimum_scr_amplitude_us": self.minimum_scr_amplitude,
            },
        )

    def _cache_path(self, signal: Signal) -> Path | None:
        if self.cache_directory is None:
            return None
        digest = hashlib.sha256()
        digest.update(np.asarray(signal.samples).tobytes())
        digest.update(np.asarray(signal.time).tobytes())
        digest.update(np.asarray([signal.sampling_frequency], dtype=float).tobytes())
        return self.cache_directory / f"{self.output_label}_{digest.hexdigest()[:16]}_decomposition.npz"


def process_eda(
    signal: Signal,
    *,
    cvxeda_path: str | Path | None = None,
    decomposition_callable: Callable | None = None,
    cache_directory: str | Path | None = None,
    output_label: str = "eda",
    minimum_scr_amplitude: float = _DEFAULT_MINIMUM_SCR_AMPLITUDE,
    show_plots: bool = False,
    save_plots: bool = False,
    output_dir: str | Path | None = None,
) -> ProcessingResult:
    """Functional wrapper for EDAProcessor.process."""

    return EDAProcessor(
        cvxeda_path=cvxeda_path,
        decomposition_callable=decomposition_callable,
        cache_directory=cache_directory,
        output_label=output_label,
        minimum_scr_amplitude=minimum_scr_amplitude,
    ).process(
        signal,
        show_plots=show_plots,
        save_plots=save_plots,
        output_dir=output_dir,
        output_label=output_label,
    )


def _resolve_cvxeda_callable(
    *,
    decomposition_callable: Callable | None = None,
    cvxeda_path: str | Path | None = None,
) -> Callable:
    if decomposition_callable is not None:
        if not callable(decomposition_callable):
            raise TypeError("decomposition_callable must be callable")
        return decomposition_callable

    if cvxeda_path is not None:
        return _load_cvxeda_from_path(cvxeda_path)

    for module_name in ("cvxeda", "cvxeda.cvxEDA", "cvxEDA"):
        try:
            module = importlib.import_module(module_name)
        except ImportError:
            continue
        cvxeda = getattr(module, "cvxEDA", None)
        if callable(cvxeda):
            return cvxeda

    raise ImportError(_MISSING_CVXEDA_MESSAGE)


def _load_cvxeda_from_path(cvxeda_path: str | Path) -> Callable:
    path = Path(cvxeda_path)
    file_path = path / "cvxEDA.py" if path.is_dir() else path
    if not file_path.is_file():
        raise ImportError(f"cvxEDA.py was not found at {file_path}")

    spec = importlib.util.spec_from_file_location("multimodalphysiokit_external_cvxEDA", file_path)
    if spec is None or spec.loader is None:
        raise ImportError(f"cvxEDA.py could not be loaded from {file_path}")
    module = importlib.util.module_from_spec(spec)
    try:
        spec.loader.exec_module(module)
    except Exception as exc:  # pragma: no cover - importlib preserves source exception detail.
        raise ImportError(f"cvxEDA.py could not be loaded from {file_path}") from exc
    return _get_cvxeda_from_module(module, file_path)


def _get_cvxeda_from_module(module: ModuleType, file_path: Path) -> Callable:
    cvxeda = getattr(module, "cvxEDA", None)
    if not callable(cvxeda):
        raise ImportError(f"cvxEDA.py at {file_path} does not define callable cvxEDA")
    return cvxeda


def _validate_signal(signal: Signal) -> None:
    if not isinstance(signal, Signal):
        raise TypeError("signal must be a Signal instance")
    ratio = signal.sampling_frequency / _FS_NEW
    if not np.isclose(ratio, round(ratio)) or ratio < 1:
        raise ValueError("EDA sampling_frequency must be an integer multiple of 100 Hz")
    if signal.samples.size < 2 or signal.time.size < 2:
        raise ValueError("EDA signal is too short")
    if float(signal.time[-1] - signal.time[0]) <= 0:
        raise ValueError("EDA signal duration must be positive")


def _cvxeda_source(decomposition_callable: Callable | None, cvxeda_path: Path | None) -> str:
    if decomposition_callable is not None:
        return "callable"
    if cvxeda_path is not None:
        return "path"
    return "import"


def _filter_eda(eda_downsampled: np.ndarray) -> np.ndarray:
    b, a = butter(5, 0.5 / (_FS_NEW / 2), btype="low")
    try:
        return filtfilt(b, a, eda_downsampled)
    except ValueError as exc:
        raise ValueError("EDA signal is too short for low-pass filtering") from exc


def _extract_eda_features(
    scl: np.ndarray,
    scr: np.ndarray,
    time: np.ndarray,
    sampling_frequency: float,
    *,
    minimum_scr_amplitude: float = _DEFAULT_MINIMUM_SCR_AMPLITUDE,
) -> dict[str, float]:
    if scl.size < 3 or scr.size < 3 or time.size < 3:
        raise ValueError("EDA signal is too short for feature extraction")
    if not (
        np.all(np.isfinite(scl))
        and np.all(np.isfinite(scr))
        and np.all(np.isfinite(time))
    ):
        raise ValueError("EDA inputs must contain finite values")
    duration = float(time[-1] - time[0])
    if duration <= 0:
        raise ValueError("EDA signal duration must be positive")

    amplitude_scr, rise_time_scr, number_scr = _scr_processing(
        scr,
        time,
        sampling_frequency,
        minimum_scr_amplitude=minimum_scr_amplitude,
    )

    if number_scr == 0:
        scr_amplitude_mean = 0.0
        scr_amplitude_std = 0.0
        scr_rise_time_mean = 0.0
        scr_rise_time_std = 0.0
    else:
        scr_amplitude_mean = float(np.mean(amplitude_scr))
        scr_amplitude_std = _matlab_std(amplitude_scr)
        scr_rise_time_mean = float(np.mean(rise_time_scr))
        scr_rise_time_std = _matlab_std(rise_time_scr)

    values = {
        "scl_mean": float(np.mean(scl)),
        "scl_std": _matlab_std(scl),
        "scl_slope": float(np.polyfit(np.arange(1, scl.size + 1), scl, 1)[0]),
        "scr_amplitude_mean": scr_amplitude_mean,
        "scr_amplitude_std": scr_amplitude_std,
        "scr_rise_time_mean": scr_rise_time_mean,
        "scr_rise_time_std": scr_rise_time_std,
        "scr_per_min": float(number_scr / duration * 60),
    }
    if not all(np.isfinite(value) for value in values.values()):
        raise ValueError("non-finite EDA feature value")
    return values


def _scr_processing(
    scr: np.ndarray,
    time: np.ndarray,
    sampling_frequency: float,
    *,
    minimum_scr_amplitude: float = _DEFAULT_MINIMUM_SCR_AMPLITUDE,
) -> tuple[np.ndarray, np.ndarray, int]:
    """Return accepted trough-to-peak SCR amplitudes (uS) and rise times."""
    del time
    if not np.isfinite(minimum_scr_amplitude) or minimum_scr_amplitude < 0:
        raise ValueError("minimum_scr_amplitude must be finite and greater than or equal to 0")
    if scr.size < 3:
        raise ValueError("EDA SCR signal is too short for SCR processing")
    b, a = butter(2, 2.0 / (sampling_frequency / 2), btype="low")
    try:
        scr_filtered = filtfilt(b, a, scr)
    except ValueError as exc:
        raise ValueError("EDA SCR signal is too short for SCR filtering") from exc

    indexes_scr_max, indexes_scr_min = _find_max_min(scr_filtered)
    paired_maxima, paired_minima = _pair_scr_events(indexes_scr_max, indexes_scr_min)
    if paired_maxima.size == 0:
        return np.asarray([], dtype=float), np.asarray([], dtype=float), 0

    amplitudes = scr_filtered[paired_maxima] - scr_filtered[paired_minima]
    rise_times = (paired_maxima - paired_minima) / sampling_frequency
    keep = amplitudes >= minimum_scr_amplitude
    amplitude_scr_checked = amplitudes[keep]
    rise_time_scr_checked = rise_times[keep]
    number_scr = int(amplitude_scr_checked.size)
    if not (np.all(np.isfinite(amplitude_scr_checked)) and np.all(np.isfinite(rise_time_scr_checked))):
        raise ValueError("SCR event features must contain finite values")
    return amplitude_scr_checked, rise_time_scr_checked, number_scr


def _pair_scr_events(
    maxima: np.ndarray,
    minima: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Pair each maximum with its immediately preceding local minimum."""
    paired_maxima: list[int] = []
    paired_minima: list[int] = []
    for maximum in np.asarray(maxima, dtype=int):
        preceding = np.asarray(minima, dtype=int)
        preceding = preceding[preceding < maximum]
        if preceding.size == 0:
            continue
        paired_maxima.append(int(maximum))
        paired_minima.append(int(preceding[-1]))
    return np.asarray(paired_maxima, dtype=int), np.asarray(paired_minima, dtype=int)


def _find_max_min(signal: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    if signal.size < 3:
        return np.asarray([], dtype=int), np.asarray([], dtype=int)
    up = signal[1] - signal[0] > 0
    down = not up
    minima: list[int] = []
    maxima: list[int] = []
    for index in range(1, signal.size - 1):
        if up and signal[index + 1] < signal[index]:
            down = True
            up = False
            maxima.append(index)
        if down and signal[index + 1] > signal[index]:
            up = True
            down = False
            minima.append(index)
    return np.asarray(maxima, dtype=int), np.asarray(minima, dtype=int)


def _matlab_std(values: np.ndarray) -> float:
    values = np.asarray(values, dtype=float).reshape(-1)
    if values.size <= 1:
        return 0.0
    return float(np.std(values, ddof=1))


def _load_decomposition_cache(cache_path: Path) -> dict[str, np.ndarray] | None:
    if not cache_path.is_file():
        return None
    with np.load(cache_path) as cache:
        required = {"scr", "scl", "residuals"}
        if not required.issubset(cache.files):
            raise ValueError(f"EDA decomposition cache is missing required arrays: {cache_path}")
        scr = np.asarray(cache["scr"], dtype=float).reshape(-1)
        scl = np.asarray(cache["scl"], dtype=float).reshape(-1)
        residuals = np.asarray(cache["residuals"], dtype=float).reshape(-1)
    if not (scr.size == scl.size == residuals.size):
        raise ValueError("cached EDA decomposition arrays must have matching lengths")
    if not (np.all(np.isfinite(scr)) and np.all(np.isfinite(scl)) and np.all(np.isfinite(residuals))):
        raise ValueError("cached EDA decomposition arrays must contain finite values")
    return {"scr": scr, "scl": scl, "residuals": residuals}


def _save_decomposition_cache(
    cache_path: Path,
    *,
    scr: np.ndarray,
    scl: np.ndarray,
    residuals: np.ndarray,
) -> None:
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    np.savez(cache_path, scr=scr, scl=scl, residuals=residuals)


def _plot_decomposition_if_requested(
    decomposition: dict[str, object],
    *,
    minimum_scr_amplitude: float,
    output_label: str,
    show_plots: bool,
    save_plots: bool,
    output_dir: str | Path | None,
) -> None:
    if not (show_plots or save_plots):
        return
    _plot_eda_decomposition(
        np.asarray(decomposition["eda_filtered"], dtype=float),
        np.asarray(decomposition["scl"], dtype=float),
        np.asarray(decomposition["scr"], dtype=float),
        np.asarray(decomposition["residuals"], dtype=float),
        np.asarray(decomposition["time"], dtype=float),
        minimum_scr_amplitude=minimum_scr_amplitude,
        output_label=output_label,
        show_plots=show_plots,
        save_plots=save_plots,
        output_dir=output_dir,
    )


def _plot_eda_decomposition(
    eda_filtered: np.ndarray,
    scl: np.ndarray,
    scr: np.ndarray,
    residuals: np.ndarray,
    time: np.ndarray,
    *,
    minimum_scr_amplitude: float = _DEFAULT_MINIMUM_SCR_AMPLITUDE,
    output_label: str,
    show_plots: bool = False,
    save_plots: bool = False,
    output_dir: str | Path | None = None,
) -> None:
    _validate_plot_options(show_plots=show_plots, save_plots=save_plots, output_dir=output_dir)
    fig, ax = plt.subplots()
    ax.plot(time, eda_filtered, linewidth=1.5, label="eda")
    ax.plot(time, scl, linewidth=1.5, label="scl")
    ax.plot(time, scr, linewidth=1.5, label="scr")
    ax.plot(time, residuals, linewidth=1.5, label="residuals")
    b, a = butter(2, 2.0 / (_FS_NEW / 2), btype="low")
    scr_filtered = filtfilt(b, a, scr)
    maxima, minima = _find_max_min(scr_filtered)
    paired_maxima, paired_minima = _pair_scr_events(maxima, minima)
    if paired_maxima.size:
        amplitudes = scr_filtered[paired_maxima] - scr_filtered[paired_minima]
        accepted_maxima = paired_maxima[amplitudes >= minimum_scr_amplitude]
        ax.scatter(
            time[accepted_maxima], scr_filtered[accepted_maxima],
            marker="x", color="black", label="accepted SCR",
        )
    ax.legend()
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("EDA (uS)")
    ax.set_title(output_label)
    ax.grid(True)

    if save_plots:
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        fig.savefig(output_path / f"{output_label}_eda_decomposition.png", dpi=150, bbox_inches="tight")
    if show_plots:
        plt.show()
    else:
        plt.close(fig)


def _validate_plot_options(
    *,
    show_plots: bool,
    save_plots: bool,
    output_dir: str | Path | None,
) -> None:
    if not isinstance(show_plots, bool):
        raise TypeError("show_plots must be a boolean")
    if not isinstance(save_plots, bool):
        raise TypeError("save_plots must be a boolean")
    if save_plots and output_dir is None:
        raise ValueError("output_dir must be provided when save_plots=True")


def _normalize_decomposition_output(output: list[object], expected_length: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    if len(output) < 6:
        raise ValueError("cvxEDA output must contain at least 6 elements")
    scr = _as_finite_1d_array(output[0], expected_length, "scr")
    scl = _as_finite_1d_array(output[2], expected_length, "scl")
    residuals = _as_finite_1d_array(output[5], expected_length, "residuals")
    return scr, scl, residuals


def _as_finite_1d_array(value: object, expected_length: int, name: str) -> np.ndarray:
    array = np.asarray(value, dtype=float).reshape(-1)
    if array.size != expected_length:
        raise ValueError(f"{name} must have the same length as filtered EDA")
    if not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must contain finite values")
    return array

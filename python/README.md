# MultimodalPhysioKit Python guide

MultimodalPhysioKit's Python processing layer is input-agnostic. Standardized
`Signal` objects flow into reusable modality processors and produce
`ProcessingResult` objects:

```text
arrays / existing Signals / mapped HDF5 / supported acquisition export / window
                                  ↓
                               Signal
                                  ↓
                         modality processor
                                  ↓
                         ProcessingResult
```

`Recording` is an optional container and offline orchestrator around this core
flow.

## Installation

MultimodalPhysioKit v0.1.0 is available from
[PyPI](https://pypi.org/project/multimodalphysiokit/):

```bash
pip install multimodalphysiokit
```

For development from a repository checkout, run from the repository root:

```bash
python -m pip install -e "python[test]"
```

Install the notebook dependencies as well with:

```bash
python -m pip install -e "python[test,examples]"
```

Python 3.10 or newer is required. See
[Dependencies](https://github.com/Gabbert97/MultimodalPhysioKit/blob/main/docs/dependencies.md)
for the separately obtained official cvxEDA v1.1.0 software and other
dependency details. EDA processing accepts a configurable
`minimum_scr_amplitude` in µS (default 0.01).

## Core architecture

### Signal

`Signal` standardizes:

- `samples`
- `time`
- `sampling_frequency`
- `label`
- `units`
- `metadata`

It may represent a complete acquisition channel, an extracted experimental
phase, an offline window, or one completed real-time buffer.

```python
import numpy as np
from multimodalphysiokit import Signal

time = np.arange(ecg_samples.size, dtype=float) / 1000.0
ecg_signal = Signal(
    samples=ecg_samples,
    time=time,
    sampling_frequency=1000.0,
    label="ecg",
    units="mV",
)
```

### Processor

A processor is the central plug-and-play component. It accepts one or more
`Signal` objects, performs modality-specific preprocessing and feature
extraction, optionally exposes intermediate outputs, and returns a
`ProcessingResult`.

```python
from multimodalphysiokit.processors import ECGProcessor

processor = ECGProcessor(
    label_frequency_analysis=1,
    return_intermediates=True,
)
result = processor.process(ecg_signal)
```

Processors analyze exactly the supplied signal or synchronized signal pair.
They are independent of its storage system and do not receive phase
definitions.

### Recording

`Recording` is an optional offline container for synchronized named signals.
It preserves modality-specific sampling frequencies, stores shared temporal
phases, maps intervals independently through each signal time vector, and
orchestrates whole-signal or phase processing.

It does not automatically synchronize independent streams, resample signals,
interpolate missing values, infer protocol phases, or apply physiological
processing.

### ProcessingResult

`ProcessingResult` contains:

- `modality`
- `processor_name`
- `features`
- `signal_labels`
- `options`
- `metadata`
- `intermediates`

`features` contains final scalar outputs. `metadata` holds concise descriptive
and scalar processing information. `intermediates` may hold optional arrays or
structured diagnostics.

## Supported input workflows

All input methods converge on `Signal` objects and the same processor API.

### 1. Existing arrays

Use `recording_from_arrays` when arrays share one sampling frequency or one
explicit time vector:

```python
from multimodalphysiokit.io import recording_from_arrays

recording = recording_from_arrays(
    {
        "ecg": ecg_samples,
        "eda": eda_samples,
    },
    sampling_frequency=1000.0,
    units={
        "ecg": "mV",
        "eda": "uS",
    },
)
```

Its signature is:

```text
recording_from_arrays(
    signals,
    *,
    sampling_frequency=None,
    time=None,
    units=None,
    labels=None,
    name="recording",
    metadata=None,
)
```

Provide exactly one of `sampling_frequency` or `time`. The helper does not
assign a different sampling frequency to every array.

For synchronized modalities with different sampling frequencies, construct
individual signals and add them:

```python
from multimodalphysiokit import Recording, Signal

ecg_signal = Signal(
    ecg_samples,
    ecg_time,
    sampling_frequency=1000.0,
    label="ecg",
    units="mV",
)
temperature_signal = Signal(
    temperature_samples,
    temperature_time,
    sampling_frequency=10.0,
    label="skin_temperature",
    units="degC",
)

recording = Recording(name="experiment")
recording.add_signal("ecg", ecg_signal)
recording.add_signal("skin_temperature", temperature_signal)
```

The user is responsible for ensuring that these signals are already
synchronized on a meaningful time axis.

### 2. Existing Signal objects

Direct processing is the most general interface:

```python
result = ECGProcessor(label_frequency_analysis=1).process(ecg_signal)
```

The signal may come from an application, acquisition adapter, phase extractor,
offline window generator, or completed external buffer.

### 3. Generic HDF5

`read_generic_hdf5` reads explicitly selected one-dimensional numeric
datasets. It does not infer arbitrary schemas.

```python
from multimodalphysiokit.io import read_generic_hdf5

recording = read_generic_hdf5(
    "recording.h5",
    signal_paths={
        "ecg": "/signals/ecg",
        "eda": "/signals/eda",
    },
    sampling_frequency_path="/metadata/sampling_frequency",
)
```

Its signature is:

```text
read_generic_hdf5(
    file_path,
    signal_paths,
    *,
    sampling_frequency=None,
    sampling_frequency_path=None,
    time_path=None,
    units=None,
    labels=None,
    name=None,
    metadata=None,
)
```

Provide exactly one timing source: a scalar `sampling_frequency`, a path to a
scalar sampling-frequency dataset, or a path to a one-dimensional time
dataset. Signal names and HDF5 dataset paths are explicit mappings.

### 4. BioSignalsPlux/OpenSignals

`read_biosignalsplux_hdf5` is a specialized convenience adapter for the
supported OpenSignals/BioSignalsPlux layout:

```python
from multimodalphysiokit.io import read_biosignalsplux_hdf5

recording = read_biosignalsplux_hdf5("recording.h5")
```

It recognizes supported sensor channels, converts raw ADC values to package
units, creates time vectors and `Signal` objects, records channel metadata, and
returns a `Recording`.

Its signature is:

```text
read_biosignalsplux_hdf5(
    file_path,
    *,
    fnirs_column_order=(0, 1),
    name="biosignalsplux",
)
```

The bundled notebooks use this adapter, but processors are not coupled to it.

## Direct processor usage

Recording is not required:

```python
from multimodalphysiokit.processors import TemperatureProcessor

processor = TemperatureProcessor()
result = processor.process(temperature_signal)
```

For a completed real-time PPG window (shorter than the 60 s needed by HRV
frequency analysis):

```python
from multimodalphysiokit.processors import PPGProcessor

processor = PPGProcessor(label_frequency_analysis=None)
result = processor.process(ppg_window)
result.features["ppg_bpm_mean"], result.metadata["n_valid_pulses"]
```

For paired inputs:

```python
from multimodalphysiokit.processors import FNIRSProcessor

processor = FNIRSProcessor(age_years=35)
result = processor.process(red_signal, infrared_signal)
```

The same processor object can be reused with another compatible signal:

```python
window_result = processor.process(red_window, infrared_window)
```

## Recording orchestration

### Inspect and manage signals

```python
recording.signal_names()
recording.has_signal("ecg")
ecg_signal = recording.get_signal("ecg")
recording.add_signal("temperature", temperature_signal)
```

Names are user-defined unless an input adapter assigns canonical names.

### Whole-signal processing

`process_signal` retrieves named signals and forwards them in order:

```python
result = recording.process_signal(
    ECGProcessor(label_frequency_analysis=1),
    signal_names="ecg",
)
```

For a paired processor:

```python
result = recording.process_signal(
    FNIRSProcessor(age_years=35),
    signal_names=("fnirs_red", "fnirs_infrared"),
    process_kwargs={
        "baseline_red_value": baseline_red_value,
        "baseline_infrared_value": baseline_infrared_value,
    },
)
```

`process_kwargs` is forwarded without modality-specific interpretation.

## Temporal phases

Phases are stored as seconds in a float matrix of shape `(2, n_phases)`.
Columns are phases, row 0 contains starts, and row 1 contains ends:

```python
import numpy as np

phase_intervals_seconds = np.array(
    [
        [0.0, 60.0],
        [60.0, 120.0],
    ],
    dtype=float,
)

recording.set_phases(
    phase_intervals_seconds,
    labels=["segment_1", "segment_2"],
)
```

This defines `segment_1 = [0, 60)` and `segment_2 = [60, 120)`. Labels are
unique, supplied order is preserved, and overlaps are allowed. Users provide
temporal intervals rather than sample-index matrices.

`phase_indices(signal_or_name, phases=None)` maps intervals through a chosen
signal's time vector. `extract_phase(signal_name, phase, reset_time=True)`
returns a new `Signal`. `clear_phases()` removes the phase definition.

### Plotting

```python
from multimodalphysiokit.utils import plot_signal_with_phases

plot_signal_with_phases(
    recording,
    signal_name="ecg",
    title="Raw ECG signal",
    show=True,
)
```

Phases must already be defined. The utility retrieves the complete stored
signal and obtains overlays through `Recording.phase_indices()`. It supports
optional saving without duplicating temporal mapping logic.

### Phase processing

```python
phase_results = recording.process_phases(
    ECGProcessor(label_frequency_analysis=1),
    signal_names="ecg",
)
```

The result is:

```python
{
    "segment_1": ProcessingResult(...),
    "segment_2": ProcessingResult(...),
}
```

Select and order phases with `phases=["segment_2", "segment_1"]`. Each phase is
mapped, extracted, and processed independently. Filtering, decomposition,
event detection, STFT, and feature extraction restart for each extracted
signal. This can produce different edge behavior from processing a complete
parent signal and slicing the processed output.

## Multimodal postprocessing

Separate plug-and-play processors can reuse one shared temporal definition:

```python
from multimodalphysiokit.processors import (
    ECGProcessor,
    RespirationProcessor,
    TemperatureProcessor,
)

results = {
    "ecg": recording.process_phases(
        ECGProcessor(label_frequency_analysis=1),
        signal_names="ecg",
    ),
    "respiration": recording.process_phases(
        RespirationProcessor(),
        signal_names="rip",
    ),
    "temperature": recording.process_phases(
        TemperatureProcessor(),
        signal_names="temp",
    ),
}
```

The `rip` and `temp` names above are assigned by the bundled input adapter;
other recordings may use different names. Every signal is mapped through its
own time vector. Processing is sequential and does not imply threading or
multiprocessing.

### Long-format phase results

```python
import pandas as pd

rows = [
    {
        "modality": modality,
        "phase": phase_label,
        "feature_name": feature_name,
        "value": feature_value,
    }
    for modality, phase_results in results.items()
    for phase_label, processing_result in phase_results.items()
    for feature_name, feature_value in processing_result.features.items()
]

results_table = pd.DataFrame(rows)
```

Phase association belongs to the outer result dictionaries. Feature names are
phase-independent. Do not flatten `intermediates` or `metadata` into the
scalar feature table. Pandas is an example dependency, not a core runtime
dependency.

## Externally managed sliding windows

Generate temporal windows outside the processing layer:

```python
from multimodalphysiokit.utils import create_sliding_phase_intervals

intervals, labels = create_sliding_phase_intervals(
    start_seconds=0.0,
    end_seconds=300.0,
    window_seconds=60.0,
    overlap_seconds=30.0,
)
recording.set_phases(intervals, labels)
```

The utility operates only in seconds, returns complete-window intervals and
labels, and does not process signals. Recording maps the intervals through
each selected signal.

## Window-based real-time integration

Acquisition and buffer management can produce a completed window:

```python
window_signal = Signal(
    samples=window_samples,
    time=window_time,
    sampling_frequency=1000.0,
    label="ecg",
    units="mV",
)

processor = ECGProcessor(label_frequency_analysis=2)
window_result = processor.process(window_signal)
```

The processor object can be reused for successive completed windows.
Acquisition, buffering, scheduling, and hard real-time guarantees remain
outside the core library. Current algorithms may use noncausal filtering or
complete-window operations; this is not sample-by-sample stateful streaming.

## Optional intermediate outputs

ECG exposes optional arrays when configured:

```python
result = ECGProcessor(
    label_frequency_analysis=1,
    return_intermediates=True,
).process(ecg_signal)
```

Keys are `filtered_signal`, `r_peak_indices`, `r_peak_times`,
`r_peak_amplitudes`, `accepted_r_peak_indices`,
`accepted_r_peak_times`, `ibi`, `ibi_times`, `bpm`, and `bpm_times`.
STFT mode also returns `stft_matrix`, `stft_power`, `stft_frequencies`, and
`stft_times`.

## fNIRS scalar baseline workflow

Compute filtered-current reference scalars once:

```python
processor = FNIRSProcessor(age_years=35)
baseline_red_value, baseline_infrared_value = (
    processor.compute_baseline_values(
        baseline_red_signal,
        baseline_infrared_signal,
    )
)
```

Then reuse them:

```python
results = recording.process_phases(
    processor,
    signal_names=("fnirs_red", "fnirs_infrared"),
    phases=["segment_2", "segment_3"],
    process_kwargs={
        "baseline_red_value": baseline_red_value,
        "baseline_infrared_value": baseline_infrared_value,
    },
)
```

Both values must be supplied together, finite, and strictly positive. Baseline
signals are not arguments to `process()`. Without supplied values, the
processor derives both scalars from the target pair after excluding ten
seconds from each temporal end. Result metadata includes `baseline_mode`,
`baseline_red_value`, and `baseline_infrared_value`.

## Modality guides

- [ECG](https://github.com/Gabbert97/MultimodalPhysioKit/blob/main/docs/modalities/ecg.md)
- [EDA](https://github.com/Gabbert97/MultimodalPhysioKit/blob/main/docs/modalities/eda.md)
- [Respiration](https://github.com/Gabbert97/MultimodalPhysioKit/blob/main/docs/modalities/respiration.md)
- [Temperature](https://github.com/Gabbert97/MultimodalPhysioKit/blob/main/docs/modalities/temperature.md)
- [PPG](https://github.com/Gabbert97/MultimodalPhysioKit/blob/main/docs/modalities/ppg.md)
- [fNIRS](https://github.com/Gabbert97/MultimodalPhysioKit/blob/main/docs/modalities/fnirs.md)

## Example notebooks

The notebooks use one bundled synthetic BioSignalsPlux/OpenSignals fixture as
a practical example source, not as an architectural requirement:

- [ECG](https://github.com/Gabbert97/MultimodalPhysioKit/blob/main/examples/python/ecg_example.ipynb)
- [EDA](https://github.com/Gabbert97/MultimodalPhysioKit/blob/main/examples/python/eda_example.ipynb)
- [Respiration](https://github.com/Gabbert97/MultimodalPhysioKit/blob/main/examples/python/respiration_example.ipynb)
- [Temperature](https://github.com/Gabbert97/MultimodalPhysioKit/blob/main/examples/python/temperature_example.ipynb)
- [Multimodal](https://github.com/Gabbert97/MultimodalPhysioKit/blob/main/examples/python/multimodal_example.ipynb)

They demonstrate the specialized reader, signal retrieval, raw plotting,
direct and phase processing, intermediate outputs, and multimodal result
tables.

## Testing and validation

Run:

```bash
cd python
python -m pytest -v
```

Unit tests, bundled demonstrations, private/local validation, and
MATLAB/Python scientific comparisons are distinct activities. The bundled
fixture contains mathematically generated demonstration data, not
scientifically validated physiology. See
[Validation](https://github.com/Gabbert97/MultimodalPhysioKit/blob/main/docs/validation.md) and
[Local validation](https://github.com/Gabbert97/MultimodalPhysioKit/blob/main/docs/local_validation.md).

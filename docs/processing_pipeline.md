# Processing pipeline

## Python workflows

Python supports direct processing and optional offline orchestration:

```text
Direct:
Signal → Processor → ProcessingResult

Whole signal:
Recording → stored Signal(s) → processor.process(...)

Phases:
Recording → extract independent Signal segment(s)
          → processor.process(...)
          → dict[phase_label, ProcessingResult]
```

## 1. Load or construct signals

Construct or receive `Signal` objects directly when an application already
manages samples and time:

```python
from multimodalphysiokit.processors import ECGProcessor

result = ECGProcessor(label_frequency_analysis=1).process(ecg_signal)
```

For in-memory arrays use `recording_from_arrays`. For HDF5, use
`read_generic_hdf5` with explicit mappings or the specialized
`read_biosignalsplux_hdf5` adapter for its supported acquisition layout. All
paths converge on the same `Signal`/processor interface.

## 2. Define temporal phases when needed

```python
import numpy as np

recording.set_phases(
    np.array([[0.0, 60.0], [60.0, 120.0]]),
    labels=["segment_1", "segment_2"],
)
```

Intervals are seconds with half-open semantics. Recording maps them through
each signal's own time vector.

## 3. Run a modality processor

Direct use:

```python
from multimodalphysiokit.processors import TemperatureProcessor

signal = recording.get_signal("temp")
result = TemperatureProcessor().process(signal)
```

Whole-signal convenience:

```python
result = recording.process_signal(
    TemperatureProcessor(),
    signal_names="temp",
)
```

Phase processing:

```python
phase_results = recording.process_phases(
    TemperatureProcessor(),
    signal_names="temp",
)
```

Every extracted phase is processed independently. Consequently, filtering,
decomposition, event detection, spectral analysis, and feature extraction
restart for each call. Processing a segment independently can differ at its
boundaries from preprocessing a parent signal and slicing afterward.

## 4. Modality processing paths

| Modality | Main path |
|---|---|
| ECG | filtering, R-peak detection, inter-beat intervals, HRV and frequency analysis |
| EDA | downsampling, filtering, cvxEDA decomposition, SCR processing |
| Respiration | downsampling, filtering, extrema and breath pairing |
| Temperature | low-pass filtering, artifact screening, derivative features |
| PPG | band-pass beat detection, pulse feet/peaks and quality check, pulse-shape, HRV and frequency analysis |
| fNIRS | paired-current downsampling/filtering, baseline scalars, concentration conversion |

In MATLAB, each processor is an object-facing adapter over its modality's
validated procedural implementation. `EDAProcessor` delegates one complete
Signal window to `matlab/eda/eda_signal_processing.m`; phase processing first
extracts each Recording interval and then runs that window independently.

Each processor validates the supplied signal or synchronized signal pair.
A successful call does not replace acquisition quality control.

## 5. Tabulate scalar phase features

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

Only final scalar `features` belong in this table. Preserve `metadata`
separately and do not flatten intermediate arrays. Pandas is used by the
notebook example, not by core processing.

## MATLAB workflow

MATLAB Recording workflows use temporal intervals in seconds with half-open
semantics. Validated procedural functions may use one-based internal indexes
after Recording extracts the requested Signal window. See the
[MATLAB guide](../matlab/README.md).

Detailed modality pages:

- [ECG](modalities/ecg.md)
- [EDA](modalities/eda.md)
- [Respiration](modalities/respiration.md)
- [Temperature](modalities/temperature.md)
- [PPG](modalities/ppg.md)
- [fNIRS](modalities/fnirs.md)

# MultimodalPhysioKit

[![Python CI](https://github.com/Gabbert97/MultimodalPhysioKit/actions/workflows/python-ci.yml/badge.svg?branch=main)](https://github.com/Gabbert97/MultimodalPhysioKit/actions/workflows/python-ci.yml)
[![MATLAB CI](https://github.com/Gabbert97/MultimodalPhysioKit/actions/workflows/matlab-ci.yml/badge.svg?branch=main)](https://github.com/Gabbert97/MultimodalPhysioKit/actions/workflows/matlab-ci.yml)
[![DOI](https://zenodo.org/badge/1330070732.svg)](https://doi.org/10.5281/zenodo.21877145)

MultimodalPhysioKit is a modular API for processing synchronized multimodal
physiological signals in offline postprocessing and window-based real-time
pipelines.

The library provides plug-and-play processing objects for ECG,
electrodermal activity, respiration, skin temperature, and fNIRS. Each
processor operates on standardized `Signal` objects and is independent of the
original data source, allowing the same processing pipeline to be used with
in-memory arrays, generic HDF5 datasets, BioSignalsPlux/OpenSignals recordings,
extracted experimental phases, or externally generated streaming windows.

## Key capabilities

- Plug-and-play physiological processors with standardized
  `ProcessingResult` outputs
- Direct processing of `Signal` objects without a `Recording`
- Optional synchronized multimodal `Recording` containers
- Modality-specific sampling frequencies and sample counts
- Shared temporal phases in seconds for offline postprocessing
- Externally generated offline or completed real-time analysis windows
- Input adapters for arrays, explicit generic HDF5 mappings, and supported
  BioSignalsPlux/OpenSignals HDF5 exports
- Python and MATLAB implementations for the supported modalities

The library consumes synchronized inputs; it does not automatically
synchronize independently acquired streams.

## Core architecture

```text
In-memory arrays
Existing Signal objects
Generic HDF5 datasets
BioSignalsPlux/OpenSignals recordings
Externally generated streaming windows
                     ↓
               Signal objects
                     ↓
       Plug-and-play modality processors
                     ↓
            ProcessingResult objects
```

`Signal` standardizes samples, time, sampling frequency, label, units, and
metadata. It may represent a complete channel, one experimental phase, an
offline window, or one completed real-time buffer.

A modality processor is the central reusable component. It accepts one or more
`Signal` objects, applies modality-specific preprocessing and feature
extraction, optionally exposes intermediate outputs, and returns a
`ProcessingResult`.

`Recording` is optional. It stores named synchronized signals while preserving
their sampling frequencies, defines shared temporal phases, maps those phases
through each signal's time vector, and orchestrates whole-signal or phase-based
offline processing.

`ProcessingResult` standardizes scalar features, metadata, options, signal
labels, and optional intermediate outputs.

## Plug-and-play processors

The same processor API applies regardless of input source:

```python
from multimodalphysiokit.processors import ECGProcessor

ecg_signal = recording.get_signal("ecg")
processor = ECGProcessor(
    label_frequency_analysis=1,
    return_intermediates=True,
)
result = processor.process(ecg_signal)

# The same reusable processor can analyze another completed Signal window.
window_result = processor.process(ecg_window)
```

Processors receive exactly the supplied signal or synchronized signal pair.
They do not receive a `Recording` or phase definitions.

## Supported input workflows

### Existing arrays

`recording_from_arrays` creates signals sharing one sampling frequency or one
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

For modalities with different sampling frequencies, construct each `Signal`
with its own time vector and add it to a `Recording`.

### Existing Signal objects

Users can create or receive a `Signal` from any acquisition pipeline and call
a processor directly:

```python
result = ECGProcessor(label_frequency_analysis=1).process(ecg_signal)
```

This is the most general processing interface.

### Generic HDF5

The generic reader uses explicit one-dimensional dataset mappings and exactly
one timing source:

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

It does not infer arbitrary HDF5 schemas.

### BioSignalsPlux/OpenSignals

The specialized convenience adapter interprets the supported
BioSignalsPlux/OpenSignals layout and converts recognized channels:

```python
from multimodalphysiokit.io import read_biosignalsplux_hdf5

recording = read_biosignalsplux_hdf5("recording.h5")
```

The bundled notebooks use this adapter as one practical example input source.

## Direct processing

Recording is not required:

```python
from multimodalphysiokit.processors import FNIRSProcessor, TemperatureProcessor

temperature_result = TemperatureProcessor().process(temperature_signal)
fnirs_result = FNIRSProcessor(age_years=35).process(
    red_signal,
    infrared_signal,
)
```

## Multimodal and phase-based postprocessing

`Recording` provides optional convenience methods:

```python
result = recording.process_signal(
    processor,
    signal_names="ecg",
)

phase_results = recording.process_phases(
    processor,
    signal_names="ecg",
)
```

Each temporal phase is mapped, extracted, and processed independently.
Filtering a phase independently can produce different boundary behavior from
filtering its complete parent signal and slicing afterward.

Separate processors can reuse one temporal phase definition:

```python
from multimodalphysiokit.processors import (
    ECGProcessor,
    RespirationProcessor,
    TemperatureProcessor,
)

results = {
    "ecg": recording.process_phases(
        ECGProcessor(label_frequency_analysis=1),
        "ecg",
    ),
    "respiration": recording.process_phases(
        RespirationProcessor(),
        "rip",
    ),
    "temperature": recording.process_phases(
        TemperatureProcessor(),
        "temp",
    ),
}
```

Here `rip` and `temp` are canonical names assigned by the bundled
BioSignalsPlux adapter. Signal names are otherwise user-defined.

## Window-based real-time integration

Completed external buffers can be wrapped as `Signal` objects:

```python
from multimodalphysiokit import Signal

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

Acquisition, buffering, window scheduling, and hard real-time guarantees remain
outside the current core library. Processors analyze completed windows and may
use noncausal filtering or complete-window operations; they are not presented
as sample-by-sample causal streaming implementations.

## Available modalities

- [ECG](docs/modalities/ecg.md)
- [EDA](docs/modalities/eda.md)
- [Respiration](docs/modalities/respiration.md)
- [Temperature](docs/modalities/temperature.md)
- [fNIRS](docs/modalities/fnirs.md)

## Examples

The [Python notebooks and MATLAB scripts](examples/README.md) demonstrate the
API using one bundled synthetic BioSignalsPlux/OpenSignals fixture. This is a
practical example input source, not a limitation of the processing
architecture. The MATLAB files are direct, editable workflows for loading
data, defining temporal phases, constructing processors, and inspecting
native results.

## Installation and dependencies

MultimodalPhysioKit v0.1.0 is not yet published to PyPI. Once available, the
planned public installation command is:

```bash
python -m pip install multimodalphysiokit
```

For development from a repository checkout, run from the repository root:

```bash
python -m pip install -e "python[test]"
```

To include dependencies used by the example notebooks:

```bash
python -m pip install -e "python[test,examples]"
```

Python 3.10 or newer is required. EDA uses external cvxEDA. See
[Dependencies](docs/dependencies.md), the detailed
[Python guide](python/README.md), and the separate
[MATLAB guide](matlab/README.md).

## Validation and scientific background

Run Python tests from `python/` with `python -m pytest -v`. See
[Validation](docs/validation.md), [local validation](docs/local_validation.md),
and [CITATION.cff](CITATION.cff). Scientific background for the supported
modalities includes the publications linked by the citation metadata and the
cvxEDA reference in the
[third-party notices](docs/legal/THIRD_PARTY_NOTICES.md). The package is a
research toolkit, not a clinical device; users remain responsible for
acquisition quality, units, processing assumptions, and study-specific
validation.

MultimodalPhysioKit is licensed under LGPL-3.0-or-later. See
[LICENSE](LICENSE), [NOTICE](docs/legal/NOTICE.md), and
[third-party notices](docs/legal/THIRD_PARTY_NOTICES.md).

Project maintenance information is available in the [changelog](CHANGELOG.md),
[contribution guide](CONTRIBUTING.md), and [security policy](SECURITY.md).

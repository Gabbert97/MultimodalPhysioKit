# MultimodalPhysioKit examples

The Python examples are self-contained Jupyter notebooks using the bundled
synthetic OpenSignals/BioSignalsPlux fixture. The mathematically generated
data demonstrate the API and are not a scientific validation dataset.

This is one practical input source, not a limitation of the processing
architecture. The same processors accept compatible `Signal` objects created
from arrays, other adapters, extracted phases, or externally managed windows.

## Bundled synthetic fixture

The file is:

```text
sample_data/synthetic_multimodal_180s.hdf5
```

Load it through the public reader:

```python
from multimodalphysiokit.io import read_biosignalsplux_hdf5

recording = read_biosignalsplux_hdf5(
    "examples/sample_data/synthetic_multimodal_180s.hdf5"
)
recording.signal_names()
```

The returned names are:

```text
ecg, eda, rip, temp, fnirs_red, fnirs_infrared
```

The reader recognizes these supported channels, converts their raw ADC values
to package units, creates time vectors and `Signal` objects, retains channel
metadata, and returns a ready-to-use `Recording`. The notebooks do not perform
manual HDF5 extraction or conversion.

## Python notebooks

From the repository root:

```bash
python -m pip install -e "python[examples]"
jupyter lab
```

Open:

- [ECG example](python/ecg_example.ipynb): loads `ecg`, plots it, runs
  `ECGProcessor`, and inspects optional ECG/STFT intermediates.
- [EDA example](python/eda_example.ipynb): loads `eda`, plots it, and runs
  `EDAProcessor`; external cvxEDA is required.
- [Respiration example](python/respiration_example.ipynb): loads `rip`, plots
  it, and runs `RespirationProcessor`.
- [Temperature example](python/temperature_example.ipynb): loads `temp`, plots
  it, and demonstrates direct and phase processing.
- [Multimodal example](python/multimodal_example.ipynb): reuses one temporal
  phase definition across all bundled modalities, processes each modality,
  and displays every scalar feature in a modality/phase table.

Every notebook independently loads the file, inspects signal names, defines
neutral temporal phases, plots raw signals with
`plot_signal_with_phases`, retrieves `Signal` objects from the `Recording`,
and calls public processor APIs.

## MATLAB scripts

The MATLAB examples are short, editable scripts:

- [Temperature](matlab/temperature_example.m)
- [Respiration](matlab/respiration_example.m)
- [EDA](matlab/eda_example.m)
- [ECG](matlab/ecg_example.m)
- [fNIRS](matlab/fnirs_example.m)
- [Multimodal](matlab/multimodal_example.m)

Each script shows the direct workflow: load the bundled synthetic fixture, define
temporal phases, construct a processor, process independently extracted
Signals, and inspect the native result. Customize a workflow by editing
`dataFile`, `phaseIntervals`, `phaseLabels`, or the processor configuration.
The EDA script requires an external cvxEDA installation and its optimization
dependency.

## Advanced cross-language validation

The [MATLAB ↔ Python validation notebook](cross_language_validation/cross_language_validation.ipynb)
is an advanced example comparing final feature extraction in both languages
from the same tracked synthetic recording. A small MATLAB script exports the
reference features; the notebook runs Python and presents the comparison.

## Temporal phases

The bundled synthetic file has no protocol annotations. The notebooks
therefore use neutral illustrative intervals in seconds:

```python
import numpy as np

phase_intervals_seconds = np.array(
    [[0.0, 90.0], [90.0, 180.0]],
    dtype=float,
)
recording.set_phases(
    phase_intervals_seconds,
    labels=["segment_1", "segment_2"],
)
```

Columns represent `[0, 90)` and `[90, 180)`. Recording maps these intervals
through each stored signal's time vector.

## Direct and offline processing

Processors can be used without Recording orchestration:

```python
from multimodalphysiokit.processors import TemperatureProcessor

temperature_signal = recording.get_signal("temp")
result = TemperatureProcessor().process(temperature_signal)
```

For offline phases:

```python
phase_results = recording.process_phases(
    TemperatureProcessor(),
    signal_names="temp",
)
```

Each dictionary value is an independently produced `ProcessingResult`.
Feature names do not include phase labels.

# Data model and temporal phases

## Python Signal

A `Signal` represents one sampled channel:

```python
Signal(
    samples,
    time,
    sampling_frequency,
    label="",
    units="",
    metadata=None,
)
```

Its public fields are `samples`, `time`, `sampling_frequency`, `label`,
`units`, and `metadata`. Samples and time have equal length, time increases
strictly, and sampling frequency is positive.

A signal may contain a complete channel, a phase extracted by `Recording`, an
external analysis window, or a completed real-time buffer window. Processors
do not need to know its origin. Current processing can include noncausal
filtering and complete-window calculations.

## Python Recording

`Recording(name="recording", signals=None, metadata=None,
phase_intervals_seconds=None, phase_labels=None)` is an optional offline
container for named signals and shared temporal phase definitions.

```python
recording.add_signal("ecg", ecg_signal)
ecg_signal = recording.get_signal("ecg")
recording.has_signal("ecg")
recording.signal_names()
```

Input adapters may return an already populated `Recording`; users may also
assemble one from existing `Signal` objects.

Recording does not synchronize independently acquired signals, resample them,
interpolate missing values, infer experimental phases, or apply
modality-specific processing.

## Python temporal phases

Python stores phase intervals once on `Recording`:

```text
recording.phase_intervals_seconds.shape == (2, n_phases)
```

Row 0 contains start times and row 1 contains end times. Columns are phases,
using half-open semantics `[start_seconds, end_seconds)`.

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

This means:

```text
segment_1 = [0, 60)
segment_2 = [60, 120)
```

Labels are nonempty and unique. Order is preserved and intervals may overlap.
Users provide seconds, not sample-index matrices. `phase_indices()` maps
selected intervals with `numpy.searchsorted` through the chosen
`Signal.time`, producing zero-based inclusive-start/exclusive-stop indexes.
Signals may have different sampling frequencies and sample counts.

`extract_phase(signal_name, phase, reset_time=True)` uses that mapping and
returns a new `Signal`. `clear_phases()` removes all temporal definitions.

## Python ProcessingResult

`ProcessingResult` contains:

- `modality`
- `processor_name`
- `features`
- `signal_labels`
- `options`
- `metadata`
- `intermediates`

`features` contains final scalar outputs. `metadata` contains concise
descriptive and scalar processing information. `intermediates` can hold
optional arrays or structured diagnostics without mixing them into features.

For phase processing, phase labels belong to the outer dictionary:

```python
{
    "segment_1": ProcessingResult(...),
    "segment_2": ProcessingResult(...),
}
```

Feature names remain phase-independent.

## Metadata boundaries

Keep metadata concise. Do not store credentials, identifying participant
information, private recording paths, or large signal arrays in metadata.

## MATLAB data model

MATLAB `mpk.Recording` stores the same 2-by-N temporal intervals in seconds
with half-open `[start,end)` semantics. It maps each interval independently
through each `mpk.Signal.Time` vector. Validated procedural MATLAB functions
may still use one-based internal indexes after a Signal has been selected.

See [Processing pipeline](processing_pipeline.md) for complete workflows.

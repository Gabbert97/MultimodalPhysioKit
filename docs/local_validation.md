# Local validation

Local validation is an optional developer workflow for recordings that cannot
be distributed with the public repository.

## Public/private boundary

- Private recordings, paths, identifiers, caches, plots, and comparison tables
  must remain outside version control.
- Public tests must run without local recordings.
- Local validation is not an implicit package dependency.
- Scientific changes still require focused public tests using distributable
  controlled inputs.

## Supported BioSignalsPlux inputs

For a local file in the supported OpenSignals/BioSignalsPlux HDF5 layout, use
the public reader with a path supplied by local configuration:

```python
from multimodalphysiokit.io import read_biosignalsplux_hdf5

recording = read_biosignalsplux_hdf5(configured_local_path)
recording.signal_names()
```

Do not publish the configured path. For another HDF5 layout, use
`read_generic_hdf5` with explicit mappings only when its documented contract
fits the file.

## Temporal phases and processing

Provide phase intervals in seconds:

```python
recording.set_phases(intervals_seconds, labels)

whole_result = processor.process(recording.get_signal(signal_name))
phase_results = recording.process_phases(
    processor,
    signal_names=signal_name,
)
```

Each phase is processed independently. Do not compare this workflow with
full-signal preprocessing followed by slicing as though boundary behavior must
be identical.

## Cross-language comparisons

Convert Python half-open temporal intervals through the relevant signal time
vector before comparing with MATLAB one-based inclusive indexes. Align units,
channel order, sampling frequencies, baseline scalars, processor options,
feature order, and tolerances.

Only publish aggregate validation summaries after removing non-public
acquisition details and confirming they may be shared.

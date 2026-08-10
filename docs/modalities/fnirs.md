# fNIRS

## Python API

```python
from multimodalphysiokit.processors import FNIRSProcessor

processor = FNIRSProcessor(age_years=35)
result = processor.process(red_signal, infrared_signal)
```

`FNIRSProcessor(*, age_years)` requires a positive finite age.
`process(red_signal, infrared_signal, *, baseline_red_value=None,
baseline_infrared_value=None, show_plots=False, save_plots=False,
output_dir=None, output_label="fnirs")` analyzes one synchronized pair.

The signals must have equal lengths, sampling frequencies, and synchronized
time vectors. Currents must be finite and strictly positive, and the sampling
rate must support downsampling to 10 Hz.

The pair may come from any compatible source or completed analysis window.
The bundled synthetic BioSignalsPlux example provides canonical
`fnirs_red` and `fnirs_infrared` signals.

## Scalar baseline workflow

Compute the validated filtered-current reference scalars once:

```python
baseline_red_value, baseline_infrared_value = (
    processor.compute_baseline_values(
        baseline_red_signal,
        baseline_infrared_signal,
    )
)
```

The helper downsamples the pair to 10 Hz, applies the canonical Butterworth
current low-pass filter (0.2 Hz passband edge and 0.4 Hz stopband edge) with
explicit, deterministic zero-phase edge handling, excludes exactly ten seconds
from each end, and returns two Python floats.

Reuse those scalars for target phases:

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

Both baseline values must be supplied together. Each must be a numeric,
finite, strictly positive scalar. Do not pass baseline signals to `process()`.

When both values are omitted, the processor self-calculates them from the
target pair using the same ten-second boundary exclusion; that target must
retain baseline samples after excluding both ends.

Metadata includes:

- `baseline_mode` (`"self"` or `"provided"`)
- `baseline_red_value`
- `baseline_infrared_value`
- `age_years`
- original/processed sampling frequencies and sample counts
- start/end times

## Processing and features

The processor downsamples and low-pass filters the target currents, applies
the selected scalar baseline to optical-density and modified Beer–Lambert
calculations, high-pass filters oxy/deoxy concentrations, and extracts 19
features for each channel.

Keys have `fnirs_oxy_` or `fnirs_deoxy_` followed by:

`min`, `max`, `mean`, `variance`, `std`, `skewness`, `kurtosis`,
`power_band`, `max_frequency`, `median_frequency`, `spectral_entropy`,
`mean_diff`, `peak_to_peak`, `slope`, `zero_crossing`, `polarity`, `entropy`,
`auc`, or `total_energy`.

Skewness and Pearson kurtosis use explicitly biased finite-sample estimators,
preserving the historical MATLAB feature convention. Peak-to-peak smoothing
uses the same explicit zero-phase SOS edge policy as the other fNIRS filters.

This produces 38 phase-independent scalar feature keys.

## Methodological references

MultimodalPhysioKit independently implements the fNIRS processing workflow
described by Luzzani et al., “ECG, Respiration, fNIRS, and Eye Tracking for
Stress and Mental Workload Monitoring in Human-Machine Interaction,” *IEEE
Access*, 2025, DOI: `10.1109/ACCESS.2025.3588384`. This publication is a
scientific reference, not a source-code origin.

## MATLAB

MATLAB uses `mpk.FNIRSProcessor` or `fNIRS_signal_processing`. Align
wavelength ordering, filtered-current baseline values, downsampling,
Beer–Lambert parameters, filtering, and feature order for comparisons.

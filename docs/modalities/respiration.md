# Respiration

## Python API

```python
from multimodalphysiokit.processors import RespirationProcessor

result = RespirationProcessor().process(respiration_signal)
```

The accepted peak-to-peak breathing-rate range is configurable:

```python
processor = RespirationProcessor(
    min_breathing_rate=3.0,
    max_breathing_rate=30.0,
)
result = processor.process(respiration_signal)
```

`min_breathing_rate` and `max_breathing_rate` are expressed in breaths per
minute and define an inclusive acceptance range. They must be finite,
`min_breathing_rate` must be greater than zero, and
`max_breathing_rate` must be greater than `min_breathing_rate`. The defaults
remain 5.0 and 30.0 breaths/minute, preserving the behavior of version 0.1.0.

`process(signal, *, show_plots=False, save_plots=False, output_dir=None)`
analyzes exactly one supplied signal.

The input can be any compatible respiration `Signal`, regardless of its
storage or acquisition source. For a concrete adapter example, the bundled
synthetic BioSignalsPlux fixture exposes its RIP channel under the name `rip`.

## Processing

Downsampling to 100 Hz, filtering, candidate and definitive extrema detection,
breath pairing, respiratory quantities, and feature extraction occur inside
the processor. Do not manually downsample or convert temporal phases to sample
indexes.

When used with `Recording.process_phases()`, every phase is downsampled and
processed independently.

## Features

The 14 keys are:

- `respiration_breath_rate_mean`
- `respiration_breath_rate_std`
- `respiration_inspiratory_time_mean`
- `respiration_inspiratory_time_std`
- `respiration_expiratory_time_mean`
- `respiration_expiratory_time_std`
- `respiration_timing_ratio_mean`
- `respiration_timing_ratio_std`
- `respiration_amplitude_mean`
- `respiration_amplitude_std`
- `respiration_minute_ventilation_mean`
- `respiration_minute_ventilation_std`
- `respiration_coefficient_of_variation`
- `respiration_autocorrelation`

Metadata includes original/processed sampling frequencies and sample counts,
plus signal start and end times.

## Assumptions

Detected extrema must correspond to physiological breath cycles. Motion,
belt shifts, apnea, sighs, and irregular breathing can affect pairing.
Relative amplitude is not absolute ventilation without calibration.

The default 5-30 breaths-per-minute acceptance range is a screening choice.
Users should review it for their population and protocol and may configure the
range through `RespirationProcessor` without changing filtering, extrema
detection, morphology validation, or feature formulas.

## Methodological references

MultimodalPhysioKit independently implements the respiration processing
workflow described by Luzzani et al., “ECG, Respiration, fNIRS, and Eye
Tracking for Stress and Mental Workload Monitoring in Human-Machine
Interaction,” *IEEE Access*, 2025, DOI: `10.1109/ACCESS.2025.3588384`. This
publication is a scientific reference, not a source-code origin.

## MATLAB

MATLAB uses `mpk.RespirationProcessor` or
`respiration_signal_processing`. Preserve the validated sampling,
filtering, extrema, pairing, and feature conventions during comparisons.

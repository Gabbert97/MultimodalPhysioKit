# Temperature

## Direct processing

```python
from multimodalphysiokit.processors import TemperatureProcessor

result = TemperatureProcessor().process(temperature_signal)
```

`TemperatureProcessor()` has no constructor parameters.
`process(signal, *, show_plots=False, save_plots=False, output_dir=None)`
filters and extracts features from exactly the supplied signal.

The signal may represent a complete temperature channel or a completed
offline/real-time window from any compatible source. The bundled
BioSignalsPlux adapter uses the concrete signal name `temp`.

## Phase processing

```python
phase_results = recording.process_phases(
    TemperatureProcessor(),
    signal_names="temp",
)
```

Every extracted phase is filtered and screened independently. This can differ
at segment edges from filtering a complete signal and slicing afterward.

## Features

The eight validated keys are:

- `temperature_delta`
- `temperature_mean`
- `temperature_std`
- `temperature_delta_over_time`
- `temperature_linear_slope`
- `temperature_derivative_mean`
- `temperature_derivative_std`
- `temperature_derivative_linear_slope`

No feature key contains a phase suffix.

Result metadata contains `source`, `sampling_frequency`, `n_samples`,
`start_time_seconds`, and `end_time_seconds`.

## Processing and assumptions

The pipeline applies the validated low-pass filter, screens short windows for
abrupt artifacts, calculates level/change/trend, and summarizes the first
derivative.

Artifact screening uses 2-second windows. A window is removed when both
(a) a local extremum in it deviates from the phase mean by more than 2.0
degrees, and (b) its mean absolute sample-to-sample change exceeds 100 times
the phase's median absolute change. The thresholds are absolute temperature
differences, so results are identical for degrees Celsius and kelvin;
convert degrees Fahrenheit to Celsius before processing. The median baseline
keeps a large artifact from masking itself, and the absolute changes make
screening symmetric for warming and cooling. Peripheral temperature depends on placement, contact,
acclimatization, environment, calibration, and phase duration.

## Methodological references

MultimodalPhysioKit independently implements the skin-temperature processing
workflow described by Luzzani et al., “EDA, PPG and Skin Temperature as
Predictive Signals for Mental Failure by a Statistical Analysis on Stress and
Mental Workload,” *IEEE Open Journal of Engineering in Medicine and Biology*,
2025, DOI: `10.1109/OJEMB.2024.3515473`. This publication is a scientific
reference, not a source-code origin.

## MATLAB

MATLAB uses `mpk.TemperatureProcessor` or
`temperature_signal_processing`. The validated eight-feature subset,
filtering, artifact screening, derivative calculation, and standard-deviation
convention must be aligned for numerical comparison.

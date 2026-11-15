# EDA

## Python API

```python
from multimodalphysiokit.processors import EDAProcessor

processor = EDAProcessor()
result = processor.process(eda_signal)
```

The constructor is:

```text
EDAProcessor(
    *,
    cvxeda_path=None,
    decomposition_callable=None,
    cache_directory=None,
    output_label="eda",
    minimum_scr_amplitude=0.01,
)
```

`process(signal, *, show_plots=False, save_plots=False, output_dir=None,
output_label=None)` analyzes exactly one supplied signal.

## Processing

The supplied signal is downsampled internally to 100 Hz, low-pass filtered,
decomposed by external cvxEDA, processed for SCR events, and summarized into
features. Input sampling frequency must meet the processor's integer
downsampling requirement.

Before downsampling, a zero-phase 4th-order Butterworth low-pass at 40 Hz
(0.4 x the 100 Hz target rate) removes content above the new Nyquist
frequency so that it cannot alias into the analysis band; every k-th
sample is then kept, so sample timing is unchanged. Inputs already at
100 Hz are not filtered by this step. A brief zero-phase boundary
transient (tens of milliseconds) can remain at each signal edge.

MultimodalPhysioKit pairs each local maximum of the 2 Hz low-pass-filtered
phasic signal with its immediately preceding local minimum. SCR amplitude is
the resulting trough-to-peak change in microSiemens, not the absolute phasic
value at the peak. A response is accepted when its amplitude is greater than
or equal to `minimum_scr_amplitude`, whose default is 0.01 µS. Rise time is
reported as `(peak_index - minimum_index) / sampling_frequency` but is not an
acceptance criterion. The former relative 20th-percentile amplitude and rise-
time screening is not part of the canonical pipeline.

If no response meets the threshold, all four SCR amplitude/rise-time summary
features and `eda_scr_per_min` are returned as zero. With one accepted event,
the amplitude and rise-time standard deviations are zero.

`Recording.process_phases()` calls this complete pipeline independently for
every extracted phase. cvxEDA therefore runs once per phase rather than once
on a parent recording.

## Features

The eight keys are:

- `eda_scl_mean`
- `eda_scl_std`
- `eda_scl_slope`
- `eda_scr_amplitude_mean`
- `eda_scr_amplitude_std`
- `eda_scr_rise_time_mean`
- `eda_scr_rise_time_std`
- `eda_scr_per_min`

Result metadata contains `source`, original/processed sampling frequencies and
sample counts, start/end times, and `cvxeda_source`. Large decomposition arrays
are not placed in metadata or scalar features.

## cvxEDA

cvxEDA is external and not bundled. Official, unmodified cvxEDA v1.1.0 is the
validated reference version in both languages, and MultimodalPhysioKit does
not clip or otherwise modify its decomposition outputs. There is no production
fallback. Install v1.1.0, provide `cvxeda_path`, or provide a compatible
`decomposition_callable`. Optional caching is keyed by the complete supplied
signal and processing dependency. See [Dependencies](../dependencies.md).

The cache key includes the decomposition source, plus the cvxEDA parameters:
in Python the module and qualified name of a `decomposition_callable`, or a
hash of the cvxEDA source file for an installed module or `cvxeda_path`; in
MATLAB `func2str` of a `DecompositionCallable`, or a hash of the resolved
`cvxEDA.m`. A cached result is therefore never reused by a different
decomposition. In Python, the decomposition source must be resolvable even
when a cached result exists.

Record cvxEDA implementation and solver versions for reproducibility.

The cvxEDA decomposition is external. MultimodalPhysioKit independently
implements the surrounding preprocessing, dependency integration, validation
of cvxEDA outputs, SCR processing and screening, feature extraction and
aggregation, metadata and result construction, plotting, and orchestration.

## Methodological references

MultimodalPhysioKit independently implements its EDA processing workflow with
reference to Luzzani et al., “EDA, PPG and Skin Temperature as Predictive
Signals for Mental Failure by a Statistical Analysis on Stress and Mental
Workload,” *IEEE Open Journal of Engineering in Medicine and Biology*, 2025,
DOI: `10.1109/OJEMB.2024.3515473`.

The external cvxEDA decomposition is described by Greco et al., “cvxEDA: A
Convex Optimization Approach to Electrodermal Activity Processing,” *IEEE
Transactions on Biomedical Engineering*, 2016, DOI:
`10.1109/TBME.2015.2474131`. cvxEDA is separately obtained software, not
MultimodalPhysioKit source code.

## MATLAB

`mpk.EDAProcessor` is the reusable object-facing adapter. It validates the
supplied `mpk.Signal`, delegates the scientific calculation to the canonical
procedural implementation under `matlab/eda`, and packages the result. Each
call processes one complete Signal; `Recording.processPhases` extracts each
phase and invokes the same pipeline independently.

MATLAB exposes the equivalent constructor option `MinimumScrAmplitude`, in µS.
It requires the official v1.1.0 external `cvxEDA.m`, a configured `CvxedaPath`,
or an injected decomposition callable. Decomposition caching is disabled unless
`CacheDirectory` is explicitly configured, and plotting is disabled unless a
processing call explicitly requests it.

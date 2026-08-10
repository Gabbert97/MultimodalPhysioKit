# MultimodalPhysioKit for MATLAB

The MATLAB implementation exposes the `mpk` package, input adapters, plotting
utilities, and processors for ECG, EDA, respiration, temperature, and fNIRS.

## Setup

From the repository root, add the MATLAB tree to the path:

```matlab
addpath(genpath('matlab'))
```

Signal Processing Toolbox and Statistics and Machine Learning Toolbox are used
by the processing routines. EDA additionally requires external cvxEDA and the
Optimization Toolbox solver used by it. See
[Dependencies](../docs/dependencies.md).

## Basic workflow

Create a recording from arrays, retrieve a signal, define phases, and call a
processor:

```matlab
fs = 100;
time = (0:(120*fs-1))' / fs;
signals = struct('respiration',sin(2*pi*0.25*time));

recording = mpk.io.recordingFromArrays( ...
    signals,'SamplingFrequency',fs,'Units',struct('respiration','a.u.'));
signal = recording.getSignal('respiration');
intervalsSeconds = [0 60;60 120];
recording.setPhases( ...
    intervalsSeconds, ...
    {'baseline','task'});

processor = mpk.RespirationProcessor();
processor.ShowPlots = false;
results = recording.processPhases(processor,'respiration');
```

Recording intervals are expressed in seconds and use half-open
`[start,end)` semantics. Procedural scientific functions may internally use
one-based indexes after Recording extracts a Signal.

## Core concepts

- `mpk.Signal` stores samples, time, sampling frequency, label, units, and
  non-sensitive metadata.
- `mpk.Recording` groups named signals and owns temporal phase intervals.
- `mpk.ProcessingResult` records the modality, processor identity, features,
  source label, options, metadata, and creation time.

MATLAB preserves the modality-specific feature cells returned by the
procedural routines. fNIRS returns named oxy and deoxy feature groups inside a
struct. See [Data model and indexing](../docs/data_model.md).

## Creating phases

Create complete sliding windows in seconds:

```matlab
[intervalsSeconds,phaseLabels] = ...
    mpk.utils.createSlidingPhaseIntervals( ...
        signal.Time(1), ...
        signal.Time(end) + 1/signal.SamplingFrequency, ...
        60,30,'LabelPrefix','window');
recording.setPhases(intervalsSeconds,phaseLabels);
```

The duration and overlap are in seconds. Incomplete trailing windows are
dropped.

Use `mpk.utils.plotSignalWithPhases(signal,intervalsSeconds, ...
'PhaseLabels',phaseLabels)` before processing to check phase placement.

## Creating recordings

Use `mpk.io.recordingFromArrays` for values already in memory. Supply either
`SamplingFrequency` or `Time`.

For files, use:

- `mpk.io.readGenericHdf5` with an explicit dataset mapping; or
- `mpk.io.readBioSignalsPluxHdf5` for supported BioSignalsPlux/OpenSignals
  exports.

fNIRS imports require the dual-wavelength column order to be specified as
`[red infrared]`. Keep file locations and identifying acquisition metadata in
local configuration, not in source or documentation. See
[Data model and indexing](../docs/data_model.md) for the shared input
containers.

## Processors

Processors accept Signals directly, while Recording optionally orchestrates
whole-signal or phase processing:

```matlab
ecgProcessor = mpk.ECGProcessor();
ecgResult = recording.processSignal(ecgProcessor,'ecg');

respProcessor = mpk.RespirationProcessor();
respResults = recording.processPhases(respProcessor,'respiration');

tempProcessor = mpk.TemperatureProcessor();
tempResults = recording.processPhases(tempProcessor,'temperature');

fnirsProcessor = mpk.FNIRSProcessor('AgeYears',35);
fnirsResults = recording.processPhases( ...
    fnirsProcessor,{'fnirs_red','fnirs_infrared'});

edaProcessor = mpk.EDAProcessor('MinimumScrAmplitude',0.01);
edaResults = recording.processPhases(edaProcessor,'eda');
```

The ECG default frequency mode uses STFT-based HRV analysis and requires
phases of at least 60 seconds. EDA requires the official, unmodified cvxEDA
v1.1.0 `cvxEDA.m` to be resolvable on the MATLAB path; its configurable
minimum trough-to-peak SCR amplitude is expressed in µS. fNIRS requires
`AgeYears` and correctly ordered wavelength
channels. Review the [processing overview](../docs/processing_pipeline.md) and
the modality guide before processing study data.

MATLAB processors are object-facing adapters over the validated procedural
implementations in the modality directories. In particular,
`mpk.EDAProcessor` delegates each complete supplied Signal to
`eda/eda_signal_processing.m`. Its decomposition cache is disabled by default,
and EDA plotting occurs only when explicitly requested for a processing call.

The [`examples/matlab`](../examples/matlab) directory contains minimal scripts
for temperature, respiration, EDA, ECG, fNIRS, and multimodal processing. Run
or edit these direct workflows; customize `dataFile`, `phaseIntervals`,
`phaseLabels`, and processor configuration in the script. The examples expose
native `ProcessingResult` objects and do not add wrapper configuration or
reporting layers.

## Testing

When MATLAB is available, run the public test suite from the repository root:

```matlab
run('matlab/tests/run_matlab_tests.m')
```

The public suite uses controlled synthetic inputs, including the bundled
synthetic fixture for reader and example smoke coverage, and
dependency-handling checks. See [Validation](../docs/validation.md) for the
distinction between public tests and local cross-language validation.

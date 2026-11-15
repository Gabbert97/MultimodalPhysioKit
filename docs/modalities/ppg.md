# PPG

## Python API

```python
from multimodalphysiokit.processors import PPGProcessor

processor = PPGProcessor(label_frequency_analysis=1, return_intermediates=True)
result = processor.process(ppg_signal, show_plots=True)
```

`PPGProcessor(*, label_frequency_analysis=1, return_intermediates=False)`
holds reusable configuration. `process(signal, *, show_plots=False,
save_plots=False, output_dir=None, output_label="ppg")` analyses exactly one
supplied `Signal`. `process_ppg(...)` is the functional wrapper, and
`detect_ppg_pulses(signal)` returns only the pulse-detection stage for
debugging.

`label_frequency_analysis` selects HRV frequency analysis: `1` STFT, `2` PSD
(both need at least 60 s of signal), or `None` to skip it. Use `None` for
short signals such as completed real-time windows.

## Processing

The supplied signal is analysed as one piece, whatever its length (a whole
recording, an extracted phase, or a completed real-time window). There is no
internal fixed-length windowing; all thresholds adapt to the pulses found.

1. **Cleaned signal**: 3rd-order zero-phase Butterworth high-pass at 0.5 Hz.
   It removes baseline drift (breathing, motion, sensor drift, often much
   larger than the pulse on wrist devices) while keeping the pulse morphology.
   **All pulse characteristics are measured on this signal.**
2. **Band-pass** 0.5-3.0 Hz: the cleaned signal low-passed at 3.0 Hz
   (3rd-order, zero-phase). It covers 30-180 bpm and is used only to locate
   beats: minima at least 1/3 s apart (180 bpm) are the beat troughs.
3. **Pulse feet (onsets)**: each trough is moved to the minimum of the
   cleaned signal within -10 %/+30 % of the median trough distance.
4. **Pulses**: pulse *k* spans foot *k* to foot *k+1*. Its systolic peak is
   the maximum between the feet.
   - amplitude = peak value minus the straight baseline between the two feet,
     evaluated at the peak time;
   - duration = foot-to-foot time; rise time = foot-to-peak time.
5. **Pulse quality**: a pulse is valid when its duration is within 0.8-1.2x
   the running median duration of up to 11 neighbouring pulses (the whole
   signal's median for short windows) and its peak lies in the first half of
   the pulse.
6. **Inter-beat intervals**: peak-to-peak intervals are kept only between two
   consecutive valid pulses with 40-180 bpm. Rejected intervals are skipped,
   never imputed. pNN50 skips pairs that span a rejected interval, and the
   HRV tachogram places each IBI on its real time and interpolates over gaps.
7. **Frequency analysis** uses the ECG HRV implementation (4 Hz tachogram,
   LF 0.04-0.15 Hz, HF 0.15-0.40 Hz, 60 s window).

Under `Recording.process_phases()`, this sequence runs independently for
every extracted phase.

## Features

Pulse and time-domain keys (always returned):

- `ppg_pulse_amplitude_mean`, `ppg_pulse_amplitude_std`, `ppg_pulse_amplitude_median`
- `ppg_pulse_duration_mean`, `ppg_pulse_duration_std`, `ppg_pulse_duration_median`
- `ppg_rise_time_mean`, `ppg_rise_time_std`, `ppg_rise_time_median`
- `ppg_bpm_mean`, `ppg_bpm_std`, `ppg_bpm_median`
- `ppg_pnn50`, `ppg_pnn50_max` (same definitions as the ECG keys)

Frequency-domain keys (only when `label_frequency_analysis` is 1 or 2):

- `ppg_plf_mean`, `ppg_plf_std`, `ppg_phf_mean`, `ppg_phf_std`,
  `ppg_lf_hf_mean`, `ppg_lf_hf_std`

Pulse statistics use valid pulses only; standard deviations use N-1.
Amplitudes are in the input signal's units, measured on the cleaned signal
(the 0.5 Hz high-pass reduces them slightly, about 0.2 % for a 75 bpm pulse). Metadata includes `n_pulses`,
`n_valid_pulses`, and `n_accepted_ibis` for quick quality checks.

## Intermediates and plots

With `return_intermediates=True`: `cleaned_signal`, `filtered_signal`, `pulse_onset_indices`
(one more than the number of pulses), `pulse_peak_indices`,
`pulse_peak_times`, `pulse_valid`, `pulse_amplitudes`, `pulse_durations`,
`pulse_rise_times`, `ibi`, `ibi_times` (time of the opening peak), and `bpm`.
STFT mode adds `stft_power`, `stft_frequencies`, and `stft_times`. Python
indexes are zero-based; MATLAB indexes are one-based.

`show_plots=True` (or `save_plots=True` with `output_dir`) draws the raw
signal, the cleaned signal with accepted pulses in blue and rejected pulses in
red and feet and peaks marked, the band-passed signal, the BPM trend, and the HRV frequency
plots when frequency analysis is enabled. Saved files: `<label>_pulses.png`,
`<label>_BPM_trend.png`, and the ECG-style `<label>_STFT_*.png` or
`<label>_PSD_LFHF.png`.

## Assumptions

The signal must be sampled above 6 Hz. Pulse detection assumes upward
systolic peaks; invert signals recorded with the opposite polarity (on the
MAUS dataset, both the Infiniti fingertip and PixArt wrist PPG are used as
recorded). The share of valid pulses (`n_valid_pulses / n_pulses`) is a useful
quality indicator: windows with less than about half valid pulses are
unreliable. Motion
artifacts, poor contact, and arrhythmias reduce the number of valid pulses;
check `n_valid_pulses` and the pulse plot.

## Methodological reference

MultimodalPhysioKit independently implements the PPG workflow described by
Luzzani et al., "EDA, PPG and Skin Temperature as Predictive Signals for
Mental Failure by a Statistical Analysis on Stress and Mental Workload,"
*IEEE Open Journal of Engineering in Medicine and Biology*, 2025, DOI:
`10.1109/OJEMB.2024.3515473`. Differences from that workflow: the beat band
is widened from 0.5-1.5 Hz to 0.5-3.0 Hz, pulse characteristics are
measured on a 0.5 Hz high-passed (drift-free) signal instead of the raw one, the signal is analysed as one piece
instead of internal 10 s windows, pulse quality uses a running median
duration, rejected intervals are skipped instead of repeating the previous
IBI, and pNN50 follows the standard definition.

## MATLAB

```matlab
processor = mpk.PPGProcessor('LabelFrequencyAnalysis',1, ...
    'ReturnIntermediates',true,'ShowPlots',true);
result = processor.process(ppgSignal);
```

Use `'LabelFrequencyAnalysis',[]` to skip frequency analysis. The procedural
implementation is `matlab/ppg/ppg_signal_processing.m`. Python and MATLAB
features agree to floating-point precision without frequency analysis and in
STFT mode. In PSD mode, MATLAB estimates the spectrum with Burg's method
(`pburg`) while Python uses a Hamming periodogram, so PSD-mode values differ
between languages (shared with ECG).

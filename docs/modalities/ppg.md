# PPG

## Python API

```python
from multimodalphysiokit.processors import PPGProcessor

processor = PPGProcessor(label_frequency_analysis=1, return_intermediates=True)
result = processor.process(ppg_signal, show_plots=True)
```

`PPGProcessor(*, label_frequency_analysis=1, return_intermediates=False,
max_peak_position=0.5)` holds reusable configuration. `process(signal, *, show_plots=False,
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

1. **Cleaned signal**: 0.5-8 Hz band-pass (3rd-order zero-phase Butterworth
   high-pass at 0.5 Hz, then low-pass at 8 Hz). The high-pass removes
   baseline drift (breathing, motion, sensor drift, often much larger than
   the pulse on wrist devices); the low-pass removes noise above the pulse
   harmonics. **All pulse characteristics are measured on this signal.**
2. **Beat band-pass** 0.5-3.0 Hz: the high-passed signal low-passed at
   3.0 Hz (3rd-order, zero-phase). It covers 30-180 bpm and is used only to
   locate beats: minima at least 1/3 s apart (180 bpm) are candidate troughs.
3. **Dicrotic notch rejection**: a candidate trough is kept only if its
   prominence is at least 30 % of the 75th percentile of the prominences of
   up to 11 neighbouring candidates. The notch between the systolic and the
   dicrotic wave is much shallower than a real foot, so it no longer splits
   a pulse in two.
4. **Pulse feet (onsets)**: each kept trough is moved to the minimum of the
   cleaned signal within -10 %/+30 % of the median trough distance.
5. **Pulses**: pulse *k* spans foot *k* to foot *k+1*. Its systolic peak is
   the **first** local maximum between the feet whose height above the
   foot-to-foot line is at least 70 % of the highest local maximum (the
   global maximum if the pulse has none). A later wave that is taller than
   the systolic one, common on the wrist, is therefore not taken as the
   peak. **Inter-beat intervals are measured between systolic peaks**.
   - amplitude = peak value minus the straight baseline between the two feet,
     evaluated at the peak time;
   - duration = foot-to-foot time; rise time = foot-to-peak time.
6. **Pulse quality**: a pulse is valid when its duration is within 0.8-1.2x
   the running median duration of up to 11 neighbouring pulses (the whole
   signal's median for short windows) and its peak lies within the first
   `max_peak_position` fraction of the pulse (default 0.5, the first half).
7. **Inter-beat intervals**: peak-to-peak intervals are kept only between two
   consecutive valid pulses with 40-180 bpm. Rejected intervals are skipped,
   never imputed. pNN50 skips pairs that span a rejected interval, and the
   HRV tachogram places each IBI on its real time and interpolates over gaps.
8. **Frequency analysis** uses the ECG HRV implementation (4 Hz tachogram,
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
(the 0.5-8 Hz band-pass reduces them slightly, about 1 % for a 75 bpm pulse). Metadata includes `n_pulses`,
`n_valid_pulses`, and `n_accepted_ibis` for quick quality checks.

## Intermediates and plots

With `return_intermediates=True`: `cleaned_signal`, `filtered_signal`, `pulse_onset_indices`
(one more than the number of pulses), `pulse_peak_indices`,
`pulse_peak_times`, `pulse_valid`, `pulse_amplitudes`, `pulse_durations`,
`pulse_rise_times`, `ibi`, `ibi_times` (systolic peak opening each IBI), and `bpm`.
STFT mode adds `stft_power`, `stft_frequencies`, and `stft_times`. Python
indexes are zero-based; MATLAB indexes are one-based.

`show_plots=True` (or `save_plots=True` with `output_dir`) draws the raw
signal, the cleaned signal with accepted pulses in blue and rejected pulses in
red, feet and peaks marked, the band-passed signal, the BPM trend, and the HRV frequency
plots when frequency analysis is enabled. Saved files: `<label>_pulses.png`,
`<label>_BPM_trend.png`, and the ECG-style `<label>_STFT_*.png` or
`<label>_PSD_LFHF.png`.

## Assumptions

The signal must be sampled above 16 Hz. Pulse detection assumes upward
systolic peaks; invert signals recorded with the opposite polarity (on the
MAUS dataset, both the Infiniti fingertip and PixArt wrist PPG are used as
recorded). The share of valid pulses (`n_valid_pulses / n_pulses`) is a useful
quality indicator: windows with less than about half valid pulses are
unreliable.

`max_peak_position` trades pulse yield against pulse quality. Fingertip PPG
peaks sit early and never approach 0.5. Rounded wrist pulses peak later, so
0.5 rejects part of them; a larger value (e.g. 0.65) keeps more wrist pulses
but also admits noisier ones. Rounded or flat-topped systolic peaks also make peak times, and
therefore beat-to-beat features (pNN50, HRV), less precise than heart rate. Motion
artifacts, poor contact, and arrhythmias reduce the number of valid pulses;
check `n_valid_pulses` and the pulse plot.

## Validation on MAUS

Reproduce with `examples/ppg_maus_validation/validate_ppg_maus.py` (needs
the public MAUS dataset and NeuroKit2; see the script header for details).

- **Data**: MAUS (22 participants; rest plus six N-back trials), fingertip
  PPG (Infiniti, 256 Hz, same device as the ECG) and wrist PPG (PixArt,
  100 Hz), cut into 459 non-overlapping 90 s windows per device. Every window
  is processed on its own with the default settings.
- **Reference**: ECG R-peaks (NeuroKit2 `ecg_peaks`, artifact correction).
  R-R intervals more than 20 % from their local median (e.g. ectopic beats)
  and their two neighbours are excluded, leaving normal-to-normal (NN)
  intervals.
- **Comparison**: NeuroKit2 0.3.0 `ppg_clean` + `ppg_peaks` on the same
  windows. NeuroKit2 rejects no beats, so all its intervals in 40-180 bpm are
  used. `bishop` is not listed: it is too slow for windowed use (about 11 s
  per 90 s window at 100 Hz).
- **Metrics**: `ppg_bpm_mean` and `ppg_pnn50` from the library; the same
  definitions for NeuroKit2 beats and for the ECG NN intervals. A PPG beat
  matches an R-peak if it lies within 150 ms of the R-peak plus the window's
  median pulse arrival delay. *NN intervals kept*: share of NN intervals
  covered by an accepted PPG interval between the two matching beats.
  *Accepted intervals not matching ECG*: accepted PPG intervals that do not
  correspond to a pair of consecutive R-peaks. *IBI error*:
  |PPG interval - R-R interval| over the kept intervals. Bias and 95 %
  limits of agreement are Bland-Altman (method - ECG). Runtimes are medians
  per 90 s window, single process, on a laptop.

### Fingertip PPG

| Method | NN intervals kept | Accepted intervals not matching ECG | IBI error, median / 95th pct (ms) | Runtime (ms) |
|---|---|---|---|---|
| MultimodalPhysioKit | 97.5 % | **0.3 %** | 3.9 / 11.7 | 2.1 |
| NK2 elgendi | 99.2 % | 0.7 % | 3.9 / 11.7 | 1.3 |
| NK2 elgendi (artifact correction) | **99.4 %** | 0.6 % | 3.9 / 11.7 | 2.8 |
| NK2 charlton | 98.8 % | 0.6 % | 3.9 / 15.6 | 22.8 |
| NK2 charlton2024 | 98.8 % | 0.6 % | 3.9 / 15.6 | 22.9 |

| Method | Heart rate: median abs. error (IQR), bpm | Bias [95 % LoA], bpm | Within 2 bpm | pNN50: median abs. error (IQR) | Bias [95 % LoA] | Within 0.05 |
|---|---|---|---|---|---|---|
| MultimodalPhysioKit | 0.05 (0.02-0.11) | **-0.02 [-0.52, +0.47]** | **100 %** | 0.014 (0.002-0.029) | **+0.007 [-0.085, +0.098]** | **90.0 %** |
| NK2 elgendi | 0.03 (0.01-0.07) | +0.12 [-1.37, +1.62] | 97.4 % | 0.014 (0.003-0.031) | +0.018 [-0.098, +0.133] | 87.6 % |
| NK2 elgendi (artifact correction) | 0.03 (0.01-0.06) | +0.09 [-1.29, +1.46] | 98.3 % | 0.014 (0.001-0.029) | +0.014 [-0.097, +0.126] | 89.3 % |
| NK2 charlton | 0.04 (0.02-0.08) | +0.14 [-1.27, +1.54] | 97.4 % | 0.016 (0.007-0.038) | +0.020 [-0.087, +0.127] | 84.5 % |
| NK2 charlton2024 | 0.04 (0.02-0.09) | +0.11 [-1.11, +1.32] | 98.5 % | 0.015 (0.004-0.038) | +0.019 [-0.086, +0.123] | 82.4 % |

On good-quality fingertip PPG all methods find the same beats; 3.9 ms is one
sample at 256 Hz. MultimodalPhysioKit keeps about 2 % fewer intervals,
because it rejects irregular pulses, and in exchange has the fewest wrong
intervals and the narrowest limits of agreement.

### Wrist PPG

| Method | NN intervals kept | Accepted intervals not matching ECG | IBI error, median / 95th pct (ms) | Runtime (ms) |
|---|---|---|---|---|
| MultimodalPhysioKit | 54.7 % | **35.7 %** | **23.8 / 123** | 1.5 |
| NK2 elgendi | 55.9 % | 42.6 % | 25.5 / 141 | 1.2 |
| NK2 elgendi (artifact correction) | **57.4 %** | 41.7 % | 26.1 / 136 | 3.2 |
| NK2 charlton | 56.7 % | 42.1 % | 25.8 / 144 | 14.4 |
| NK2 charlton2024 | 56.8 % | 41.9 % | 25.8 / 143 | 14.4 |

| Method | Heart rate: median abs. error (IQR), bpm | Bias [95 % LoA], bpm | Within 2 bpm | pNN50: median abs. error (IQR) | Bias [95 % LoA] | Within 0.05 |
|---|---|---|---|---|---|---|
| MultimodalPhysioKit | **0.56 (0.20-1.10)** | **-0.17 [-2.59, +2.25]** | **91.1 %** | **0.272 (0.123-0.465)** | +0.301 [-0.120, +0.722] | 12.0 % |
| NK2 elgendi | 0.98 (0.40-1.98) | +0.53 [-3.54, +4.61] | 75.4 % | 0.365 (0.191-0.564) | +0.378 [-0.089, +0.844] | 8.5 % |
| NK2 elgendi (artifact correction) | 0.89 (0.31-1.88) | +0.47 [-3.33, +4.27] | 77.8 % | 0.341 (0.147-0.536) | +0.349 [-0.118, +0.816] | 10.0 % |
| NK2 charlton | 1.01 (0.43-1.98) | +0.59 [-2.83, +4.01] | 75.2 % | 0.359 (0.182-0.567) | +0.371 [-0.094, +0.836] | 8.7 % |
| NK2 charlton2024 | 1.01 (0.39-1.93) | +0.52 [-2.92, +3.97] | 76.9 % | 0.357 (0.182-0.568) | +0.370 [-0.094, +0.834] | 8.7 % |

Wrist heart rate is usable, but **wrist pNN50 is not reliable with any
method** (it is overestimated by about 0.3). The wrist PPG and the ECG come
from different devices whose alignment within a window was not verified, so
the absolute wrist errors may be overstated; the ranking of the methods is
the informative part.

## Methodological reference

MultimodalPhysioKit independently implements the PPG workflow described by
Luzzani et al., "EDA, PPG and Skin Temperature as Predictive Signals for
Mental Failure by a Statistical Analysis on Stress and Mental Workload,"
*IEEE Open Journal of Engineering in Medicine and Biology*, 2025, DOI:
`10.1109/OJEMB.2024.3515473`. Differences from that workflow: the beat band
is widened from 0.5-1.5 Hz to 0.5-3.0 Hz, shallow troughs (dicrotic notches)
are not used as feet, the systolic peak is the first sufficiently tall local
maximum instead of the global maximum, pulse characteristics are measured on
a 0.5-8 Hz band-passed signal instead of the raw one, the signal is analysed as one piece
instead of internal 10 s windows, pulse quality uses a running median
duration, rejected intervals are skipped instead of repeating the previous
IBI, and pNN50 follows the standard definition.

## MATLAB

```matlab
processor = mpk.PPGProcessor('LabelFrequencyAnalysis',1, ...
    'ReturnIntermediates',true,'ShowPlots',true,'MaxPeakPosition',0.5);
result = processor.process(ppgSignal);
```

Use `'LabelFrequencyAnalysis',[]` to skip frequency analysis. The procedural
implementation is `matlab/ppg/ppg_signal_processing.m`. Python and MATLAB
features agree to floating-point precision without frequency analysis and in
STFT mode. In PSD mode, MATLAB estimates the spectrum with Burg's method
(`pburg`) while Python uses a Hamming periodogram, so PSD-mode values differ
between languages (shared with ECG).

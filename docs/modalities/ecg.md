# ECG

## Python API

```python
from multimodalphysiokit.processors import ECGProcessor

processor = ECGProcessor(
    label_frequency_analysis=1,
    return_intermediates=True,
)
result = processor.process(ecg_signal)
```

`ECGProcessor(*, label_frequency_analysis=1, return_intermediates=False)`
contains reusable configuration. `process(signal, *, show_plots=False,
save_plots=False, output_dir=None, output_label="ecg")` analyzes exactly one
supplied `Signal`.

Frequency mode `1` uses the validated STFT path; mode `2` uses the validated
PSD path. STFT duration requirements apply to each supplied signal.

## Processing

The processor filters ECG, detects and screens R-peaks, calculates inter-beat
intervals and BPM, extracts time-domain HRV features, and performs the selected
frequency analysis. Under `Recording.process_phases()`, this complete sequence
restarts independently for every extracted phase, including STFT generation.

## Features

Time-domain keys:

- `ecg_bpm_mean`
- `ecg_bpm_std`
- `ecg_bpm_median`
- `ecg_pnn50`
- `ecg_pnn50_max`

Frequency-domain keys:

- `ecg_plf_mean`
- `ecg_plf_std`
- `ecg_phf_mean`
- `ecg_phf_std`
- `ecg_lf_hf_mean`
- `ecg_lf_hf_std`

Feature names never include a phase label.

`ecg_pnn50` is the share of successive differences between adjacent accepted
inter-beat intervals that exceed 50 ms (Task Force of the ESC/NASPE, 1996),
reported as a fraction. A difference is evaluated only when both neighbouring
intervals were accepted by IBI screening; pairs spanning a rejected interval
are skipped. `ecg_pnn50_max` is the longest run of consecutive differences
exceeding 50 ms, divided by the expected number of IBIs in two minutes
(`120 * n_ibi / sum(ibi)`). Both are zero when no difference qualifies.

In PSD mode (`label_frequency_analysis=2`), every accepted IBI is placed at
the real time of the R-peak that closes it before resampling the tachogram
at 4 Hz. Rejected intervals therefore leave a gap that is interpolated over,
rather than being removed from the timeline (which would compress time and
shift spectral content upward).

Result metadata contains `source`, `sampling_frequency`, `n_samples`,
`start_time_seconds`, `end_time_seconds`, `frequency_analysis_method`, and
`return_intermediates`.

## Optional intermediates

With `return_intermediates=True`, results include:

- `filtered_signal`
- `r_peak_indices`
- `r_peak_times`
- `r_peak_amplitudes`
- `accepted_r_peak_indices`
- `accepted_r_peak_times`
- `ibi`
- `ibi_times`
- `bpm`
- `bpm_times`

STFT mode also includes `stft_matrix`, `stft_power`, `stft_frequencies`, and
`stft_times`. These arrays belong in `ProcessingResult.intermediates`, not in
scalar feature tables.

## Assumptions

R-peak quality controls downstream results. Motion artifacts, ectopic or
missing beats, morphology, signal polarity, acquisition rate, and units
require study-specific review. Independently processing extracted phases may
produce different edge behavior from filtering a complete parent signal.

## Methodological references

MultimodalPhysioKit independently implements the ECG processing workflow
described by Luzzani et al., “ECG, Respiration, fNIRS, and Eye Tracking for
Stress and Mental Workload Monitoring in Human-Machine Interaction,” *IEEE
Access*, 2025, DOI: `10.1109/ACCESS.2025.3588384`.

The independently implemented QRS detector follows the methodological lineage
of J. Pan and W. J. Tompkins, “A Real-Time QRS Detection Algorithm,” *IEEE
Transactions on Biomedical Engineering*, BME-32(3), 230-236, 1985, DOI:
`10.1109/TBME.1985.325532`. These publications are scientific references, not
source-code origins.

## MATLAB

MATLAB uses `mpk.ECGProcessor` or the validated procedural ECG functions.
Public `mpk.Recording` phases are temporal half-open intervals. The validated
procedural ECG entry point internally retains one-based inclusive sample
indexes; map temporal intervals through the Signal time vector before direct
procedural parity comparisons.

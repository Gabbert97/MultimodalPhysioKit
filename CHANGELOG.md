# Changelog

Notable changes to MultimodalPhysioKit are documented in this file.

## [Unreleased]

All changes below apply to both Python and MATLAB unless stated otherwise.

### Fixed

- ECG `ecg_pnn50` / `ecg_pnn50_max`: after a rejected inter-beat interval,
  successive IBI differences were taken between non-adjacent intervals (the
  IBI counter fell one position behind per rejected interval). Differences
  are now only taken between truly adjacent accepted IBIs.
- BioSignalsPlux NTC temperature conversion subtracted 273.5 instead of
  273.15 when converting kelvin to degrees Celsius, reporting every
  temperature 0.35 °C too low. **Feature values change:** `temperature_mean`
  increases by 0.35 °C; differences, slopes, derivatives, and standard
  deviations are unchanged. The synthetic fixture generator uses the same
  constant, and the bundled sample file was regenerated (only its raw
  temperature channel changed; the encoded physical temperatures are
  identical).
- EDA decomposition cache keys now include the decomposition source (Python:
  callable module/qualified name or a hash of the cvxEDA source file;
  MATLAB: `func2str` of the callable or a hash of `cvxEDA.m`) and the cvxEDA
  parameters. Previously a different decomposition silently reused another
  source's cached result. Existing cache files are no longer matched and are
  recomputed. In Python, the decomposition source must now be resolvable even
  when a cached result exists.

### Changed

- **ECG pNN50 definition (feature values change).** `ecg_pnn50` is now the
  standard share of successive differences > 50 ms among evaluated adjacent
  pairs (Task Force, 1996). Previously only runs of two or more consecutive
  differences >= 50 ms were counted, divided by the number of IBIs, so
  isolated large differences were ignored. `ecg_pnn50_max` now counts runs of
  any length (including 1) of differences > 50 ms. MATLAB computes both in
  the new `matlab/ecg/ecg_pnn50_features.m`.
- **EDA and respiration anti-alias downsampling (feature values change for
  inputs above 100 Hz).** A zero-phase 4th-order Butterworth low-pass at
  40 Hz is applied before keeping every k-th sample. Previously samples were
  decimated directly, so content above 50 Hz (e.g. mains harmonics near
  100/200 Hz) aliased into the analysis band. Sample selection and timing
  are unchanged; 100 Hz inputs are unaffected.
- **ECG PSD mode (`label_frequency_analysis=2`) tachogram timing (feature
  values change when IBI screening rejects intervals).** Each accepted IBI is
  placed at the real time of its closing R-peak, so rejected intervals are
  interpolated over instead of being removed from the timeline, which
  compressed time and shifted spectral content upward. Gap-free signals are
  unaffected. In MATLAB, each phase's PSD now uses its own real-time
  tachogram. STFT mode is unchanged.
- **Temperature artifact screening (feature values change when screening
  removes windows).** The amplitude criterion is now an absolute deviation of
  more than 2.0 degrees from the phase mean (previously ±6 % of the mean,
  which depended on the unit and never triggered in kelvin). The rate
  criterion compares a window's mean absolute change with 100 times the
  phase's *median* absolute change (previously the signed mean change, which
  made the criterion always true while cooling and so wrongly removed edges of
  smooth cooling trends).

## [0.2.0]

### Added

- Configurable minimum and maximum accepted breathing rates for the Python
  `RespirationProcessor` and `process_respiration` APIs. Defaults remain 5 and
  30 breaths/minute to preserve version 0.1.0 behavior.

## [0.1.0]

### Added

- Standardized `Signal` and `Recording` data models.
- Standardized `ProcessingResult` outputs.
- ECG processing.
- EDA processing with optional external cvxEDA integration.
- Respiration processing.
- Skin-temperature processing.
- fNIRS processing.
- Generic HDF5 loading.
- BioSignalsPlux/OpenSignals HDF5 adapter.
- MATLAB implementation.
- Python implementation.
- Deterministic synthetic multimodal sample dataset.
- Python and MATLAB test suites.
- Python packaging.
- Python CI workflow.
- MATLAB CI workflow.

# Changelog

Notable changes to MultimodalPhysioKit are documented in this file.

## [Unreleased]

_No unreleased changes yet._

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

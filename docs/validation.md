# Validation

MultimodalPhysioKit separates unit testing, bundled demonstrations, local
validation, and cross-language numerical comparisons.

## Python unit tests

The public pytest suite uses controlled synthetic inputs and bundled fixtures
to test data validation, temporal phase mapping, I/O adapters, dependency
handling, processor APIs, intermediate outputs, and numerical behavior.

```bash
cd python
python -m pytest -v
```

## Bundled demonstration notebooks

The notebooks under `examples/python/` use
`read_biosignalsplux_hdf5`, public plotting utilities, and public processor
APIs. The bundled HDF5 fixture contains only mathematically generated signals
and demonstrates the workflow; it is not a scientific validation dataset.

## Local validation

Private recordings and derived comparisons remain outside the public
repository. They are not required by package installation or public tests.
See [Local validation](local_validation.md).

## MATLAB tests

When MATLAB and required toolboxes are available:

```matlab
run('matlab/tests/run_matlab_tests.m')
```

MATLAB dependency tests may report unavailable toolboxes or external cvxEDA.

## MATLAB/Python numerical comparisons

Cross-language comparisons align:

- samples, units, sampling frequencies, and channel meanings;
- shared Python and MATLAB temporal intervals in seconds;
- preprocessing parameters and fNIRS baseline values;
- stable feature names, feature order, and phase order; and
- floating-point values within documented tolerances.

Python associates phases through the outer dictionary returned by
`Recording.process_phases()` and keeps feature names phase-independent.
MATLAB retains its modality-specific feature containers.

For fNIRS, compare using the same filtered-current scalar baseline values and
wavelength ordering. For temperature, the validated Python feature set is:

- `temperature_delta`
- `temperature_mean`
- `temperature_std`
- `temperature_delta_over_time`
- `temperature_linear_slope`
- `temperature_derivative_mean`
- `temperature_derivative_std`
- `temperature_derivative_linear_slope`

Parity checks increase confidence in equivalent implementations; they do not
establish clinical validity or suitability for every protocol.

## Reproducible reporting

Record software and dependency versions, processor options, sampling
frequencies, temporal phase definitions, MATLAB index conversions, and
comparison tolerances. Share only data and outputs approved for distribution.

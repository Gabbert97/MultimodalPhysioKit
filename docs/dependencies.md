# Dependencies

## Core Python dependencies

The package requires Python 3.10 or newer and declares:

- NumPy
- SciPy
- h5py
- Matplotlib

Install from the repository root:

```bash
python -m pip install -e python
```

Matplotlib supports package plotting utilities and processor diagnostics.

## Notebook and example dependencies

Jupyter is needed to run the example notebooks. The multimodal notebook uses
pandas to display a long-format feature table. Jupyter and pandas are example
dependencies and are not declared as core package requirements.

```bash
python -m pip install jupyter pandas
```

## Testing dependencies

The declared `test` extra installs pytest:

```bash
python -m pip install -e "python[test]"
```

The project currently declares no separate development extra. Contributors
install any additional formatting, documentation, or notebook tools explicitly
for the task at hand.

## Python cvxEDA dependency

EDA requires official, unmodified cvxEDA v1.1.0 as its validated external
reference implementation. There is no silent
fallback decomposition. `EDAProcessor` can resolve an installed implementation,
accept `cvxeda_path`, or accept a compatible `decomposition_callable`.

The processor's missing-dependency error recommends:

```bash
python -m pip install "cvxeda==1.1.0"
```

The upstream project and scientific reference are:
<https://github.com/lciti/cvxEDA>

> A. Greco, G. Valenza, A. Lanata, E. P. Scilingo, and L. Citi, “cvxEDA: A
> Convex Optimization Approach to Electrodermal Activity Processing,” *IEEE
> Transactions on Biomedical Engineering*, 63(4), 797–804, 2016.
> <https://doi.org/10.1109/TBME.2015.2474131>

Review its independent license and citation requirements. See
[Third-Party Notices](legal/THIRD_PARTY_NOTICES.md).

## MATLAB dependencies

MATLAB processing can require:

- Signal Processing Toolbox
- Statistics and Machine Learning Toolbox
- Optimization Toolbox for the default MATLAB cvxEDA solver

For MATLAB EDA, obtain the official, unmodified `cvxEDA.m` from the cvxEDA
v1.1.0 release and place it on the MATLAB path. Do not use a locally modified
copy that clips negative phasic values. MultimodalPhysioKit supplies the
scientific parameters explicitly (`tau0=2.0`, `tau1=0.7`, `delta_knot=10.0`,
`alpha=8e-4`, and `gamma=1e-2`) while leaving the returned decomposition
unchanged. Requirements vary by modality; see the MATLAB tests and modality
guides.

## Data are not dependencies

The bundled HDF5 file is a mathematically generated demonstration fixture.
Private/local validation recordings are not package dependencies and are not
required by public tests.

# Release checklist

Use this checklist for every release `vX.Y.Z`. Publishing a GitHub Release
triggers `.github/workflows/release.yml`, which rebuilds the distributions
from the release tag and uploads them to
[PyPI](https://pypi.org/project/multimodalphysiokit/). A version can be
uploaded to PyPI only once.

## Safety rules

- [ ] Keep `external/cvxEDA.m`, `validation-output/`, `local_validation/`,
      and `CLAUDE.md` untracked.
- [ ] Confirm no private or raw physiological recordings are tracked. The
      deterministic synthetic HDF5 fixture is the only tracked recording.
- [ ] Search tracked content for secrets, credentials, identifying data,
      local paths, private dataset references, and machine-specific
      configuration.
- [ ] Push only `main` and the new release tag; never use
      `git push --mirror`, `git push --all`, or `git push --tags`.

## Before the release

- [ ] All intended changes are merged into `main` through a pull request,
      and the working tree is clean.
- [ ] Python CI and MATLAB CI are green on `main`.
- [ ] The complete local Python suite passes:
      `python -m pytest python/tests -q`.
- [ ] The complete MATLAB suite passes locally
      (`run('matlab/tests/run_matlab_tests.m')`), allowing only documented
      optional-dependency skips.
- [ ] The MATLAB/Python cross-language notebook
      (`examples/cross_language_validation/`) passes for ECG, EDA,
      respiration, temperature, PPG, and fNIRS (needs official cvxEDA v1.1.0
      in both languages).
- [ ] If PPG processing changed, rerun
      `examples/ppg_maus_validation/validate_ppg_maus.py` on MAUS and update
      the tables in `docs/modalities/ppg.md`.
- [ ] `python tools/generate_synthetic_sample.py --check` passes.
- [ ] `python -m build python --outdir dist` and `python -m twine check dist/*`
      pass.
- [ ] The wheel installs in a fresh environment and the Python tests pass
      against the installed distribution.
- [ ] The three `python/` license copies match their repository-root
      counterparts byte-for-byte.

## Version and notes

- [ ] `python/pyproject.toml` and `CITATION.cff` declare version `X.Y.Z`.
- [ ] `CHANGELOG.md` has a `## [X.Y.Z]` section that accurately describes the
      release, including every change that alters feature values.
- [ ] Documentation (READMEs, `docs/`, examples) matches the released
      behavior.

## Publish

- [ ] Create the annotated tag on the verified `main` commit and push only
      that tag:

      ```bash
      git tag -a vX.Y.Z -m "MultimodalPhysioKit vX.Y.Z"
      git push origin vX.Y.Z
      ```

- [ ] Create the GitHub Release from `vX.Y.Z`, using the changelog section as
      the basis for concise release notes, and publish it.
- [ ] Confirm the "Publish Python package" workflow succeeds and the new
      version appears on PyPI.

## After the release

- [ ] Confirm Zenodo archived the release, and update the DOI in
      `CITATION.cff` if appropriate.
- [ ] Install `multimodalphysiokit==X.Y.Z` from PyPI in a fresh environment
      and import the package.
- [ ] Verify the public README, PyPI page, release links, citation metadata,
      and CI badges.

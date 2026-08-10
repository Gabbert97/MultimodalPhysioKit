# Release checklist

This checklist records the clean-snapshot release process for
MultimodalPhysioKit v0.1.0. The existing development repository remains a
private archive; its Git history is never transferred to the public repository.

## Safety rules

- [ ] **Never make the old private development repository public.**
- [ ] **Never copy its `.git` directory into the public snapshot.**
- [ ] **Never push old branches, tags, pull-request refs, or history to the new
      public repository.**
- [ ] **Never use `git push --mirror`, `git push --all`, or `git push --tags`
      during migration.**
- [ ] Keep `external/cvxEDA.m` local and untracked.
- [ ] Keep `validation-output/` and `local_validation/` untracked.
- [ ] Confirm no private or raw physiological recordings are included. The
      tracked deterministic synthetic HDF5 fixture is allowed.

## Pre-release validation

- [ ] Finalize `public-release-v0.1` in the private development repository.
- [ ] Confirm the working tree is clean and all intended changes are committed.
- [ ] Confirm Python CI and the complete local Python test suite are green.
- [ ] Confirm the complete MATLAB test suite is green locally, allowing only
      documented optional-dependency assumptions.
- [ ] Confirm the MATLAB/Python cross-language validation passes for ECG, EDA,
      respiration, temperature, and fNIRS.
- [ ] Run `python tools/generate_synthetic_sample.py --check` successfully.
- [ ] Build the Python sdist and wheel with
      `python -m build python --outdir dist`.
- [ ] Run `python -m twine check dist/*` successfully.
- [ ] Install the wheel in a fresh environment and run the Python tests against
      the installed distribution.
- [ ] Confirm the three `python/` license copies match their repository-root
      counterparts byte-for-byte.
- [ ] Confirm `CITATION.cff` and `python/pyproject.toml` declare version `0.1.0`.
- [ ] Confirm `CHANGELOG.md` accurately describes v0.1.0.
- [ ] Search tracked content for secrets, credentials, identifying data, local
      paths, private dataset references, and machine-specific configuration.

## Final private release commit

- [ ] Merge or otherwise finalize the release branch internally if desired.
- [ ] Select and record the exact final release commit SHA.
- [ ] Confirm that commit contains the complete intended v0.1.0 tracked tree.

## Clean snapshot creation

- [ ] Export the selected commit with `git archive` or an equivalent tracked-
      file-only method.
- [ ] Confirm the snapshot excludes `.git`, ignored local files, external
      dependencies, validation outputs, and all old history.
- [ ] Compare the snapshot manifest with the tracked files in the selected
      commit and confirm they match exactly.

## New Git repository initialization

- [ ] Run `git init --initial-branch=main` inside the clean snapshot.
- [ ] Create one commit named `Initial public release snapshot for v0.1.0`.
- [ ] Confirm the new repository has exactly one commit.
- [ ] Confirm its branch is `main`, it has no tags or historical refs, and no
      remote is configured yet.

## New GitHub repository creation

- [ ] Rename the existing GitHub repository to a private archive name such as
      `MultimodalPhysioKit-private-archive` and confirm it remains private.
- [ ] Create `Gabbert97/MultimodalPhysioKit` as a new private repository.
- [ ] Do not initialize the new GitHub repository with a README, license,
      `.gitignore`, or template.
- [ ] Add the new repository as the clean snapshot's remote.
- [ ] Push only the new `main` branch; do not push any other ref.

## Pre-public verification

- [ ] Fresh-clone the new GitHub repository into a separate directory.
- [ ] Confirm the clone has exactly one commit, only the `main` branch, no
      tags, no pull-request refs, and a clean working tree.
- [ ] Confirm the fresh clone contains exactly the expected release files.
- [ ] From the fresh clone, rerun the Python tests, synthetic-fixture check,
      package build, and `twine check`.
- [ ] Confirm the old private archive still exists and remains private.

## Public visibility

- [ ] Make only the new `Gabbert97/MultimodalPhysioKit` repository public.
- [ ] Confirm again that the old development/archive repository remains private.
- [ ] Confirm GitHub Python CI is green on public `main`.
- [ ] Confirm MATLAB CI is green using MathWorks public-repository licensing.
- [ ] Confirm both README CI badges are green.

## v0.1.0 release

- [ ] Only after public CI succeeds, create the annotated tag `v0.1.0` from the
      verified public release commit.
- [ ] Push only the `v0.1.0` tag.
- [ ] Create the GitHub Release from `v0.1.0`, using the changelog entry as the
      basis for concise release notes.

## Post-release publication

- [x] Publish the validated Python artifacts to
      [PyPI](https://pypi.org/project/multimodalphysiokit/).
- [x] Update README installation wording to reflect PyPI availability.
- [ ] Connect the public repository to Zenodo and archive the release.
- [ ] Record the resulting DOI and update `CITATION.cff` when appropriate.
- [ ] Verify the public README, release links, citation metadata, and CI badges.

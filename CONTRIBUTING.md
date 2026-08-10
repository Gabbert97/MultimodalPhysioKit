# Contributing to MultimodalPhysioKit

Contributions should be general-purpose, reproducible, and suitable for a
public physiological-signal-processing repository.

## Before starting

- Check existing issues and pull requests.
- Discuss large features, API changes, numerical-method changes, or new
  dependencies before implementation.
- Contribute only material you have the right to distribute.
- Do not include recordings, identifying information, credentials, restricted
  dataset details, or machine-specific configuration.

## Development expectations

- Keep changes focused and follow the surrounding Python or MATLAB style.
- Add or update documentation for user-visible behavior.
- Add focused synthetic tests for features and bug fixes.
- Document dependencies, their purpose, and their license implications.
- Keep local data, local validation outputs, figures, caches, and generated
  files out of git.

Run Python tests from `python/`:

```bash
python -m pytest -v
```

When MATLAB is available, run its public tests from the repository root:

```matlab
run('matlab/tests/run_matlab_tests.m')
```

Public tests must use synthetic data only. Local validation may use authorized
non-public recordings, but those inputs and outputs must remain local. See
[Validation](docs/validation.md) and
[Local validation](docs/local_validation.md).

## Pull requests

Describe the purpose, key decisions, limitations, and validation performed.
Include tests, or explain why tests do not apply. Remove sensitive or local
information from logs and screenshots before sharing them.

Contributions to MultimodalPhysioKit are provided under the repository's
LGPL-3.0-or-later license unless explicitly stated otherwise.

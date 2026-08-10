"""Python package for MultimodalPhysioKit."""

from importlib.metadata import PackageNotFoundError, version

from .core import ProcessingResult, Recording, Signal

try:
    __version__ = version("multimodalphysiokit")
except PackageNotFoundError:  # Direct import from an uninstalled source tree.
    __version__ = "unknown"

__all__ = [
    "__version__",
    "Signal",
    "Recording",
    "ProcessingResult",
]

"""Python processor entry points."""

from .ecg import ECGProcessor, detect_ecg_r_peaks, process_ecg
from .eda import EDAProcessor, process_eda
from .fnirs import FNIRSProcessor, process_fnirs
from .respiration import RespirationProcessor, process_respiration
from .temperature import TemperatureProcessor, process_temperature

__all__ = [
    "ECGProcessor",
    "EDAProcessor",
    "FNIRSProcessor",
    "RespirationProcessor",
    "TemperatureProcessor",
    "detect_ecg_r_peaks",
    "process_ecg",
    "process_eda",
    "process_fnirs",
    "process_respiration",
    "process_temperature",
]

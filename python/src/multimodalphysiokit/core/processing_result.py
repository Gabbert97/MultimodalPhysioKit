from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field
from typing import Any


@dataclass
class ProcessingResult:
    """Lightweight container for processor outputs and metadata."""

    modality: str
    processor_name: str
    features: dict[str, Any] | None = None
    signal_labels: list[str] | tuple[str, ...] | str | None = None
    options: dict[str, Any] | None = None
    metadata: dict[str, Any] | None = None
    intermediates: dict[str, object] = field(default_factory=dict, repr=False)

    def __post_init__(self) -> None:
        self.modality = str(self.modality)
        self.processor_name = str(self.processor_name)
        self.features = {} if self.features is None else dict(self.features)
        self.options = {} if self.options is None else dict(self.options)
        self.metadata = {} if self.metadata is None else dict(self.metadata)
        if not isinstance(self.intermediates, Mapping):
            raise TypeError("intermediates must be a mapping")
        self.intermediates = dict(self.intermediates)

        if self.signal_labels is None:
            self.signal_labels = []
        elif isinstance(self.signal_labels, str):
            self.signal_labels = [self.signal_labels]
        else:
            self.signal_labels = [str(label) for label in self.signal_labels]

    def feature_names(self) -> list[str]:
        """Return feature names in insertion order."""

        return list(self.features.keys())

    def has_feature(self, name: str) -> bool:
        """Return whether a feature exists."""

        return str(name) in self.features

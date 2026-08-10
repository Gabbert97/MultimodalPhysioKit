from __future__ import annotations

import importlib
from types import SimpleNamespace

import numpy as np
import pytest

import multimodalphysiokit.processors.eda as eda_module
from multimodalphysiokit.core import Signal
from multimodalphysiokit.processors import EDAProcessor
from multimodalphysiokit.processors.eda import _filter_eda, _resolve_cvxeda_callable


def make_signal(*, offset=0.0, duration=70.0, sampling_frequency=1000.0):
    time = offset + np.arange(0, duration, 1 / sampling_frequency)
    samples = 2 + 0.2 * np.sin(2 * np.pi * 0.03 * time) + 0.03 * np.sin(2 * np.pi * 0.8 * time)
    return Signal(samples, time, sampling_frequency, label="eda", units="uS")


def fake_cvxeda(eda, delta):
    n = len(eda)
    zeros = np.zeros(n)
    return np.full(n, 0.1), zeros, np.asarray(eda) - 0.1, zeros, zeros, zeros, 0.0


def make_processor(**kwargs):
    return EDAProcessor(decomposition_callable=fake_cvxeda, **kwargs)


def test_configuration_is_stored_without_data():
    processor = EDAProcessor(
        decomposition_callable=fake_cvxeda, cache_directory="cache", output_label="subject"
    )
    assert processor.decomposition_callable is fake_cvxeda
    assert str(processor.cache_directory) == "cache"
    assert processor.output_label == "subject"
    assert not hasattr(processor, "signal")


@pytest.mark.parametrize("name,value", [
    ("cvxeda_path", object()), ("cache_directory", object()),
    ("decomposition_callable", "bad"), ("output_label", 1),
])
def test_invalid_configuration(name, value):
    with pytest.raises(TypeError, match=name):
        EDAProcessor(**{name: value})


def test_downsampling_values_and_time_are_exact_slices():
    signal = make_signal(offset=4.2)
    result = make_processor().downsample_signal(signal)
    np.testing.assert_array_equal(result["eda"], signal.samples[::10])
    np.testing.assert_array_equal(result["time"], signal.time[::10])
    assert result["sampling_frequency"] == 100


def test_filtering_matches_existing_helper_exactly():
    signal = make_signal()
    downsampled = signal.samples[::10]
    result = make_processor().filter_signal(signal)
    np.testing.assert_allclose(result["eda_filtered"], _filter_eda(downsampled), rtol=0, atol=0)


def test_cvxeda_receives_filtered_signal_and_existing_delta():
    captured = {}
    def capture(eda, delta):
        captured["eda"] = np.asarray(eda).copy()
        captured["delta"] = delta
        return fake_cvxeda(eda, delta)
    signal = make_signal()
    processor = EDAProcessor(decomposition_callable=capture)
    result = processor.decompose(signal)
    np.testing.assert_array_equal(captured["eda"], processor.filter_signal(signal)["eda_filtered"])
    assert captured["delta"] == pytest.approx(0.01)
    np.testing.assert_array_equal(result["scr"], np.full(result["scr"].size, 0.1))
    np.testing.assert_array_equal(result["scl"], result["eda_filtered"] - 0.1)
    np.testing.assert_array_equal(result["residuals"], 0.0)


def test_external_cvxeda_receives_explicit_scientific_parameters(monkeypatch):
    captured = {}
    def capture(eda, delta, **kwargs):
        captured["delta"] = delta
        captured["parameters"] = kwargs
        return fake_cvxeda(eda, delta)
    monkeypatch.setattr(eda_module, "_resolve_cvxeda_callable", lambda **kwargs: capture)
    EDAProcessor().decompose(make_signal())
    assert captured["delta"] == pytest.approx(0.01)
    assert captured["parameters"] == {
        "tau0": 2.0, "tau1": 0.7, "delta_knot": 10.0,
        "alpha": 8e-4, "gamma": 1e-2,
    }


def test_generator_callable_is_supported():
    def generator(eda, delta):
        return (value for value in fake_cvxeda(eda, delta))
    assert EDAProcessor(decomposition_callable=generator).decompose(make_signal())["scr"].size == 7000


def test_missing_dependency_error_is_clear(monkeypatch):
    monkeypatch.setattr(importlib, "import_module", lambda name: (_ for _ in ()).throw(ImportError(name)))
    with pytest.raises(ImportError, match="pip install cvxeda"):
        EDAProcessor().decompose(make_signal())


def test_dependency_resolution_order(monkeypatch):
    calls = []
    def importer(name):
        calls.append(name)
        if name == "cvxeda.cvxEDA":
            return SimpleNamespace(cvxEDA=fake_cvxeda)
        raise ImportError(name)
    monkeypatch.setattr(importlib, "import_module", importer)
    assert _resolve_cvxeda_callable() is fake_cvxeda
    assert calls == ["cvxeda", "cvxeda.cvxEDA"]


def test_cvxeda_file_and_directory_paths_are_supported(tmp_path):
    source = tmp_path / "cvxEDA.py"
    source.write_text(
        "import numpy as np\n"
        "def cvxEDA(eda, delta, **kwargs):\n"
        " n=len(eda); z=np.zeros(n); return np.ones(n),z,np.asarray(eda),z,z,z,0.0\n"
    )
    for path in (source, tmp_path):
        result = EDAProcessor(cvxeda_path=path).decompose(make_signal())
        np.testing.assert_array_equal(result["scr"], 1.0)


def test_cache_reuses_identical_signal_and_invalidates_changed_signal(tmp_path):
    calls = 0
    def counting(eda, delta):
        nonlocal calls
        calls += 1
        return fake_cvxeda(eda, delta)
    processor = EDAProcessor(decomposition_callable=counting, cache_directory=tmp_path)
    signal = make_signal()
    processor.decompose(signal)
    processor.decompose(signal)
    assert calls == 1
    changed = Signal(signal.samples + 0.01, signal.time, signal.sampling_frequency, label="eda")
    processor.decompose(changed)
    assert calls == 2
    assert len(list(tmp_path.glob("eda_*_decomposition.npz"))) == 2


def test_cache_can_be_read_by_another_processor(tmp_path):
    signal = make_signal()
    EDAProcessor(decomposition_callable=fake_cvxeda, cache_directory=tmp_path).decompose(signal)
    def fail(eda, delta):
        raise AssertionError("cache should be used")
    result = EDAProcessor(decomposition_callable=fail, cache_directory=tmp_path).decompose(signal)
    assert result["scr"].size == 7000


def test_too_short_filtering_error_is_preserved():
    signal = make_signal(duration=0.02)
    with pytest.raises(ValueError, match="too short for low-pass filtering"):
        make_processor().decompose(signal)

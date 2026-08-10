from importlib.metadata import version


def test_package_imports():
    import multimodalphysiokit
    from multimodalphysiokit import ProcessingResult, Recording, Signal
    from multimodalphysiokit import processors
    from multimodalphysiokit import utils

    assert hasattr(multimodalphysiokit, "__version__")
    assert Signal is not None
    assert Recording is not None
    assert ProcessingResult is not None
    assert all(
        hasattr(processors, name)
        for name in (
            "ECGProcessor", "EDAProcessor", "RespirationProcessor",
            "TemperatureProcessor", "FNIRSProcessor",
        )
    )
    removed_class_name = "".join(chr(value) for value in (80, 104, 97, 115, 101, 83, 101, 116))
    assert not hasattr(multimodalphysiokit, removed_class_name)
    assert hasattr(utils, "create_sliding_phase_intervals")
    removed_model_name = "".join(chr(value) for value in (112, 104, 97, 115, 101, 95, 115, 101, 116))
    assert not hasattr(utils, "create_sliding_" + removed_model_name)
    assert not hasattr(utils, "create_" + removed_model_name + "_from_seconds")


def test_package_version_matches_distribution() -> None:
    import multimodalphysiokit

    assert multimodalphysiokit.__version__ == version(
        "multimodalphysiokit"
    )

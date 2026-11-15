import matplotlib

matplotlib.use("Agg")

import numpy as np
import pytest

from multimodalphysiokit import ProcessingResult, Recording, Signal
from multimodalphysiokit.processors import PPGProcessor, process_ppg

FS = 100.0
AMPLITUDE = 1.0
RISE_TIME = 0.15

PULSE_AND_TIME_DOMAIN_NAMES = [
    "ppg_pulse_amplitude_mean", "ppg_pulse_amplitude_std", "ppg_pulse_amplitude_median",
    "ppg_pulse_duration_mean", "ppg_pulse_duration_std", "ppg_pulse_duration_median",
    "ppg_rise_time_mean", "ppg_rise_time_std", "ppg_rise_time_median",
    "ppg_bpm_mean", "ppg_bpm_std", "ppg_bpm_median", "ppg_pnn50", "ppg_pnn50_max",
]
FREQUENCY_NAMES = [
    "ppg_plf_mean", "ppg_plf_std", "ppg_phf_mean", "ppg_phf_std",
    "ppg_lf_hf_mean", "ppg_lf_hf_std",
]


def _pulse(n_samples, rise_samples):
    """Foot at 0, smooth rise to AMPLITUDE at rise_samples, decay back to 0."""
    tau = np.arange(n_samples, dtype=float)
    rise = AMPLITUDE * (0.5 - 0.5 * np.cos(np.pi * tau / rise_samples))
    decay = AMPLITUDE * (1 - (tau - rise_samples) / (n_samples - rise_samples)) ** 2
    return np.where(tau < rise_samples, rise, decay)


def make_ppg(ibis_seconds, *, rise_times=None, label="ppg"):
    """Synthetic PPG with known feet, peaks, amplitudes, and inter-beat intervals."""
    lengths = np.round(np.asarray(ibis_seconds) * FS).astype(int)
    rises = np.full(lengths.size, RISE_TIME) if rise_times is None else np.asarray(rise_times)
    samples = np.concatenate(
        [_pulse(n, int(round(r * FS))) for n, r in zip(lengths, rises)] + [np.zeros(1)]
    )
    time = np.arange(samples.size) / FS
    return Signal(samples, time, FS, label=label, units="a.u.")


def rsa_ibis(duration_seconds, rsa_frequency=0.25):
    ibis, elapsed = [], 0.0
    while elapsed < duration_seconds:
        ibis.append(0.9 + 0.05 * np.sin(2 * np.pi * rsa_frequency * elapsed))
        elapsed += ibis[-1]
    return ibis


def test_constant_rate_recovers_known_pulse_shape_and_heart_rate():
    signal = make_ppg([0.8] * 150)  # 75 bpm, 120 s

    features = process_ppg(signal, label_frequency_analysis=None).features

    assert features["ppg_bpm_mean"] == pytest.approx(75.0, abs=1e-6)
    assert features["ppg_bpm_std"] == pytest.approx(0.0, abs=1e-6)
    # The 0.5-8 Hz cleaning trims amplitude slightly (~1 %) and rounds the sharp
    # synthetic foot, which moves it up to 2 samples earlier.
    assert features["ppg_pulse_amplitude_mean"] == pytest.approx(AMPLITUDE, rel=0.02)
    assert features["ppg_rise_time_mean"] == pytest.approx(RISE_TIME, abs=2 / FS)
    assert features["ppg_pulse_duration_mean"] == pytest.approx(0.8, abs=1 / FS)
    assert features["ppg_pnn50"] == 0.0


def test_detects_heart_rates_above_90_bpm():
    # 150 bpm (2.5 Hz) lies outside a 0.5-1.5 Hz detection band (max 90 bpm),
    # where the attenuated beat loses to small in-band interference (here a
    # 1 Hz component at 10 % of the pulse amplitude, e.g. motion).
    clean = make_ppg([0.4] * 300)
    signal = Signal(
        clean.samples + 0.1 * AMPLITUDE * np.sin(2 * np.pi * 1.0 * clean.time),
        clean.time,
        FS,
    )

    features = process_ppg(signal, label_frequency_analysis=None).features

    assert features["ppg_bpm_mean"] == pytest.approx(150.0, abs=1e-6)


def test_pnn50_of_alternating_intervals_is_one():
    # Every successive difference is 100 ms > 50 ms.
    signal = make_ppg([0.8, 0.9] * 70)

    features = process_ppg(signal, label_frequency_analysis=None).features

    assert features["ppg_pnn50"] == pytest.approx(1.0)


def test_malformed_pulse_is_rejected_and_its_intervals_are_skipped():
    ibis = [0.8] * 100
    rises = [RISE_TIME] * 100
    rises[50] = 0.6  # systolic peak in the second half: malformed shape
    signal = make_ppg(ibis, rise_times=rises)

    result = process_ppg(signal, label_frequency_analysis=None, return_intermediates=True)

    valid = result.intermediates["pulse_valid"]
    onset_times = signal.time[result.intermediates["pulse_onset_indices"][:-1]]
    rejected_onsets = onset_times[~valid]
    np.testing.assert_allclose(rejected_onsets, [50 * 0.8], atol=0.05)
    # The two intervals touching the rejected pulse are skipped, not imputed.
    assert result.metadata["n_accepted_ibis"] == valid.size - 1 - 2
    # Neighbouring beats may shift by one sample (band-pass spreads the bad pulse).
    assert result.features["ppg_bpm_std"] < 0.5
    # The 8 Hz low-pass rounds the sharp synthetic foot by up to 2 samples.
    assert result.features["ppg_rise_time_mean"] == pytest.approx(RISE_TIME, abs=2 / FS)


def test_short_real_time_window_without_frequency_analysis():
    signal = make_ppg([0.8] * 13)  # ~10 s window

    result = process_ppg(signal, label_frequency_analysis=None)

    assert list(result.features) == PULSE_AND_TIME_DOMAIN_NAMES
    with pytest.raises(ValueError, match="label_frequency_analysis=None"):
        process_ppg(signal, label_frequency_analysis=1)


@pytest.mark.parametrize("mode", [1, 2])
def test_respiratory_sinus_arrhythmia_is_high_frequency_power(mode):
    # IBIs modulated at 0.25 Hz (breathing) must show up in the HF band.
    signal = make_ppg(rsa_ibis(300))

    features = process_ppg(signal, label_frequency_analysis=mode).features

    assert list(features) == PULSE_AND_TIME_DOMAIN_NAMES + FREQUENCY_NAMES
    assert features["ppg_phf_mean"] > 0.8
    assert features["ppg_lf_hf_mean"] < 0.25


def test_result_metadata_and_intermediates():
    signal = make_ppg([0.8] * 100, label="")

    result = PPGProcessor(label_frequency_analysis=None, return_intermediates=True).process(signal)

    assert isinstance(result, ProcessingResult)
    assert result.modality == "ppg"
    assert result.signal_labels == ["ppg"]
    assert result.metadata["n_pulses"] == result.intermediates["pulse_valid"].size
    n_pulses = result.metadata["n_pulses"]
    assert result.intermediates["pulse_onset_indices"].size == n_pulses + 1
    for key in ("pulse_peak_indices", "pulse_amplitudes", "pulse_durations", "pulse_rise_times"):
        assert result.intermediates[key].size == n_pulses
    assert result.intermediates["filtered_signal"].size == signal.n_samples
    assert result.intermediates["ibi"].size == result.intermediates["bpm"].size
    assert process_ppg(signal, label_frequency_analysis=None).intermediates == {}


def test_save_plots_writes_diagnostic_figures(tmp_path):
    signal = make_ppg(rsa_ibis(70))

    process_ppg(signal, label_frequency_analysis=1, save_plots=True, output_dir=tmp_path)

    for name in ("ppg_pulses.png", "ppg_BPM_trend.png", "ppg_STFT_HRV.png", "ppg_STFT_LFHF.png"):
        assert (tmp_path / name).is_file(), name


def test_recording_phase_processing_uses_the_processor():
    signal = make_ppg([0.8] * 150)
    recording = Recording(signals={"ppg": signal}, phase_intervals_seconds=[[0, 60], [60, 119]])

    results = recording.process_phases(PPGProcessor(label_frequency_analysis=None), "ppg")

    assert set(results) == {"phase_1", "phase_2"}
    for result in results.values():
        assert result.features["ppg_bpm_mean"] == pytest.approx(75.0, abs=1e-6)


def test_invalid_inputs_raise_clear_errors():
    with pytest.raises(ValueError, match="None, 1, or 2"):
        PPGProcessor(label_frequency_analysis=3)
    with pytest.raises(TypeError, match="Signal"):
        process_ppg(np.zeros(100))
    with pytest.raises(ValueError, match="output_dir"):
        process_ppg(make_ppg([0.8] * 20), label_frequency_analysis=None, save_plots=True)
    slow = Signal(np.zeros(100), np.arange(100) / 5.0, 5.0)
    with pytest.raises(ValueError, match="16 Hz"):
        process_ppg(slow, label_frequency_analysis=None)
    flat = Signal(np.zeros(2000), np.arange(2000) / FS, FS)
    with pytest.raises(ValueError, match="too few"):
        process_ppg(flat, label_frequency_analysis=None)


def test_pulse_shape_is_measured_on_drift_free_cleaned_signal():
    # Wrist PPG often rides on a slow baseline drift much larger than the
    # pulse (here 0.05 Hz, 8x the pulse amplitude). Pulse characteristics
    # must still reflect the pulse itself: all pulses valid, rise time and
    # amplitude close to their true values.
    clean = make_ppg([0.8] * 150)
    drift = 8 * AMPLITUDE * np.sin(2 * np.pi * 0.05 * clean.time)
    signal = Signal(clean.samples + drift, clean.time, FS)

    result = process_ppg(signal, label_frequency_analysis=None, return_intermediates=True)

    assert result.metadata["n_valid_pulses"] >= result.metadata["n_pulses"] - 2
    assert result.features["ppg_rise_time_mean"] == pytest.approx(RISE_TIME, abs=0.02)
    assert result.features["ppg_pulse_amplitude_mean"] == pytest.approx(AMPLITUDE, rel=0.05)
    assert result.intermediates["cleaned_signal"].size == signal.n_samples


def test_ibis_are_measured_between_systolic_peaks():
    # Feet every 0.8 s, rise time alternating 0.15 / 0.25 s: the systolic peaks
    # are therefore 0.9 / 0.7 s apart. IBIs must be peak-to-peak intervals.
    ibis = [0.8] * 150
    rises = [0.15, 0.25] * 75
    signal = make_ppg(ibis, rise_times=rises)

    result = process_ppg(signal, label_frequency_analysis=None, return_intermediates=True)

    peak_times = signal.time[result.intermediates["pulse_peak_indices"]]
    expected = np.diff(peak_times)[: result.intermediates["ibi"].size]
    np.testing.assert_allclose(result.intermediates["ibi"], expected, atol=1e-9)
    # The 8 Hz low-pass can move these asymmetric peaks by one sample (10 ms).
    ibi = result.intermediates["ibi"]
    np.testing.assert_allclose(ibi[::2], ibi[0], atol=1e-9)
    np.testing.assert_allclose(sorted([ibi[0], ibi[1]]), [0.7, 0.9], atol=1 / FS + 1e-9)
    assert result.features["ppg_pnn50"] == pytest.approx(1.0)   # every difference is 200 ms


def _symmetric_ppg(n_beats=150, ibi=0.8, peak_fraction=0.55):
    """Rounded, well-formed pulses whose peak sits at 55 % of the pulse."""
    n = int(ibi * FS)
    r = int(peak_fraction * n)
    tau = np.arange(n, dtype=float)
    pulse = np.where(tau < r, 0.5 - 0.5 * np.cos(np.pi * tau / r),
                     0.5 + 0.5 * np.cos(np.pi * (tau - r) / (n - r)))
    samples = np.concatenate([pulse] * n_beats + [np.zeros(1)])
    return Signal(samples, np.arange(samples.size) / FS, FS)


def test_max_peak_position_controls_the_late_peak_rule():
    signal = _symmetric_ppg()

    with pytest.raises(ValueError, match="too few valid"):
        process_ppg(signal, label_frequency_analysis=None)  # default 0.5 rejects them

    result = process_ppg(signal, label_frequency_analysis=None, max_peak_position=0.65)
    assert result.metadata["n_valid_pulses"] == result.metadata["n_pulses"]
    # Rounded peaks can land one sample apart, so IBIs are exact to one sample.
    assert result.features["ppg_bpm_mean"] == pytest.approx(75.0, abs=0.1)


def _two_hump_ppg(first, second, first_at, second_at, n_beats=150, ibi=0.8):
    """Pulses made of two Gaussian humps (width 50 ms) at known times after the beat start."""
    tau = np.arange(int(ibi * FS)) / FS
    pulse = (first * np.exp(-0.5 * ((tau - first_at) / 0.05) ** 2)
             + second * np.exp(-0.5 * ((tau - second_at) / 0.05) ** 2))
    samples = np.concatenate([pulse] * n_beats + [np.zeros(1)])
    return Signal(samples, np.arange(samples.size) / FS, FS)


def test_dicrotic_wave_does_not_split_the_pulse():
    # Systolic wave at 0.15 s and a dicrotic wave half as tall at 0.45 s, one
    # beat every 0.8 s (75 bpm). The notch between them is not a pulse foot:
    # every beat must be one accepted pulse and every IBI 0.8 s.
    signal = _two_hump_ppg(1.0, 0.5, 0.15, 0.45)

    result = process_ppg(signal, label_frequency_analysis=None, return_intermediates=True)

    assert result.metadata["n_pulses"] >= 145
    assert result.metadata["n_valid_pulses"] >= result.metadata["n_pulses"] - 2
    np.testing.assert_allclose(result.intermediates["ibi"], 0.8, atol=0.011)


def test_systolic_peak_is_the_first_hump_even_when_a_later_one_is_taller():
    # Wrist-like pulse: systolic hump (0.85) at 0.15 s, taller late hump (1.0)
    # at 0.32 s. The systolic peak is the first hump, not the global maximum.
    signal = _two_hump_ppg(0.85, 1.0, 0.15, 0.32)

    result = process_ppg(signal, label_frequency_analysis=None, return_intermediates=True)

    peak_times = signal.time[result.intermediates["pulse_peak_indices"]]
    np.testing.assert_allclose(np.mod(peak_times, 0.8), 0.15, atol=0.02)


def test_cleaning_removes_content_above_8_hz():
    # A 30 Hz interference must not reach the cleaned signal on which pulse
    # characteristics are measured.
    clean = make_ppg([0.8] * 150)
    noisy = Signal(clean.samples + 0.05 * np.sin(2 * np.pi * 30 * clean.time), clean.time, FS)

    reference = process_ppg(clean, label_frequency_analysis=None, return_intermediates=True)
    result = process_ppg(noisy, label_frequency_analysis=None, return_intermediates=True)

    difference = result.intermediates["cleaned_signal"] - reference.intermediates["cleaned_signal"]
    assert np.max(np.abs(difference[200:-200])) < 0.005


@pytest.mark.parametrize("value", [0.0, 1.5, "0.5", float("nan")])
def test_invalid_max_peak_position_is_rejected(value):
    with pytest.raises((ValueError, TypeError), match="max_peak_position"):
        PPGProcessor(max_peak_position=value)

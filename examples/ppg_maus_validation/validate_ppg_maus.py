"""Validate the PPG processor on the public MAUS dataset against ECG and NeuroKit2.

MAUS (Beh et al., 2021) records ECG and fingertip PPG (Infiniti, 256 Hz, same
device and clock) plus wrist PPG (PixArt, 100 Hz) for 22 participants: one
rest segment and six N-back trials each. This script does not ship any data;
download MAUS and pass the folder that holds the participant directories.

Every segment is cut into non-overlapping 90 s windows; each window is
processed on its own, as a completed real-time window would be.

Reference: ECG R-peaks from NeuroKit2 (``ecg_clean`` + ``ecg_peaks``, method
"neurokit", artifact correction). R-R intervals more than 20 % away from the
median of their ~10 neighbours are irregular (e.g. ectopic beats); they and
their two neighbouring intervals are excluded, leaving normal-to-normal (NN)
intervals in 40-180 bpm, as HRV standards require.

Methods: MultimodalPhysioKit ``PPGProcessor`` (defaults) and NeuroKit2
``ppg_clean`` + ``ppg_peaks`` (elgendi, elgendi with artifact correction,
charlton, charlton2024; bishop with ``--bishop``, very slow). NeuroKit2 does
not reject beats, so all of its intervals in 40-180 bpm are used.

Per window and method, with the definitions of the library features:
  bpm_mean  mean of 60 / IBI over accepted IBIs        (ppg_bpm_mean)
  pnn50     share of |successive IBI difference| > 50 ms
            over pairs of adjacent accepted IBIs          (ppg_pnn50)
and the same values on the ECG NN intervals. Beat level: a PPG beat matches
an R-peak when it lies within 150 ms of the R-peak plus the window's median
pulse arrival delay.

Wrist caveat: the PixArt and the ECG are separate devices; their alignment
within a window is assumed, not verified, so wrist results rank methods but
may overstate absolute errors.

Requirements: multimodalphysiokit, neurokit2 (tested with 0.3.0), pandas.
Run:
    python examples/ppg_maus_validation/validate_ppg_maus.py /path/to/MAUS/Data/Raw_data
Outputs (default validation-output/ppg_maus/, git-ignored):
    windows.csv     one row per window x device x method
    summary.md      the tables reported in docs/modalities/ppg.md
"""
import argparse
import time
import warnings
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import neurokit2 as nk
import numpy as np
import pandas as pd

from multimodalphysiokit import Signal
from multimodalphysiokit.processors import PPGProcessor

warnings.filterwarnings("ignore")
WINDOW_S = 90.0
MATCH_S = 0.15
ECG_FS = 256.0
DEVICES = {  # device: (sampling rate, trial file, rest file, rest column)
    "fingertip": (256.0, "inf_ppg.csv", "inf_resting.csv", "Resting_PPG"),
    "wrist": (100.0, "pixart.csv", "pixart_resting.csv", "Resting"),
}
NK_METHODS = {"NK2 elgendi": ("elgendi", False), "NK2 elgendi_corrected": ("elgendi", True),
              "NK2 charlton": ("charlton", False), "NK2 charlton2024": ("charlton2024", False)}


# ------------------------------------------------------------------ loading
def trimmed(x):
    x = np.asarray(x, float)
    finite = np.flatnonzero(np.isfinite(x))
    if finite.size == 0:
        return None
    x = x[: finite[-1] + 1]
    return x if np.all(np.isfinite(x)) else None


def segments(folder: Path):
    """Yield (segment, ecg, {device: ppg}) for the rest segment and each trial."""
    ecg_trials = pd.read_csv(folder / "inf_ecg.csv")
    rest = {"ecg": pd.read_csv(folder / "inf_resting.csv")["Resting_ECG"]}
    trials = {}
    for device, (_, trial_file, rest_file, rest_column) in DEVICES.items():
        rest[device] = pd.read_csv(folder / rest_file)[rest_column]
        trials[device] = pd.read_csv(folder / trial_file)
    yield "rest", trimmed(rest["ecg"]), {d: trimmed(rest[d]) for d in DEVICES}
    for i in range(ecg_trials.shape[1]):
        yield (f"trial{i + 1}", trimmed(ecg_trials.iloc[:, i]),
               {d: trimmed(trials[d].iloc[:, i]) for d in DEVICES})


# --------------------------------------------------------------- reference
def ecg_reference(ecg):
    cleaned = nk.ecg_clean(ecg, sampling_rate=int(ECG_FS))
    _, info = nk.ecg_peaks(cleaned, sampling_rate=int(ECG_FS), method="neurokit", correct_artifacts=True)
    r = np.asarray(info["ECG_R_Peaks"], float) / ECG_FS
    rr = np.diff(r)
    local = np.array([np.median(rr[max(0, i - 5):i + 6]) for i in range(rr.size)])
    irregular = np.abs(rr / local - 1) > 0.2
    excluded = irregular | np.r_[False, irregular[:-1]] | np.r_[irregular[1:], False]
    return r, ~excluded & (rr >= 60 / 180) & (rr <= 60 / 40)


# ------------------------------------------------------------------ methods
def mpk_beats(x, fs):
    signal = Signal(x, np.arange(x.size) / fs, fs)
    result = PPGProcessor(label_frequency_analysis=None, return_intermediates=True).process(signal)
    beats = signal.time[result.intermediates["pulse_peak_indices"]]
    valid = result.intermediates["pulse_valid"]
    ibi = np.diff(beats)
    usable = valid[:-1] & valid[1:] & (ibi >= 60 / 180) & (ibi <= 60 / 40)
    return beats, usable, result.features


def nk_beats(x, fs, method, correct):
    fs = int(round(fs))
    cleaned = nk.ppg_clean(x, sampling_rate=fs, method="elgendi")
    _, info = nk.ppg_peaks(cleaned, sampling_rate=fs, method=method, correct_artifacts=correct)
    beats = np.asarray(info["PPG_Peaks"], float) / fs
    ibi = np.diff(beats)
    return beats, (ibi >= 60 / 180) & (ibi <= 60 / 40), None


# ------------------------------------------------------------------ metrics
def features(beats, usable):
    """bpm_mean and pnn50 with the definitions of ppg_bpm_mean and ppg_pnn50."""
    ibi = np.diff(beats)
    if usable.sum() < 3:
        return np.nan, np.nan
    pairs = usable[:-1] & usable[1:]
    differences = np.abs(np.diff(ibi))[pairs]
    pnn50 = float(np.mean(differences > 0.05)) if differences.size else 0.0   # as ppg_pnn50
    return float(np.mean(60 / ibi[usable])), pnn50


def beat_level(r, nn, beats, usable):
    """NN intervals kept, accepted intervals not matching the ECG, |IBI - RR| (ms)."""
    if beats.size < 3:
        return 0, 0, []
    previous = np.searchsorted(r, beats) - 1
    delay = np.median(beats[previous >= 0] - r[previous[previous >= 0]])
    match = np.full(r.size, -1)
    for i, expected in enumerate(r + delay):
        j = int(np.argmin(np.abs(beats - expected)))
        if abs(beats[j] - expected) <= MATCH_S:
            match[i] = j
    true_pairs, kept, errors = set(), 0, []
    for i in range(r.size - 1):
        a, b = match[i], match[i + 1]
        if a >= 0 and b == a + 1:
            true_pairs.add(a)
            if nn[i] and usable[a]:
                kept += 1
                errors.append(1000 * abs((beats[b] - beats[a]) - (r[i + 1] - r[i])))
    wrong = sum(1 for a in np.flatnonzero(usable) if a not in true_pairs)
    return kept, wrong, errors


def run_participant(args):
    folder, use_bishop = args
    methods = {"mpk": None, **NK_METHODS}
    if use_bishop:
        methods["NK2 bishop"] = ("bishop", False)
    rows = []
    for segment, ecg, ppg in segments(folder):
        if ecg is None:
            continue
        for w in range(int(ecg.size / ECG_FS // WINDOW_S)):
            lo = w * WINDOW_S
            try:
                r, nn = ecg_reference(ecg[int(lo * ECG_FS):int((lo + WINDOW_S) * ECG_FS)])
            except Exception:
                continue
            ref_bpm, ref_pnn50 = features(r, nn)
            for device, (fs, *_) in DEVICES.items():
                x = ppg[device]
                if x is None or x.size < (lo + WINDOW_S) * fs:
                    continue
                x = x[int(lo * fs):int((lo + WINDOW_S) * fs)]
                base = dict(participant=folder.name, segment=segment, window=w + 1, device=device,
                            ecg_bpm=ref_bpm, ecg_pnn50=ref_pnn50, n_nn=int(nn.sum()))
                for name, spec in methods.items():
                    start = time.perf_counter()
                    try:
                        beats, usable, lib = mpk_beats(x, fs) if spec is None else nk_beats(x, fs, *spec)
                    except Exception as error:
                        rows.append({**base, "method": name, "failed": True, "error": str(error)[:80]})
                        continue
                    runtime = 1000 * (time.perf_counter() - start)
                    bpm, pnn50 = features(beats, usable)
                    if lib is not None:  # the library's own features must match these definitions
                        assert np.isclose(bpm, lib["ppg_bpm_mean"], equal_nan=True)
                        assert np.isclose(pnn50, lib["ppg_pnn50"], equal_nan=True)
                    kept, wrong, errors = beat_level(r, nn, beats, usable)
                    rows.append({**base, "method": name, "failed": False, "runtime_ms": runtime,
                                 "bpm": bpm, "pnn50": pnn50, "nn_kept": kept,
                                 "n_accepted": int(usable.sum()), "n_wrong": wrong,
                                 "ibi_errors_ms": " ".join(f"{e:.2f}" for e in errors)})
    return rows


# ------------------------------------------------------------------ summary
def bland_altman(d):
    d = d.dropna()
    return d.mean(), d.mean() - 1.96 * d.std(), d.mean() + 1.96 * d.std()


def summarize(table: pd.DataFrame) -> str:
    out = ["# MAUS PPG validation\n"]
    for device, g in table.groupby("device", sort=False):
        ok = g[~g.failed]
        n_windows = g.groupby("method").size().max()
        out.append(f"\n## {device} ({n_windows} windows of {WINDOW_S:.0f} s)\n")
        beat, hr, pnn = [], [], []
        for method, m in ok.groupby("method", sort=False):
            errors = np.array([float(e) for s in m.ibi_errors_ms.dropna() for e in str(s).split()])
            failed = int(g[(g.method == method)].failed.sum())
            beat.append(f"| {method} | {failed} | {100 * m.nn_kept.sum() / m.n_nn.sum():.1f} % | "
                        f"{100 * m.n_wrong.sum() / max(m.n_accepted.sum(), 1):.1f} % | "
                        f"{np.median(errors):.1f} | {np.percentile(errors, 95):.1f} | "
                        f"{m.runtime_ms.median():.1f} |")
            dh = m.bpm - m.ecg_bpm
            bias, low, high = bland_altman(dh)
            a = dh.abs().dropna()
            hr.append(f"| {method} | {a.median():.2f} ({a.quantile(.25):.2f}-{a.quantile(.75):.2f}) | "
                      f"{bias:+.2f} [{low:+.2f}, {high:+.2f}] | {100 * (a <= 2).mean():.1f} % |")
            dp = m.pnn50 - m.ecg_pnn50
            bias, low, high = bland_altman(dp)
            a = dp.abs().dropna()
            pnn.append(f"| {method} | {a.median():.3f} ({a.quantile(.25):.3f}-{a.quantile(.75):.3f}) | "
                       f"{bias:+.3f} [{low:+.3f}, {high:+.3f}] | {100 * (a <= 0.05).mean():.1f} % |")
        out += ["### Beats\n",
                "| Method | Failed windows | NN intervals kept | Accepted intervals not matching ECG | "
                "Median IBI error (ms) | 95th pct IBI error (ms) | Runtime per window (ms) |",
                "|---|---|---|---|---|---|---|", *beat,
                "\n### Heart rate (`ppg_bpm_mean` vs ECG NN)\n",
                "| Method | Median abs. error, bpm (IQR) | Bias [95 % limits of agreement], bpm | Windows within 2 bpm |",
                "|---|---|---|---|", *hr,
                "\n### pNN50 (`ppg_pnn50` vs ECG NN)\n",
                "| Method | Median abs. error (IQR) | Bias [95 % limits of agreement] | Windows within 0.05 |",
                "|---|---|---|---|", *pnn]
    return "\n".join(out) + "\n"


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("maus_root", type=Path, help="folder containing the MAUS participant folders")
    parser.add_argument("--output-dir", type=Path, default=Path("validation-output/ppg_maus"))
    parser.add_argument("--jobs", type=int, default=1, help="parallel participants (runtimes assume 1)")
    parser.add_argument("--bishop", action="store_true", help="also run NeuroKit2 bishop (very slow)")
    args = parser.parse_args()
    folders = sorted(p for p in args.maus_root.iterdir() if (p / "inf_ecg.csv").exists())
    with ProcessPoolExecutor(args.jobs) as pool:
        rows = [row for part in pool.map(run_participant, [(f, args.bishop) for f in folders]) for row in part]
    table = pd.DataFrame(rows)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    table.to_csv(args.output_dir / "windows.csv", index=False)
    report = summarize(table)
    (args.output_dir / "summary.md").write_text(report)
    print(report)


if __name__ == "__main__":
    main()

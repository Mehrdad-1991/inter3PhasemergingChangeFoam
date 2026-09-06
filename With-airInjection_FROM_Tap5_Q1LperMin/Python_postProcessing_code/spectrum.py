#!/usr/bin/env python3
################################################################################
# MIT License
# 
# Copyright (c) 2026 - Universität Rostock 
#                       Fakultät für Maschinenbau und Schiffstechnik
#                       Lehrstuhl für Modellierung und Simulation
#                       
# Author: Mehrdad Kazemi
# Date: Dec 2025 (Improved low-frequency quality)
# Address: Albert-Einstein-Str. 2, 18059 Rostock, Deutschland
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
################################################################################

"""
Author : Mehrdad
Date   : Dec 2025 (Improved low-frequency quality)

Key improvements vs. original:
1) Uses ONLY NumPy FFT (no SciPy -> avoids your SciPy/NumPy mismatch).
2) Selects a *steady* time window (discard initial transient).
3) Uses Welch-style averaging (segmenting + overlap) for smoother low-frequency spectra.
4) Keeps frequency resolution controlled by segment length: df = 1/(Nseg*dt).
5) Exports comparison CSV in the SAME paired format: [Frequency, case1, Frequency, case2, ...]
6) Produces comparison PDFs per tap as before.

Notes:
- True low-frequency bin spacing is df = 1/Tseg. To improve df you must increase Tseg.
- Welch reduces variance and makes low-frequency peaks more stable.
"""

import csv
import shutil
from pathlib import Path
from datetime import datetime

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rc('figure', max_open_warning=0)

# -----------------------------
# Global config
# -----------------------------
FIXED_DT = 5e-05                  # sampling interval (s)
FS = 1.0 / FIXED_DT               # sampling rate (Hz) = 20000 Hz
FREQ_MAX = 10000                  # Nyquist = 10000 Hz for dt=5e-05

# --- Low-frequency improvement knobs ---
DISCARD_FRACTION = 0.3            # discard first 30% of samples as transient
REMOVE_MEAN = True                # remove mean (DC)
DETREND_LINEAR = False            # optional linear detrend (slow drift)

# Welch parameters (choose segment length by desired df)
DESIRED_DF_HZ = 0.25              # e.g. 0.25 Hz -> Tseg = 4 s
OVERLAP_FRACTION = 0.5            # 50% overlap
WINDOW = "hann"                   # Hann window

# --- Plot styling ---
PLOT_FMIN = 10.0
PLOT_FMAX_ALL = 10000.0

# === USER CONFIG ===
BASE_DIR = Path("..")
LOG_FILE = BASE_DIR / "spectrum_log.txt"
COMPARISON_DIR = BASE_DIR / "Pressure_Spectrum_Comparison"

POINTS = range(1, 12)
LINEWIDTH = [0.6]
ALPHA = [0.95]


# -----------------------------
# Logging and I/O
# -----------------------------
def log(msg: str) -> None:
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    with open(LOG_FILE, 'a') as f:
        f.write(f"[{timestamp}] {msg}\n")


def clean_output_dirs() -> None:
    for subdir in ["pressurespectrumsFFT", "pressurespectrumsPSD",
                   "velocityspectrumsFFT", "velocityspectrumsPSD"]:
        full_path = BASE_DIR / "postProcessing" / subdir
        shutil.rmtree(full_path, ignore_errors=True)

    shutil.rmtree(COMPARISON_DIR, ignore_errors=True)
    COMPARISON_DIR.mkdir(parents=True, exist_ok=True)


def save_csv(rows, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, 'w', newline='') as f:
        csv.writer(f).writerows(rows)


def read_time_series(pt: int):
    pp = BASE_DIR / f"postProcessing/tap{pt}/tap{pt}pressure.dat"
    log(f"Reading: {pp}")
    if not pp.is_file():
        log(f"[WARN] Missing pressure data at tap{pt}")
        return None

    t_list, p_list = [], []
    with open(pp) as f:
        for L in f:
            if L.startswith('#') or L.startswith('['):
                continue
            parts = L.replace('(', ' ').replace(')', ' ').split()
            if len(parts) < 2:
                continue
            t = float(parts[0])
            p = float(parts[1])

            # keep your conversion exactly
            p = p * 1e-5 - 1

            t_list.append(t)
            p_list.append(p)

    t = np.asarray(t_list, dtype=float)
    y = np.asarray(p_list, dtype=float)

    if t.size < 10:
        return None
    return t, y


# -----------------------------
# Signal processing
# -----------------------------
def detrend_linear(y: np.ndarray) -> np.ndarray:
    n = np.arange(y.size, dtype=float)
    A = np.vstack([n, np.ones_like(n)]).T
    coeff, _, _, _ = np.linalg.lstsq(A, y, rcond=None)
    return y - (A @ coeff)


def get_window(N: int) -> np.ndarray:
    if WINDOW.lower() == "hann":
        return np.hanning(N)
    raise ValueError(f"Unsupported window: {WINDOW}")


def select_steady_window(t: np.ndarray, y: np.ndarray):
    """Discard initial transient samples by fraction."""
    N = y.size
    i0 = int(np.floor(DISCARD_FRACTION * N))
    i0 = min(max(i0, 0), N - 2)
    return t[i0:], y[i0:]


def welch_average_spectrum(y: np.ndarray, dt: float):
    """
    Welch-like averaged *magnitude spectrum* (not PSD):
    - Split into segments of length Nseg (from DESIRED_DF_HZ)
    - Apply Hann window
    - Compute rFFT magnitude, normalize by (Nseg * coherent_gain)
    - Average magnitudes across segments

    Returns:
      f (Hz), Pavg (averaged magnitude)
    """
    N = y.size
    # Segment length from desired df: df = 1/(Nseg*dt) => Nseg = 1/(df*dt)
    Nseg = int(round(1.0 / (DESIRED_DF_HZ * dt)))
    # Make Nseg even-ish and not larger than available data
    Nseg = max(256, min(Nseg, N))

    step = int(round(Nseg * (1.0 - OVERLAP_FRACTION)))
    step = max(1, step)

    if Nseg < 8 or N < Nseg:
        # fallback: single segment FFT
        w = np.hanning(N)
        cg = np.mean(w)
        F = np.fft.rfft(y * w)
        f = np.fft.rfftfreq(N, d=dt)
        P = np.abs(F) / (N * max(cg, 1e-12))
        return f, P

    w = get_window(Nseg)
    cg = max(np.mean(w), 1e-12)

    specs = []
    nseg_count = 0
    for start in range(0, N - Nseg + 1, step):
        seg = y[start:start + Nseg]

        if REMOVE_MEAN:
            seg = seg - np.mean(seg)
        if DETREND_LINEAR:
            seg = detrend_linear(seg)

        F = np.fft.rfft(seg * w)
        P = np.abs(F) / (Nseg * cg)
        specs.append(P)
        nseg_count += 1

    f = np.fft.rfftfreq(Nseg, d=dt)
    Pavg = np.mean(np.vstack(specs), axis=0)

    return f, Pavg


def compute_spectrum(t: np.ndarray, y: np.ndarray, dt: float = FIXED_DT):
    # check dt quality (informational)
    d = np.diff(t)
    dt_mean = float(np.mean(d))
    dt_std = float(np.std(d))
    dt_maxdev = float(np.max(np.abs(d - dt)))
    log(f"[INFO] dt_mean={dt_mean:.8e}, dt_std={dt_std:.3e}, max|dt-dt0|={dt_maxdev:.3e}")

    # select steady window
    t2, y2 = select_steady_window(t, y)
    N2 = y2.size
    T2 = N2 * dt
    log(f"[INFO] After discard: N={N2}, T={T2:.6f} s")

    # Welch-averaged magnitude spectrum
    f, P = welch_average_spectrum(y2, dt=dt)

    # small floor for log plots
    P[P <= 0] = np.finfo(float).tiny
    return f, P


# -----------------------------
# Pipeline
# -----------------------------
def process_individual() -> None:
    out_p = BASE_DIR / "postProcessing/pressurespectrumsFFT"
    out_p.mkdir(parents=True, exist_ok=True)

    for pt in POINTS:
        data = read_time_series(pt)
        if data is None:
            log(f"[WARN] No usable time series for tap{pt}")
            continue

        t, y = data
        f, P = compute_spectrum(t, y, dt=FIXED_DT)

        if f.size == 0 or P.size == 0:
            log(f"[SKIP] Empty spectrum for tap{pt}")
            continue

        # Save spectrum to .dat (Frequency, Magnitude)
        save_csv([[float(fv), float(pv)] for fv, pv in zip(f, P)],
                 out_p / f"Spectrum_tap{pt}.dat")

        # Plot (log-log)
        fig, ax = plt.subplots(figsize=(12, 8))
        ax.plot(f, P, '-k', lw=0.8)
        ax.set_title(f"tap{pt} (Welch-avg FFT magnitude)")
        ax.set_xlabel("Frequency [Hz]")
        ax.set_ylabel("FFT magnitude (normalized)")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.grid(which="both", ls=":")
        ax.set_xlim(PLOT_FMIN, min(FREQ_MAX, float(f[-1])))
        fig.savefig(out_p / f"pressure{pt}.pdf")
        plt.close(fig)


def load_spectrum_file(pt: int):
    path = BASE_DIR / "postProcessing/pressurespectrumsFFT" / f"Spectrum_tap{pt}.dat"
    if not path.is_file():
        return [], []
    freqs, specs = [], []
    with open(path) as f:
        for L in f:
            if L.startswith('#') or L.startswith('['):
                continue
            parts = L.replace('(', ' ').replace(')', ' ').replace(',', ' ').split()
            if len(parts) < 2:
                continue
            freqs.append(float(parts[0]))
            specs.append(float(parts[1]))
    return freqs, specs


def plot_single(pt: int) -> None:
    f, P = load_spectrum_file(pt)
    if not f:
        log(f"[INFO] No data for tap{pt}")
        return

    # build rows: Frequency, Magnitude
    rows = [["Frequency", "Magnitude"]]
    for fv, pv in zip(f, P):
        rows.append([fv, pv])

    fig, ax = plt.subplots(figsize=(12, 8))
    ax.plot(f, P, lw=0.6, alpha=0.95)

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Frequency [Hz]")
    ax.set_ylabel("FFT magnitude (normalized)")
    ax.set_title(f"tap{pt} (Welch-avg magnitude)")
    ax.set_xlim(PLOT_FMIN, min(FREQ_MAX, float(f[-1])))
    ax.grid(which="both", ls=":")

    fig.savefig(COMPARISON_DIR / f"FFT_tap{pt}.pdf")
    save_csv(rows, COMPARISON_DIR / f"FFT_tap{pt}.csv")
    plt.close(fig)


def main():
    LOG_FILE.unlink(missing_ok=True)
    log("=== Spectrum Processing Started (Welch + steady-window) ===")

    clean_output_dirs()

    # Log expected segment info
    Tseg = 1.0 / DESIRED_DF_HZ
    Nseg = int(round(Tseg / FIXED_DT))
    log(f"[CONFIG] dt={FIXED_DT}, fs={FS}, desired_df={DESIRED_DF_HZ} Hz -> Tseg={Tseg:.3f} s, Nseg={Nseg}")

    process_individual()

    for pt in POINTS:
        plot_single(pt)

    log("=== Spectrum Processing Finished ===")


if __name__ == "__main__":
    main()


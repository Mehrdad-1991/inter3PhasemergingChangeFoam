#!/usr/bin/env python3
# -*- coding: utf-8 -*-
################################################################################
# MIT License
# 
# Copyright (c) 2026 - Universität Rostock 
#                       Fakultät für Maschinenbau und Schiffstechnik
#                       Lehrstuhl für Modellierung und Simulation
#                       
# Author: Mehrdad Kazemi
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

import re
import shutil
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# -----------------------------
# Settings
# -----------------------------
INPUT_DIR = Path("../Pressure_Spectrum_Comparison/")   # <-- results folder backward
OUT_DIR = INPUT_DIR / "figures_fft"
SUMMARY_CSV = INPUT_DIR / "fft_localmax_summary_10_20_10_30_20_30Hz.csv"

FREQ_MIN_ALL = 10
FREQ_MAX_ALL = 10000
FREQ_MAX_COMP = 2000

plt.rcParams.update({
    "figure.dpi": 150,
    "savefig.dpi": 600,
    "font.size": 11,
    "axes.labelsize": 12,
    "axes.titlesize": 12,
    "legend.fontsize": 10,
    "xtick.labelsize": 11,
    "ytick.labelsize": 11,
    "axes.grid": True,
    "grid.alpha": 0.25,
    "grid.linestyle": "-",
    "axes.linewidth": 1.0,
})


# -----------------------------
# Helpers
# -----------------------------
def sanitize_filename(s: str) -> str:
    s = re.sub(r"[^\w\-. ]+", "_", s)
    return s.strip().replace(" ", "_")[:200]


def classify_case(case_name: str) -> str:
    """
    User rule:
      1_*  -> without injection
      2_*  -> with injection
    """
    s = str(case_name).strip().lower()
    m = re.match(r"^(\d+)_", s)
    if m:
        return "without" if m.group(1) == "1" else "with" if m.group(1) == "2" else "without"

    # Fallbacks (rare)
    if "no_injection" in s or "without_injection" in s:
        return "without"
    if "with_injection" in s:
        return "with"

    return "without"


def extract_mesh_size(case_name: str):
    """
    Extract something like: mehsSize_0.49_million (typo included in your names)
    Returns float or None.
    """
    s = str(case_name)
    m = re.search(r"(?:meshsize|mehsSize|meshSize)_([0-9]*\.?[0-9]+)_million", s, flags=re.IGNORECASE)
    return float(m.group(1)) if m else None


def extract_pair_columns(df: pd.DataFrame):
    """
    Your CSV format: [Frequency, CaseA, Frequency, CaseB, ...]
    Returns list of dicts: {"case":..., "freq_col":..., "p_col":...}
    """
    cols = list(df.columns)
    pairs = []
    i = 0
    while i < len(cols) - 1:
        c_freq = cols[i]
        c_p = cols[i + 1]
        if str(c_freq).strip().lower().startswith("frequency"):
            pairs.append({"case": str(c_p), "freq_col": c_freq, "p_col": c_p})
            i += 2
        else:
            i += 1

    if not pairs:
        raise ValueError("Could not detect (Frequency, Case) column pairs in the CSV.")
    return pairs


def prepare_xy(df, freq_col, p_col, fmin, fmax):
    x = pd.to_numeric(df[freq_col], errors="coerce").to_numpy()
    y = pd.to_numeric(df[p_col], errors="coerce").to_numpy()
    m = np.isfinite(x) & np.isfinite(y) & (x > 0) & (y > 0) & (x >= fmin) & (x <= fmax)
    x, y = x[m], y[m]
    if x.size > 1:
        idx = np.argsort(x)
        x, y = x[idx], y[idx]
    return x, y


def local_max_in_band(x, y, f_lo, f_hi):
    m = (x >= f_lo) & (x <= f_hi)
    if not np.any(m):
        return np.nan, np.nan
    xx, yy = x[m], y[m]
    j = int(np.nanargmax(yy))
    return float(xx[j]), float(yy[j])


def plot_group(tap_name, title, curves, xlim, outpath: Path):
    fig = plt.figure(figsize=(7.2, 4.6))
    ax = fig.add_subplot(1, 1, 1)

    for c in curves:
        ax.loglog(c["x"], c["y"], linewidth=1.6, label=c["label"])

    ax.set_xlim(xlim[0], xlim[1])
    ax.set_xlabel("Frequency [Hz]")
    ax.set_ylabel("Pressure FFT [Pa]")
    ax.set_title(f"{tap_name}: {title}")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(loc="best", frameon=True, framealpha=0.9)

    fig.tight_layout()
    outpath.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(outpath, bbox_inches="tight")
    plt.close(fig)


# -----------------------------
# Main
# -----------------------------
def main():
    # --- Clean previous outputs (if any) ---
    if OUT_DIR.exists() and OUT_DIR.is_dir():
        shutil.rmtree(OUT_DIR)
    if SUMMARY_CSV.exists() and SUMMARY_CSV.is_file():
        SUMMARY_CSV.unlink()

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    csv_files = sorted(INPUT_DIR.glob("FFT_tap*.csv"))
    if not csv_files:
        raise FileNotFoundError(f"No files found: {INPUT_DIR / 'FFT_tap*.csv'}")

    for csv_path in csv_files:
        tap_name = csv_path.stem.replace("FFT_", "")
        print(f"Processing: {csv_path.name}")

        df = pd.read_csv(csv_path)
        
        # For single simulation: columns are [Frequency, Magnitude]
        if len(df.columns) >= 2:
            freq_col = df.columns[0]
            mag_col = df.columns[1]
            
            x, y = prepare_xy(df, freq_col, mag_col, FREQ_MIN_ALL, FREQ_MAX_ALL)
            
            if x.size > 1:
                # Plot spectrum
                fig, ax = plt.subplots(figsize=(12, 8))
                ax.loglog(x, y, linewidth=1.6, label="FFT Magnitude")
                ax.set_xlim(FREQ_MIN_ALL, FREQ_MAX_ALL)
                ax.set_xlabel("Frequency [Hz]")
                ax.set_ylabel("Pressure FFT [Pa]")
                ax.set_title(f"{tap_name}: FFT Spectrum (10–10000 Hz)")
                ax.grid(True, which="both", alpha=0.25)
                ax.legend(loc="best", frameon=True, framealpha=0.9)
                
                fig.tight_layout()
                out_file = OUT_DIR / f"{sanitize_filename(tap_name)}_spectrum_loglog_10_10000.png"
                fig.savefig(out_file, bbox_inches="tight")
                plt.close(fig)
                
                print(f"Saved plot to: {out_file}")
            else:
                print(f"[SKIP] No valid data for {tap_name}")
        else:
            print(f"[SKIP] Unexpected CSV format in {csv_path.name}")

    print(f"\nSaved figures to: {OUT_DIR.resolve()}")
    print(f"Saved summary CSV to: {SUMMARY_CSV.resolve()}")


if __name__ == "__main__":
    main()


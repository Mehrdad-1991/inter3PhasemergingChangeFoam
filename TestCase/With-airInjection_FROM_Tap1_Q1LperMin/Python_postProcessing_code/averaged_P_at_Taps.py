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

import numpy as np
from pathlib import Path
from datetime import datetime
import csv

# === USER CONFIG ===
BASE_DIR = Path("..")
LOG_FILE = BASE_DIR / "avg_pressure_per_tap_log.txt"
OUTPUT_FILE = BASE_DIR / "averaged_pressure_per_tap.csv"

POINTS = range(1, 12)

# Averaging window in seconds
AVERAGING_WINDOW = (0.06, 0.37)

def log(msg):
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    with open(LOG_FILE, 'a') as f:
        f.write(f"[{timestamp}] {msg}\n")

def read_and_average_tap(pt, t_start, t_end):
    pp = BASE_DIR / f"postProcessing/tap{pt}/tap{pt}pressure.dat"
    if not pp.is_file():
        log(f"[WARN] Missing: {pp}")
        return None
    ts_pp, spp = [], []
    with open(pp) as f:
        for L in f:
            if L.startswith('#') or L.startswith('['): continue
            parts = L.replace('(', ' ').replace(')', ' ').split()
            try:
                t, p = float(parts[0]), float(parts[1])
                ts_pp.append(t)
                spp.append(p)
            except:
                continue
    t = np.array(ts_pp)
    p = np.array(spp)
    mask = (t >= t_start) & (t <= t_end)
    if not np.any(mask):
        log(f"[WARN] Tap{pt} has no data in window.")
        return None
    return np.mean(p[mask])

def main():
    LOG_FILE.unlink(missing_ok=True)
    log("=== Per-Tap Averaging Started ===")

    rows = [["Tap", "AveragePressure"]]

    for pt in POINTS:
        avg = read_and_average_tap(pt, *AVERAGING_WINDOW)
        if avg is not None:
            rows.append([f"tap{pt}", avg])
            log(f"[OK] tap{pt} → avg = {avg:.6f}")
        else:
            rows.append([f"tap{pt}", "NaN"])

    with open(OUTPUT_FILE, 'w', newline='') as f:
        csv.writer(f).writerows(rows)

    log(f"[DONE] Wrote tap-level averages to {OUTPUT_FILE}")
    log("=== Finished ===")

if __name__ == "__main__":
    main()


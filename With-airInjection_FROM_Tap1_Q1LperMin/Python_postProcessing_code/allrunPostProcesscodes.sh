#!/bin/bash
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

# Script to run Python post-processing codes in order
# Order: PU.py -> averaged_P_at_Taps.py -> spectrum.py -> Cd_cl_comparison.py -> plot_fft_all_taps.py

echo "Starting post-processing pipeline..."

# 1. Run PU.py to merge U and p data into .dat files
echo "Running PU.py..."
python3 PU.py
if [ $? -ne 0 ]; then
    echo "Error in PU.py"
    exit 1
fi

# 2. Run averaged_P_at_Taps.py to average pressures
echo "Running averaged_P_at_Taps.py..."
python3 averaged_P_at_Taps.py
if [ $? -ne 0 ]; then
    echo "Error in averaged_P_at_Taps.py"
    exit 1
fi

# 3. Run spectrum.py to compute FFT spectra
echo "Running spectrum.py..."
python3 spectrum.py
if [ $? -ne 0 ]; then
    echo "Error in spectrum.py"
    exit 1
fi

# 4. Run Cd_cl_comparison.py for force coefficients
echo "Running Cd_cl_comparison.py..."
python3 Cd_cl_comparison.py
if [ $? -ne 0 ]; then
    echo "Error in Cd_cl_comparison.py"
    exit 1
fi

# 5. Run plot_fft_all_taps.py to plot spectra
echo "Running plot_fft_all_taps.py..."
python3 plot_fft_all_taps.py
if [ $? -ne 0 ]; then
    echo "Error in plot_fft_all_taps.py"
    exit 1
fi

echo "All post-processing scripts completed successfully."

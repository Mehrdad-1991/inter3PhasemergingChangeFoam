################################################################################
# MIT License
# 
# Copyright (c) 2026 - Universität Rostock 
#                       Fakultät für Maschinenbau und Schiffstechnik
#                       Lehrstuhl für Modellierung und Simulation
#                       
# Author: Mehrdad Kazemi
# Date: April 2025 (validated)
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
import os
import csv

def sample_records_by_time(records, time_step=0.001, tol=1e-5):
    """
    Given a sorted list of records (each record is a list [Time, Cd, Cl], with Time as a string),
    this function returns a new list in which only one record per time step (time_step) is kept.
    
    It selects the first record, then for each subsequent record it selects it only if its time 
    is at least (previous_time + time_step) (within a tolerance).
    """
    if not records:
        return []
    sampled = []
    last_time = None
    for rec in records:
        try:
            t_val = float(rec[0])
        except ValueError:
            continue
        if last_time is None:
            sampled.append(rec)
            last_time = t_val
        else:
            if t_val >= last_time + time_step - tol:
                sampled.append(rec)
                last_time = t_val
    return sampled

def read_force_coeffs_data(base_folder):
    """
    This function:
      - Navigates to <base_folder>/postProcessing/forceCoeffs/
      - Lists all subdirectories with numeric names and sorts them in ascending order.
      - In each subdirectory, reads the file 'coefficient.dat' (ignoring lines starting with "#").
      - To avoid duplication, only imports lines with a time greater than the last time from the previous subdirectory.
      - Extracts Time (col 0), Cd (col 1), and Cl (col 4) from each valid line.
      - After aggregating all records, the records are "sampled" so that only one record per 0.001 time step is kept.
      
    Returns the aggregated (and sampled) data records (each record is [Time, Cd, Cl]).
    """
    sim_force_dir = os.path.join(base_folder, "postProcessing", "forceCoeffs")
    if not os.path.isdir(sim_force_dir):
        print(f"Directory '{sim_force_dir}' does not exist.")
        return []

    # List numeric subdirectories.
    subdirs = []
    for name in os.listdir(sim_force_dir):
        full_path = os.path.join(sim_force_dir, name)
        if os.path.isdir(full_path):
            try:
                numeric_val = float(name)
                subdirs.append((numeric_val, name))
            except ValueError:
                continue

    if not subdirs:
        print(f"No numeric subdirectories found in '{sim_force_dir}'.")
        return []

    # Sort subdirectories in ascending order.
    sorted_subdirs = sorted(subdirs, key=lambda x: x[0])
    data_records = []
    prev_last_time = None  # Track the latest time from the previous subdirectory.
    for _, subdir in sorted_subdirs:
        coeff_file = os.path.join(sim_force_dir, subdir, "coefficient.dat")
        if not os.path.exists(coeff_file):
            print(f"Coefficient file '{coeff_file}' does not exist.")
            continue

        with open(coeff_file, "r") as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split()
                if len(parts) < 5:
                    continue
                try:
                    t_val = float(parts[0])
                except ValueError:
                    continue
                # Only include data that is strictly later than the last time from previous subdirectory.
                if prev_last_time is not None and t_val <= prev_last_time:
                    continue
                # Append the record: [Time, Cd, Cl] (adjust indices if necessary).
                data_records.append([parts[0], parts[1], parts[4]])
        # Update prev_last_time if any records were read from this subdirectory.
        if data_records:
            try:
                prev_last_time = float(data_records[-1][0])
            except ValueError:
                pass
    # Sample the aggregated records so that only one record per 0.001 time step is kept.
    sampled_records = sample_records_by_time(data_records, time_step=0.001, tol=1e-5)
    print(f"Collected {len(sampled_records)} records")
    return sampled_records

def compute_cumulative_averages(records):
    """
    Given a list of records (each record is [Time, Cd, Cl] with Cd and Cl as strings),
    return a new list where each record is augmented with two additional columns:
    [Time, Cd, Cl, Cd_averaged, Cl_averaged]. The cumulative average is computed from the start of the data.
    """
    cum_sum_cd = 0.0
    cum_sum_cl = 0.0
    new_records = []
    for i, rec in enumerate(records):
        try:
            cd_val = float(rec[1])
            cl_val = float(rec[2])
        except ValueError:
            cd_val = 0.0
            cl_val = 0.0
        cum_sum_cd += cd_val
        cum_sum_cl += cl_val
        avg_cd = cum_sum_cd / (i + 1)
        avg_cl = cum_sum_cl / (i + 1)
        # Format the averages to 6 decimal places (adjust formatting as needed).
        new_rec = [rec[0], rec[1], rec[2], f"{avg_cd:.6f}", f"{avg_cl:.6f}"]
        new_records.append(new_rec)
    return new_records

def write_force_coeffs_file(base_folder, data_records):
    """
    Writes the force coefficients data into a CSV file ('forceCoeffs.csv')
    with columns: Time, Cd, Cl, Cd_averaged, Cl_averaged
    """
    output_file = os.path.join(base_folder, "forceCoeffs.csv")
    
    # Compute the cumulative averages.
    data_with_averages = compute_cumulative_averages(data_records)
    
    # Write to CSV.
    with open(output_file, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["Time", "Cd", "Cl", "Cd_averaged", "Cl_averaged"])
        for row in data_with_averages:
            writer.writerow(row)
    
    print(f"Output file '{output_file}' has been created with force coefficients.")

def main():
    base_folder = ".."  # Results folder backward of current folder.
    
    data_records = read_force_coeffs_data(base_folder)
    if not data_records:
        print("No force coefficients data was read. Please check the folder names and file paths.")
    else:
        write_force_coeffs_file(base_folder, data_records)

if __name__ == "__main__":
    main()

#!/usr/bin/env python3
# Filename:         bin_aligned_reads.py
# Author:           Gregory Wickham
# Created:          2024-11-14
# Version:          1.1
# Date modified:    2024-11-22
# Description:      Bin reads into 1 Kb bins from tabular alignment files in .parquet format

import os
import argparse
import numpy as np
import pandas as pd

def read_file(file_path):
    """
    Reads a file in Parquet format and converts it to a numpy array.
    """
    if file_path.endswith('.parquet'):
        data = pd.read_parquet(file_path)
    else:
        raise ValueError(f"Unsupported file format: {file_path}")
    
    # Convert to numpy array and ensure required columns are present
    return data[['Region Start', 'Region End', 'Length', 'Midpoint']].to_numpy()

def bin_aligned_reads(datafiles):
    """
    Processes the input data, bins reads by their midpoints, and saves the result to a CSV file.
    """
    longest_read = 0
    for name, file_path in datafiles.items():
        # Load the data for each sample
        data = read_file(file_path)
        
        # Filter out reads > 1000 or < -1000
        data = data[np.abs(data[:, 2]) <= 1000]
        datafiles[name] = data 
        longest_read = max(longest_read, len(data))
        
    # Populate raw matrix with midpoints for each sample
    raw_data = np.full((longest_read, len(datafiles)), np.nan)  # Initialize data matrix
    for i, (name, data) in enumerate(datafiles.items()):
        midpoints = data[:, 3]
        raw_data[:len(midpoints), i] = midpoints
        
    # Define bin size and create bins
    maxbin = np.ceil(np.nanmax(raw_data) / 1000) * 1000
    minbin = np.floor(np.nanmin(raw_data) / 1000) * 1000
    
    bins = np.arange(
        minbin,
        maxbin, 
        1000
    )
    numbins = len(bins)
    
    # Bin the data for each sample
    histtable = np.zeros((numbins - 1, len(datafiles)))     
    for i in range(len(datafiles)):
        hist_counts, _ = np.histogram(raw_data[:, i], bins=bins)
        histtable[:, i] = hist_counts
    
    # Prepare header and data for export
    header = ['Bins'] + list(datafiles.keys())
    binned_data = np.column_stack((bins[:-1], histtable))  # Exclude last bin edge for histtable

    # Convert data to DataFrame
    df = pd.DataFrame(binned_data, columns=header)
    df.to_csv("binned_reads.csv", index=False)

    print("Data has been successfully exported to binned_reads.csv")
    
def main():
    # Set up command-line argument parsing
    parser = argparse.ArgumentParser(
        description='Bin reads data from multiple Parquet files and export to CSV.'
    )
    # Positional argument for input Parquet files
    parser.add_argument(
        'parquet_files', 
        type=str, 
        nargs='+',  # Accept multiple arguments
        help='Paths to one or more Parquet files to process.'
    )

    args = parser.parse_args()
    
    # Ensure that only valid Parquet files are passed
    datafiles = {}
    for file_path in args.parquet_files:
        if file_path.endswith('.parquet'):
            sample_name = os.path.splitext(os.path.basename(file_path))[0]  # Get name without extension
            datafiles[sample_name] = file_path
        else:
            raise ValueError(f"Invalid file format: {file_path}. Only Parquet files are allowed.")
    
    # Call the function to bin aligned reads
    bin_aligned_reads(datafiles)

if __name__ == "__main__":
    main()

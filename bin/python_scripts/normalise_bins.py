#!/usr/bin/env python3
# Filename:         normalise_bins.py
# Author:           Gregory Wickham
# Created:          2024-11-25
# Version:          1.0
# Date modified:    2024-11-25
# Description:      Normalise bins and transform to log2 fold change

import argparse
import numpy as np
import pandas as pd

def normalise_bins(binned_reads):
    """
    Takes binned_reads.csv and normalises by mean and transforms to log2 fold change.
    """
    # Read the input CSV file
    df = pd.read_csv(binned_reads)

    bins = df['Bins']                    # Save the 'bin' column for later use
    data = df.drop(columns=['Bins'])     # Remove the 'bin' column
  
    means = data.mean()                         # Compute the mean for each column
    log2_fold_change = np.log2(data / means)    # Calculate the log2 fold change from the mean for each value
    log2_fold_change.insert(0, 'Bins', bins)     # Add the 'Bins' column back to the transformed data

    # Save the transformed data to a new CSV file
    log2_fold_change.to_csv("log2FC_binned_reads.csv", index=False)
    print("Normalised data saved to log2FC_binned_reads.csv")

def main():
    # Set up command-line argument parsing
    parser = argparse.ArgumentParser(
        description='Normalise bins and transform to log2FC.'
    )
    # Positional argument for binned_reads.csv
    parser.add_argument(
        'binned_reads', 
        type=str, 
        help='Path to the binned_reads.csv file.'
    )

    args = parser.parse_args()

    # Ensure the input file ends with binned_reads.csv
    if not args.binned_reads.endswith('binned_reads.csv'):
        raise ValueError(f"Invalid file: {args.binned_reads}. The file must be named 'binned_reads.csv'.")

    # Call the function to normalise bins
    normalise_bins(args.binned_reads)

if __name__ == "__main__":
    main()

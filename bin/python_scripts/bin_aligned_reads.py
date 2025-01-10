#!/usr/bin/env python3
# Filename:         bin_aligned_reads.py
# Author:           Gregory Wickham
# Created:          2024-11-14
# Version:          1.1
# Date modified:    2024-11-22
# Description:      Bin reads into 1 Kb bins from tabular alignment files in .parquet format and normalises data to log2 fold change

import os
import argparse
import numpy as np
import pandas as pd

class AlignmentProcessor:
    def __init__(self, datafiles):
        self.datafiles = {os.path.splitext(os.path.basename(file))[0]: file for file in datafiles}
        self.longest_read = 0
        self.raw_data = None
        self.bins = None


    def read_file(self, file_path):
        """
        Reads a file in Parquet format and converts it to a numpy array.
        """
        if file_path.endswith('.parquet'):
            data = pd.read_parquet(file_path)
        else:
            raise ValueError(f"Unsupported file format: {file_path}")
        return data[['Region Start', 'Region End', 'Length', 'Midpoint']].to_numpy()


    def bin_aligned_reads(self):
        """
        Processes the input data, bins reads by their midpoints, and saves the result to a CSV file.
        """
        for name, file_path in self.datafiles.items():
            data = self.read_file(file_path)
            
            # Filter out reads > 1000 or < -1000
            data = data[np.abs(data[:, 2]) <= 1000]
            self.datafiles[name] = data 
            self.longest_read = max(self.longest_read, len(data))
            
        # Populate raw matrix with midpoints for each sample
        self.raw_data = np.full((self.longest_read, len(self.datafiles)), np.nan)  # Initialize data matrix
        for i, (name, data) in enumerate(self.datafiles.items()):
            midpoints = data[:, 3]
            self.raw_data[:len(midpoints), i] = midpoints
            
        # Define binwidth and create bins
        maxbin = np.ceil(np.nanmax(self.raw_data) / 1000) * 1000
        minbin = np.floor(np.nanmin(self.raw_data) / 1000) * 1000
        
        self.bins = np.arange(
            minbin,
            maxbin, 
            1000
        )
        numbins = len(self.bins)
        
        # Bin the data for each sample
        histtable = np.zeros((numbins - 1, len(self.datafiles)))     
        for i in range(len(self.datafiles)):
            hist_counts, _ = np.histogram(self.raw_data[:, i], bins=self.bins)
            histtable[:, i] = hist_counts
        
        header = ['Bins'] + list(self.datafiles.keys())
        binned_data = np.column_stack((self.bins[:-1], histtable))

        df = pd.DataFrame(binned_data, columns=header)
        df.to_csv("raw_binned_reads.csv", index=False)

        print("Data has been successfully exported to raw_binned_reads.csv")
        

    def normalise_bins(self, binned_reads):
        """
        Takes binned_reads.csv and normalises by mean and transforms to log2 fold change.
        """
        # Read the input CSV file
        df = pd.read_csv(binned_reads)

        bins = df['Bins']              
        data = df.drop(columns=['Bins']) 
      
        means = data.mean()                         
        log2_fold_change = np.log2(data / means)   
        log2_fold_change.insert(0, 'Bins', bins) 

        # Save the transformed data to a new CSV file
        log2_fold_change.to_csv("log2FC_binned_reads.csv", index=False)
        print("Normalised data saved to log2FC_binned_reads.csv")

    def process(self):
        """
        Main function to process BAM files, bin reads, and normalise the data.
        """
        # Step 1: Bin reads from the given datafiles
        self.bin_aligned_reads()
        
        # Step 2: Normalize the binned reads
        self.normalise_bins("raw_binned_reads.csv")


def main():
    """
    Main function to parse arguments and initiate processing.
    """
    parser = argparse.ArgumentParser(description="Process BAM files, bin reads, and normalise data.")
    
    # Command-line arguments for BAM files and parameters
    parser.add_argument("datafiles", nargs="+", help="Paths to Parquet files to process.")
    
    args = parser.parse_args()
    
    processor = AlignmentProcessor(args.datafiles)
    processor.process()

if __name__ == "__main__":
    main()

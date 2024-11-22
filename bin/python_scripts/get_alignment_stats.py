#!/usr/bin/env python3
# Filename:         get_alignment_stats.py
# Author:           Gregory Wickham
# Created:          2024-11-14
# Version:          1.1
# Date modified:    2024-11-22
# Description:      Tabulate read lengths and positions of sorted indexed .bam files into .parquet format

import sys
import csv
import pysam
import numpy as np
import pandas as pd
import os

def process_bam_file(bam_path):
    """
    Processes a single BAM file to compute alignment stats and writes output to a Parquet file.
    """
    # Read BAM file
    bamfile = pysam.AlignmentFile(bam_path, 'rb')
    data = []
    
    # Fetch reads between pos 1 and 4641652
    for read in bamfile.fetch("U00096.3", 0, 4641652):
        if read.template_length != 0:                       # Get origin and length
            origin = int(read.reference_start) + 1
            read_length = int(read.template_length)
            # Get read start/end pos
            if read_length < 0:                             # Pos if read length is < 0 (3')
                read_start_pos = origin + read_length
                read_end_pos = origin
            else:                                           # Pos if read length is > 0 (5')
                read_start_pos = origin
                read_end_pos = origin + read_length

            # Make origin centric
            if read_start_pos >= 3925875:                   # Right arm ori to '1'
                read_start_pos -= 3925875
                read_end_pos -= 3925875
            elif read_start_pos <= 1590764:                 # Right arm past '1' to dif
                read_start_pos += 715777
                read_end_pos += 715777
            else:                                           # Left arm
                read_start_pos -= 3925875
                read_end_pos -= 3925875     

            # Append data: [first, last, bases, midpoint]
            midpoint = round((read_start_pos + read_end_pos) / 2)
            data.append([read_start_pos, read_end_pos, abs(read_length), midpoint])

    bamfile.close()
    
    # Convert to numpy array and sort by midpoint
    data = np.array(data)
    data = data[data[:, 3].argsort()]
    
    # Remove rows where the first position is zero
    data = data[data[:, 0] != 0].astype(float)
    
    # Write to Parquet
    output_parquet = os.path.splitext(bam_path)[0] + "_readlist.parquet"
    df = pd.DataFrame(data, columns=['Region Start', 'Region End', 'Length', 'Midpoint'])
    df.to_parquet(output_parquet, engine='pyarrow', index=False)
    print(f"Processed {bam_path} -> {output_parquet}")


def main():
    """
    Main function to process multiple BAM files provided as arguments.
    """
    # Check if at least one BAM file is provided
    if len(sys.argv) < 2:
        print("Usage: python get_alignment_stats.py <file1.bam> <file2.bam> ...")
        sys.exit(1)

    # Process each BAM file
    for bam_file in sys.argv[1:]:
        process_bam_file(bam_file)


if __name__ == "__main__":
    main()

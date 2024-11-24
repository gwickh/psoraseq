#!/usr/bin/env python3
# Filename:         get_alignment_stats.py
# Author:           Gregory Wickham
# Created:          2024-11-14
# Version:          1.3
# Date modified:    2024-11-24
# Description:      Tabulate read lengths and positions of sorted indexed .bam files into .parquet format

import sys
import pysam
import numpy as np
import pandas as pd
import os
import argparse

def process_bam_file(bam_path, ori_centric_offset):
    """
    Processes a single BAM file to compute alignment stats and writes output to a Parquet file.
    """
    data = []
    
    # Read BAM file
    bamfile = pysam.AlignmentFile(bam_path, 'rb')
    
    # Fetch reads from the first reference
    for read in bamfile.fetch(bamfile.references[0], 0, bamfile.lengths[0]):
        if read.template_length != 0:  # Get origin and length
            alignment_start = int(read.reference_start) + 1
            read_length = int(read.template_length)
            
            # Get read start/end pos
            if read_length < 0:                             # Pos if read length is < 0 (3')
                print(read_length) 
                read_start_pos = alignment_start + read_length
                read_end_pos = alignment_start
            else:                                           # Pos if read length is > 0 (5')
                read_start_pos = alignment_start
                read_end_pos = alignment_start + read_length

            # Make origin centric
            if ori_centric_offset:
                read_start_pos = ((read_start_pos - ori_centric_offset + (bamfile.lengths[0]/2)) % bamfile.lengths[0]) - (bamfile.lengths[0]/2)
                read_end_pos = ((read_start_pos - ori_centric_offset + (bamfile.lengths[0]/2)) % bamfile.lengths[0]) - (bamfile.lengths[0]/2)

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
    print(f"Processed {bam_path} written to {output_parquet}")


def main():
    """
    Main function to process multiple BAM files provided as arguments.
    """
    # Parse arguments
    parser = argparse.ArgumentParser(description="Tabulate read lengths and positions of sorted indexed BAM files into Parquet format.")
    parser.add_argument(
        "bam_files", 
        nargs="+", 
        help="Paths to one or more BAM files to process."
    )
    parser.add_argument(
        "--make_ori_centric",
        type=int,
        default=0,
        help="Numeric offset to make coordinates origin-centric. Default: 0 (disabled)."
    )
    args = parser.parse_args()

    # Process each BAM file
    for bam_file in args.bam_files:
        process_bam_file(bam_file, args.make_ori_centric)


if __name__ == "__main__":
    main()
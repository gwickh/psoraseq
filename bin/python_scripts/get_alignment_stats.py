import sys
import csv
import pysam
import numpy as np 
import pandas as pd

def get_alignment_stats(name, output_csv):
    #read BAM file as bamfile var
    bamfile = pysam.AlignmentFile(name, 'rb')
    data = []
    
    #fetch reads between pos 1 and 4641652
    for read in bamfile.fetch("U00096.3", 0, 4641652):
        if read.template_length != 0:                       #get origin and length
                origin = int(read.reference_start) + 1
                read_length = int(read.template_length)
                #get read start/end pos
                if read_length < 0:                         #pos if read length is < 0 (3')
                    read_start_pos = origin + read_length
                    read_end_pos = origin
                else:                                       #pos if read length is > 0 (5')
                    read_start_pos = origin
                    read_end_pos = origin + read_length

        #Make origin centric
        if read_start_pos >= 3925875:       # Right arm ori to '1'
            read_start_pos -= 3925875
            read_end_pos -= 3925875
        elif read_start_pos <= 1590764:     # Right arm past '1' to dif
            read_start_pos += 715777
            read_end_pos += 715777
        else:                               #Left arm
            read_start_pos -= 3925875
            read_end_pos -= 3925875     
        
        
        # Append data: [first, last, bases, midpoint]
        midpoint = round((read_start_pos + read_end_pos) / 2)
        data.append(
            [read_start_pos, read_end_pos, abs(read_length), midpoint]
        )
    bamfile.close()
    
    # Convert to np array and sort by midpoint
    data = np.array(data)
    data = data[data[:, 3].argsort()]
    
    # Remove rows where the first position is zero
    data = data[data[:, 0] != 0]
    
    # Convert data to double precision
    data = data.astype(float)
    
    # Write to CSV
    with open(output_csv, 'w', newline='') as csvfile:
        csv_writer = csv.writer(csvfile)
        csv_writer.writerow(['Region Start', 'Region End', 'Length', 'Midpoint'])  # Header row
        csv_writer.writerows(data)
        
if __name__ == "__main__":
    import sys
    get_alignment_stats(sys.argv[1], sys.argv[2])
# psoraseq
Pipeline for psoralen sequencing analysis
```
Installation:
    sudo apt update
    sudo apt-get install -y openjdk-11-jdk build-essential gcc zlib1g-dev make

    wget -qO- https://get.nextflow.io | bash
    sudo mv nextflow /usr/local/bin/
    nextflow -version
    
    git clone https://github.com/gwickh/psoraseq.git

Usage:
    Run the pipeline with alignment:
        nextflow run psoraseq/main.nf --reads_dir <directory_with_reads> --bowtie2_ref <path_to_index> --output_dir <output_directory> --ori_centric_offset <ori_position>

    Skip alignment and run downstream processes:
        nextflow run psoraseq/main.nf --skip_alignment --output_dir <output_directory>

Options:
    --reads_dir             Path to directory containing paired reads.
    --bowtie2_ref           Path to Bowtie2 index or reference .fasta file to make new index
    --output_dir            Path to output directory.
    --ori_centric_offset    (OPTIONAL) Numeric offset to make coordinates origin-centric (default: 0).
    --smoothing_window      (OPTIONAL) Smooth trendline with binsized sliding window (default: 1)
    --skip_alignment        (OPTIONAL) Skip alignment step and use precomputed BAM/BAI files in output directory.
    --help                  Display help message.

Example:
    nextflow run psoraseq/main.nf --reads_dir reads/ --bowtie2_ref e_coli_K12 --output_dir output_test --ori_centric_offset 3925744 --smoothing_window 100
```

#!/usr/bin/env nextflow

// Display help message
if (params.help) {
    println """
    Usage:
        Run the pipeline with alignment:
            nextflow run psoraseq/main.nf --reads_dir <directory_with_reads> --bowtie2_ref <path_to_index> --output_dir <output_directory>

        Skip alignment and run downstream processes:
            nextflow run psoraseq/main.nf --skip_alignment --output_dir <output_directory>
    Options:
        --reads_dir             Path to directory containing paired reads.
        --bowtie2_ref           Path to Bowtie2 index or reference .fasta file to make new index
        --output_dir            Path to output directory.
        --ori_centric_offset    (OPTIONAL) Numeric offset to make coordinates origin-centric (default: 0).
        --skip_alignment        (OPTIONAL) Skip alignment step and use precomputed BAM/BAI files in output directory.
        --help                  Display help message.
    
    Example:
        nextflow run psoraseq/main.nf --reads_dir reads/ --bowtie2_ref e_coli_K12 --output_dir output_test
    """
    exit 0
}

// bowtie2
process BOWTIE2 {
    conda "${projectDir}/environment.yml"

    input:
    path reads_dir
    path bowtie2_ref
    path output_dir

    output:
    path "${bowtie2_ref.simpleName}.*.bt2", emit: index_files
    path "*.bam", emit: bam_file
    path "*.bai", emit: bam_index_file

    publishDir file(params.bowtie2_ref).parent.toString(), 
        mode: 'copy', 
        pattern: "*.bt2"
    publishDir file(params.output_dir).toString(), 
        mode: 'copy',
        pattern: "*.bam"
    publishDir file(params.output_dir).toString(), 
        mode: 'copy',
        pattern: "*.bai"

    script:
    """
    bash ${projectDir}/bin/shell_scripts/bowtie2_alignment.sh \
        $reads_dir $bowtie2_ref $output_dir
    """
}

// tabulate aligned reads
process ALIGNMENT_STATS {
    conda "${projectDir}/environment.yml"
    
    input:
    path bam_file
    path bam_index_file

    output:
    path "*.parquet", emit: readlist 

    script:
    """
    python3 ${projectDir}/bin/python_scripts/get_alignment_stats.py $bam_file --make_ori_centric $params.ori_centric_offset
    """
}

// bin reads
process BIN_READS {
    conda "${projectDir}/environment.yml"

    input:
    path readlist

    output:
    path "*.csv", emit: binned

    publishDir file(params.output_dir).toString(), 
        mode: 'copy',
        pattern: "*.csv"

    script:
    """
    python3 ${projectDir}/bin/python_scripts/bin_aligned_reads.py $readlist
    """
}

// normalise reads
process NORMALISE_READS {
    conda "${projectDir}/environment.yml"

    input:
    path binned

    output:
    path "*.csv", emit: binned_norm

    publishDir file(params.output_dir).toString(), 
        mode: 'copy',
        pattern: "*.csv"

    script:
    """
    python3 ${projectDir}/bin/python_scripts/normalise_bins.py $binned
    """
}

// define workflow
workflow {
    if (!params.skip_alignment) {
        bam_output = BOWTIE2(
            file(params.reads_dir), 
            file(params.bowtie2_ref), 
            file(params.output_dir)
        )

        tabular_alignment = ALIGNMENT_STATS(
            bam_output.bam_file,
            bam_output.bam_index_file
        )
    } else {
        // Assume BAM and BAI files are already in output_dir
        bam_files = file("${params.output_dir}/*.bam")
        bai_files = file("${params.output_dir}/*.bai")
        
        tabular_alignment = ALIGNMENT_STATS(
            bam_files,
            bai_files
        )
    }

    binning = BIN_READS(
        tabular_alignment.readlist
    )

    NORMALISE_READS(
        binning.binned
    )

}
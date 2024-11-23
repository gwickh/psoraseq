#!/usr/bin/env nextflow

// params
params.reads_dir = null
params.bowtie2_ref = null
params.output_dir = null
params.skip_alignment = false  // Default: run alignment

// bowtie2
process BOWTIE2 {
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
    input:
    path bam_file
    path bam_index_file

    output:
    path "*.parquet", emit: readlist 

    script:
    """
    python3 ${projectDir}/bin/python_scripts/get_alignment_stats.py $bam_file
    """
}

// bin reads
process BIN_READS {
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

    BIN_READS(
        tabular_alignment.readlist
    )
}
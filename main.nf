#!/usr/bin/env nextflow

// params
params.reads_dir = null
params.bowtie2_ref = null
params.output_dir = null
params.output_csv = 'binned_data.csv'

// bowtie2
process BOWTIE2 {
    input:
    path reads_dir
    path bowtie2_ref
    path output_dir

    output:
    path "${bowtie2_ref.simpleName}.*.bt2", emit: index_files
    path "*.bam*", emit: bam_file

    publishDir file(params.bowtie2_ref).parent.toString(), 
        mode: 'copy', 
        pattern: "*.bt2"
    publishDir file(params.output_dir).toString(), 
        mode: 'copy',
        pattern: "*.bam*"

    script:
    """
    bash ~/webber_group/Gregory_Wickham/psoraseq/psoraseq/bin/shell_scripts/bowtie2_alignment.sh \
        $reads_dir $bowtie2_ref $output_dir
    """
}

// bin aligned reads
process ALIGNMENT_STATS {
    input:
    path bam_file
    path output_dir

    output:
    path "${params.output_csv}"

    publishDir file(params.output_dir).toString(), 
        mode: 'copy',
        pattern: "*.csv"

    script:
    """
    python3 ~/webber_group/Gregory_Wickham/psoraseq/psoraseq/bin/python_scripts/get_alignment_stats.py \
        $bam_file ${params.output_csv}
    """
}

workflow {
    bam_output = BOWTIE2(file(params.reads_dir), file(params.bowtie2_ref), file(params.output_dir))
    ALIGNMENT_STATS(bam_output.bam_file, file(params.output_dir))
}
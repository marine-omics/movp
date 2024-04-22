
process name_by_sample {

    input:
    tuple val(sample), path(bam)

    output:
    tuple val(sample), path ("${sample}.bam"), emit: merged_bam

    script:
    """
    cp "${bam}" "${sample}.bam"
    """
}
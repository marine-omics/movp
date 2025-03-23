

process faidx {

    input:
    path fasta

    output:
    path ("*.fai"), emit: fai

    script:
    """
    samtools faidx $fasta
    """
}

process sidx {

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.bai")

    script:
    """
    samtools index $bam
    """
}


process flagstat {

    publishDir "$params.outdir/bamstats", mode: 'copy'  

    input:
    tuple val(sample), path(bam), path(bai)

    output:
    path ("*.flagstat"), emit: flagstat

    script:

    def outfile = "${sample}.flagstat"

    """
    samtools flagstat $bam > $outfile
    """
}

process stat {

    publishDir "$params.outdir/bamstats", mode: 'copy'

    input:
    tuple val(sample), path(bam), path(bai)

    output:
    path ("*.stat"), emit: stat

    script:

    def outfile = "${sample}.stat"

    """
    samtools stat $bam > $outfile
    """
}

process idxstat {

    publishDir "$params.outdir/bamstats", mode: 'copy'  

    input:
    tuple val(sample), path(bam), path(bai)

    output:
    path ("*.idxstat"), emit: idxstat

    script:

    def outfile = "${sample}.idxstat"

    """
    samtools idxstats $bam > $outfile
    """
}

process samtools_merge {

    input:
    tuple val(sample), path(bam)

    output:
    tuple val(sample), path ("${sample}.bam"), emit: merged_bam

    script:

    def outfile = "${sample}.bam"

    """
    samtools merge -O BAM $outfile $bam 
    """
}
